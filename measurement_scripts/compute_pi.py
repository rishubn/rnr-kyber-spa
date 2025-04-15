import numpy as np
from pathlib import Path
from argparse import ArgumentParser
from scalib.metrics import SNR
from tqdm import tqdm
from scalib.preprocessing import Quantizer
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor
from scalib.tools import ContextExecutor
import time
import json
import matplotlib.pyplot as plt
import zarr
from numcodecs import Blosc
import os
from tqdm.contrib import itertools
from loguru import logger
import pickle
from profiling import LDAModel
import matplotlib.pyplot as plt
from nttwrapper import KYBERQ


def parse_args():
    argparse = ArgumentParser()
    argparse.add_argument("--train-set", type=Path)
    argparse.add_argument("--test-set", type=Path)
    argparse.add_argument("--lda", type=Path)
    argparse.add_argument("--nc", type=int)
    argparse.add_argument("--outfile", type=Path)
    argparse.add_argument("--coordinate", type=str)
    return argparse.parse_args()


def get_information(trace_set, model, ent, coordinate):
    with zarr.open(trace_set, mode="r") as trace_set, open(model, "rb") as f:
        model = pickle.load(f)
        lda = model.lda
        pois = model.pois
        i, j = map(int, coordinate.split("_"))
        chunk_size = trace_set.qtraces.chunks[0]
        num_chunks = trace_set.qtraces.shape[0] // chunk_size
        pbar = tqdm(total=num_chunks)
        sumx = 0.0
        sumsq = 0.0

        for k in range(num_chunks):
            pr_X_L = np.log2(
                lda.predict_proba(
                    trace_set.qtraces[chunk_size * k : (k + 1) * chunk_size, pois]
                )[
                    np.arange(chunk_size),
                    trace_set.labels[chunk_size * k : (k + 1) * chunk_size, i, j]
                    % KYBERQ,
                ]
            )
            if np.any(
                pr_X_L == -np.inf
                ):
                continue
            sumx += pr_X_L.sum()
            sumsq += (pr_X_L**2).sum()
            current_num = (k+1)*chunk_size
            it = ent + (1.0 / current_num) * sumx
            dev = (sumsq / current_num) - (
                sumx / current_num
            ) ** 2
            ci = 2 * dev / np.sqrt(current_num)
            logger.info(f"mean {it} std {dev} ci {ci}")
            if ci < 0.1: 
                return {"mean": it, "std_dev": dev, "ci": ci}

            pbar.update(1)
        it = ent + (1.0 / trace_set.qtraces.shape[0]) * sumx
        dev = (sumsq / trace_set.qtraces.shape[0]) - (
            sumx / trace_set.qtraces.shape[0]
        ) ** 2
        ci = 2 * dev / np.sqrt(trace_set.qtraces.shape[0])
        return {"mean": it, "std_dev": dev, "ci": ci}


def main():

    args = parse_args()
    H_X = np.log2(args.nc)
    pi = get_information(args.test_set, args.lda, H_X, args.coordinate)
    logger.info(f'PI: {pi["mean"]} +/- {pi["ci"]}, std_dev: {pi["std_dev"]}')
    ti = get_information(args.train_set, args.lda, H_X, args.coordinate)
    logger.info(f'TI: {ti["mean"]} +/- {ti["ci"]}, std_dev: {ti["std_dev"]}')
    with open(args.outfile, "w") as f:
        json.dump(
            {"PI": pi, "TI": ti},
            f,
        )


if __name__ == "__main__":
    main()
