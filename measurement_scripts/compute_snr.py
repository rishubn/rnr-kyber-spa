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
from nttwrapper import KYBERQ


def parse_args():
    argparse = ArgumentParser()
    argparse.add_argument("--trace-set", type=Path)
    argparse.add_argument("--snr", type=Path)
    argparse.add_argument("--nc", type=int)
    argparse.add_argument("--coordinate", type=str)
    argparse.add_argument("--plot", action="store_true")
    return argparse.parse_args()

    # for npz in [np.load(f) for f in trace_path.glob("*.npz")]:
    #     traces = npz["traces"]
    #     labels = npz["labels"]
    #     qtraces = npz["qtraces"]
    #     snr.fit_u(
    #         qtraces, labels[:, var_num[0], var_num[1]][:, np.newaxis].astype(np.uint16)
    #     )
    # return var_num, snr.get_snr()


def main():
    args = parse_args()
    i, j = map(int, args.coordinate.split("_"))
    compressor = Blosc(cname="lz4")
    logger.add(f"{args.snr}/snr.log", level="INFO")
    # snrs = [[SNR(args.nc) for _ in range(9)] for _ in range(256)]
    with zarr.open(args.trace_set, mode="r") as trace_set, zarr.open(
        args.snr / "snr.zarr", mode="a"
    ) as snr_set:
        sync = zarr.ProcessSynchronizer(f"{args.snr}/snr.sync")
        if "SNR" not in snr_set:
            snrs_results = snr_set.zeros(
                "SNR",
                shape=(256, 9, trace_set.qtraces.shape[1]),
                chunks=(1, 1, None),
                dtype=np.float64,
                compressor=compressor,
                overwrite=True,
                synchronizer=sync,
            )
        else:
            snrs_results = snr_set["SNR"]

        logger.info("Starting")
        snr = SNR(nc=args.nc)
        chunk_size = trace_set.qtraces.chunks[0]
        num_chunks = trace_set.qtraces.shape[0] // chunk_size

        for k in tqdm(range(num_chunks)):
            snr.fit_u(
                trace_set.qtraces[chunk_size * k : (k + 1) * chunk_size],
                trace_set.labels[chunk_size * k : (k + 1) * chunk_size, i, j][
                    :, np.newaxis
                ]
                % KYBERQ,
            )
        snrs_results[i, j] = snr.get_snr()[0]

        logger.info(f"Done {i},{j}")
        Path(args.snr / f"{args.coordinate}.done").touch()
        if args.plot:
            plt.figure()
            plt.plot(snrs_results[i, j])
            plt.savefig(args.snr / f"plots/snr_{i}_{j}.png")
            plt.close()
    # with zarr.open(args.trace_set, mode="r+") as trace_set, ContextExecutor(
    #     max_workers=1
    # ) as pool:
    #     snrs_results = trace_set.zeros(
    #         "SNRs",
    #         shape=(256, 9, trace_set.traces.shape[1]),
    #         chunks=(1, 1, None),
    #         dtype=np.float64,
    #         compressor=compressor,
    #         overwrite=True,
    #     )

    #     pbar = tqdm(total=9 * 256)
    #     for i in range(256):
    #         for j in range(9):
    #             futures.append(pool.submit(compute_snr, snrs[i][j], (i, j), trace_set))
    #     for f in as_completed(futures):
    #         i, j = f.result()

    #         snrs_results[i, j] = snrs[i][j].get_snr()[0]
    #         plt.figure()
    #         plt.plot(snrs_results[i, j])
    #         plt.savefig(args.snr_path / f"snr_{i}_{j}.png")
    #         del snrs[i][j]
    #         pbar.update(1)
    # futures.append(pool.submit())
    # with ContextExecutor(max_workers=64) as pool:
    #     for i in range(1):
    #         for j in range(1):
    #             futures.append(pool.submit(compute_snr, (i, 6), args.trace_path))

    # for f in as_completed(futures):
    #     pbar.update(1)
    #     coor, snrv = f.result()
    #     np.savez_compressed(
    #         args.snr_path / f"snr_{coor[0]}_{coor[1]}.npz", coordinate=coor, snr=snrv
    #     )
    #     plt.figure()
    #     plt.plot(snrv[0])
    #     plt.savefig(args.snr_path / f"snr_{coor[0]}_{coor[1]}.png")

    # snr = SNR(nc=3329)

    # for npz in tqdm([np.load(f) for f in args.trace_path.glob("*.npz")]):
    #     traces = npz["traces"]
    #     labels = npz["labels"]
    #     qtraces = Quantizer.fit(traces).quantize(traces)
    #     snr.fit_u(qtraces, labels[:, 0:1, 0].astype(np.uint16))
    # s = snr.get_snr()[0]
    # print(s)


if __name__ == "__main__":

    main()
