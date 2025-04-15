from argparse import ArgumentParser
from scalib.metrics import SNR
from scalib.preprocessing import Quantizer
import numpy as np
import matplotlib.pyplot as plt
from scalib.modeling import RLDAClassifier, LDAClassifier
from scalib.metrics import RLDAInformationEstimator
from pathlib import Path
import sys
from tqdm import tqdm
from scalib.tools import ContextExecutor
from concurrent.futures import as_completed
import pickle
import json
import time
import zarr
from math import ceil, log2
from loguru import logger
from dataclasses import dataclass
import numpy.typing as npt
from nttwrapper import KYBERQ


@dataclass
class LDAModel:
    coordinate: str
    pois: npt.ArrayLike
    p: int
    lda: LDAClassifier


def parse_args():
    argparse = ArgumentParser()
    argparse.add_argument("--train-set", type=Path)
    argparse.add_argument("--snr", type=Path)
    argparse.add_argument("--lda-path", type=Path)
    argparse.add_argument("--p")
    argparse.add_argument("--best-model", type=Path)
    argparse.add_argument("--nc", type=int)
    argparse.add_argument("--snr-threshold", type=float)
    argparse.add_argument("--coordinate", type=str)
    argparse.add_argument("--outfile", type=Path)
    return argparse.parse_args()


def get_pois(snr, t):
    x = np.where(snr >= t)[0]
    y = snr[x]
    s = np.argsort(y)[::-1]
    pois = x[s]
    if len(pois) > 1000:
        return pois[:1000]
    else:
        return pois


def main():
    args = parse_args()
    with zarr.open(args.train_set, mode="r") as trace_set, zarr.open(
        args.snr, mode="r"
    ) as snr_set:
        logger.add(f"{args.outfile.parent / args.outfile.stem}.log", level="INFO")
        i, j = map(int, args.coordinate.split("_"))
        pois = get_pois(snr_set.SNR[i, j], args.snr_threshold)
        print(f"POIs {pois}\nNum:{len(pois)}")
        if args.best_model:
            with open(args.best_model, "r") as bmf:
                bm = json.load(bmf)
            p = bm[f"{j}"]
        else:
            p = int(args.p)
        if len(pois) < p:
            p = len(pois)
        logger.info(f"Dims: {p}")
        lda = LDAClassifier(args.nc, p)
        chunk_size = trace_set.qtraces.chunks[0]
        num_chunks = trace_set.qtraces.shape[0] // chunk_size
        logger.info("Starting")
        pbar = tqdm(total=num_chunks)
        arr_chunk = np.arange(num_chunks)
        np.random.shuffle(arr_chunk)
        for k in arr_chunk:
            assert (
                trace_set.qtraces[chunk_size * k : (k + 1) * chunk_size, pois] > -32768
            ).all(), f"{k}\n{np.where(trace_set.qtraces[chunk_size * k : (k + 1) * chunk_size, pois] <= -32768)}"
            assert (
                trace_set.qtraces[chunk_size * k : (k + 1) * chunk_size, pois] < 32768
            ).all(), f"{k}\n{np.where(trace_set.qtraces[chunk_size * k : (k + 1) * chunk_size, pois] >= 32768)}"
            lda.fit_u(
                trace_set.qtraces[chunk_size * k : (k + 1) * chunk_size, pois],
                trace_set.labels[
                    chunk_size * k : (k + 1) * chunk_size,
                    i,
                    j,
                ]
                % KYBERQ,
            )
            pbar.update(1)
        lda.solve(done=True)
    with open(args.outfile, "wb") as f:
        pickle.dump(
            LDAModel(
                **{
                    "coordinate": args.coordinate,
                    "pois": pois,
                    "p": p,
                    "lda": lda,
                }
            ),
            f,
        )
        logger.info("Done")
    # args = parse_args()
    # futures = []
    # pbar = tqdm(total=9 * 256)
    # with ContextExecutor(max_workers=1) as pool:
    #     for i in range(1):
    #         for j in range(1):
    #             snrnpz = np.load(args.snr_path / f"snr_{i}_6.npz")
    #             var_num = snrnpz["coordinate"]
    #             snr = snrnpz["snr"][0]
    #             futures.append(
    #                 pool.submit(
    #                     compute_rlda,
    #                     var_num,
    #                     args.p,
    #                     args.num_pois,
    #                     snr,
    #                     args.trace_path,
    #                 )
    #             )
    # for f in as_completed(futures):
    #     pbar.update(1)
    #     coordinate, pois, rlda = f.result()
    #     with open(
    #         args.rlda_path / f"rlda_{coordinate[0]}_{coordinate[1]}.pkl", "wb"
    #     ) as f:
    #         pickle.dump({"coordinate": coordinate, "pois": pois, "rlda": rlda}, f)

    # with open(args.rlda_path / "rlda_stats.json", "w") as f:
    #     json.dump(
    #         {
    #             "elapsed_time": time.strftime(
    #                 "%H:%M:%S", time.gmtime(pbar.format_dict["elapsed"])
    #             ),
    #             "num_pois": args.num_pois,
    #             "num_dims": args.p,
    #         },
    #         f,
    #         indent=4,
    #     )
    # for npz in tqdm([np.load(f) for f in args.ptrace_path.glob("*.npz")]):
    #     traces = npz["traces"]
    #     labels = npz["labels"]
    #     qtraces = Quantizer.fit(traces).quantize(traces)
    #     l = labels[:, 0:1, 6].astype(np.uint16)
    #     snr.fit_u(qtraces, l)

    # s = snr.get_snr()[0]
    # pois = s.argsort()[-250:][::-1][1:]
    # plt.plot(s)
    # plt.plot(pois, s[pois], "ro")
    # plt.savefig("snr.png", dpi=1200)
    # for npz in tqdm([np.load(f) for f in args.ptrace_path.glob("*.npz")]):
    #     traces = npz["traces"]
    #     labels = npz["labels"]
    #     qtraces = Quantizer.fit(traces).quantize(traces)
    #     l = labels[:, 0:1, 6].astype(np.uint16)
    #     rlda.fit_u(qtraces[:, pois], l.astype(np.uint64))

    # rlda.solve()

    # cl = rlda.get_clustered_model(0, 0.01, 10000000, False)
    # it = RLDAInformationEstimator(cl, 0)
    # for npz in tqdm([np.load(f) for f in args.atrace_path.glob("*.npz")]):
    #     traces = npz["traces"]
    #     labels = npz["labels"]
    #     qtraces = Quantizer.fit(traces).quantize(traces)
    #     al = labels[:, 0:1, 6].astype(np.uint64)
    #     it.fit_u(qtraces[:, pois], al[:, 0])
    # # attack_npz = np.load(args.attack_traces)
    # # atraces = attack_npz["traces"]
    # # alabels = attack_npz["labels"]
    # # quantizer = Quantizer.fit(atraces)
    # # qatraces = quantizer.quantize(atraces)

    # # al = alabels[:, 0:1, 0].astype(np.uint64)
    # # it.fit_u(qatraces[:, pois], al[:, 0])

    # pi_l, pi_u = it.get_information()
    # print(f"{pi_l} < PI < {pi_u}")


if __name__ == "__main__":
    main()
