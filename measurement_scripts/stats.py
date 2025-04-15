import zarr
from argparse import ArgumentParser
from pathlib import Path
import pandas as pd
import numpy as np
from tqdm import tqdm

def parse_args():
    argparse = ArgumentParser()
    argparse.add_argument("--attack-set-prefix", type=Path)
    argparse.add_argument("--result-path", type=Path)
    argparse.add_argument("--outfile", type=Path)
    # argparse.add_argument("--num-attacks", type=int)
    return argparse.parse_args()


def GE(distri, label):
    return np.log2(np.where((label % 3329) == np.argsort(distri)[::-1])[0] + 1)[0]


def main():
    args = parse_args()
    pbar = tqdm(total=25*256*8)
    stats = {"num_traces":[], "attack_num": [], "var": [], "layer":[], "prior_ge": [], "posterior_ge": []}
    for p in args.result_path.glob("**/results.zarr"):
        num_traces = int(p.parent.stem)
        with zarr.open(p, mode="r") as result_set:
            for j in range(result_set.posteriors.shape[0]):  # num attacks

                with zarr.open(
                    f"{args.attack_set_prefix}{j}/trace_set.zarr", mode="r"
                ) as attack_set:

                    for v in range(result_set.posteriors.shape[1]):  # num vars
                        for l in range(result_set.posteriors.shape[2]):
                            stats["num_traces"].append(num_traces)
                            stats["attack_num"].append(j)
                            stats["var"].append(v)
                            stats["layer"].append(l)
                            stats["prior_ge"].append(
                                GE(result_set.priors[j, v,l], attack_set.labels[0, v, l])
                            )
                            stats["posterior_ge"].append(
                                GE(result_set.posteriors[j, v,l], attack_set.labels[0, v, l])
                            )
                            pbar.update(1)
        df = pd.DataFrame(stats)
        df.to_csv(args.outfile, index=False)


if __name__ == "__main__":
    main()
