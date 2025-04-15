from argparse import ArgumentParser
from pathlib import Path
from nttwrapper import INTTImpl, INTT, KYBERQ, WIDE_S_KYBERQ, WIDE_U_KYBERQ, ZETAS
from scalib.attacks import FactorGraph, GenFactor, BPState
from os.path import exists
import numpy as np
import zarr
from profiling import LDAModel
import pickle
import matplotlib.pyplot as plt
from tqdm import tqdm
from loguru import logger
import json

def parse_args():
    argparse = ArgumentParser()
    argparse.add_argument("--attack-set", type=Path)
    argparse.add_argument("--lda-path", type=Path)
    argparse.add_argument("--nc", type=int)
    argparse.add_argument("--fg", type=Path)
    argparse.add_argument("--num-attack-traces", type=int, default=1000)
    argparse.add_argument(
        "--intt-type",
        type=lambda c: INTTImpl[c.upper()],
        choices=list(INTTImpl),
        required=True,
    )
    argparse.add_argument("--plot-priors", action="store_true")
    argparse.add_argument("--outpath", type=Path)
    argparse.add_argument("--it-path", type=Path)
    argparse.add_argument("--attack-num", type=int)
    argparse.add_argument("--pi-threshold", type=float, default=3.0)
    return argparse.parse_args()


def _inv_butterfly_s(nc):
    if not exists(f"./factor_cache/inv_butterfly_s.npy"):
        bff = []
        for a in range(nc):
            for b in range(nc):
                bff.append([a, b, ((a + b) % KYBERQ), ((b - a) % KYBERQ)])
        x = np.array(bff, dtype=np.uint32)
        np.save(f"./factor_cache/inv_butterfly_s.npy", x)
        return GenFactor.sparse_functional(x)
    else:
        x = np.load(f"./factor_cache/inv_butterfly_s.npy")
        return GenFactor.sparse_functional(x)


def plot_distri(arr, domain, outpath, title="", axvline=None):

    plt.figure()
    plt.title(title)
    plt.stem(arr)
    if axvline:
        plt.axvline(axvline, color="r")
    plt.savefig(outpath / f"{title}.png", dpi=300)
    plt.close()


def get_prior(lda_model, traces):
    p = lda_model.lda.predict_proba(traces[:, lda_model.pois])
    p = np.where(p == 0, 1e-40, p)
    log_prior = np.sum(np.log2(p), axis=0)
    prior = (log_prior - np.min(log_prior)) / (np.max(log_prior) - np.min(log_prior))
    prior = prior/prior.sum()
    return np.where(prior == 0, 1e-40, prior)
def get_prior2(lda_model, traces):
    n = traces.shape[0]
    distris = lda_model.lda.predict_proba(traces[:, lda_model.pois])
    prior = distris[0]
    i = 1
    while i < n:
        prior *= distris[i]
        prior = prior/prior.sum()
        i += 1
    return np.where(prior == 0, 1e-40, prior)


def compute_ge(distri, labels):
    ge = np.log2(np.where(labels == np.argsort(distri)[::-1])[0] + 1)
    return ge


def get_zetas():
    z = ZETAS.copy()
    z[ZETAS < 0] = ZETAS[ZETAS < 0] + KYBERQ
    return z


def attack(priors, fg, num_iter=20):
    zetas = get_zetas()
    factor_graph = FactorGraph(fg.read_text())
    bp = BPState(
        factor_graph,
        1,
        public_values={f"zeta{i}": zetas[i] * 169 for i in range(1, 256 // 2)},
        gen_factors={"_INV_BUTTERFLY": _inv_butterfly_s(nc=KYBERQ)},
    )
    for i in range(256):
        for j in range(8):

            bp.set_evidence(f"v{i}_{j}", priors[i, j])
    bp.bp_loopy(num_iter, initialize_states=True)
    results = np.zeros((256, 8, KYBERQ))
    for i in range(256):
        for j in range(8):
            results[
                i,
                j,
            ] = bp.get_distribution(f"v{i}_{j}")
    return results


def main():
    args = parse_args()
    priors = np.zeros((256, 8, args.nc))
    pbar = tqdm(total=256 * 8)
    with zarr.open(args.attack_set, mode="r") as attack_set, zarr.open(
        args.outpath / "results.zarr", mode="a"
    ) as result_set:
        for i in range(256):
            for j in range(8):
                if exists(args.lda_path / f"{i}_{j}.pkl"):
                    with open(args.it_path / f"{i}_{j}.json") as f:
                        d = json.load(f)
                        if d["PI"]["mean"] < 0.0:
                            p = np.ones((KYBERQ,))
                            priors[i,j] = p/p.sum()
                        else:
                            model = pickle.load(open(args.lda_path / f"{i}_{j}.pkl", "rb"))
                            priors[i, j] = get_prior2(
                                model, attack_set.qtraces[0 : args.num_attack_traces]
                            )
                    logger.info(priors[i,j])
                else:
                    priors[i, j, 0] = 1.0

                if args.plot_priors:
                    plot_distri(
                        priors[i, j],
                        np.arange(args.nc),
                        outpath=args.outpath / "plots",
                        title=f"v{i}_{j}_prior",
                        axvline=attack_set.labels[0, i, j] % 3329,
                    )
                pbar.update(1)
        prior_res = result_set["priors"]
        for j in range(32):
            ge = compute_ge(priors[j, 0], labels=attack_set.labels[0, j, 0] % 3329)
            logger.info(ge)
        prior_res[args.attack_num] = priors
        results = attack(priors, args.fg)
        post_res = result_set["posteriors"]
        post_res[args.attack_num] = results
        logger.info("DONE")
        for j in range(32):
            post_ge = compute_ge(results[j, 0], labels=attack_set.labels[0, j, 0] % 3329)
            logger.info(post_ge)
    # with zarr.open(args.attack_set, mode="r") as attack_set, open(args.rlda_path, mode='rb') as f:

    #     model = pickle.load(f)

    # plot_distri(
    #                 pr[0],
    #                 np.arange(3329),
    #                 outpath="pr.png",
    #                 title=f"v0_0_prior",
    #                 axvline=attack_set.labels[0, 0, 4],
    #             )

    # for i in range(256):
    #     for j in range(7):
    #         if exists(args.rlda_path / f"{i}_{j}.pkl"):
    #             model = pickle.load(
    #                 open(
    #                     args.rlda_path / f"{i}_{j}.pkl",
    #                     "rb",
    #                 )
    #             )
    #             pr = model.rlda.predict_proba(
    #                 attack_set.qtraces[0 : args.num_attack_traces, model.pois]
    #             )
    #             print(pr)
    #             # pr = pr[:, 0 : args.nc]
    #         else:
    #             pr = np.zeros((1, args.nc))
    #             pr[:, 0] = 1.0
    #         if args.plot_priors:
    #             plot_distri(
    #                 pr[0],
    #                 np.arange(4096),
    #                 outpath=args.outpath / "plots",
    #                 title=f"v{i}_{j}_prior",
    #                 axvline=attack_set.labels[0, i, j],
    #             )
    #             print(attack_set.labels[0, i, j])


if __name__ == "__main__":
    main()
