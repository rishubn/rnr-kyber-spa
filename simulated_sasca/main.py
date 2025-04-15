from __future__ import annotations
import sys
from argparse import ArgumentParser
from pathlib import Path
from loguru import logger
from intt import INTT, SparseInputGenerator, INTTImpl
from consts import Domain, REDUNDANCY
from attack import BPAttack
from leakage_models import LeakageModel
import numpy as np
from tqdm import tqdm
from stats import Stats
from plots import plot_distri
import numpy.typing as npt
from time import time
import random
from utils import Artifact


def get_args():
    parser = ArgumentParser()
    parser.add_argument("--fg", type=Path, help="Factor graph file")
    parser.add_argument("--iterations", type=int, help="BP Iterations")
    parser.add_argument(
        "--num-expr", type=int, default=1, help="Number of simul. experiments (threads)"
    )
    parser.add_argument("--sigma", type=float, default=0.0001, help="std.dev. of AWGN")
    parser.add_argument(
        "--load-artifact",
        type=Artifact.load_artifact,
        default=None,
        help="Load artifact file",
    )
    parser.add_argument(
        "--save-artifact", type=str, default="default", help="Save path to artifact"
    )
    parser.add_argument(
        "--non-zero",
        type=int,
        default=32,
        help="Number of non-zero coeffs in INTT input",
    )
    parser.add_argument(
        "--word-size",
        type=int,
        default=16,
        help="Size of register (for leakage generation)",
    )
    parser.add_argument(
        "--gen-graph", action="store_true", help="Generate SASCA graph for INTT"
    )
    parser.add_argument("--height", type=int, default=256)
    parser.add_argument("--layers", type=int, default=7)
    parser.add_argument("--seed", type=int, default=int(time()))
    parser.add_argument(
        "--intt-impl",
        type=lambda c: INTTImpl[c.upper()],
        choices=list(INTTImpl),
    )
    return parser.parse_args()


def attack(
    fg: str, num_expr: int, priors: npt.NDArray, num_iter, height=256
) -> npt.NDArray:
    attack = BPAttack(fg, height)
    progress = tqdm(range(num_expr))
    return attack.attack_batch(priors, progress=progress, num_iter=num_iter)


def main():
    args = get_args()
    random.seed(args.seed)
    np.random.seed(args.seed)
    logger.debug(args)
    if args.intt_impl == INTTImpl.SIGNED or args.intt_impl == INTTImpl.UNSIGNED:
        assert (
            REDUNDANCY == 1
        ), "REDUNDANCY must be set to 1 for select INTT implementation"
    elif args.intt_impl == INTTImpl.WIDE_UNSIGNED:
        assert REDUNDANCY == 5, "REDUNDANCY must equal 5 for WIDE_UNSIGNED"
    elif args.intt_impl == INTTImpl.WIDE_SIGNED:
        assert REDUNDANCY == 9, "REDUNDANCY must equal 9 for WIDE_SIGNED"
    if args.gen_graph:
        print(INTT.make_graph(height=args.height, layers=args.layers))
        sys.exit(0)
    if args.load_artifact is None:
        intt = INTT(
            SparseInputGenerator(
                domain=args.intt_impl.to_domain(), non_zero=args.non_zero
            ),
            inttimpl=args.intt_impl,
            height=args.height,
            layers=args.layers,
        )

        model = LeakageModel(
            args.intt_impl.to_leakage_function(), word_size=args.word_size
        )

        values = intt.simulate(args.num_expr)
        leakage = model.generate_leakage(
            values,
            domain=args.intt_impl.to_domain(),
            sigma=args.sigma,
            non_zero_sparse=args.non_zero,
        )
        artifact = args.load_artifact

        fg = args.fg.read_text()

        results = attack(fg, args.num_expr, leakage, args.iterations, args.height)
        artifact = Artifact(
            values,
            args.sigma,
            args.num_expr,
            args.iterations,
            args.non_zero,
            args.word_size,
            leakage,
            results,
            args.intt_impl,
        )
        Artifact.save_artifact(args.save_artifact, artifact)
        plot_distri(
            artifact.attack_results[0, 0],
            Domain.KYBER_0_Q().to_range(),
            axvline=artifact.intt_values[0, 0, 0],
        )
        stats = Stats(args.intt_impl, args.non_zero)
        print(
            stats.avg_rank(
                artifact.attack_results,
                artifact.intt_values[:, :, 0],
            )
        )
        print(
            stats.success_rate(
                artifact.attack_results,
                artifact.intt_values[:, :, 0],
            )
        )
    else:
        artifact = args.load_artifact
        print(artifact)
        stats = Stats(artifact.intt_impl, args.non_zero)
        plot_distri(
            artifact.attack_results[0, 0],
            Domain.KYBER_0_Q().to_range(),
            axvline=artifact.intt_values[0, 0, 0],
        )
        print(
            stats.avg_rank(
                artifact.attack_results,
                artifact.intt_values[:, :, 0],
            )
        )


if __name__ == "__main__":
    main()
    # small_attack()
#  results = attack.attack_batch(priors, progress)
#  np.save("results", results)
#    for i in range(KYBERQ):
#        if i < KYBERQ // 2:
#            assert np.allclose(
#                priors[0, 0, 0, i], priors2[0, 0, 0, i + KYBERQ]
#            ), f"{i}, {priors[0,0,0,i]}, {priors2[0,0,0,i+KYBERQ]}"
#        else:
#            assert np.allclose(
#                priors[0, 0, 0, i], priors2[0, 0, 0, i - KYBERQ]
#            ), f"{i}, {priors[0,0,0,i]}, {priors2[0,0,0,i-KYBERQ]}"

# attack = BPAttack(fg)
#    results = attack.attack_batch(priors, progress)
#   np.save("results", results)
# results = np.load("results.npy")
# print(rank(results, None))
