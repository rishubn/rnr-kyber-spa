from scalib.attacks import FactorGraph, GenFactor, BPState
from scalib.tools import ContextExecutor
import numpy as np
import numpy.typing as npt
from consts import KYBERQ, K_KYBERQ
import intt
from os.path import exists
from os import cpu_count
from concurrent.futures import as_completed


class BPAttack:
    def __init__(self, fg: str, height: int = 256):
        self.factor_graph = FactorGraph(fg)
        poszetas = BPAttack.map_zetas(intt.zetas)
        self.pubs = {f"zeta{i}": 169 * poszetas[i] for i in range(1, height // 2)}

    @classmethod
    def map_zetas(cls, zetas: npt.NDArray) -> npt.NDArray:
        out = zetas.copy()
        out[zetas < 0] = zetas[zetas < 0] + KYBERQ
        return out

    @classmethod
    def _inv_butterfly(cls, nc=KYBERQ) -> GenFactor:
        if not exists(f"./factor_cache/inv_butterfly_{nc}.npy"):
            bff = []
            for a in range(nc):
                for b in range(nc):
                    bff.append([a, b, ((a + b) % nc), ((b - a) % nc)])
            x = np.array(bff, dtype=np.uint32)
            np.save(f"./factor_cache/inv_butterfly_{nc}.npy", x)
            return GenFactor.sparse_functional(x)
        else:
            x = np.load(f"./factor_cache/inv_butterfly_{nc}.npy")
            return GenFactor.sparse_functional(x)

    @classmethod
    def _inv_butterfly_wide_u(cls, nc=KYBERQ) -> GenFactor:
        if not exists(f"./factor_cache/inv_butterfly_wide_u{nc}.npy"):
            bff = []
            for a in range(nc):
                for b in range(nc):
                    bff.append([a, b, ((a + b) % nc), ((2 * K_KYBERQ + (b - a)) % nc)])
            x = np.array(bff, dtype=np.uint32)
            np.save(f"./factor_cache/inv_butterfly_wide_u{nc}.npy", x)
            return GenFactor.sparse_functional(x)
        else:
            x = np.load(f"./factor_cache/inv_butterfly_wide_u{nc}.npy")
            return GenFactor.sparse_functional(x)

    def attack_batch(self, priors: npt.NDArray, num_iter=20, progress=None):
        num_expr = priors.shape[0]
        queue = []
        results = np.zeros((num_expr, priors.shape[1], KYBERQ))
        with ContextExecutor(max_workers=cpu_count()) as pool:
            for expr in range(num_expr):
                future = pool.submit(self.attack, prior=priors[expr], num_iter=num_iter)
                queue.append(future)
            for i, future in enumerate(queue):
                result = future.result()
                results[i,] = result
                progress.update()
        return results

    def attack(self, prior: npt.NDArray, num_iter: int = 20):
        bp = BPState(
            self.factor_graph,
            1,
            public_values=self.pubs,
            gen_factors={"_INV_BUTTERFLY": BPAttack._inv_butterfly(nc=KYBERQ)},
        )
        for i in range(prior.shape[1] - 1):
            for j in range(prior.shape[0]):
                bp.set_evidence(f"v{i}_{j}", prior[j, i])
        bp.bp_loopy(num_iter, initialize_states=True)  # TODO temp set to 1
        results = np.zeros((prior.shape[0], KYBERQ))
        for j in range(prior.shape[0]):
            results[j,] = bp.get_distribution(f"v0_{j}")
        return results
