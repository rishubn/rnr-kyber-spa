import numpy as np
import numpy.typing as npt
from consts import Domain, KYBERQ


def normalize(arr):
    return arr / np.sum(arr)


def entropy(arr: npt.NDArray) -> npt.NDArray:
    return -1 * np.nansum(arr * np.log2(arr), axis=-1)


class Stats:
    def __init__(self, intt_impl, non_zero_sparse):
        self.intt_impl = intt_impl
        self.key_domain = intt_impl.to_domain()
        self.non_zero_sparse = non_zero_sparse

    def rank(self, distri, key):
        key = self.key_domain.map_to_Q(key)
        ranking = np.argsort(distri)[..., ::-1]
        rank = np.zeros((key.shape[0], self.non_zero_sparse))
        for i in range(key.shape[0]):
            for j in range(self.non_zero_sparse):
                rank[i, j] = (np.where(ranking[i, j] == key[i, j])[0][0]) + 1

        rank = np.log2(rank)
        return rank

    def avg_rank(self, distri, key):
        key = self.key_domain.map_to_Q(key)
        ranki = self.rank(distri, key)
        return ranki.mean(axis=1).mean()

    def success_rate(self, distri: npt.NDArray, key: npt.NDArray) -> npt.NDArray:
        key = self.key_domain.map_to_Q(key)
        ranking = np.argsort(distri)[..., ::-1]
        return (ranking[..., 0] == key).all(axis=-1).sum() / (key.shape[0])
