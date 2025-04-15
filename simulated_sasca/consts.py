from __future__ import annotations
from enum import Enum
import numpy as np
import numpy.typing as npt
from typing import Final
from dataclasses import dataclass

KYBERQ: Final[int] = 3329
REDUNDANCY: int = 1
K_KYBERQ: Final[int] = REDUNDANCY * KYBERQ


@dataclass(frozen=True)
class Domain:
    start: int
    end: int

    def to_range(self) -> npt.NDArray:
        return np.arange(self.start, self.end)

    def map_to_Q(self, arr):
        if self == Domain.KYBER_0_Q():
            return arr
        elif self == Domain.KYBER_Q2_Q2():
            arr[arr < 0] += KYBERQ
            return arr
        elif self == Domain.KYBER_0_KQ():
            arr %= KYBERQ
            return arr
        elif self == Domain.KYBER_KQ2_KQ2():
            arr %= KYBERQ
            return arr

    @classmethod
    def KYBER_0_Q(cls):
        return Domain(start=0, end=KYBERQ)

    @classmethod
    def KYBER_Q2_Q2(cls):
        return Domain(start=-KYBERQ // 2 + 1, end=KYBERQ // 2 + 1)

    @classmethod
    def KYBER_0_KQ(cls):
        return Domain(start=0, end=K_KYBERQ)

    @classmethod
    def KYBER_KQ2_KQ2(cls):
        return Domain(start=-K_KYBERQ // 2 + 1, end=K_KYBERQ // 2 + 1)


# class Domain(Enum):
#    KYBER_Q2_Q2 = [-KYBERQ // 2 + 1, KYBERQ // 2 + 1]
#    KYBER_0_Q = [0, KYBERQ]
#    KYBER_0_KQ = [0, KYBERQ_ENC]
#
#    def to_range(self) -> npt.NDArray:
#        return np.arange(self.value[0], self.value[1])
#
#    @classmethod
#    def map_domain_values(cls, a: Domain, b: Domain, arr: npt.NDArray) -> npt.NDArray:
#        if a == b:
#            return arr
#        match (a, b):
#            case (Domain.KYBER_Q2_Q2, Domain.KYBER_0_Q):
#                arr[arr < 0] += KYBERQ
#                return arr
#
#    @classmethod
#    def map_domain(cls, a: Domain, b: Domain, arr: npt.NDArray) -> npt.NDArray:
#        match (a, b):
#            case (Domain.KYBER_Q2_Q2, Domain.KYBER_0_Q):
#                return arr[
#                    ...,
#                    np.concatenate(
#                        [
#                            np.where(np.arange(a.value[0], a.value[1]) >= 0)[0],
#                            np.where(np.arange(a.value[0], a.value[1]) < 0)[0],
#                        ],
#                        axis=-1,
#                    ),
#                ]
#            case (Domain.KYBER_0_Q, Domain.KYBER_Q2_Q2):
#                return arr[
#                    ...,
#                    np.concatenate(
#                        [
#                            np.where(np.arange(b.value[0], b.value[1]) > 0)[0],
#                            np.where(np.arange(b.value[0], b.value[1]) <= 0)[0],
#                        ],
#                        axis=-1,
#                    ),
#                ]
#

# class DomainMap(Enum):
#    def KYBER_INTT(arr):
#        return arr[
#            np.concatenate([np.where(arr >= 0)[0], np.where(arr < 0)[0]], axis=-1)
#        ]
