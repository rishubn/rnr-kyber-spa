import ctypes
from enum import Enum, auto
import numpy as np
import numpy.typing as npt
from dataclasses import dataclass
from .consts import KYBERQ, WIDE_S_KYBER_Q, WIDE_U_KYBER_Q
import os

kyber_funcs = ctypes.CDLL("libpqcrystals_ntt_ref.so")
kyber_funcs.pqcrystals_kyber768_ref_invntt.restype = None
kyber_funcs.pqcrystals_kyber768_ref_invntt.argtype = [ctypes.POINTER(ctypes.c_short)]

kyber_funcs.pqcrystals_kyber768_ref_invntt_full.argtypes = [
    ctypes.POINTER(ctypes.c_short),
    ctypes.POINTER((ctypes.c_short * 9) * 256),
]
kyber_funcs.pqcrystals_kyber768_ref_invntt.restype = None


class INTTImpl(Enum):
    SIGNED = auto()


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
    def from_impl(cls, impl):
        match impl:
            case INTTImpl.SIGNED:
                return Domain.KYBER_Q2_Q2()

    @classmethod
    def KYBER_0_Q(cls):
        return Domain(start=0, end=KYBERQ)

    @classmethod
    def KYBER_Q2_Q2(cls):
        return Domain(start=-KYBERQ // 2 + 1, end=KYBERQ // 2 + 1)

    @classmethod
    def KYBER_WSQ2_WSQ2(cls):
        return Domain(start=-WIDE_S_KYBER_Q // 2 + 1, end=WIDE_S_KYBER_Q // 2 + 1)

    @classmethod
    def KYBER_0_WUQ(cls):
        return Domain(start=0, end=WIDE_U_KYBER_Q)


class INTT:
    def __init__(self, impl: INTTImpl):
        self.impl = impl
        self.domain = Domain.from_impl(impl)

    def compute(self, r):
        pass

    def compute_full(self, _r, rlen):
        r = [0] * 256
        arr = np.zeros((256, 9), dtype=np.int16)

        for i in range(rlen):
            r[i] = _r[i]
        arr = np.zeros((256, 9), dtype=np.int16)
        r_c = (ctypes.c_short * 256)(*r)
        arr_c = np.ctypeslib.as_ctypes(arr)
        match self.impl:
            case INTTImpl.SIGNED:
                kyber_funcs.pqcrystals_kyber768_ref_invntt_full(r_c, arr_c)
        return np.ctypeslib.as_array(arr_c)
