import ctypes
from enum import Enum, auto
import numpy as np
import numpy.typing as npt
from dataclasses import dataclass
from .consts import KYBERQ, WIDE_S_KYBERQ, WIDE_U_KYBERQ
import os


lib_path = os.path.join(os.path.dirname(__file__), "libpqcrystals_ntt_ref.so")
kyber_funcs = ctypes.CDLL(lib_path)
kyber_funcs.pqcrystals_kyber768_ref_invntt.restype = None
kyber_funcs.pqcrystals_kyber768_ref_invntt.argtype = [ctypes.POINTER(ctypes.c_short)]

# SIGNED
kyber_funcs.pqcrystals_kyber768_ref_invntt_full.argtypes = [
    ctypes.POINTER(ctypes.c_short),
    ctypes.POINTER((ctypes.c_short * 9) * 256),
]
kyber_funcs.pqcrystals_kyber768_ref_invntt.restype = None

# WIDE SIGNED
kyber_funcs.pqcrystals_kyber768_ref_invntt_wide_s_full.argtypes = [
    ctypes.POINTER(ctypes.c_short),
    ctypes.POINTER((ctypes.c_short * 9) * 256),
]
kyber_funcs.pqcrystals_kyber768_ref_invntt_wide_s_full.restype = None

# UNSIGNED
kyber_funcs.pqcrystals_kyber768_ref_invntt_u_full.argtypes = [
    ctypes.POINTER(ctypes.c_ushort),
    ctypes.POINTER((ctypes.c_ushort * 9) * 256),
]
kyber_funcs.pqcrystals_kyber768_ref_invntt_u_full.restype = None

# WIDE UNSIGNED
kyber_funcs.pqcrystals_kyber768_ref_invntt_wide_u_full.argtypes = [
    ctypes.POINTER(ctypes.c_ushort),
    ctypes.POINTER((ctypes.c_ushort * 9) * 256),
]
kyber_funcs.pqcrystals_kyber768_ref_invntt_wide_u_full.restype = None


class INTTImpl(Enum):
    SIGNED = auto()
    WIDE_SIGNED = auto()
    UNSIGNED = auto()
    WIDE_UNSIGNED = auto()


@dataclass(frozen=True)
class Domain:
    start: int
    end: int

    def to_range(self) -> npt.NDArray:
        return np.arange(self.start, self.end)

    # def map_to_Q(self, arr):
    #     if self == Domain.KYBER_0_Q():
    #         return arr
    #     elif self == Domain.KYBER_Q2_Q2():
    #         arr[arr < 0] += KYBERQ
    #         return arr
    #     elif self == Domain.KYBER_0_WUQ():
    #         arr %= KYBERQ
    #         return arr
    #     elif self == Domain.KYBER_KQ2_KQ2():
    #         arr %= KYBERQ
    #         return arr

    @classmethod
    def from_impl(cls, impl):
        match impl:
            case INTTImpl.SIGNED:
                return Domain.KYBER_Q2_Q2()
            case INTTImpl.WIDE_SIGNED:
                return Domain.KYBER_WSQ2_WSQ2()
            case INTTImpl.UNSIGNED:
                return Domain.KYBER_0_Q()
            case INTTImpl.WIDE_UNSIGNED:
                return Domain.KYBER_0_WUQ()
            case _:
                raise NotImplementedError()

    @classmethod
    def KYBER_0_Q(cls):
        return Domain(start=0, end=KYBERQ)

    @classmethod
    def KYBER_Q2_Q2(cls):
        return Domain(start=-KYBERQ // 2 + 1, end=KYBERQ // 2 + 1)

    @classmethod
    def KYBER_WSQ2_WSQ2(cls):
        return Domain(start=-WIDE_S_KYBERQ // 2 + 1, end=WIDE_S_KYBERQ // 2 + 1)

    @classmethod
    def KYBER_0_WUQ(cls):
        return Domain(start=0, end=WIDE_U_KYBERQ)


class INTT:
    def __init__(self, impl: INTTImpl):
        self.impl = impl
        self.domain = Domain.from_impl(impl)

    def compute(self, r):
        pass

    def compute_full(self, _r):
        match self.impl:
            case INTTImpl.SIGNED | INTTImpl.WIDE_SIGNED:
                return self.compute_full_s(_r)
            case INTTImpl.UNSIGNED | INTTImpl.WIDE_UNSIGNED:
                return self.compute_full_u(_r)

    def compute_full_s(self, _r):
        r = [0] * 256
        arr = np.zeros((256, 9), dtype=np.int16)

        for i in range(len(_r)):
            r[i] = _r[i]

        r_c = (ctypes.c_short * 256)(*r)
        arr_c = np.ctypeslib.as_ctypes(arr)
        match self.impl:
            case INTTImpl.SIGNED:
                kyber_funcs.pqcrystals_kyber768_ref_invntt_full(r_c, arr_c)
            case INTTImpl.WIDE_SIGNED:
                kyber_funcs.pqcrystals_kyber768_ref_invntt_wide_s_full(r_c, arr_c)
            case _:
                raise NotImplementedError()
        return np.ctypeslib.as_array(arr_c)

    def compute_full_u(self, _r):
        r = [0] * 256
        arr = np.zeros((256, 9), dtype=np.uint16)
        for i in range(len(_r)):
            r[i] = _r[i]

        r_c = (ctypes.c_ushort * 256)(*r)
        arr_c = np.ctypeslib.as_ctypes(arr)
        match self.impl:
            case INTTImpl.UNSIGNED:
                kyber_funcs.pqcrystals_kyber768_ref_invntt_u_full(r_c, arr_c)
            case INTTImpl.WIDE_UNSIGNED:
                kyber_funcs.pqcrystals_kyber768_ref_invntt_wide_u_full(r_c, arr_c)
            case _:
                raise NotImplementedError()
        return np.ctypeslib.as_array(arr_c)
