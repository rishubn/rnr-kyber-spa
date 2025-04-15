import ctypes
import numpy as np
from enum import Enum

kyber_funcs = ctypes.CDLL("../kyber-unsigned/ref/lib/libpqcrystals_ntt_ref.so")
kyber_funcs.pqcrystals_kyber768_ref_invntt.restype = None
kyber_funcs.pqcrystals_kyber768_ref_invntt.argtype = [ctypes.POINTER(ctypes.c_short)]

kyber_funcs.pqcrystals_kyber768_ref_invntt_full.argtypes = [
    ctypes.POINTER(ctypes.c_short),
    ctypes.POINTER((ctypes.c_short * 9) * 256),
]
kyber_funcs.pqcrystals_kyber768_ref_invntt.restype = None


class INTTImpl(Enum):
    SIGNED = "p"
    WIDE_SIGNED = "w"


class INVNTT:
    @property
    def key(self):
        return [0]

    @property
    def output_len(self):
        return 64

    @property
    def dtype(self):
        return np.uint8

    @property
    def sdtype(self):
        return "B"

    @property
    def pt_len(self):
        return 32

    def __init__(self, impl, fixed=False):
        self.intt_impl = impl
        self.fixed = fixed
        if self.fixed:
            self.fixed_msg = self.gen_msg()

    def gen_msg(self, rand=True):
        match self.intt_impl:
            case INTTImpl.SIGNED:
                return self.gen_msg_s(rand=rand)
            case INTTImpl.WIDE_SIGNED:
                return self.gen_msg_ws(rand=rand)

    def gen_msg_ws(self, rand=True):
        pt = (
            np.random.randint(low=-1664, high=1664, size=(32,), dtype=np.int16)
            if rand
            else np.zeros((32,), dtype=np.np.int16)
        )

    def gen_msg_s(self, rand=True):
        pt = (
            np.random.randint(low=-1664, high=1664, size=(32,), dtype=np.int16)
            if rand
            else np.zeros((32,), dtype=np.np.int16)
        )
        return pt

    def invntt_full(self, msg):
        r = [0] * 256
        arr = np.zeros((256, 9), dtype=np.int16)

        for i in range(32):
            r[i] = msg[i]
        arr = np.zeros((256, 9), dtype=np.int16)
        r_c = (ctypes.c_short * 256)(*r)
        arr_c = np.ctypeslib.as_ctypes(arr)
        kyber_funcs.pqcrystals_kyber768_ref_invntt_full(r_c, arr_c)
        return np.ctypeslib.as_array(arr_c)

    def gen_labels(self, msg):
        labels = self.invntt_full(msg)
        labels[labels < 0] += 3329
        return labels
