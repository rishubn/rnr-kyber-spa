import numpy as np
import ctypes

kyber_funcs = ctypes.CDLL(
    "../kyber-unsigned/py/nttwrapper/src/nttwrapper/libpqcrystals_ntt_ref.so"
)
kyber_funcs.pqcrystals_kyber768_ref_barrett_reduce.restype = ctypes.c_short
kyber_funcs.pqcrystals_kyber768_ref_barrett_reduce_u.restype = ctypes.c_ushort
kyber_funcs.pqcrystals_kyber768_ref_barrett_reduce_wide_u.restype = ctypes.c_ushort
kyber_funcs.pqcrystals_kyber768_ref_barrett_reduce_wide_s.restype = ctypes.c_short
kyber_funcs.pqcrystals_kyber768_ref_fqmul.restype = ctypes.c_short
kyber_funcs.pqcrystals_kyber768_ref_fqmul_u.restype = ctypes.c_ushort
kyber_funcs.pqcrystals_kyber768_ref_fqmul_wide_u.restype = ctypes.c_ushort
kyber_funcs.pqcrystals_kyber768_ref_fqmul_wide_s.restype = ctypes.c_short

KYBER_Q = np.int16(3329)
QINV = np.int32(62209)
BARRET_V = np.int16(((1 << 26) + KYBER_Q // 2) // KYBER_Q)
MONT = np.int32(-1044)  # 2^16 mod q

QINV_U = np.int32(3327)
MONT_SQ_U = np.int32(1353)
INTT_CONST_U = np.int32(1441)


def barrett_reduce(a):
    return kyber_funcs.pqcrystals_kyber768_ref_barrett_reduce(ctypes.c_short(a))


def barrett_reduce_u(a):
    return kyber_funcs.pqcrystals_kyber768_ref_barrett_reduce_u(ctypes.c_ushort(a))


def barrett_reduce_wide_u(a):
    return kyber_funcs.pqcrystals_kyber768_ref_barrett_reduce_wide_u(ctypes.c_ushort(a))


def barrett_reduce_wide_s(a):
    return kyber_funcs.pqcrystals_kyber768_ref_barrett_reduce_wide_s(ctypes.c_short(a))


def fqmul(a, b):
    return kyber_funcs.pqcrystals_kyber768_ref_fqmul(
        ctypes.c_short(a), ctypes.c_short(b)
    )


def fqmul_u(a, b):
    return kyber_funcs.pqcrystals_kyber768_ref_fqmul_u(
        ctypes.c_ushort(a), ctypes.c_ushort(b)
    )


def fqmul_wide_u(a, b):
    return kyber_funcs.pqcrystals_kyber768_ref_fqmul_wide_u(
        ctypes.c_ushort(a), ctypes.c_ushort(b)
    )


def fqmul_wide_s(a, b):
    return kyber_funcs.pqcrystals_kyber768_ref_fqmul_wide_s(
        ctypes.c_short(a), ctypes.c_short(b)
    )


# def fqmul(a: np.int16, b: np.int16) -> np.int16:
#     return montgomery_reduce(np.int32(a) * b)
#

# /*************************************************
# * Name:        montgomery_reduce
# *
# * Description: Montgomery reduction; given a 32-bit integer a, computes
# *              16-bit integer congruent to a * R^-1 mod q, where R=2^16
# *
# * Arguments:   - int32_t a: input integer to be reduced;
# *                           has to be in {-q2^15,...,q2^15-1}
# *
# * Returns:     integer in {-q+1,...,q-1} congruent to a * R^-1 modulo q.
# **************************************************/


def montgomery_reduce(a: np.int32) -> np.int16:
    # dropping upper bits on purpose
    u = np.int16(np.int64(a) * np.int64(QINV))
    t = np.int32(u) * np.int32(KYBER_Q)
    t = a - t
    t >>= 16
    return np.int16(t)


montgomery_reduce_outrange = (-KYBER_Q + 1, KYBER_Q - 1)


# def fqmul_u(a: np.uint16, b: np.uint16) -> np.uint16:
#     return montgomery_reduce_u(np.uint32(a) * b)


def montgomery_reduce_u(a: np.uint32) -> np.uint16:
    assert a < (KYBER_Q << 16)
    m = np.uint16(np.uint16(a) * np.uint64(QINV_U))
    t = a + (np.uint32(m) * KYBER_Q) >> 16
    assert t < (2 * KYBER_Q)
    return t


# /*************************************************
# * Name:        barrett_reduce
# *
# * Description: Barrett reduction; given a 16-bit integer a, computes
# *              centered representative congruent to a mod q in {-(q-1)/2,...,(q-1)/2}
# *
# * Arguments:   - int16_t a: input integer to be reduced
# *
# * Returns:     integer in {-(q-1)/2,...,(q-1)/2} congruent to a modulo q.
# **************************************************/

# def barrett_reduce(a: np.int16) -> np.int16:
#     t = np.int16((np.int32(BARRET_V) * np.int32(a) + (1 << 25)) >> 26)
#     t = np.int16(np.int32(t) * np.int32(KYBER_Q))
#     return np.int16(np.int32(a) - np.int32(t))
#
#
# barret_reduce_outrange = (-(KYBER_Q - 1) // 2, (KYBER_Q - 1) // 2)
#
# BARRETT_K = 26
# BARRETT_CONST = np.uint16(((1 << BARRETT_K) + KYBER_Q // 2) // KYBER_Q)
#
#
# def barrett_reduce_u(a: np.uint16) -> np.uint16:
#     t = (np.uint32(a) * BARRETT_CONST) >> BARRETT_K
#     r = a - (t * KYBER_Q)
#     assert r < KYBER_Q
#     return r
#


# /*************************************************
# * Name:        csubq
# *
# * Description: Conditionallly subtract q
# *
# * Arguments:   - int16_t x: input integer
# *
# * Returns:     a - q if a >= q, else a
# **************************************************/
def csubq(a: np.int16) -> np.int16:
    a -= KYBER_Q
    a += (a >> 15) & KYBER_Q
    return a


if __name__ == "__main__":
    import argparse

    cli = argparse.ArgumentParser()
    cli.add_argument("--lista", type=str)
    for x in range(1 << 16):
        print(barrett_reduce_u(x))
