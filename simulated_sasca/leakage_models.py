from __future__ import annotations
import numpy as np
from enum import Enum
import numpy.typing as npt
from scipy.stats import norm
from loguru import logger
from consts import KYBERQ, REDUNDANCY, Domain


class LeakageFunction(Enum):
    def HW(x, word_size=None):
        return bin(x).count("1")

    def HW_signed(x, word_size=16):
        if x < 0:
            return word_size - LeakageFunction.HW(-x - 1, word_size)
        else:
            return LeakageFunction.HW(x)


class LeakageModel:
    def __init__(self, leakage_function: LeakageFunction, word_size: int) -> None:
        self.word_size = word_size
        self.leakage_function = leakage_function
        self.leakage_function_vec = np.vectorize(self.leakage_function)

    def generate_leakage(
        self,
        values: npt.NDArray,
        domain: Domain,
        sigma: float,
        non_zero_sparse: int = None,
    ) -> npt.NDArray:
        joint_distri = np.histogram2d(
            domain.map_to_Q(domain.to_range()),
            [self.leakage_function(i, self.word_size) for i in domain.to_range()],
            bins=[np.arange(0, KYBERQ + 1), np.arange(0, self.word_size + 2)],
            density=False,
        )[0]
        hwx = np.array(
            [self.leakage_function(x, self.word_size) for x in values.reshape(-1)],
            dtype=values.dtype,
        )
        # hwx = self.leakage_function_vec(values, self.word_size).reshape(-1)
        # Norm.pdf(HW(x) - H, 0, \sigma)
        z = norm.pdf(
            hwx - np.arange(self.word_size + 1)[..., np.newaxis], loc=0.0, scale=sigma
        )
        # Pr[H] = \sum_X Pr[X,H]
        H = joint_distri.sum(axis=0)
        # Pr[H + z]
        H_z = np.dot(H / (KYBERQ * REDUNDANCY), z)
        # Pr[H=hwx | X = x]
        H_hwx_X_x = (joint_distri).dot(z)
        i = 0
        for a in H_z:
            i += 1
            assert a != 0, f"{i} {a}"
        X_H_hwx = (H_hwx_X_x / REDUNDANCY) / (H_z * KYBERQ)

        out = X_H_hwx.T.reshape(*values.shape, KYBERQ)
        if non_zero_sparse is not None:
            hw0 = np.zeros((KYBERQ))
            hw0[0] = 1.0
            out[:, -(values.shape[1] - non_zero_sparse), 0] = hw0
        return out


def vedad(x, Q, R, D, lf, sigma=0.0001):
    hw_hist = np.zeros((D + 1))
    hw_hist_joint = np.zeros((Q, D + 1))
    A_H_h = np.zeros((Q))
    hwx = lf(x)
    QR = Q * R
    for ai in range(QR):
        hwa = lf(ai)
        hw_hist[hwa] += 1
        real_a = ai % Q
        hw_hist_joint[real_a][hwa] += 1

    for a in range(Q):
        H_h = 0
        H_h_A_a = 0
        for hw in range(D):
            Z_hw = norm.pdf(hwx - hw, loc=0.0, scale=sigma)
            H_h += (hw_hist[hw] / QR) * Z_hw
            H_h_A_a += (hw_hist_joint[a][hw] / R) * Z_hw
        A_H_h[a] = H_h_A_a / (H_h * Q)

    print(A_H_h)
    return A_H_h


if __name__ == "__main__":
    # x = 16  # np.random.randint(KYBERQ_ENC)

    # hist = np.zeros((KYBERQ, 16))
    # for i in range(KYBERQ):
    #     hist[i,] = np.histogram(
    #         LeakageFunction.HW_enc(i, word_size=16), bins=np.arange(17)
    #     )[0]
    # hw_hist = np.histogram(
    #     [LeakageFunction.HW(i, word_size=16) for i in Domain.KYBER_0_KQ.to_range()],
    #     bins=np.arange(17),
    # )
    sigma = 0.0001
    # x = np.random.randint(0, KYBERQ, (2, 256, 9))
    lm = LeakageModel(LeakageFunction.HW_signed, word_size=4)
    llm = lm.generate_leakage(np.array([-1]), sigma=sigma, domain=Domain.KYBER_Q2_Q2)
    print(llm)
    # vedad(0, 13, 1, 4, LeakageFunction.HW)
    # vv = np.array([vedad(LeakageFunction.HW(x), sigma=sigma) for x in [0]])
    # print(vv)
    # assert np.allclose(ll.T, vv)
    # for x in range(13):
    #     ll = leak(x, sigma=sigma)
    #     print("=======================================")
    #     vv = vedad(LeakageFunction.HW(x), sigma=sigma)
    #     assert np.allclose(ll, vv)
    print("DONE")
# q = 7
# x = np.arange(q)

# print("-----------x-----------------")
# print(x)
# x[x > q // 2] -= q
# print("-----------x-----------------")
# print(x)
# hw = LeakageModel(LeakageFunction.HW_signed, word_size=4)
# leak = hw.generate_leakage(x, sigma=0.001)
# print("-----------leak-----------------")
# print(leak)
# prior = hw.leakage_to_prior(leak, np.arange(-q // 2 + 1, q // 2 + 1))
# print("------------prior---------------")
# for i, j in zip(np.arange(-q // 2 + 1, q // 2 + 1), prior[6]):
#     print(f"{i} -> {j}: {LeakageFunction.HW_signed(i, 4)}")
