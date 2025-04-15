import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--R", type=int, default=2**16)
    parser.add_argument("--Q", type=int, default=3329)
    parser.add_argument("--K", type=int, default=1)

    args = parser.parse_args()
    R = args.R
    Q = args.Q
    K = args.K
    MONT = R % Q
    MINV = pow(R, -1, Q)
    QINV = pow(-Q, -1, R)
    WIDE_Q = K * Q
    WIDE_MONT = R % WIDE_Q
    WIDE_QINV = pow(-WIDE_Q, -1, R)
    WIDE_MINV = pow(R, -1, WIDE_Q)

    print(f"""R: {R}
Q: {Q}
K: {K}
MONT: {MONT}
MINV: {MINV}
QINV: {QINV}

WIDE_Q: {WIDE_Q}
WIDE_MONT: {WIDE_MONT}
WIDE_MINV: {WIDE_MINV}
WIDE_QINV: {WIDE_QINV}""")
