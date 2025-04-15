import reduce
import numpy as np
from math import log2
from abc import ABC, abstractmethod
import numpy.typing as npt
from loguru import logger
from string import Template
from consts import Domain, KYBERQ, K_KYBERQ, REDUNDANCY
from enum import Enum, auto
from leakage_models import LeakageFunction

# fmt: off
zetas = np.array([-1044, -758,  -359,  -1517, 1493,  1422,  287,   202,  -171,  622,   1577,
182,   962,   -1202, -1474, 1468,  573,   -1325, 264,  383,   -829,  1458,
-1602, -130,  -681,  1017,  732,   608,   -1542, 411,  -205,  -1571, 1223,
652,   -552,  1015,  -1293, 1491,  -282,  -1544, 516,  -8,    -320,  -666,
-1618, -1162, 126,   1469,  -853,  -90,   -271,  830,  107,   -1421, -247,
-951,  -398,  961,   -1508, -725,  448,   -1065, 677,  -1275, -1103, 430,
555,   843,   -1251, 871,   1550,  105,   422,   587,  177,   -235,  -291,
-460,  1574,  1653,  -246,  778,   1159,  -147,  -777, 1483,  -602,  1119,
-1590, 644,   -872,  349,   418,   329,   -156,  -75,  817,   1097,  603,
610,   1322,  -1285, -1465, 384,   -1215, -136,  1218, -1335, -874,  220,
-1187, -1659, -1185, -1530, -1278, 794,   -1510, -854, -870,  478,   -108,
-308,  996,   991,   958,   -1460, 1522,  1628], dtype=np.int16)
# fmt: on


class INTTImpl(Enum):
    SIGNED = auto()
    UNSIGNED = auto()
    WIDE_UNSIGNED = auto()
    WIDE_SIGNED = auto()

    def to_domain(self):
        match self:
            case INTTImpl.SIGNED:
                return Domain.KYBER_Q2_Q2()
            case INTTImpl.UNSIGNED:
                return Domain.KYBER_0_Q()
            case INTTImpl.WIDE_UNSIGNED:
                return Domain.KYBER_0_KQ()
            case INTTImpl.WIDE_SIGNED:
                return Domain.KYBER_KQ2_KQ2()

    def to_leakage_function(self):
        match self:
            case INTTImpl.SIGNED | INTTImpl.WIDE_SIGNED:
                return LeakageFunction.HW_signed
            case INTTImpl.UNSIGNED | INTTImpl.WIDE_UNSIGNED:
                return LeakageFunction.HW

    def to_string(self):
        match self:
            case INTTImpl.SIGNED:
                return "signed"
            case INTTImpl.UNSIGNED:
                return "unsigned"
            case INTTImpl.WIDE_UNSIGNED:
                return f"wide unsigned K={REDUNDANCY}"
            case INTTImpl.WIDE_SIGNED:
                return f"wide signed K={REDUNDANCY}"


class INTTInputGenerator(ABC):
    @abstractmethod
    def generate(num_expr: int):
        pass


class RandomInputGenerator(INTTInputGenerator):
    def __init__(self, domain: Domain = Domain.KYBER_Q2_Q2()) -> None:
        self.domain = domain

    def generate(self, num_expr: int, size: int = 256):
        return np.random.randint(
            low=self.domain.start,
            high=self.domain.end,
            size=(num_expr, size),
            dtype=(
                np.int16
                if self.domain == Domain.KYBER_Q2_Q2()
                or self.domain == Domain.KYBER_KQ2_KQ2()
                else np.uint16
            ),
        )


class SparseInputGenerator(RandomInputGenerator):
    def __init__(
        self, domain: Domain = Domain.KYBER_Q2_Q2(), non_zero: int = 32
    ) -> None:
        super().__init__(domain)
        self.non_zero = non_zero

    def generate(self, num_expr: int, size=256):
        arr = super().generate(num_expr, size)
        # arr = np.ones((num_expr, size), dtype=np.int16)
        arr[:, -(size - self.non_zero) :] = 0
        return arr


class INTT:
    def __init__(
        self,
        input_gen: INTTInputGenerator,
        inttimpl: INTTImpl,
        height: int = 256,
        layers: int = 7,
    ):
        self.input_gen = input_gen
        self.height = height
        self.layers = layers
        self.inttimpl = inttimpl

    def simulate(self, num_expr: int) -> npt.NDArray:
        input = self.input_gen.generate(num_expr, size=self.height)
        logger.debug("INTT input:\n {}", input)
        output = np.zeros((num_expr, self.height, self.layers + 2), dtype=input.dtype)
        for i in range(num_expr):
            output[i] = self.compute(input[i])
        logger.debug("INTT output:\n {}", output)
        return output

    @classmethod
    def make_graph(cls, height=256, layers=7) -> str:
        template = Template("NC 3329\n$vars\n$consts\n$bff\n$factors")
        properties = []
        lvars = []
        consts = []
        k = height // 2 - 1
        dist = 2
        layer = 0
        i = 0
        while dist <= (height // 2):
            start = 0
            while start < height:
                j = start
                while j < (start + dist):
                    a = f"v{j}_{layer}"
                    b = f"v{j+dist}_{layer}"
                    c = f"v{j}_{layer+1}"
                    d = f"v{j+dist}_{layer+1}"
                    dsub = f"v{j+dist}_{layer+1}_sub"
                    lvars.append(f"VAR SINGLE {a}")
                    lvars.append(f"VAR SINGLE {b}")
                    lvars.append(f"VAR SINGLE {c}")
                    lvars.append(f"VAR SINGLE {dsub}")
                    lvars.append(f"VAR SINGLE {d}")
                    properties.append(
                        f"PROPERTY F{i}: _INV_BUTTERFLY({a}, {b}, {c}, {dsub})"
                    )
                    properties.append(f"PROPERTY F{i+1}: {d} = zeta{k} * {dsub}")
                    i += 2
                    j += 1
                consts.append(f"PUB SINGLE zeta{k}")
                k -= 1
                start = j + dist
            dist <<= 1
            layer += 1

        return template.substitute(
            {
                "factors": "\n".join(properties),
                "vars": "\n".join(sorted(list(set(lvars)))),
                "consts": "\n".join(consts),
                "bff": "GENERIC SINGLE _INV_BUTTERFLY",
            }
        )

    def fqmul(self, a, b):
        match self.inttimpl:
            case INTTImpl.SIGNED:
                return reduce.fqmul(a, b)
            case INTTImpl.UNSIGNED:
                return reduce.fqmul_u(a, b)
            case INTTImpl.WIDE_UNSIGNED:
                return reduce.fqmul_wide_u(a, b)
            case INTTImpl.WIDE_SIGNED:
                return reduce.fqmul_wide_s(a, b)

    def barrett_reduce(self, a):
        match self.inttimpl:
            case INTTImpl.SIGNED:
                return reduce.barrett_reduce(a)
            case INTTImpl.UNSIGNED:
                return reduce.barrett_reduce_u(a)
            case INTTImpl.WIDE_UNSIGNED:
                return reduce.barrett_reduce_wide_u(a)
            case INTTImpl.WIDE_SIGNED:
                return reduce.barrett_reduce_wide_s(a)

    def inv_butterfly(self, a, b, zeta):
        match self.inttimpl:
            case INTTImpl.SIGNED:
                t = a
                c = self.barrett_reduce(t + b)
                d = b - t
                d = self.fqmul(zeta, d)
                return c, d
            case INTTImpl.UNSIGNED:
                t = a
                c = self.barrett_reduce(t + b)
                d = self.barrett_reduce(2 * KYBERQ + b - t)
                d = self.fqmul(zeta, d)
                return c, d
            case INTTImpl.WIDE_UNSIGNED:
                t = a
                c = self.barrett_reduce(t + b)
                d = self.barrett_reduce(2 * K_KYBERQ + b - t)
                d = self.fqmul(zeta, d)
                return c, d
            case INTTImpl.WIDE_SIGNED:
                t = a
                c = self.barrett_reduce(t + b)
                d = b - t
                d = self.fqmul(zeta, d)
                return c, d

    def get_zeta(self, k):
        match self.inttimpl:
            case INTTImpl.SIGNED | INTTImpl.WIDE_SIGNED:
                return zetas[k]
            case INTTImpl.UNSIGNED | INTTImpl.WIDE_UNSIGNED:
                return zetas[k] + reduce.KYBER_Q if zetas[k] < 0 else zetas[k]

    def compute(self, r: npt.NDArray) -> npt.NDArray:
        assert log2(self.height).is_integer(), "Num vars is not a power of 2"
        intermediates = np.zeros((self.height, self.layers + 2), dtype=np.int16)
        intermediates[:, 0] = r[:]
        k = self.height // 2 - 1
        dist = 2  # len
        f: np.int16 = 1441
        layer = 1
        while dist <= (self.height // 2):
            start = 0
            while start < self.height:
                zeta = self.get_zeta(k)
                k -= 1
                j = start
                while j < (start + dist):
                    c, d = self.inv_butterfly(r[j], r[j + dist], zeta)
                    r[j] = c
                    r[j + dist] = d
                    j += 1
                start = j + dist
            intermediates[:, layer] = r[:]
            dist <<= 1
            layer += 1
        for j in range(self.height):
            r[j] = self.fqmul(r[j], f)
        intermediates[:, layer] = r[:]
        return intermediates


if __name__ == "__main__":
    print(INTT.make_graph())
