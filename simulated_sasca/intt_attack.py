from os.path import exists
import sys
from scipy.stats import norm
import re
from pathlib import Path
import numpy as np
from argparse import ArgumentParser
from math import log2
from reduce import barrett_reduce, fqmul
import toml
from scalib.attacks import BPState, FactorGraph, GenFactor
from dataclasses import dataclass
zetas = [
    -1044,
    -758,
    -359,
    -1517,
    1493,
    1422,
    287,
    202,
    -171,
    622,
    1577,
    182,
    962,
    -1202,
    -1474,
    1468,
    573,
    -1325,
    264,
    383,
    -829,
    1458,
    -1602,
    -130,
    -681,
    1017,
    732,
    608,
    -1542,
    411,
    -205,
    -1571,
    1223,
    652,
    -552,
    1015,
    -1293,
    1491,
    -282,
    -1544,
    516,
    -8,
    -320,
    -666,
    -1618,
    -1162,
    126,
    1469,
    -853,
    -90,
    -271,
    830,
    107,
    -1421,
    -247,
    -951,
    -398,
    961,
    -1508,
    -725,
    448,
    -1065,
    677,
    -1275,
    -1103,
    430,
    555,
    843,
    -1251,
    871,
    1550,
    105,
    422,
    587,
    177,
    -235,
    -291,
    -460,
    1574,
    1653,
    -246,
    778,
    1159,
    -147,
    -777,
    1483,
    -602,
    1119,
    -1590,
    644,
    -872,
    349,
    418,
    329,
    -156,
    -75,
    817,
    1097,
    603,
    610,
    1322,
    -1285,
    -1465,
    384,
    -1215,
    -136,
    1218,
    -1335,
    -874,
    220,
    -1187,
    -1659,
    -1185,
    -1530,
    -1278,
    794,
    -1510,
    -854,
    -870,
    478,
    -108,
    -308,
    996,
    991,
    958,
    -1460,
    1522,
    1628,
]
@dataclass
class RandomInput:
    pass

@dataclass
class SparseInput:
    sparse_factor: int

InputType = RandomInput | SparseInput

def HW(x, wordsize):
    return bin(x).count("1")

def HW_arr(func, x, wordsize):
    print(x)
    return np.array([func(n, wordsize) for n in x])

def HW_signed(x, word_size):
    if x < 0:
        return word_size - HW(-x -1, word_size)
    else:
        return HW(x, word_size)

def parse_circuit(circuit_path):
    with open(Path(circuit_path)) as f:
        cf = toml.loads(f.read())
    fg = FactorGraph(cf["factor_graph"])
    leaking_vars = cf["params"]["leaking_vars"]
    recover_vars = cf["params"]["recover_vars"]
    consts = cf["params"]["consts"]
    return fg, leaking_vars, recover_vars, consts
    

def invntt_model(r: np.ndarray, height: int = 256, layers: int = 7):
    assert log2(height).is_integer(), "Num vars is not a power of 2"
    # assert layers >= log2(height), "Too many layers for num vars"
    intermediates = np.zeros((height, layers + 2), dtype=np.int16)
    intermediates[:, 0] = r[:]
    k = height // 2 - 1
    dist = 2  # len
    f: np.int16 = 1441
    layer = 1
    while dist <= (height // 2):
        start = 0
        while start < height:
            zeta = zetas[k]
            k -= 1
            j = start
            while j < start + dist:
                t = r[j]
                r[j] = barrett_reduce(t + r[j + dist])
                r[j + dist] = r[j + dist] - t
                r[j + dist] = fqmul(zeta, r[j + dist])
                j += 1
            start = j + dist
        intermediates[:, layer] = r[:]
        dist <<= 1
        layer += 1
    #for j in range(height):
    #    r[j] = fqmul(r[j], f)
    intermediates[:, layer] = r[:]
    return intermediates


def leakage_to_prior(leakage, word_size=16):
    x = np.zeros((3329,))
    for v in range(-1664, 1665):
        if v < 0:
            x[v+3329] = leakage[HW_signed(v, word_size)]
        else:
            x[v] = leakage[HW_signed(v, word_size)]
    return x

def generate_leakage(intt_model, leaking_vars, sigma=1.0, word_size=16):
    leakage = {}
    for var in leaking_vars:
        var_loc = list(map(int, re.findall(r"\d+", var)))
        val = intt_model[var_loc[1], var_loc[0]]
        hw = HW_signed(val, word_size)
        leakage[var] = norm.pdf(hw + np.random.normal(0, sigma, 1), range(0, word_size+1), sigma)
        leakage[var] = leakage_to_prior(leakage[var])
    return leakage
def generate_intt_input(input_type: InputType):
    match input_type:
        case RandomInput:
            return np.random.randint(-1664, 1665, (256,), dtype=np.int16)
    

def create_butterfly_factor():
    nc = 3329
    if not exists(f"./factor_cache/inv_butterfly_{nc}.npy"):
        bff = []
        for a in range(nc):
            for b in range(nc):
                bff.append([a, b, ((a + b) % nc), ((b - a) % nc)])
        x = np.array(bff, dtype=np.uint32)
        np.save(f"./factor_cache/inv_butterfly_{nc}.npy", x)
    else:
        x = np.load(f"./factor_cache/inv_butterfly_{nc}.npy")
        return x

def do_bp(fg, priors, recover_vars, consts, iterations=1):
    gf = GenFactor.sparse_functional(create_butterfly_factor())
    bp = BPState(fg, 1, public_values=consts, gen_factors={"_INV_BUTTERFLY":gf})
    for p, dist in priors.items():
        bp.set_evidence(p, dist)
    bp.bp_loopy(iterations, True)
    
    out = {}
    for v in recover_vars:
        out[v] = bp.get_distribution(v)
    return out

def test_intt(input, intt):
    inputarr = np.array(input, dtype=np.int16)
    outputarr = np.array(intt, dtype=np.int16)
    res = invntt_model(inputarr)
    if not np.allclose(res[:, -1], outputarr):
        print("INV NTTs do not match")
        print(res[:, -1])
        print(outputarr)
        sys.exit(1)
    sys.exit(0)


def compute_entropy(posteriors):
    return {
            k: -1*np.nansum(
                v * np.log2(v))
        for k, v in posteriors.items()
    }

def get_args():
    argparse = ArgumentParser()
    argparse.add_argument("--intt-file")
    argparse.add_argument("--test")
    argparse.add_argument("--noise-std", type=float)
    argparse.add_argument("--input", nargs="+", type=int)
    argparse.add_argument("--intt", nargs="+", type=int)
    return argparse.parse_args()


if __name__ == "__main__":
    args = get_args()
    if args.test:
        test_intt(args.input, args.intt)
    fg, leaking_vars, recover_vars, consts = parse_circuit(args.intt_file)
    intt_input = generate_intt_input(RandomInput())
    intt_model = invntt_model(intt_input)
    priors = generate_leakage(intt_model, leaking_vars)
    posteriors = do_bp(fg, priors, recover_vars, consts)
    entropy = compute_entropy(posteriors)
