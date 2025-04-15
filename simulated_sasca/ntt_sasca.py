import argparse
import numpy as np
from scalib.attacks import BPState, FactorGraph
from scalib.attacks.factor_graph import GenFactor
from scalib.postprocessing.rankestimation import rank_accuracy
from scalib.tools import ContextExecutor
from argparse import ArgumentParser
from pathlib import Path
import re
import toml
from scipy.stats import norm
import matplotlib.pyplot as plt
import tqdm
from concurrent.futures import as_completed
from multiprocessing import cpu_count


def get_argparse():
    def check_file(file: str):
        p = Path(file)
        if p.exists():
            return p
        else:
            raise FileNotFoundError(f"{file} does not exist")

    parser = ArgumentParser()
    parser.add_argument("--nc", type=int, required=True)
    parser.add_argument("--factor_graph", type=check_file)
    parser.add_argument("--noise-std", type=float)
    parser.add_argument("--num-experiments", type=int, default=100)
    parser.add_argument("--v-thresh", type=float, default=1)
    parser.add_argument("--outfile", type=str, default="entropies.txt")
    return parser

def variance(distribution, recenter=True):
    if recenter == True:
        center = int(np.round(np.nansum(np.arange(len(distribution)) * distribution)))
        p = np.roll(distribution, (len(distribution)//2) - center)
    else:
        p = distribution
    x = np.arange(len(p))
    prod = x * p
    return np.nansum(x * prod) - (np.nansum(prod))**2

def compute_stats(posteriors, recover_vars):
    entropy = 0
    entropies = []
    variances = []
    costs = []
    key = []
    for var in recover_vars.keys():
        entropy += -1 * np.nansum(posteriors[var] * np.log2(posteriors[var]))
        entropies.append(-1 * np.nansum(posteriors[var] * np.log2(posteriors[var])))
        variances.append(variance(posteriors[var]))
        costs.append(posteriors[var])
        key.append(recover_vars[var])

    # costs = np.array(costs, dtype=np.float64)
    # rmin, r, rmax = rank_accuracy(costs, key)
    return entropy, 1, entropies, variances  # , np.log2(r)
    return entropy, 1  # , np.log2(r)


def check_pin(user_input_pin):
    correct_pin = [1, 2, 3, 4]
    sum = 0
    for i in range(4):  # count from 0 to 3
        # add 1 if pins are equal else 0
        sum += user_input_pin[i] == correct_pin[i]
    return sum == 4


def run_sasca(priors, fg: FactorGraph, public_values=None, iterations=20):
    bp = BPState(fg, 1, public_values=public_values)
    posteriors = {}
    # set prior evidence
    for var, data in priors.items():
        bp.set_evidence(var, data)
    # do BP
    bp.bp_loopy(iterations, initialize_states=False)

    # get posteriors
    for var in fg.vars():
        posteriors[var] = bp.get_distribution(var)

    return posteriors


## trace/prior generation


def gen_priors(model, leaking_vars, noise_std, nc, rng):
    traces = {}
    for var in leaking_vars:
        data = model[var]
        hw = HW(data)
        lx = rng.normal(hw, noise_std)
        p = np.zeros((nc,))
        for i in range(nc):
            p[i] = norm.pdf(lx, HW(i), noise_std)
        p = np.round(p / np.sum(p), 10)
        traces[var] = p

    # zero out half of ntt
    nv = len(leaking_vars)
    a = np.zeros((nv,))
    a[:nv] = 1
    # np.random.shuffle(a)
    for ftr, var in zip(a, sorted(leaking_vars)):
        zero_pdf = np.zeros((nc,))
        zero_pdf[0] = 1
        if ftr:
            traces[var] = zero_pdf

    return traces


def HW(x):
    return bin(x).count("1")


## Factorgraph parsing


def evaluate_fg(circ, IV):
    # initialize model with constants and initial values
    evaluated_model = {x[1]: 0 for x in circ["vars"]}
    for k, v in circ["pubs"].items():
        evaluated_model[k] = v
    for k, v in IV.items():
        evaluated_model[k] = v

    # compute values of intermediates from IV + constants
    for prop in circ["properties"]:
        res = prop[1]  # store result
        if ("+" in prop) or ("*" in prop):
            operands = {
                prop[2]: evaluated_model[prop[2]],
                prop[4]: evaluated_model[prop[4]],
            }
            expr = "(" + " ".join(prop[2:]) + f') % {circ["nc"]}'
            evaluated_model[res] = eval(expr, operands)
        elif "invert" in prop:
            evaluated_model[res] = circ["tables"]["invert"][evaluated_model[prop[3]]]

    return evaluated_model


def parse_factorgraph(fg, consts=None, tables=None, nc=None):
    circ = {
        "nc": nc if nc else int(fg.split("\n")[0].split()[1]),
        "vars": [],
        "properties": [],
        "pubs": {},
        "tables": {},
        "generics": [],
    }

    for line in [s.strip() for s in fg.split("\n") if s]:
        if line.startswith("PROPERTY"):
            circ["properties"].append(
                re.findall(r"[\w]+|\^|\||\&|\!|\+|\*", line.split("PROPERTY")[1])
            )
        elif line.startswith("VAR"):
            circ["vars"].append(re.findall(r"[\w]+", line)[1:])
        elif line.startswith("GENERIC"):
            circ["generics"].append(re.findall(r"[\w]+", line)[1:])

    if consts:
        for k, v in consts.items():
            circ["pubs"][k] = v
    if tables:
        for k, v in tables.items():
            circ["tables"][k] = v
    return circ


def invert_table(nc):
    return np.array([(-x) % nc for x in range(nc)], np.uint32)


def do_experiment(model, priors, params, rng_state):
    fg = FactorGraph(df["factor_graph"], tables={"invert": invert_table(args.nc)})
    posteriors = run_sasca(priors, fg, public_values=params["consts"])

    # prior_e, prior_r = compute_stats(priors, {x: model[x] for x in params["inputs"]})
    post_e, post_r, post_es, post_vs = compute_stats(posteriors, {x: model[x] for x in params["inputs"]})

    return (post_e, post_r, post_es, post_vs, rng_state)


if __name__ == "__main__":
    import cProfile, sys

    args = get_argparse().parse_args()
    with open(args.factor_graph, "r") as f:
        df = toml.loads(f.read())

    params = df["params"]
    entropy = np.zeros((args.num_experiments,))
    data = np.zeros((2,args.num_experiments,2*len(params["consts"])))
    rng = np.random.default_rng()

    fg = parse_factorgraph(
        df["factor_graph"],
        consts=params["consts"],
        tables={"invert": invert_table(args.nc)},
        nc=args.nc,
    )

    i = 0
    progress = tqdm.tqdm()
    with open(args.outfile, "w") as f, ContextExecutor(max_workers=cpu_count()) as pool:
        try:
            # while True:
            for index in range(args.num_experiments):
                queue = []
                for _ in range(1): #2 * cpu_count()):
                    state = rng.bit_generator.state
                    iv = {x: rng.integers(args.nc) for x in params["inputs"]}
                    print("iv:", iv)
                    model = evaluate_fg(fg, iv)
                    priors = gen_priors(
                        model, params["leaking_vars"], args.noise_std, args.nc, rng
                    )
                    params["consts"] = {k: v % args.nc for k, v in params["consts"].items()}
                    future = pool.submit(
                        do_experiment,
                        model=model,
                        priors=priors,
                        params=params,
                        rng_state=state,
                    )
                    queue.append(future)
                for future in as_completed(queue):
                    post_e, post_r, post_es, post_vs, rng_state = future.result()
                    f.write(f"{post_e},{post_vs},{rng_state}\n")
                    np.copyto(data[0][i], post_es)
                    np.copyto(data[1][i], post_vs)
                    i += 1
                    progress.update()
            np.save(args.outfile[:args.outfile.index('.')], data)
            f.close()
            print("closing")
            sys.exit(0)
        except KeyboardInterrupt:
            f.close()
            print("closing")
            sys.exit(0)
    for i in tqdm.tqdm(range(args.num_experiments)):
        with cProfile.Profile() as pr:
            prior_stats, post_stats = do_experiment(
                parse_factorgraph(
                    df["factor_graph"],
                    consts=params["consts"],
                    tables={"invert": invert_table(args.nc)},
                    nc=args.nc,
                ),
                df["factor_graph"],
                params,
            )
            entropy[i] = post_stats[0]
            pr.print_stats(sort="time")
        sys.exit(0)
    plt.figure()
    plt.boxplot(entropy, positions=[args.noise_std])
    plt.savefig("out.png")

