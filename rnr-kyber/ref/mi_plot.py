import csv
import matplotlib.pyplot as plt
# import tikzplotlib

def get_data(filename):
    H, CH, I = [], [], []
    with open(filename, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            for D, d in zip((H, CH, I), row[1:]):
                D.append(float(d))
    return (H, CH, I)


signed = get_data("../data/mi_ring_size_s.csv")
unsigned = get_data("../data/mi_ring_size_u.csv")
wide_u = get_data("../data/mi_ring_size_w_u.csv")
wide_s = get_data("../data/mi_ring_size_w_s.csv")


# for D, l in zip((H, CH, I), ("H(X)", "H(X|hw(X))", "I(X;hw(X))")):
#     plt.plot(range(1, len(D)+1), D, label=l)
for (H, CH, I), name in zip((signed, unsigned, wide_u, wide_s), ("signed", "unsigned", "wide_u", "wide_s")):
    # plt.plot(range(2, len(H)+1), [i/h for i,h in zip(I[1:], H[1:])], label=("I(X;hw(X)) / H(X) %s" % name))
    plt.plot(range(1, len(I)+1), I, label=("I(X;hw(X)) %s" % name))
plt.ylabel("Value")
plt.xlabel("Ring Size")
plt.xscale("log", base=2)
plt.legend()
plt.show()
# plt.savefig("../tikz/mi_ring_size.tikz", backend="tikz", bbox_inches="tight", pad_inches=0.0)