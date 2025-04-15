import matplotlib.pyplot as plt
import numpy as np
from consts import Domain
from utils import Artifact
from argparse import ArgumentParser


def plot_distri(arr, domain, title="", axvline=None):
    print(arr.shape)
    plt.figure()
    plt.title(title)
    plt.stem(domain, arr)
    if axvline:
        plt.axvline(axvline, color="r")
    plt.show()


def main():
    pass


if __name__ == "__main__":
    main()
