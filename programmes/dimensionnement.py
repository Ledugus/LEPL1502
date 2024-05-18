from numpy import *
import matplotlib.pyplot as plt


w = linspace(0, 4000, 1000)
R2 = 12
L1 = 0.39 * 10 ** (-3)
L2 = 0.79 * 10 ** (-3)
M = 0.25 * 10 ** (-3)

# Capacités possibles
capas = array([10**5, 10**4, 10**3, 10**1]) * 10 * 10 ** (-12)
resistances = array([68, 200])


def opti_M(M, R1, R2):
    return sqrt(R1 * R2) / (M)


def combine(poss_list):
    half = poss_list / 2
    double = poss_list * 2
    sums = [i + j for i in poss_list for j in poss_list]
    return concatenate((half, double, sums, poss_list))


def capas_from_freq(L, freq):
    return 1 / (L * freq**2)


def plot_freq():
    R1 = combine(resistances)
    # create 2 subplots, one with a line and one with a scatter
    plt.subplot(211)
    plt.scatter(R1, opti_M(M, R1, R2) / 1000)
    plt.xlabel("Résistance R1")
    plt.xlim(0, max(R1 + 10))
    plt.ylabel("Fréquence optimale (kHz)")

    plt.subplot(212)
    C1 = combine(capas)
    for x, c in enumerate(C1):
        plt.plot(R1, ones_like(R1) * c, color=(x / len(C1), 0, 0), linestyle="--")
    plt.scatter(R1, capas_from_freq(L1, opti_M(M, R1, R2)))
    # set a log scale on plot
    plt.xlabel("Résistance R1")
    plt.xlim(0, max(R1 + 10))
    plt.yscale("log")
    plt.show()


plot_freq()
