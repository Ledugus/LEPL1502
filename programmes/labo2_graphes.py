import matplotlib.pyplot as plt
import numpy as np


@np.vectorize
def constant_function(x, t):
    return t


def graphe_v_in_v_dd():
    t1 = np.arange(0, 400, 1)
    plt.plot(t1, constant_function(t1, 5), label=r"$V_{dd}$")
    plt.plot(t1, constant_function(t1, 5 * 0.3334), label=r"$V_{dd}$")
    plt.ylim((0, 5.2))
    plt.xlabel(r"Temps ($\mu$s)")
    plt.ylabel(r"Tension (V)")
    plt.savefig("graphe_v_in_v_dd.png")


def triangle2(length, amplitude):
    section = length // 8
    x = np.linspace(0, amplitude, section + 1)
    return np.r_[x, x[-2:1:-1], x, x[-2::-1]]


def comparateur(input, value, top):
    x = np.array([top if i > value else 0 for i in input])
    return x


def graphe_v_in_v_out():
    t1 = np.arange(0, 200, 1)
    v_inp = triangle2(400, 2)
    v_out = comparateur(v_inp, 5 * 0.3334, 5)
    plt.plot(t1, v_inp, label=r"$V_{in+}$")
    plt.plot(t1, constant_function(t1, 5 * 0.3334), label=r"$V_{in-}$")
    plt.plot(t1, v_out, label=r"$V_{out}$")
    plt.ylim((0, 5.2))
    plt.legend()
    plt.xlabel(r"Temps ($\mu$s)")
    plt.ylabel(r"Tension (V)")
    plt.savefig("graphe_v_in_v_out.png")


graphe_v_in_v_out()
