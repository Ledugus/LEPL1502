import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


def square_signal(t, v_max):
    return (signal.square(np.pi * t) + 1) * v_max / 2


def v_c(t, RC, v_max):
    t1 = t[: len(t) // 4]

    v_c1 = v_max * (1 - np.exp(-t1 / RC))
    v_c2 = v_c1[-1] * (np.exp(-t1 / RC))
    v_c3 = v_max + (v_c2[-1] - v_max) * np.exp(-t1 / RC)
    v_c4 = v_c3[-1] * (np.exp(-t1 / RC))

    return np.concatenate((v_c1, v_c2, v_c3, v_c4))


def rc_plot():
    v_max = 2
    x = np.linspace(0, 4, 1000)
    square = square_signal(x, v_max)
    v_c_ = v_c(x, 1, v_max)
    plt.axes().set_aspect("equal")
    plt.plot(x, square, label="V(t)")
    plt.plot(x, v_c_, label="$V_C(t)$")

    plt.xlabel("Temps")
    plt.ylabel("Tension")
    plt.xticks([], [])
    plt.yticks([], [])
    plt.legend()

    plt.savefig("rc.png")
    plt.show()


rc_plot()
