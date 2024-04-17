import matplotlib.pyplot as plt
import numpy as np

R = 200
C = 10 ** (-6)
V_0 = 3
L = np.array([10 ** (-2), 10 ** (-1), 10 ** (0)])
f = np.linspace(0, 4000, 1000)
w = 2 * np.pi * f

colors = ["r", "g", "b"]


def i_max(w, l, r):
    return V_0 / np.sqrt(r**2 + (1 / (C * w) - l * w) ** 2)


for l, color in zip(L, colors):
    plt.plot(w, i_max(w, l, R), color=color, label=f"L = {l*10**3} $mH$")
    plt.axvline(x=1 / np.sqrt(l * C), color=color, linestyle="--")
plt.legend()
plt.xlabel("Fréquence (Hz)")
plt.xscale("log")
plt.ylabel("Courant max (A)")
plt.show()
