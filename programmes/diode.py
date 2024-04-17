import matplotlib.pyplot as plt
import numpy as np

r_interne = 0
i_s = 25 * (10 ** (-9))
R = np.array([68, 200, 10000]) + r_interne
V_D_mesure = np.array([0.883, 0.795, 0.618])

t = np.linspace(0, 6, 400)
courbe_diode = i_s * (np.exp(t / (2 * 0.026)) - 1)
I_calcul = (6 - V_D_mesure) / R


for r in R:
    plt.plot(t, (6 - t) / r, "-k")
plt.plot(V_D_mesure, I_calcul, ".b", markersize=20)
plt.plot(t, courbe_diode, "-g")


plt.xlim([0, 2])
plt.ylim([-(10 ** (-4)), 10 ** (-1)])
plt.show()
