from numpy import *
import numpy as np
from scipy.special import roots_legendre
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------------
#
# Calcul des composantes radiales et verticales
# du champ magnétique généré
# par une série de boucles circulaires de courants !
#
#    X,Z : tableaux numpy des coordonnées x,z des points en [m]
#    Rsource : liste des rayons des boucles en [m]
#    Zsource : liste de hauteurs des boucles en [m]
#    Isource : liste des courants des boucles en [A]
#    data : structure contenant les paramètres matériels
#
#  La fonction renvoie une liste [Bx,Bz] où Bx,Bz sont des tableaux/variables
#  du même type que X,Z.   Le résultat est exprimé en [Tesla] ou [kg /s2 A]
#


def inductanceMegneticField(X, Z, Rsource, Zsource, Isource, data):

    mu0 = data.mu0
    Xcircle = data.Xcircle
    Ycircle = data.Ycircle
    n = len(Rsource)
    dtheta = (2 * pi) / len(Xcircle)

    Bx = zeros(shape(X))
    Bz = zeros(shape(X))
    for k in range(n):
        for l in range(len(Xcircle)):
            dx = -Rsource[k] * Ycircle[l] * dtheta
            dy = Rsource[k] * Xcircle[l] * dtheta
            rx = X - Rsource[k] * Xcircle[l]
            ry = -Rsource[k] * Ycircle[l]
            rz = Z - Zsource[k]
            r = (rx * rx + ry * ry + rz * rz) ** (3 / 2)
            Bx += Isource[k] * (dy * rz) / r
            Bz += Isource[k] * (dx * ry - dy * rx) / r
    coeff = (mu0) / (4 * pi)
    Bx *= coeff
    Bz *= coeff

    return [Bx, Bz]


# ------------------------------------------------------------------------------------
#
# Calcul des points et points d'intégration
# pour une intégration double de Gauss-Legendre
# afin de calculer les flux du champs magnétique
#
#    X0,Xf : intervalle radial d'intégration
#    Z0,Zf = intervalle vertical de moyenne du flux
#    nX : nombre de sous-intervalles en x
#    nZ : nombre de sous-intervalles en z
#    nGL : nombre de points de Gauss-Legendre par intervalle (en x et en z)
#
#  La fonction renvoie une liste [X,Z,W] contenant les abscisses et les poids.
#  Ce seront des tableaux unidimensionnels de taille nX*nGL*nZ*nGL
#  Si nX ou nZ sont nuls, on fera une intégration unidimensionnelle en utilisant
#  X0 ou Z0 et les tableaux auront une taille nZ*nGL, nX*nGL ou nGL respectivement
#


def inductanceGaussLegendre(X0, Xf, Z0, Zf, n, m, nGaussLegendre):

    xi, we = roots_legendre(nGaussLegendre)

    X = X0 * ones(nGaussLegendre)
    Z = Z0 * ones(nGaussLegendre)
    W = ones(nGaussLegendre) / nGaussLegendre

    if n > 0:
        X = zeros(nGaussLegendre * n)
        Z = Z0 * ones(nGaussLegendre * n)
        W = zeros(nGaussLegendre * n)
        Xnode = linspace(X0, Xf, n + 1)
        h = (Xf - X0) / n
        for i in range(n):
            Xlocal = Xnode[i] + h / 2 + xi * h / 2
            map = range(i * nGaussLegendre, (i + 1) * nGaussLegendre)
            X[map] = Xlocal
            W[map] = Xlocal * we * pi * h

    if m > 0:
        V = W
        Z = zeros(nGaussLegendre * m)
        W = zeros(nGaussLegendre * m)
        Znode = linspace(Z0, Zf, m + 1)
        h = (Zf - Z0) / m
        for i in range(m):
            Zlocal = Znode[i] + h / 2 + xi * h / 2
            map = range(i * nGaussLegendre, (i + 1) * nGaussLegendre)
            Z[map] = Zlocal
            W[map] = we / (2 * m)
        if n == 0:
            X = repeat(X, m)
        else:
            Z = tile(Z, nGaussLegendre * max(n, 1))
            X = repeat(X, nGaussLegendre * m)
            W = repeat(V, nGaussLegendre * m) * tile(W, nGaussLegendre * max(n, 1))

    return [X, Z, W]


# ------------------------------------------------------------------------------------
#
# Calcul des points et points d'intégration
# pour une intégration double de Simpson
# afin de calculer les flux du champs magnétique
#
#    X0,Xf : intervalle radial d'intégration
#    Z0,Zf = intervalle vertical de moyenne du flux
#    nX : nombre de sous-intervalles en x
#    nZ : nombre de sous-intervalles en z
#
#  La fonction renvoie une liste [X,Z,W] contenant les abscisses et les poids.
#  Ce seront des tableaux unidimensionnels de taille (2nX+1)*(2nZ+1)
#  Si nX ou nZ sont nuls, on fera une intégration unidimensionnelle en utilisant
#  X0 ou Z0 et les tableaux auront une taille (2nZ+1), (2nX/1) ou 1 respectivement
#


def inductanceSimpson(X0, Xf, Z0, Zf, n, m):

    xi = array([-1.0, 0, 1.0])
    we = array([1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0])

    size = (2 * n + 1) * (2 * m + 1)
    Xnode = linspace(X0, Xf, 2 * n + 1)
    Znode = linspace(Z0, Zf, 2 * m + 1)
    weX = ones(2 * n + 1)
    weZ = ones(2 * m + 1)
    weX[:] = 1.0 / 3.0
    weX[1::2] = 4.0 / 3.0
    weX[2:-2:2] = 2.0 / 3.0
    weZ[:] = 1.0 / 3.0
    weZ[1::2] = 4.0 / 3.0
    weZ[2:-2:2] = 2.0 / 3.0

    X = X0 * ones(size)
    Z = Z0 * ones(size)
    W = ones(size)

    if n > 0 and m > 0:
        h = (Xf - X0) / n
        for i in range(2 * n + 1):
            for j in range(2 * m + 1):
                Xlocal = X0 + (Xf - X0) * i / (2 * n)
                Zlocal = Z0 + (Zf - Z0) * j / (2 * m)
                Wlocal = Xlocal * pi
                X[(2 * n + 1) * j + i] = Xlocal
                Z[(2 * n + 1) * j + i] = Zlocal
                W[(2 * n + 1) * j + i] = Wlocal * weX[i] * h * weZ[j] / (2 * m)
    elif n > 0:
        h = (Xf - X0) / n
        for i in range(2 * n + 1):
            Xlocal = X0 + (Xf - X0) * i / (2 * n)
            Zlocal = Z0
            Wlocal = Xlocal * pi
            X[i] = Xlocal
            Z[i] = Zlocal
            W[i] = Wlocal * weX[i] * h
    elif m > 0:
        for j in range(2 * m + 1):
            Xlocal = X0
            Zlocal = Z0 + (Zf - Z0) * j / (2 * m)
            Wlocal = Xlocal * pi
            X[j] = Xlocal
            Z[j] = Zlocal
            W[j] = Wlocal * weZ[j] / (2 * m)
    return [X, Z, W]


class MutualInductanceProject:
    mu0 = 4e-7 * pi  # permeabilité du vide en [H/m]

    Rcoil = 1.0e-2  # rayon des bobines [m]
    Hcoil = 1.0e-2  # épaisseur des bobines [m]
    Zcoil = -Hcoil / 2  # position verticale de la bobine primaire [m]
    nSpires = 190  # nombre de spires
    I = 1  # courant dans la bobine [A]
    Hsecond = 2e-2  # position relative de la bobine secondaire [m]

    nZcoil = 100  # discrétisation vertical de la bobine
    nTheta = 50  # discrétisation du cercle
    Theta = linspace(0, 2 * pi, nTheta + 1)[:-1] + pi / nTheta
    Xcircle = cos(Theta)
    Ycircle = sin(Theta)

    def streamPlot(self, X, Z, Bx, Bz, colorStream, colorIplate):
        plt.streamplot(X, Z, Bx, Bz, density=1.0, linewidth=None, color=colorStream)
        x = array([-self.Rcoil, self.Rcoil])
        z = repeat([linspace(0, -self.Hcoil, self.nZcoil)], 2, axis=0)
        plt.plot(x, z, "o-r", linewidth=2)
        plt.plot(x, z + self.Hsecond, "o-b", linewidth=2)
        plt.xlim((-0.03, 0.03))
        plt.ylim((-0.03, 0.03))


def main():

    # ------------------------------------------------------------------------------------
    #
    # Script de test
    #
    #
    # -0- Initialisation du projet
    #
    # ------------------------------------------------------------------------------------

    p = MutualInductanceProject()

    # ------------------------------------------------------------------------------------
    #
    # -1- Calcul de l'inductance, inductance mutuelle, inductance induite
    #     et inductance effective pour Mr. Oestges !
    #
    # ------------------------------------------------------------------------------------

    Zcoil = -linspace(0, p.Hcoil, p.nZcoil)
    Rcoil = p.Rcoil * ones_like(Zcoil)
    Icoil = p.I * ones_like(Zcoil) * p.nSpires / p.nZcoil

    #
    # -1.1- Flux dans la bobine primaire du champs induit par la bobine primaire
    #

    nX = 20
    nZ = 5
    nGL = 3
    X, Z, W = inductanceGaussLegendre(0, p.Rcoil, 0, -p.Hcoil / 2, nX, nZ, nGL)
    Bz = inductanceMegneticField(X, Z, Rcoil, Zcoil, Icoil, p)[1]
    flux = W @ Bz
    L = p.nSpires * (flux / p.I)
    print("Using Gauss-Legendre integration rule (nX = %d nZ = %d):-)" % (nX, nZ))
    print("Flux across primary        = %14.7e [Tesla m2]" % flux)
    print("Inductance                 = %14.7e [Henry]" % L)

    #
    # -1.2- Flux dans la bobine secondaire du champs induit par la bobine primaire
    #

    nX = 20
    nZ = 5
    nGL = 3
    X, Z, W = inductanceGaussLegendre(
        0, p.Rcoil, p.Hsecond, -p.Hcoil + p.Hsecond, nX, nZ, nGL
    )
    Bz = inductanceMegneticField(X, Z, Rcoil, Zcoil, Icoil, p)[1]
    flux = W @ Bz
    M = p.nSpires * (flux / p.I)
    print("Flux across secondary      = %14.7e [Tesla m2]" % flux)
    print("Mutual inductance          = %14.7e [Henry]" % M)

    # ------------------------------------------------------------------------------------
    #
    # -6- La jolie courbe pour Claude Oestges et Nicolas Roisin
    #
    # ------------------------------------------------------------------------------------

    nD = 21
    nX = 20
    nZ = 5
    nGL = 3
    distances = linspace(0, 1e-1, nD)
    M = zeros(nD)
    for i in range(nD):
        h = distances[i] + 0.01
        X, Z, W = inductanceGaussLegendre(0, p.Rcoil, h, h - p.Hcoil, nX, nZ, nGL)
        Bz = inductanceMegneticField(X, Z, Rcoil, Zcoil, Icoil, p)[1]
        flux = W @ Bz
        M[i] = p.nSpires * (flux / p.I)

    def m_sans_bord(distance):
        m = 0
        for x in range(p.nSpires):
            for j in range(p.nSpires):
                m += 1 / (
                    2
                    * (p.Rcoil**2 + (distance + (j + x) * p.Hcoil / p.nSpires) ** 2)
                    ** (3 / 2)
                )
        return p.mu0 * pi * m * p.Rcoil**4

    distances_exp = [0, 1, 2, 3, 4, 5, 6, 10]
    mutuelles_exp = array([1.89, 1.68, 1.58, 1.56, 1.55, 1.54, 1.535, 1.53])
    m = m_sans_bord(distances)
    plt.figure("Inductuance mutuelle en fonction de la distance")
    plt.plot(distances * 1e2, M * 1e3, "ob-", label="Méthode numérique")
    plt.plot(distances * 1e2, m * 1e3, "or-", label="M(d) : Effets de bord négligés")
    plt.plot(
        distances_exp,
        (mutuelles_exp - (0.77 + 0.75)) / 2,
        "og-",
        label="Valeurs expérimentales",
    )
    ax = plt.gca()
    ax.set_ylabel("Inductance mutuelle : M [mH]")
    ax.set_xlabel("Distance bobines : d [cm]")
    ax.grid()
    plt.legend()

    plt.show()


main()
