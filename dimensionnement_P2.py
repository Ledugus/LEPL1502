from numpy import *


class Resistance:
    def __init__(self, resistance: float) -> None:
        self.resistance = resistance


def get_values(parameters):
    values = {}
    print("Entrez les valeurs sous le format float :")
    for param_name in parameters:
        while True:
            try:
                val = float(input(f"{param_name} : "))
            except ValueError:
                continue
            break
        values[param_name] = val
    return values


def get_optimal_frequency(mutual: float, r1: Resistance, r2: Resistance):
    """Retourne la fréquence optimale pour le transfert de puissance"""
    return sqrt(r1.resistance * r2.resistance) / mutual


def main():
    parameters = ["R1", "R2", "M", "L1", "L2"]
    R1, R2, M, L1, L2 = get_values(parameters).values()

    goal_f = get_optimal_frequency(M, R1, R2)


if __name__ == "__main__":
    main()
