import json

import scipy.constants as sc


# constants
A0 = 5.29177210e-09  # cm
C = sc.c * 1e2  # cm s-1
ELECTRON_MASS = (sc.electron_mass * sc.c**2) / sc.elementary_charge  # eV
EPS0 = 5.52634940e05  # e2 eV-1 cm-1
H = sc.h / sc.elementary_charge  # eV s
PI = sc.pi

# units conversion
EV_TO_ERG = 1.602176e-12


with open("data/atomic/ionization.json", "r") as f:
    IONIZATION_ENERGY = json.load(f)
