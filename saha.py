import os
import numpy as np

from fac_utils.read import read_lev

import constants as cn


ATOMIC_PATH = "data/atomic"


def get_de_broglie_thermal_length(t_elec: float) -> float:
    return (cn.H * cn.C) / np.sqrt(2 * cn.PI * cn.ELECTRON_MASS * t_elec)


def compute_start_ion(elem: str, elem_z: int, ipd_list: list) -> int:
    elem_ionization_energy = cn.IONIZATION_ENERGY[elem]

    for nion, ipd in enumerate(ipd_list):
        if ipd < elem_ionization_energy[nion]:
            return nion

    return elem_z


def get_internal_partition_function(
    elem: str,
    elem_z: int,
    nion: int,
    t_elec: float,
    ipd_elem_list: list[float] | None = None,
) -> float:
    if not ipd_elem_list:
        ipd = 0.0
    else:
        ipd = ipd_elem_list[nion]

    if nion == elem_z:
        return 1

    nelec = elem_z - nion
    ionization_energy = cn.IONIZATION_ENERGY[elem][nion]

    path_to_lev = os.path.join(ATOMIC_PATH, f"{elem}_{nelec:02d}_DCA.lev")
    tot_lev_df = read_lev(path_to_lev)
    lev_df = tot_lev_df[tot_lev_df["E"] < (ionization_energy - ipd)]

    return lev_df.apply(
        lambda row: (row["2J"] + 1) * np.exp(-row["E"] / t_elec), axis=1
    ).sum()


def solve_saha(
    elem: str,
    elem_z: int,
    t_elec: float,
    d_elec: float,
    ipd_elem_list: list[float] | None = None,
) -> list[float]:
    elem_ionization_energy = cn.IONIZATION_ENERGY[elem]

    de_broglie_length = get_de_broglie_thermal_length(t_elec)

    if ipd_elem_list:
        ion_start = compute_start_ion(elem, elem_z, ipd_elem_list)
    else:
        ion_start = 0

    pop_ion_list = [
        1.0 if nion >= ion_start else 0.0 for nion in range(elem_z + 1)
    ]

    for nion in range(ion_start + 1, elem_z + 1):
        if ipd_elem_list:
            prev_ipd = ipd_elem_list[nion - 1]
        else:
            prev_ipd = 0.0

        ionization_energy = elem_ionization_energy[nion - 1]

        z_up = get_internal_partition_function(
            elem, elem_z, nion, t_elec, ipd_elem_list
        )
        z_down = get_internal_partition_function(
            elem, elem_z, nion - 1, t_elec, ipd_elem_list
        )

        pop_ion_list[nion] = (
            2
            * (1 / de_broglie_length) ** 3
            * (1 / d_elec)
            * (z_up / z_down)
            * np.exp(-(ionization_energy - prev_ipd) / t_elec)
            * pop_ion_list[nion - 1]
        )
    return np.asarray([pop / sum(pop_ion_list) for pop in pop_ion_list])


def get_saha_rel_lev_pop(
    elem: str, elem_z: int, t_elec: float, d_elec: float, ipd_correction: float = 0.0
):
    pop_ion_list = solve_saha(elem, elem_z, t_elec, d_elec, ipd_correction)

    pop_lev_list = []

    for nion in range(0, elem_z + 1):
        nelec = elem_z - nion

        path_to_lev = os.path.join(ATOMIC_PATH, f"{elem}_{nelec:02d}_DCA.lev")
        lev_df = read_lev(path_to_lev)

        pop_ion = pop_ion_list[nion]
        boltzman_factor_list = lev_df.apply(
            lambda row: (row["2J"] + 1) * np.exp(-row["E"] / t_elec), axis=1
        ).values

        pop_lev_list.append(boltzman_factor_list * pop_ion / sum(boltzman_factor_list))

    return pop_lev_list


def get_mean_ionization(pop_ion_list: list[float]) -> float:
    return np.sum([zeta * pop for zeta, pop in enumerate(pop_ion_list)])


def get_mean_square_ionization(pop_ion_list: list[float]) -> float:
    return np.sum([zeta**2 * pop for zeta, pop in enumerate(pop_ion_list)])
