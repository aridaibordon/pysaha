import pysaha.constants as cn
import numpy as np

from pysaha.saha import get_mean_ionization, get_mean_square_ionization


def get_debye_length(ion_pop_list: list[float], t_elec: float, d_elec: float) -> float:
    mean_ion = get_mean_ionization(ion_pop_list)
    mean_ion2 = get_mean_square_ionization(ion_pop_list)

    return np.sqrt(
        (t_elec * cn.EPS0) / (4 * cn.PI * (d_elec / mean_ion) * (mean_ion + mean_ion2))
    )


def get_ionic_radius(nion, d_elec):
    return ((3 * (nion + 1)) / (4 * cn.PI * d_elec)) ** (1 / 3)


def sp_ion(
    elem: str, nion: int, ion_pop_list: list[float], t_elec: float, d_elec: float
) -> float:
    h_ionization_energy = cn.IONIZATION_ENERGY["H"][0]

    debye_length = get_debye_length(ion_pop_list, t_elec, d_elec)
    R_zeta = get_ionic_radius(nion, d_elec)

    ratio_debye_R = debye_length / R_zeta

    return (
        3
        * h_ionization_energy
        * (cn.A0 / R_zeta)
        * (nion + 1)
        * ((1 + ratio_debye_R**3) ** (2 / 3) - ratio_debye_R**2)
    )


def sp(
    elem: str,
    elem_z: int,
    t_elec: float,
    d_elec: float,
    ion_pop_list: list[float],
) -> list[float]:
    return [
        (sp_ion(elem, nion, ion_pop_list, t_elec, d_elec) if nion < elem_z else 0.0)
        for nion in range(elem_z + 1)
    ]


def sp_average(
    elem: str,
    elem_z: int,
    t_elec: float,
    d_elec: float,
    ion_pop_list: list[float],
):
    h_ionization_energy = cn.IONIZATION_ENERGY["H"][0]
    mean_ion = get_mean_ionization(ion_pop_list)

    debye_length = get_debye_length(ion_pop_list, t_elec, d_elec)
    mean_R = get_ionic_radius(mean_ion, d_elec)
    ratio_debye_R = debye_length / mean_R

    sp_average_ipd = (
        3
        * h_ionization_energy
        * (cn.A0 / mean_R)
        * (mean_ion + 1)
        * ((1 + ratio_debye_R**3) ** (2 / 3) - ratio_debye_R**2)
    )

    return [sp_average_ipd for nion in range(elem_z + 1)]
