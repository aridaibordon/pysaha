import constants as cn
import numpy as np

from saha import get_ionization_n_moment


def get_debye_length(ion_pop_list: list[float], t_elec: float, d_elec: float) -> float:
    mean_ion = get_ionization_n_moment(1, ion_pop_list)
    mean_ion2 = get_ionization_n_moment(2, ion_pop_list)

    return np.sqrt(
        (t_elec * cn.EPS0) / (4 * cn.PI * (d_elec / mean_ion) * (mean_ion + mean_ion2))
    )


def get_ionic_radius(nion, d_elec):
    return ((3 * (nion + 1)) / (4 * cn.PI * d_elec)) ** (1 / 3)


def get_sp_ipd_nion(
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


def get_sp_ipd_elem(
    elem: str,
    elem_z: int,
    t_elec: float,
    d_elec: float,
    ion_pop_list: list[float],
) -> list[float]:
    return [
        (
            get_sp_ipd_nion(elem, nion, ion_pop_list, t_elec, d_elec)
            if nion < elem_z
            else 0.0
        )
        for nion in range(elem_z + 1)
    ]