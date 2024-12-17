# PySaha: a Saha ionization equation solver for Python.

### Requirements
This repository requieres a Python version >= 12.0 and the [fac-utils](https://github.com/aridaibordon/fac-utils) library.

## Description
PySaha is a python implementation of the Saha ionization equation to solve ionic populations for monocomponent homogeneous plasmas in local thermal equilibrium. Atomic energy structure is obtained using the [FACOF](https://github.com/aridaibordon/facof) interface of the [Flexible Atomic Code](https://github.com/flexible-atomic-code/fac).

Corrections due to plasma environment to the isolated atomic levels are included as a static correction to the ionization energy.

## Usage
The main functionality is included in the `solve_saha` function that takes as parameters the element symbol, its atomic number and the electronic temperature and density (expressed in eV and particles per cubic centimeter).

Different standard models for the ionization potential depression are included in the ipd.py file.

**Code example**
```python
from pysaha import solve_saha, ipd


POP_THRESHOLD = 1e-3

# define element parameters and plasma conditions
elem, elem_z = "C", 6
t_elec, d_elec = 50, 1e23

# run saha solver without ionization potential depression
prev_pop = solve_saha(elem, elem_z, t_elec, d_elec)

# run saha solver until populations converge
while True:
    pop_ipd = ipd.sp(elem, elem_z, t_elec, d_elec, prev_pop)
    pop = solve_saha(elem, elem_z, t_elec, d_elec, pop_ipd)

    if abs(pop - prev_pop).max() < POP_THRESHOLD:
        break

    prev_pop = pop
```