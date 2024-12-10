# PySaha: a Saha ionization equation solver for Python.

### Requirements
This repository requieres a Python version >= 12.0 and the [fac-utils](https://github.com/aridaibordon/fac-utils) library.

## Description
PySaha is a python implementation of the Saha ionization equation to solve ionic populations for monocomponent homogeneous plasmas in local thermal equilibrium. Atomic energy structure is obtained using the [FACOF](https://github.com/aridaibordon/facof) interface of the [Flexible Atomic Code](https://github.com/flexible-atomic-code/fac).

Corrections due to plasma environment to the isolated atomic levels are included as a static correction to the ionization energy. All models available are included in the ipd.py file.

## Usage
The main functionality is included in the `solve_saha` function that takes as parameters the element symbol, its atomic number and the electronic temperature and density (expressed in eV and particles per cubic centimeter).
