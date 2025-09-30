n2p2 - A neural network potential package
=========================================

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1344446.svg)](https://doi.org/10.5281/zenodo.1344446)
[![GitHub release](https://img.shields.io/github/release/CompPhysVienna/n2p2.svg)](https://GitHub.com/CompPhysVienna/n2p2/releases/)
[![Build Status](https://app.travis-ci.com/CompPhysVienna/n2p2.svg?branch=master)](https://app.travis-ci.com/CompPhysVienna/n2p2)
[![Coverage](https://codecov.io/gh/CompPhysVienna/n2p2/branch/master/graph/badge.svg)](https://codecov.io/gh/CompPhysVienna/n2p2)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This repository provides ready-to-use software for high-dimensional neural
network potentials in computational physics and chemistry.

# [Documentation](http://compphysvienna.github.io/n2p2)

This package uses automatic documentation generation via
[Sphinx](http://www.sphinx-doc.org),
[Breathe](https://breathe.readthedocs.io/en/latest/#) and
[doxygen](http://www.doxygen.nl/). An online version of the documentation which
is automatically updated with the master branch of the repository can be found
[__here__](http://compphysvienna.github.io/n2p2).

## LAMMPS Q-HDNNP Charges

The `pair_style hdnnp` interface now forwards per-atom charges predicted by
Q-HDNNP models directly into LAMMPS. Ensure `atom_style charge` (or a derived
style) is active so the `q` attribute is allocated, then inspect the values via
standard tooling, for example:

- `compute myCharge all property/atom q`
- `dump 1 all custom 1 charges.dump id type q`
- `thermo_style custom step pe etotal` (for checking energy consistency)

The charges are updated every time `pair_hdnnp` invokes the network evaluation,
which allows downstream polarization or post-processing workflows to consume
self-consistent NN charge predictions without any additional coupling code.

# Authors

See [AUTHORS.rst](https://github.com/CompPhysVienna/n2p2/blob/master/AUTHORS.rst) for a list of contributions.

# License

This software is licensed under the [GNU General Public License version 3 or any later version (GPL-3.0-or-later)](https://www.gnu.org/licenses/gpl.txt).
