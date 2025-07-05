# -*- coding: utf-8 -*-
"""
AutoPoly: Automated Polymer Generation and Simulation Package

This package provides tools for generating polymer structures and preparing them
for molecular dynamics simulations using LAMMPS and Moltemplate.

Main Components:
- Polymer: Class for defining polymer properties and sequences
- Polymerization: Core class for generating polymer structures using Moltemplate
- BeadSpringPolymer: Simplified bead-spring polymer model generator
- System: Utility class for managing file paths and system operations

External Dependencies:
- Moltemplate: For generating LAMMPS data files from molecular templates
- LAMMPS: Molecular dynamics simulation engine
- OPLS-AA: Force field parameters for atomistic simulations

Created on Fri Dec 21 12:19:08 2018
@author: zwu
"""

__author__ = "Zhenghao Wu"
__license__ = "BSD License"
__version__ = "1.0.0"
__description__ = "Automated Polymer Generation and Simulation Package"

from .system import *
from .polymer import *
from .polymerization import *
from .extern.rdlt import *