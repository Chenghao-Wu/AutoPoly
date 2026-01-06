# -*- coding: utf-8 -*-
"""
AutoPoly: Automated Polymer Generation and Simulation Package

This package provides tools for generating polymer structures and preparing them
for molecular dynamics simulations using LAMMPS and Moltemplate.

Main Components:
- System: Utility class for managing file paths and system operations
- Polymer: Class for defining polymer properties and sequences
- Polymerization: Core class for generating polymer structures using Moltemplate
- BeadSpringPolymer: Simplified bead-spring polymer model generator

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

# Explicit exports
__all__ = [
    "System",
    "Polymer",
    "Polymerization",
    "BeadSpringPolymer",
    "MonomerGenerator",
    "PolymerizationMechanism",
    "ConnectionPointModifier",
    "detect_mechanism",
    "get_pattern",
    "get_all_mechanisms",
    "__version__",
    "__author__",
]

from .system import System
from .polymer import Polymer
from .polymerization import Polymerization
from .bead_spring import BeadSpringPolymer
from .monomer_generator import MonomerGenerator
from .polymerization_mechanism import PolymerizationMechanism, detect_mechanism
from .connection_point import ConnectionPointModifier
from .polymerization_patterns import get_pattern, get_all_mechanisms