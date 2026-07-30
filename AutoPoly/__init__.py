# -*- coding: utf-8 -*-
"""
AutoPoly: Automated Polymer Generation and Simulation Package

This package provides tools for generating polymer structures and preparing them
for molecular dynamics simulations using LAMMPS and Moltemplate.

Main Components:
- System: Utility class for managing file paths and system operations
- Polymer: Class for defining polymer properties and sequences
- Polymerization: Core class for generating polymer structures using Moltemplate
- Molecule: Class for defining small molecule structures (water, benzene, etc.)
- BeadSpringPolymer: Simplified bead-spring polymer model generator
- mc: Monte Carlo placement module for chain growth and molecular placement

Three-Stage Pipeline (Geometry -> Typing -> Packing):
- GeometryBuilder: Stage 1, force-field-agnostic coordinates -> geometry.json
- UnitTyper: Stage 2, force-field assignment -> build/<ff>/ + units.json
- BoxPacker: Stage 3, box packing + moltemplate -> system.data
- UnitLibrary: units.json manifest contract between stages 2 and 3
- PlacementStrategy / register_strategy: pluggable packing strategies

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
    "Molecule",
    "BeadSpringPolymer",
    "BeadType",
    "AngleType",
    "MCConfig",
    "SAWConfig",
    "MonomerGenerator",
    "GeometryBuilder",
    "GeometryConfig",
    "UnitTyper",
    "BoxPacker",
    "UnitLibrary",
    "UnitSpec",
    "PlacementStrategy",
    "register_strategy",
    "get_strategy",
    "mc",
    "agent",
    "__version__",
    "__author__",
]

# Always available imports (no rdkit dependency)
from .system import System
from .bead_spring import BeadSpringPolymer, BeadType, AngleType, MCConfig, SAWConfig
from . import agent

# Optional imports that require rdkit
try:
    from .polymer import Polymer
    from .polymerization import Polymerization
    from .molecule import Molecule
    from .monomer_generator import MonomerGenerator
    from .geometry import GeometryBuilder, GeometryConfig
    from .typing import UnitTyper
    from .packer import BoxPacker
    from .units import UnitLibrary, UnitSpec
    from .packing import PlacementStrategy, register_strategy, get_strategy
    from . import mc
except ImportError as e:
    import warnings
    warnings.warn(
        f"Some AutoPoly components require rdkit: {e}. "
        "Install rdkit for full functionality. "
        "BeadSpringPolymer and System are still available."
    )
    Polymer = None
    Polymerization = None
    Molecule = None
    MonomerGenerator = None
    mc = None