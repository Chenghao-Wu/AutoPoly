#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Monte Carlo Placement Module for AutoPoly

This module provides Monte Carlo methods for polymer chain growth and
molecular placement in simulation boxes with collision detection.

Key Features:
- Self-Avoiding Random Walk (SAW) for chain growth
- Random molecular placement with collision avoidance
- Efficient spatial hashing via cell-linked lists
- Moltemplate command generation

Main Classes:
- CollisionDetector: Efficient 3D collision detection using cell-linked lists
- ChainGrowthMC: Self-avoiding random walk for building polymer chains
- MolecularPlacementMC: Random placement of polymers/molecules in simulation box

Example Usage:
    >>> from AutoPoly.mc import CollisionDetector, ChainGrowthMC, MolecularPlacementMC
    >>> from AutoPoly.mc import calculate_box_size
    >>>
    >>> # Calculate box size for 100 monomers
    >>> box_size = calculate_box_size(100, monomer_density=0.1)
    >>> bounds = ((-box_size/2, box_size/2),) * 3
    >>>
    >>> # Initialize collision detector
    >>> detector = CollisionDetector(bounds, cell_size=5.0)
    >>>
    >>> # Chain growth for a polymer
    >>> chain_mc = ChainGrowthMC(detector, max_attempts=1000)
    >>> placements = chain_mc.grow_chain(monomer_lt_files, chain_id=0)
    >>> commands = chain_mc.generate_lt_commands(placements)
    >>>
    >>> # Random placement of polymers
    >>> placer = MolecularPlacementMC(bounds, detector, max_attempts=10000)
    >>> polymer_placements = placer.place_all_polymers(polymer_specs)
    >>> system_commands = placer.generate_polymer_lt_commands(polymer_placements)

Created on 2026-01-29
@author: zwu
"""

from .collision import (
    MonomerSphere,
    CollisionDetector,
    calculate_box_size,
)

from .chain_growth import (
    AtomData,
    MonomerTemplate,
    MonomerPlacement,
    ChainGrowthMC,
    parse_lt_file,
    rotation_matrix_from_axis_angle,
    rotation_matrix_align_vectors,
    rotation_matrix_to_axis_angle,
    random_rotation_matrix,
)

from .placement import (
    PolymerPlacement,
    MoleculePlacement,
    MolecularPlacementMC,
)

__all__ = [
    # Collision detection
    "MonomerSphere",
    "CollisionDetector",
    "calculate_box_size",
    # Chain growth
    "AtomData",
    "MonomerTemplate",
    "MonomerPlacement",
    "ChainGrowthMC",
    "parse_lt_file",
    "rotation_matrix_from_axis_angle",
    "rotation_matrix_align_vectors",
    "rotation_matrix_to_axis_angle",
    "random_rotation_matrix",
    # Placement
    "PolymerPlacement",
    "MoleculePlacement",
    "MolecularPlacementMC",
]
