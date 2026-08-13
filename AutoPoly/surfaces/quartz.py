#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Crystalline alpha-quartz (0001) substrate slab builder.

The orthorhombic surface cell (a x a*sqrt(3) x c) is derived from the
alpha-quartz structure of Tucker, Dove & Keen (Mineralogical Magazine
2001, 65, 489; COD 1526860), space group P3121 (No. 152):
a = 4.9019 A, c = 5.3988 A.  Each Si is tetrahedrally coordinated
(Si-O = 1.607-1.616 A) and every oxygen bridges two silicons.

The (0001) cleavage leaves two dangling oxygens on each surface silicon
(geminal Q2 termination, ~9.6 silanols/nm^2); the surface silicons are
~4.9 A apart, so Si-O-Si bridging cannot reduce this density — the
builder warns when a lower oh_density target cannot be met.

Created on 2026-08-12
@author: zwu
"""
import math

from .base import SlabBuilder, SurfaceSlab

#: Backwards-compatible alias.
QuartzSlab = SurfaceSlab

#: Lattice constant a of the orthorhombic quartz(0001) surface cell (A).
QUARTZ_A = 4.9019
#: Lattice constant c (slab z direction) (A).
QUARTZ_C = 5.3988
#: Lateral b of the orthorhombic surface cell (A) = a*sqrt(3).
QUARTZ_B = QUARTZ_A * math.sqrt(3.0)

#: Fractional coordinates of the 18-atom orthorhombic quartz(0001) cell
#: (derived from COD 1526860; verified: 4 Si-O bonds per Si at
#: 1.607-1.616 A, O-O >= 2.62 A, Si-Si >= 3.05 A).
ORTHO_CELL_FRACTIONAL = (
    ("Si", (0.46730000, 0.00000000, 0.33330000)),
    ("Si", (0.96730000, 0.50000000, 0.33330000)),
    ("Si", (0.76635000, 0.23365000, 0.66663333)),
    ("Si", (0.26635000, 0.73365000, 0.66663333)),
    ("Si", (0.26635000, 0.26635000, 0.99996667)),
    ("Si", (0.76635000, 0.76635000, 0.99996667)),
    ("O",  (0.15795000, 0.42905000, 0.11613333)),
    ("O",  (0.65795000, 0.92905000, 0.11613333)),
    ("O",  (0.27745000, 0.13555000, 0.21720000)),
    ("O",  (0.77745000, 0.63555000, 0.21720000)),
    ("O",  (0.77745000, 0.36445000, 0.44946667)),
    ("O",  (0.27745000, 0.86445000, 0.44946667)),
    ("O",  (0.65795000, 0.07095000, 0.55053333)),
    ("O",  (0.15795000, 0.57095000, 0.55053333)),
    ("O",  (0.06460000, 0.20650000, 0.78280000)),
    ("O",  (0.56460000, 0.70650000, 0.78280000)),
    ("O",  (0.56460000, 0.29350000, 0.88386667)),
    ("O",  (0.06460000, 0.79350000, 0.88386667)),
)


class QuartzBuilder(SlabBuilder):
    """Hydroxylated alpha-quartz (0001) slab; see SlabBuilder for args."""

    CELL_FRACTIONAL = ORTHO_CELL_FRACTIONAL
    CELL_A = QUARTZ_A
    CELL_B = QUARTZ_B
    CELL_C = QUARTZ_C
    BUILDER_NAME = "alpha_quartz"
    ATOM_PREFIX = "qz"
    DEFAULT_CLASS_NAME = "QuartzSlab"
    DEFAULT_LT_FILENAME = "alpha_quartz_slab.lt"
