#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pluggable box-packing strategies for AutoPoly.

This package contains the stage-3 placement strategies that arrange typed
units (from a UnitLibrary manifest) inside the simulation box:

- base.py:     BoxSpec, PackingContext, PlacementResult, PlacementStrategy ABC
- registry.py: strategy registry (built-ins + user-defined at runtime)
- grid.py:     deterministic planar grid placement
- random_mc.py: Monte Carlo random placement with collision detection

Created on 2026-07-30
@author: zwu
"""
from .base import (
    BoxSpec,
    PackingContext,
    PlacementRecord,
    PlacementResult,
    PlacementStrategy,
)
from .registry import (
    get_strategy,
    register_strategy,
    registered_strategies,
    RESERVED_STRATEGY_NAMES,
)
from .grid import GridStrategy
from .random_mc import RandomMCStrategy

# Register built-in strategies.
register_strategy(GridStrategy.name, GridStrategy)
register_strategy(RandomMCStrategy.name, RandomMCStrategy)

__all__ = [
    "BoxSpec",
    "PackingContext",
    "PlacementRecord",
    "PlacementResult",
    "PlacementStrategy",
    "GridStrategy",
    "RandomMCStrategy",
    "get_strategy",
    "register_strategy",
    "registered_strategies",
    "RESERVED_STRATEGY_NAMES",
]
