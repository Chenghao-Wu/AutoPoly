#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Base classes for AutoPoly box-packing strategies.

A placement strategy takes a PackingContext (typed units + box spec) and
returns a PlacementResult: one lt instantiation command per placed instance
plus the box bounds actually used. Strategies read unit metadata only —
they never parse .lt files.

Created on 2026-07-30
@author: zwu
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import ClassVar, List, Optional, Tuple

from ..units import UnitLibrary

# ((xmin, xmax), (ymin, ymax), (zmin, zmax))
BoxBounds = Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]]


def symmetric_bounds(box_size: float) -> BoxBounds:
    """Build box bounds centered on the origin from a cubic box size."""
    half = box_size / 2.0
    return ((-half, half), (-half, half), (-half, half))


def bounds_size(bounds: BoxBounds) -> float:
    """Side length of cubic bounds."""
    return bounds[0][1] - bounds[0][0]


@dataclass
class BoxSpec:
    """
    Box specification handed to a strategy.

    Attributes:
        requested_box_size: Explicit cubic box side in Angstrom, or None
                            to let the strategy auto-size the box.
    """
    requested_box_size: Optional[float] = None


@dataclass
class PackingContext:
    """
    Everything a placement strategy needs.

    Attributes:
        units: Typed unit manifest from stage 2.
        box: Box specification (explicit size or auto).
        mc_max_attempts: Maximum random placement attempts per instance.
        rng_seed: Optional seed for reproducible stochastic strategies.
        offset: Nominal monomer-monomer spacing in Angstrom (grid layout,
                chain-length box estimates).
        monomer_density: Target monomer density (monomers/A^3) for box sizing.
    """
    units: UnitLibrary
    box: BoxSpec = field(default_factory=BoxSpec)
    mc_max_attempts: int = 10000
    rng_seed: Optional[int] = None
    offset: float = 4.0
    monomer_density: float = 0.085


@dataclass
class PlacementRecord:
    """
    One placed instance.

    Attributes:
        unit_id: ID of the UnitSpec this instance belongs to.
        instance_name: Moltemplate instance name ("polymer_1", "molecule_3").
        lt_command: Full moltemplate instantiation command line.
    """
    unit_id: str
    instance_name: str
    lt_command: str


@dataclass
class PlacementResult:
    """
    Result of a placement strategy run.

    Attributes:
        records: Placement commands in system.lt emission order.
        box_bounds: Box bounds actually used (after auto-sizing).
    """
    records: List[PlacementRecord]
    box_bounds: BoxBounds


class PlacementStrategy(ABC):
    """
    Abstract base class for box-packing strategies.

    Subclasses set the class attribute `name` and implement `place()`.
    Register custom strategies at runtime with
    ``AutoPoly.packing.register_strategy(name, cls)``.

    Reserved built-in names for future sub-projects: ``grafted_surface``,
    ``nanopore``.
    """

    name: ClassVar[str] = ""

    @abstractmethod
    def place(self, ctx: PackingContext) -> PlacementResult:
        """
        Compute placements for all units in the context.

        Args:
            ctx: Packing context with units, box spec, and MC config.

        Returns:
            PlacementResult with per-instance lt commands and box bounds.
        """
        raise NotImplementedError
