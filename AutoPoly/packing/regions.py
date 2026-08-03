#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Carve regions for subtract-style packing post-processing.

A CarveRegion describes a volume of the simulation box. After placement,
any instance whose center falls inside a region is removed whole (no
covalent bonds are ever cut — removal granularity is one chain/molecule,
because moltemplate instantiates whole .lt classes). With
``conservative=True`` the region is effectively inflated by the instance's
bounding radius, so surviving instances are guaranteed not to intersect
the carved volume at all.

Regions carry an ``apply_to`` selector ("film", "substrate", or "all") so
e.g. a film-thickness cut never touches the substrate slab.

Strategies that support subtraction collect PlacedItem views of every
placement and run them through :func:`apply_carve_regions` before
emitting moltemplate commands.

Created on 2026-07-30
@author: zwu
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, List, Optional, Sequence, Tuple

import numpy as np

from ..core.exceptions import ValidationError
from ..core.system import logger

#: Valid carve targets (which placed instances a region may remove).
APPLY_TO_CHOICES = ("film", "substrate", "all")

Range = Optional[Tuple[float, float]]


class CarveRegion(ABC):
    """
    A box volume that removes placed instances (whole-instance subtract).

    Args:
        apply_to: "film" (default), "substrate", or "all".
        conservative: If True, an instance is removed when its bounding
                      sphere *touches* the region (center distance test
                      inflated by the instance radius). If False (default),
                      only the center point is tested.
    """

    def __init__(self, apply_to: str = "film", conservative: bool = False) -> None:
        if apply_to not in APPLY_TO_CHOICES:
            raise ValidationError(
                f"CarveRegion apply_to must be one of "
                f"{list(APPLY_TO_CHOICES)}, got '{apply_to}'"
            )
        self.apply_to = apply_to
        self.conservative = bool(conservative)

    def matches_role(self, role: str) -> bool:
        """True if this region applies to an instance with the given role."""
        return self.apply_to == "all" or self.apply_to == role

    @abstractmethod
    def contains(self, position: Sequence[float], radius: float = 0.0) -> bool:
        """
        Test whether an instance falls inside the carved volume.

        Args:
            position: Instance center (x, y, z) in Angstrom.
            radius: Instance bounding radius (used only when conservative).

        Returns:
            True if the instance should be removed.
        """
        raise NotImplementedError

    def _margin(self, radius: float) -> float:
        return radius if self.conservative else 0.0


class CutAbove(CarveRegion):
    """Remove everything with center z above a plane (film thickness cut)."""

    def __init__(self, z: float, apply_to: str = "film",
                 conservative: bool = False) -> None:
        super().__init__(apply_to=apply_to, conservative=conservative)
        self.z = float(z)

    def contains(self, position: Sequence[float], radius: float = 0.0) -> bool:
        return position[2] + self._margin(radius) > self.z

    def __repr__(self) -> str:
        return f"CutAbove(z={self.z}, apply_to='{self.apply_to}')"


class CutBelow(CarveRegion):
    """Remove everything with center z below a plane."""

    def __init__(self, z: float, apply_to: str = "film",
                 conservative: bool = False) -> None:
        super().__init__(apply_to=apply_to, conservative=conservative)
        self.z = float(z)

    def contains(self, position: Sequence[float], radius: float = 0.0) -> bool:
        return position[2] - self._margin(radius) < self.z

    def __repr__(self) -> str:
        return f"CutBelow(z={self.z}, apply_to='{self.apply_to}')"


class Cylinder(CarveRegion):
    """
    Remove everything inside an infinite cylinder (hole / nanopore).

    Args:
        axis: Cylinder axis: "x", "y", or "z".
        center: Cylinder center in the plane perpendicular to `axis`,
                e.g. (x, y) for axis="z". Defaults to the box center.
        radius: Cylinder radius in Angstrom (must be > 0).
    """

    _PLANAR_AXES = {"x": (1, 2), "y": (0, 2), "z": (0, 1)}

    def __init__(self, axis: str = "z",
                 center: Tuple[float, float] = (0.0, 0.0),
                 radius: float = 5.0,
                 apply_to: str = "film",
                 conservative: bool = False) -> None:
        super().__init__(apply_to=apply_to, conservative=conservative)
        if axis not in self._PLANAR_AXES:
            raise ValidationError(
                f"Cylinder axis must be 'x', 'y', or 'z', got '{axis}'"
            )
        if radius <= 0:
            raise ValidationError(
                f"Cylinder radius must be > 0, got {radius}"
            )
        self.axis = axis
        self.center = (float(center[0]), float(center[1]))
        self.radius = float(radius)

    def contains(self, position: Sequence[float], radius: float = 0.0) -> bool:
        i, j = self._PLANAR_AXES[self.axis]
        dist = np.hypot(position[i] - self.center[0],
                        position[j] - self.center[1])
        return dist - self._margin(radius) < self.radius

    def __repr__(self) -> str:
        return (
            f"Cylinder(axis='{self.axis}', center={self.center}, "
            f"radius={self.radius}, apply_to='{self.apply_to}')"
        )


class BoxRegion(CarveRegion):
    """
    Remove everything inside an axis-aligned box (trench / pattern).

    Args:
        x, y, z: (lo, hi) ranges in Angstrom, or None for "full span".
                 At least one axis must be bounded.
    """

    def __init__(self, x: Range = None, y: Range = None, z: Range = None,
                 apply_to: str = "film",
                 conservative: bool = False) -> None:
        super().__init__(apply_to=apply_to, conservative=conservative)
        self.ranges = []
        for axis_name, rng in (("x", x), ("y", y), ("z", z)):
            if rng is None:
                self.ranges.append(None)
                continue
            lo, hi = float(rng[0]), float(rng[1])
            if lo >= hi:
                raise ValidationError(
                    f"BoxRegion {axis_name} range must satisfy lo < hi, "
                    f"got ({lo}, {hi})"
                )
            self.ranges.append((lo, hi))
        if all(r is None for r in self.ranges):
            raise ValidationError(
                "BoxRegion requires at least one bounded axis "
                "(an unbounded box would carve the entire system)"
            )

    def contains(self, position: Sequence[float], radius: float = 0.0) -> bool:
        margin = self._margin(radius)
        for axis, rng in enumerate(self.ranges):
            if rng is None:
                continue
            lo, hi = rng
            if not (lo - margin < position[axis] < hi + margin):
                return False
        return True

    def __repr__(self) -> str:
        x, y, z = self.ranges
        return (
            f"BoxRegion(x={x}, y={y}, z={z}, apply_to='{self.apply_to}')"
        )


@dataclass
class PlacedItem:
    """
    A placed instance as seen by the subtract pass.

    Attributes:
        unit_id: ID of the UnitSpec this instance belongs to (the literal
                 string "external" for an external substrate slab).
        instance_name: Moltemplate instance name.
        role: "film" or "substrate".
        position: Instance center (x, y, z) in Angstrom.
        radius: Bounding radius in Angstrom.
        payload: ("polymer", PolymerPlacement) | ("molecule",
                 MoleculePlacement) | ("command", lt_command_string) —
                 whatever the strategy needs to emit the lt command after
                 filtering.
    """
    unit_id: str
    instance_name: str
    role: str
    position: Sequence[float]
    radius: float
    payload: Tuple[str, Any]


def apply_carve_regions(
    items: List[PlacedItem],
    regions: Optional[List[CarveRegion]],
) -> List[PlacedItem]:
    """
    Remove items falling inside any carve region (whole-instance subtract).

    Args:
        items: Placed instances (film and substrate).
        regions: Carve regions applied in order; None/empty = no-op.

    Returns:
        The kept items, in original order. Removals are logged per region
        (never silent).
    """
    if not regions:
        return list(items)

    kept: List[PlacedItem] = []
    removal_counts = {id(region): 0 for region in regions}
    for item in items:
        hit = None
        for region in regions:
            if region.matches_role(item.role) and \
                    region.contains(item.position, item.radius):
                hit = region
                break
        if hit is None:
            kept.append(item)
        else:
            removal_counts[id(hit)] += 1

    for region in regions:
        count = removal_counts[id(region)]
        if count:
            logger.info(
                f"Subtract: removed {count} instance(s) via {region!r}"
            )
    removed_total = len(items) - len(kept)
    if removed_total:
        logger.info(
            f"Subtract complete: {removed_total}/{len(items)} instances "
            f"removed, {len(kept)} kept"
        )
    return kept
