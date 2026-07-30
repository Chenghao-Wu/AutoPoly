#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Workflow Convenience Function for AutoPoly Package

This module provides :func:`generate`, a thin one-shot convenience function
that composes the three pipeline stages:

    GeometryBuilder  (stage 1: FF-agnostic coordinates -> geometry/)
    UnitTyper        (stage 2: force-field typing      -> build/<ff>/)
    BoxPacker        (stage 3: box packing + moltemplate -> moltemplate/)

All stage logic lives in the stage modules (geometry.py, typing.py,
packer.py, packing/). Use the stage classes directly when you need control
between stages (e.g. typing one geometry under multiple force fields); use
generate() when a single force field one-shot run is all you need.

Created on 2026-01-06
@author: zwu
"""
from typing import List, Optional

from ..core.system import logger
from ..core.exceptions import WorkflowError
from ..packing import PlacementResult
from .geometry import GeometryBuilder, GeometryConfig
from .typing import UnitTyper
from .packer import BoxPacker

__all__ = ["generate"]


def generate(
    system: object,
    name: str,
    models: List[object],
    force_field: str = "oplsaa",
    *,
    strategy: str = "mc_random",
    geometry_config: Optional[GeometryConfig] = None,
    box_size: Optional[float] = None,
    mc_max_attempts: int = 10000,
    monomer_density: float = 0.085,
    rng_seed: Optional[int] = None,
    run_moltemplate: bool = True,
) -> PlacementResult:
    """
    Generate a complete LAMMPS system in one call (three-stage pipeline).

    Runs GeometryBuilder -> UnitTyper -> BoxPacker, producing:

    - ``<out>/<name>/geometry/geometry.json``
    - ``<out>/<name>/build/<ff>/*.lt`` + ``units.json``
    - ``<out>/<name>/moltemplate/system.lt`` -> moltemplate -> ``system.data`` etc.

    Args:
        system: System object providing get_folder_path().
        name: Project name; outputs go to ``<out>/<name>/``.
        models: List of Polymer and/or Molecule objects.
        force_field: Force field name (see FORCE_FIELD_REGISTRY).
            Defaults to "oplsaa".
        strategy: Box packing strategy ("mc_random", "grid", ...).
            Defaults to "mc_random".
        geometry_config: Optional GeometryConfig for stage 1 (MC chain
            growth parameters). Defaults to GeometryConfig().
        box_size: Optional explicit box edge length (Angstrom) for stage 3.
        mc_max_attempts: Maximum placement attempts for MC placement.
            Defaults to 10000.
        monomer_density: Target monomer density (monomers/A^3) for MC box
            sizing. Defaults to 0.085.
        rng_seed: Optional RNG seed for reproducible placement.
        run_moltemplate: Run moltemplate at the end of stage 3.
            Defaults to True.

    Returns:
        PlacementResult from the packing stage.

    Raises:
        WorkflowError: If any critical step fails or required files are missing

    Example:
        >>> from AutoPoly import System, Polymer, generate
        >>> system = System(out="pmma_out")
        >>> polymer = Polymer(chain_num=4, sequence=sequence,
        ...                   topology="linear", tacticity="atactic")
        >>> generate(system, "pmma", [polymer], force_field="oplsaa")
    """
    try:
        logger.info("Starting LAMMPS data file generation using Moltemplate")

        # Stage 1: force-field-agnostic geometry
        geometry_result = GeometryBuilder(system, name, config=geometry_config).build(models)

        # Stage 2: force-field typing on the stored geometry
        units = UnitTyper(geometry_result.dir, force_field).type()

        # Stage 3: box packing + moltemplate
        result = BoxPacker(
            system,
            name,
            strategy=strategy,
            box_size=box_size,
            mc_max_attempts=mc_max_attempts,
            monomer_density=monomer_density,
            rng_seed=rng_seed,
            run_moltemplate=run_moltemplate,
        ).pack(units)

        logger.info("Successfully completed polymer generation")
        return result

    except Exception as e:
        raise WorkflowError(
            f"Error in generate: {str(e)}"
        ) from e
