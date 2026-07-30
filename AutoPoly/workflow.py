#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Workflow Orchestration Module for AutoPoly Package

This module provides the WorkflowManager class, a thin orchestrator that
composes the three pipeline stages:

    GeometryBuilder  (stage 1: FF-agnostic coordinates -> geometry/)
    UnitTyper        (stage 2: force-field typing      -> build/<ff>/)
    BoxPacker        (stage 3: box packing + moltemplate -> moltemplate/)

All stage logic lives in the stage modules (geometry.py, typing.py,
packer.py, packing/). This class only wires them together from the
Polymerization facade's configuration.

Created on 2026-01-06
@author: zwu
"""
from .system import logger
from .exceptions import WorkflowError
from .geometry import GeometryBuilder, GeometryConfig
from .typing import UnitTyper
from .packer import BoxPacker


class WorkflowManager:
    """
    Thin orchestrator composing the three pipeline stages.

    Attributes:
        poly (Polymerization): Reference to the main Polymerization instance
    """

    def __init__(self, polymerization_instance):
        """
        Initialize the WorkflowManager with a reference to Polymerization.

        Args:
            polymerization_instance (Polymerization): The main Polymerization object
        """
        self.poly = polymerization_instance

    def make_lmp_data_file_by_moltemplate(self) -> None:
        """
        Generate the LAMMPS data file via the three-stage pipeline.

        Runs GeometryBuilder -> UnitTyper -> BoxPacker, producing:
        - <out>/<name>/geometry/geometry.json
        - <out>/<name>/build/<ff>/*.lt + units.json
        - <out>/<name>/moltemplate/system.lt -> moltemplate -> system.data etc.

        Raises:
            WorkflowError: If any critical step fails or required files are missing
        """
        try:
            poly = self.poly
            logger.info("Starting LAMMPS data file generation using Moltemplate")

            # Stage 1: force-field-agnostic geometry
            geometry_config = GeometryConfig(
                use_mc_chain_growth=poly.use_mc_chain_growth,
                mc_max_attempts=poly.mc_max_attempts,
                mc_bond_angle_min=poly.mc_bond_angle_min,
                mc_bond_angle_max=poly.mc_bond_angle_max,
                mc_intrachain_exclude_neighbors=poly.mc_intrachain_exclude_neighbors,
                offset=poly.offset,
                offset_spacing=poly.offset_spacing,
                rotate=poly.rotate,
            )
            geometry = GeometryBuilder(poly.system, poly.name, geometry_config)
            geometry_result = geometry.build(poly.model)

            # Stage 2: force-field typing on the stored geometry
            units = UnitTyper(geometry_result.dir, poly.force_field).type()

            # Stage 3: box packing + moltemplate
            BoxPacker(
                poly.system,
                poly.name,
                strategy=getattr(poly, 'placement_method', 'mc_random'),
                mc_max_attempts=poly.mc_max_attempts,
                monomer_density=poly.mc_monomer_density,
            ).pack(units)

            logger.info("Successfully completed polymer generation")

        except Exception as e:
            raise WorkflowError(
                f"Error in make_lmp_data_file_by_moltemplate: {str(e)}"
            ) from e
