#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Box Packer — Stage 3 of the AutoPoly three-stage pipeline.

Packs typed units into a simulation box and runs moltemplate:

    build/<ff>/ (units.json + typed .lt) ──▶ moltemplate/ ──▶ system.data ...

The BoxPacker:
1. Accepts a UnitLibrary or a build-directory path (loads + validates units.json)
2. Copies the typed unit/monomer .lt files into <name>/moltemplate/ (keeps the
   self-contained layout; build/ remains the reusable source)
3. Generates the force-field files from the manifest's monomer files
   (no model objects needed — subsetting reads file paths)
4. Runs the selected packing strategy -> system.lt
5. Runs moltemplate (-nocheck), validates output, post-processes
   (get_rid_of_lj_cut_coul_long, mv_files)

Created on 2026-07-30
@author: zwu
"""
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Union

import numpy as np

from ..core.system import logger
from ..core.conf import FORCE_FIELD_REGISTRY
from ..core.exceptions import WorkflowError
from ..core.file_management import get_rid_of_lj_cut_coul_long, mv_files
from ..forcefields.force_field import ForceFieldManager
from ..packing import (
    BoxSpec,
    PackingContext,
    PlacementResult,
    get_strategy,
)
from .units import UnitLibrary

MOLTEMPLATE_DIRNAME = "moltemplate"


class BoxPacker:
    """
    Stage 3: pack typed units into a box and produce LAMMPS input files.

    Example:
        >>> packer = BoxPacker(system, "peo", strategy="mc_random")
        >>> packer.pack("peo_run/peo/build/gaff2")
    """

    def __init__(
        self,
        system: object,
        name: str,
        strategy: str = "mc_random",
        box_size: Optional[float] = None,
        mc_max_attempts: int = 10000,
        monomer_density: float = 0.085,
        rng_seed: Optional[int] = None,
        run_moltemplate: bool = True,
    ) -> None:
        """
        Args:
            system: System object providing get_folder_path().
            name: Project name; output goes to <out>/<name>/moltemplate/.
            strategy: Registered packing strategy name ("mc_random", "grid",
                      or a user-registered strategy).
            box_size: Explicit cubic box side (Angstrom); None = auto-size.
            mc_max_attempts: Max random placement attempts per instance.
            monomer_density: Target monomer density (monomers/A^3) for box sizing.
            rng_seed: Optional seed for reproducible stochastic strategies.
            run_moltemplate: Run moltemplate + post-processing (False = stop
                             after writing system.lt; useful for testing).

        Raises:
            ValidationError: If the strategy name is not registered.
        """
        self.system = system
        self.name = name
        self.project_dir = Path(system.get_folder_path()) / name
        self.moltemplate_dir = self.project_dir / MOLTEMPLATE_DIRNAME
        # Validate the strategy name eagerly (raises for unknown names)
        self.strategy = get_strategy(strategy)
        self.strategy_name = strategy
        self.box_size = box_size
        self.mc_max_attempts = mc_max_attempts
        self.monomer_density = monomer_density
        self.rng_seed = rng_seed
        self.run_moltemplate = run_moltemplate

        self.path_master = str(Path(__file__).parent.parent.resolve() / "extern/")
        self.path_moltemplatesrc = str(
            Path(self.path_master) / "moltemplate" / "scripts/"
        )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def pack(self, build: Union[str, Path, UnitLibrary]) -> PlacementResult:
        """
        Pack the units from a stage-2 build into the simulation box.

        Args:
            build: UnitLibrary, or a build directory / units.json path.

        Returns:
            PlacementResult from the strategy (instances + box bounds).

        Raises:
            WorkflowError: moltemplate failures, missing system.lt, missing
                           output files (same messages as the legacy workflow).
        """
        units = self._load_units(build)
        build_dir = self._resolve_build_dir(build)

        # 1. Fresh moltemplate dir + copy typed .lt files from build/<ff>/
        self._prepare_moltemplate_dir(units, build_dir)

        # 2. Force-field files, driven by the manifest's monomer files
        self._generate_force_field_files(units)

        # 3. Strategy -> placements -> system.lt
        result = self._run_strategy(units)
        self._write_system_lt(units, result)

        # 4. moltemplate + validation + post-processing
        if self.run_moltemplate:
            self._run_moltemplate_and_validate(units.force_field)
            logger.info("Processing output files")
            get_rid_of_lj_cut_coul_long(self.moltemplate_dir)
            mv_files(self.moltemplate_dir)

        logger.info("Successfully completed polymer generation")
        return result

    # ------------------------------------------------------------------
    # Build dir / manifest loading
    # ------------------------------------------------------------------
    @staticmethod
    def _load_units(build: Union[str, Path, UnitLibrary]) -> UnitLibrary:
        if isinstance(build, UnitLibrary):
            return build
        return UnitLibrary.load(build)

    @staticmethod
    def _resolve_build_dir(build: Union[str, Path, UnitLibrary]) -> Optional[Path]:
        """Directory holding the typed .lt files (None if unknown)."""
        if isinstance(build, UnitLibrary):
            return Path(build.source_dir) if build.source_dir else None
        build = Path(build)
        return build if build.is_dir() else build.parent

    # ------------------------------------------------------------------
    # Step 1: moltemplate dir + file copying
    # ------------------------------------------------------------------
    def _prepare_moltemplate_dir(
        self, units: UnitLibrary, build_dir: Optional[Path]
    ) -> None:
        """
        Create the moltemplate dir and copy typed .lt files into it.

        Any previous moltemplate/ contents are removed first: everything in
        it is regenerated from build/<ff>/, and stale files (moved output of
        an earlier pack) would otherwise collide with the new run.
        """
        if self.moltemplate_dir.exists():
            shutil.rmtree(self.moltemplate_dir)
        self.moltemplate_dir.mkdir(parents=True, exist_ok=True)

        if build_dir is None:
            # In-memory library: .lt files must already be in place
            return

        units.validate(build_dir)
        for rel in units.monomer_files:
            shutil.copy2(build_dir / rel, self.moltemplate_dir / rel)
        for unit in units.units:
            src = build_dir / unit.lt_file
            if src.is_file():
                shutil.copy2(src, self.moltemplate_dir / unit.lt_file)

    # ------------------------------------------------------------------
    # Step 2: force-field files
    # ------------------------------------------------------------------
    def _generate_force_field_files(self, units: UnitLibrary) -> None:
        """Generate force field files from the manifest's monomer files."""
        force_field = units.force_field
        logger.info(f"Generating {force_field}.lt")

        monomer_paths = [
            str(self.moltemplate_dir / rel) for rel in units.monomer_files
        ]
        ff_manager = ForceFieldManager(
            path_cwd=str(self.moltemplate_dir),
            path_master=self.path_master,
            path_moltemplatesrc=self.path_moltemplatesrc,
            force_field=force_field,
            ff_modify_dihedral=np.array(
                [0.6446926386, -0.2143420172, 0.1782194073, 0.0]
            ),
            monomer_files=monomer_paths,
        )
        ff_manager.make_force_field_lt()

        # Modify alkyl dihedral coefficients if needed (oplsaa only; LOPLS has
        # optimized alkyl dihedrals built in, other FFs use their own)
        if force_field == "oplsaa":
            ff_manager.FFmodify_alkyl_dihedral_oplsaa()

    # ------------------------------------------------------------------
    # Step 3: strategy + system.lt
    # ------------------------------------------------------------------
    def _run_strategy(self, units: UnitLibrary) -> PlacementResult:
        """Run the configured packing strategy."""
        offset = float(units.build_config.get("offset", 4.0))
        ctx = PackingContext(
            units=units,
            box=BoxSpec(requested_box_size=self.box_size),
            mc_max_attempts=self.mc_max_attempts,
            rng_seed=self.rng_seed,
            offset=offset,
            monomer_density=self.monomer_density,
        )
        logger.info(f"Placing units with strategy '{self.strategy_name}'")
        return self.strategy.place(ctx)

    def _write_system_lt(
        self, units: UnitLibrary, result: PlacementResult
    ) -> Path:
        """Write system.lt: imports + placements + Data Boundary."""
        ff_lt = FORCE_FIELD_REGISTRY[units.force_field]["lt_file"]
        output = self.moltemplate_dir / "system.lt"

        with open(output, "w") as f:
            f.write(f'import "{ff_lt}"\n\n')

            # Unit imports (polymer chain files import their own monomers)
            seen = set()
            for unit in units.units:
                if unit.lt_file in seen:
                    continue
                seen.add(unit.lt_file)
                f.write(f'import "{unit.lt_file}"\n')
            f.write("\n")

            for record in result.records:
                f.write(record.lt_command + "\n")
            f.write("\n")

            (xmin, xmax), (ymin, ymax), (zmin, zmax) = result.box_bounds
            f.write('write_once("Data Boundary") {\n')
            f.write(f"   {xmin:.4f}  {xmax:.4f}  xlo xhi\n")
            f.write(f"   {ymin:.4f}  {ymax:.4f}  ylo yhi\n")
            f.write(f"   {zmin:.4f}  {zmax:.4f}  zlo zhi\n")
            f.write("}\n")

        return output

    # ------------------------------------------------------------------
    # Step 4: moltemplate + validation + post-processing
    # ------------------------------------------------------------------
    def _run_moltemplate_and_validate(self, force_field: str) -> None:
        """Run moltemplate and validate the output."""
        logger.info("Running moltemplate")
        self.invoke_moltemplate()
        self._validate_moltemplate_output()
        self._check_required_files()
        self._log_gaff_charges_warning(force_field)

    def invoke_moltemplate(self) -> None:
        """
        Invoke moltemplate on system.lt (with -nocheck).

        Raises:
            WorkflowError: If system.lt is missing or moltemplate fails.
        """
        try:
            system_lt = self.moltemplate_dir / "system.lt"
            if not system_lt.exists():
                raise WorkflowError(
                    f"system.lt not found in {self.moltemplate_dir}"
                )

            moltemplate_sh = Path(self.path_moltemplatesrc) / "moltemplate.sh"
            process = subprocess.run(
                ["bash", str(moltemplate_sh), "-nocheck", "./system.lt"],
                cwd=self.moltemplate_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            if process.returncode != 0:
                logger.error("Moltemplate execution failed with the following error:")
                logger.error(process.stderr)
                stdout_lines = process.stdout.splitlines()
                if stdout_lines:
                    logger.error("Last output lines:")
                    for line in stdout_lines[-5:]:
                        logger.error(line)
                raise WorkflowError("Moltemplate execution failed")

        except Exception as e:
            raise WorkflowError(f"Error running Moltemplate: {str(e)}") from e

    def _validate_moltemplate_output(self) -> None:
        """Validate that moltemplate created the required output files."""
        output_ttree_dir = self.moltemplate_dir / "output_ttree"
        if not output_ttree_dir.exists():
            raise WorkflowError(
                f"Moltemplate failed to create output_ttree directory at "
                f"{output_ttree_dir}"
            )

        required_data_files = ["Data Atoms", "Data Bond List"]
        missing_data_files = []

        for data_file in required_data_files:
            data_file_path = output_ttree_dir / data_file
            if not data_file_path.exists():
                missing_data_files.append(data_file)
            elif data_file_path.stat().st_size == 0:
                logger.error(f"Data file exists but is empty: {data_file}")
                missing_data_files.append(f"{data_file} (empty)")

        if missing_data_files:
            logger.error(
                f"Moltemplate failed to generate required data files: "
                f"{missing_data_files}"
            )
            raise WorkflowError(
                f"Critical data files missing in {output_ttree_dir}. "
                "Ensure write(\"Data Atoms\") sections are present in monomer files"
            )

    def _check_required_files(self) -> None:
        """Check that required LAMMPS input files were generated."""
        required_files = ["system.in.settings", "system.data"]
        missing_files = [
            f for f in required_files
            if not (self.moltemplate_dir / f).exists()
        ]

        if missing_files:
            logger.error(
                f"Moltemplate failed to generate required files: "
                f"{', '.join(missing_files)}"
            )
            raise WorkflowError(
                "Required LAMMPS input files not generated. "
                "Check that polymer .lt files and system.lt are properly formatted"
            )

    def _log_gaff_charges_warning(self, force_field: str) -> None:
        """Log the legacy warning about missing charges for GAFF."""
        if force_field == "gaff":
            charges_file = self.moltemplate_dir / "system.in.charges"
            if not charges_file.exists():
                logger.warning("Note: system.in.charges not generated for GAFF")
                logger.warning(
                    "GAFF requires manual charge calculation using AM1-BCC or RESP"
                )
                logger.warning("All atomic charges are currently set to 0.00")
