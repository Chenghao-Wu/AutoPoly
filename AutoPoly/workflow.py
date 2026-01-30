#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Workflow Orchestration Module for AutoPoly Package

This module provides the WorkflowManager class that handles the main workflow
orchestration for polymer structure generation using Moltemplate. It manages
the complete process from monomer template preparation to LAMMPS data file creation.

Key Features:
- Main polymerization workflow orchestration
- Moltemplate invocation and execution
- System.lt file generation
- Polymer .lt file creation for both linear and ring topologies
- Integration with force field management
- Monte Carlo placement methods (grid, mc_random)
- Self-avoiding random walk for chain growth

Created on 2026-01-06
@author: zwu
"""
import subprocess
from pathlib import Path
import numpy as np
from typing import List, Optional, Dict, Any, Tuple
from .system import logger
from .exceptions import WorkflowError
from .monomer_processing import read_lt_end_atoms, evaluate_offset


class WorkflowManager:
    """
    Manages the workflow orchestration for polymer structure generation.

    This class handles the main workflow methods extracted from the Polymerization
    class, providing a clean separation of concerns between workflow orchestration
    and other polymerization functionality.

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
        Generate the LAMMPS data file using Moltemplate.

        This is the main method that orchestrates the complete polymer generation
        process. It performs the following steps:

        1. Validates polymer models and monomer availability
        2. Copies monomer templates to working directory
        3. Generates polymer .lt files for each chain
        4. Creates force field parameter files (OPLS-AA)
        5. Generates system.lt file
        6. Runs Moltemplate to create LAMMPS data files
        7. Processes and organizes output files

        Raises:
            WorkflowError: If any critical step fails or required files are missing
        """
        try:
            logger.info("Starting LAMMPS data file generation using Moltemplate")

            # Process all models (polymers and molecules)
            poly_index = self._process_all_models()

            # Generate force field and system files
            self._generate_force_field_files()
            self._generate_system_file()

            # Run moltemplate and validate output
            self._run_moltemplate_and_validate()

            # Finalize output
            logger.info("Processing output files")
            self.poly.get_rid_of_lj_cut_coul_long()
            self.poly.mv_files()
            logger.info("Successfully completed polymer generation")

        except Exception as e:
            raise WorkflowError(
                f"Error in make_lmp_data_file_by_moltemplate: {str(e)}"
            ) from e

    def _process_all_models(self) -> int:
        """
        Process all models in the system.

        Returns:
            int: The final polymer index
        """
        poly_index = 0
        for modelii in self.poly.model:
            if hasattr(modelii, '_is_molecule') and modelii._is_molecule:
                self._generate_single_molecule(modelii)
            else:
                poly_index = self._process_polymer_model(modelii, poly_index)
        return poly_index

    def _process_polymer_model(self, model: object, poly_index: int) -> int:
        """
        Process a single polymer model using complement SMILES.

        Args:
            model: The polymer model to process
            poly_index: Current polymer index

        Returns:
            int: Updated polymer index
        """
        topology = getattr(model, 'topology', 'linear')
        dop = model.dop

        logger.info(f"Generating sequence variants for {len(model.sequence)} complement SMILES, DOP={dop}, topology={topology}")

        variant_mapping = self.poly.generate_sequence_variants_for_polymer(
            smiles_list=model.sequence,
            topology=topology,
            base_name_prefix="monomer"
        )

        logger.info(f"Generated {len(variant_mapping)} unique variant types: {list(variant_mapping.keys())}")

        # Build sequenceSet with appropriate .lt files for each position
        self._build_sequence_set(model, variant_mapping, topology, dop)
        logger.info(f"Built {len(model.sequenceSet)} chain(s) with sequence variants")

        # Generate polymer .lt files
        use_mc_chain_growth = getattr(self.poly, 'use_mc_chain_growth', True)

        for chain_idx in range(len(model.sequenceSet)):
            if model.dop > 1:
                logger.info(f"Creating poly_{poly_index+1}.lt")
                if use_mc_chain_growth:
                    self.make_poly_lt_mc(poly_index, model.sequenceSet[chain_idx], model)
                else:
                    self.make_poly_lt(poly_index, model.sequenceSet[chain_idx], model)
                poly_index += 1

        return poly_index

    def _build_sequence_set(
        self,
        model: object,
        variant_mapping: dict,
        topology: str,
        dop: int
    ) -> None:
        """
        Build the sequenceSet with appropriate .lt files for each position.

        Args:
            model: The polymer model
            variant_mapping: Mapping of variant types to filenames
            topology: Polymer topology
            dop: Degree of polymerization
        """
        for chain_idx, chain_smiles in enumerate(model.sequenceSet):
            monomer_names = []
            for pos_idx, smiles in enumerate(chain_smiles):
                variant_type = self._determine_variant_type(topology, pos_idx, dop)
                has_t1 = "_T1" in smiles
                mapping_key = f"{variant_type}_T1" if has_t1 else variant_type

                if mapping_key not in variant_mapping:
                    raise WorkflowError(
                        f"Variant key '{mapping_key}' not found in mapping. "
                        f"Available variants: {list(variant_mapping.keys())}"
                    )
                monomer_names.append(variant_mapping[mapping_key])

            model.sequenceSet[chain_idx] = monomer_names

    def _determine_variant_type(self, topology: str, pos_idx: int, dop: int) -> str:
        """
        Determine the variant type based on position and topology.

        Args:
            topology: Polymer topology (linear or ring)
            pos_idx: Position in the chain
            dop: Degree of polymerization

        Returns:
            str: Variant type (first, middle, last, or ring)
        """
        if topology == "ring":
            return "ring"
        if pos_idx == 0:
            return "first"
        if pos_idx == dop - 1:
            return "last"
        return "middle"

    def _generate_force_field_files(self) -> None:
        """Generate force field parameter files."""
        logger.info(f"Generating {self.poly.force_field}.lt")
        self.poly.ff_manager.model = self.poly.model
        self.poly.ff_manager.make_force_field_lt()

        # Modify alkyl dihedral coefficients if needed
        # Only for oplsaa - LOPLS already has optimized alkyl dihedrals built-in
        # GAFF/GAFF2/DREIDING/COMPASS use their own dihedral parameters
        if self.poly.force_field == "oplsaa":
            self.poly.ff_manager.FFmodify_alkyl_dihedral_oplsaa()

    def _generate_system_file(self) -> None:
        """Generate the system.lt file."""
        placement_method = getattr(self.poly, 'placement_method', 'mc_random')
        logger.info(f"Creating system.lt (placement_method={placement_method})")

        if placement_method == "mc_random":
            self.make_system_lt_mc(placement_method)
        else:
            self.make_system_lt()

    def _run_moltemplate_and_validate(self) -> None:
        """Run moltemplate and validate the output."""
        logger.info("Running moltemplate")
        self.invoke_moltemplate()
        self._validate_moltemplate_output()
        self._check_required_files()
        self._log_gaff_charges_warning()

    def _validate_moltemplate_output(self) -> None:
        """Validate that moltemplate created the required output files."""
        output_ttree_dir = Path(self.poly.path_cwd) / "output_ttree"
        if not output_ttree_dir.exists():
            raise WorkflowError(
                f"Moltemplate failed to create output_ttree directory at {output_ttree_dir}"
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
            logger.error(f"Moltemplate failed to generate required data files: {missing_data_files}")
            raise WorkflowError(
                f"Critical data files missing in {output_ttree_dir}. "
                "Ensure write(\"Data Atoms\") sections are present in monomer files"
            )

    def _check_required_files(self) -> None:
        """Check that required LAMMPS input files were generated."""
        required_files = ['system.in.settings', 'system.data']
        missing_files = [
            f for f in required_files
            if not (Path(self.poly.path_cwd) / f).exists()
        ]

        if missing_files:
            logger.error(f"Moltemplate failed to generate required files: {', '.join(missing_files)}")
            raise WorkflowError(
                "Required LAMMPS input files not generated. "
                "Check that polymer .lt files and system.lt are properly formatted"
            )

    def _log_gaff_charges_warning(self) -> None:
        """Log warning about missing charges file for GAFF."""
        if self.poly.force_field == "gaff":
            charges_file = Path(self.poly.path_cwd) / "system.in.charges"
            if not charges_file.exists():
                logger.warning("Note: system.in.charges not generated for GAFF")
                logger.warning("GAFF requires manual charge calculation using AM1-BCC or RESP")
                logger.warning("All atomic charges are currently set to 0.00")

    def invoke_moltemplate(self) -> None:
        """
        Invokes Moltemplate to generate the LAMMPS data file.

        This method runs the moltemplate.sh script with the system.lt file
        as input, capturing and handling any errors that occur during execution.

        Raises:
            SystemExit: If moltemplate execution fails or system.lt is not found
        """
        try:
            # First check if system.lt exists
            system_lt = Path(self.poly.path_cwd) / "system.lt"
            if not system_lt.exists():
                raise WorkflowError(f"system.lt not found in {self.poly.path_cwd}")

            # Run moltemplate with output capture for error cases
            # Use cwd parameter and explicit bash path for security
            # Use -nocheck flag to bypass post-processing validation that may fail on first run
            moltemplate_sh = Path(self.poly.path_moltemplatesrc) / "moltemplate.sh"
            process = subprocess.run(
                ["bash", str(moltemplate_sh), "-nocheck", "./system.lt"],
                cwd=self.poly.path_cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            if process.returncode != 0:
                logger.error("Moltemplate execution failed with the following error:")
                logger.error(process.stderr)
                # Also log the last few lines of stdout which might contain useful info
                stdout_lines = process.stdout.splitlines()
                if stdout_lines:
                    logger.error("Last output lines:")
                    for line in stdout_lines[-5:]:
                        logger.error(line)
                raise WorkflowError("Moltemplate execution failed")

        except Exception as e:
            raise WorkflowError(f"Error running Moltemplate: {str(e)}") from e

    def _generate_single_molecule(self, molecule) -> None:
        """
        Generate .lt file for a single molecule.

        This method generates the .lt template file for a non-polymer molecule
        (e.g., water, benzene, ethanol) using the MonomerGenerator infrastructure.

        Args:
            molecule: Molecule object containing Smiles, Count, and molecule_name

        Raises:
            SystemExit: If molecule generation fails
        """
        try:
            # Generate molecule .lt file using Polymerization method
            filename, counter = self.poly.generate_molecule_from_smiles(
                smiles=molecule.Smiles,
                molecule_name=molecule.molecule_name
            )

            # Update the molecule's sequenceSet to reference the generated .lt file
            # This is needed for system.lt generation
            for i in range(len(molecule.sequenceSet)):
                molecule.sequenceSet[i] = [filename]
                molecule.sequenceName[i] = [filename]

            logger.info(f"Generated molecule .lt file: {filename}")

        except Exception as e:
            raise WorkflowError(
                f"Error generating molecule {molecule.molecule_name}: {str(e)}"
            ) from e

    def make_system_lt(self) -> None:
        """
        Creates the system.lt file for the polymerization.

        This method generates the main system.lt file that imports all necessary
        polymer files, force field parameters, and defines the simulation box
        with proper spacing and positioning for all polymer chains and molecules.
        """
        output = Path(self.poly.path_cwd) / "system.lt"

        # Determine force field import based on force_field type
        if self.poly.force_field == "gaff":
            ff_import = 'import "gaff.lt"\n\n'
        elif self.poly.force_field == "gaff2":
            ff_import = 'import "gaff2.lt"\n\n'
        elif self.poly.force_field == "lopls":
            ff_import = 'import "loplsaa.lt"\n\n'
        elif self.poly.force_field == "dreiding":
            ff_import = 'import "dreiding.lt"\n\n'
        elif self.poly.force_field == "compass":
            ff_import = 'import "compass_published.lt"\n\n'
        else:
            ff_import = 'import "oplsaa.lt"\n\n'

        with open(output, "w") as write_f:
            # Write force field import at the top
            write_f.write(ff_import)

            polyindex = 0
            for modelii in self.poly.model:
                # Check if this is a Molecule or Polymer
                is_molecule = hasattr(modelii, '_is_molecule') and modelii._is_molecule

                if is_molecule:
                    # For molecules, import the molecule .lt file (already generated)
                    write_f.write(f"import \"{modelii.sequenceSet[0][0]}\"\n")
                    write_f.write("\n")
                elif modelii.dop > 1:
                    # For polymers with DOP>1, import poly_N.lt files
                    n_poly = len(modelii.sequenceSet)
                    for indexi in range(n_poly):
                        write_f.write(f"import \"poly_{polyindex+1}.lt\"\n")
                        polyindex += 1
                    write_f.write("\n")
                else:
                    # For polymers with DOP=1 (single monomers)
                    if len(modelii.merSet) > 1:
                        raise WorkflowError(
                            f"sequenceLen = {modelii.dop}, merSet should only have one mer type!"
                        )
                    # Import constituent monomer.lt's
                    unique_Sequence = [i[0] for i in modelii.sequenceSet]
                    print(unique_Sequence)
                    for sequenceii in range(len(unique_Sequence)):
                        write_f.write("import \""+unique_Sequence[sequenceii]+"\"\n")
                    write_f.write("\n")

            polyindex = 0
            index = 0

            # Calculate spacing based on polymer type and size
            for modelii in self.poly.model:
                # Check if this is a Molecule or Polymer
                is_molecule = hasattr(modelii, '_is_molecule') and modelii._is_molecule

                if is_molecule:
                    # For molecules, use a fixed spacing based on molecular size
                    # Use a default spacing of 5.0 Angstroms for molecules
                    spacing = 5.0
                    n_instances = modelii.Count
                else:
                    # For polymers, calculate spacing based on polymer type
                    n_poly = len(modelii.sequenceSet)
                    n_instances = n_poly
                    is_ring = hasattr(modelii, 'topology') and modelii.topology == "ring"

                    if is_ring:
                        # For ring polymers, calculate radius and use it for spacing
                        radius = self.poly.offset * len(modelii.sequenceSet[0]) / (2 * np.pi)
                        spacing = radius * 2.5  # Use 2.5x the ring radius for good separation
                    else:
                        spacing = self.poly.offset * (modelii.dop + 2)

                # Calculate grid arrangement
                grid_size = int(np.ceil(np.sqrt(n_instances)))  # Arrange in a square grid

                for moleii in range(n_instances):
                    # Calculate grid position
                    grid_x = moleii % grid_size
                    grid_y = moleii // grid_size

                    # Calculate actual position with spacing
                    pos_x = grid_x * spacing
                    pos_y = grid_y * spacing
                    pos_z = 0.0  # Keep all polymers/molecules in the same plane initially

                    if is_molecule:
                        # For molecules, instantiate the molecule type
                        write_f.write(f"molecule_{index+1} = new {modelii.molecule_name}")
                        write_f.write(f".move({pos_x:.4f},{pos_y:.4f},{pos_z:.4f})\n")
                    elif modelii.dop > 1:
                        # For polymers, instantiate the polymer type
                        write_f.write(f"polymer_{index+1} = new poly_{polyindex+1}")
                        write_f.write(f".move({pos_x:.4f},{pos_y:.4f},{pos_z:.4f})\n")
                        polyindex += 1
                    else:
                        # For single monomers (DOP=1)
                        write_f.write(f"molecule_{index+1} = new {modelii.merSet[0]}")
                        write_f.write(f".move({pos_x:.4f},{pos_y:.4f},{pos_z:.4f})\n")

                    index += 1

                write_f.write("\n")

            # Adjust box size based on total system size
            max_coord = max([
                grid_size * spacing for modelii in self.poly.model
                if len(modelii.sequenceSet) > 0
            ])
            box_size = max_coord * 1.2  # Add 20% padding

            # Write box boundaries
            write_f.write("write_once(\"Data Boundary\") {\n")
            write_f.write(f"   -{box_size/2:.4f}  {box_size/2:.4f}  xlo xhi\n")
            write_f.write(f"   -{box_size/2:.4f}  {box_size/2:.4f}  ylo yhi\n")
            write_f.write(f"   -{box_size/2:.4f}  {box_size/2:.4f}  zlo zhi\n")
            write_f.write("}\n")

    def make_poly_lt(self, poly_index: int, monomer_set: list, model: object) -> None:
        """
        Creates a poly.lt file for the specified polymer.

        This method generates the polymer .lt file that defines how monomers
        are connected to form the polymer chain. It supports both linear and
        ring polymer topologies.

        Args:
            poly_index (int): The index of the polymer.
            monomer_set (list): The list of monomers in the polymer.
            model (object): The polymer model object containing topology information.
        """
        output = Path(self.poly.path_cwd) / f"poly_{poly_index+1}.lt"

        # Determine force field import and inheritance based on force_field type
        if self.poly.force_field == "gaff":
            ff_import = 'import "gaff.lt"\n'
            ff_inherits = "GAFF"
        elif self.poly.force_field == "gaff2":
            ff_import = 'import "gaff2.lt"\n'
            ff_inherits = "GAFF2"
        elif self.poly.force_field == "lopls":
            ff_import = 'import "loplsaa.lt"\n'
            ff_inherits = "OPLSAA"  # LOPLS extends OPLSAA
        elif self.poly.force_field == "dreiding":
            ff_import = 'import "dreiding.lt"\n'
            ff_inherits = "DREIDING"
        elif self.poly.force_field == "compass":
            ff_import = 'import "compass_published.lt"\n'
            ff_inherits = "COMPASS"
        else:
            ff_import = 'import "oplsaa.lt"\n'
            ff_inherits = "OPLSAA"

        with open(output, "w") as write_f:
            write_f.write(ff_import)

            # Import unique monomers - ensure proper .lt extension
            unique_monomers = list(dict.fromkeys(monomer_set))
            for monomer in unique_monomers:
                # Remove .lt if it exists, then add it back
                base_name = monomer[:-3] if monomer.endswith('.lt') else monomer
                write_f.write(f"import \"{base_name}.lt\"\n")

            write_f.write("\n")

            # Define combined molecule (ex.polymer)
            write_f.write(f"poly_{poly_index+1} inherits {ff_inherits} {{\n\n")
            write_f.write("    create_var {$mol}\n\n")

            # Check if this is a ring polymer
            is_ring = hasattr(model, 'topology') and model.topology == "ring"

            if is_ring:
                # Calculate radius based on number of monomers and offset
                n_monomers = len(monomer_set)
                radius = self.poly.offset * n_monomers / (2 * np.pi)

                for i in range(n_monomers):
                    # Calculate position on the ring
                    angle = 2 * np.pi * i / n_monomers
                    x = radius * np.cos(angle)
                    y = radius * np.sin(angle)

                    # Calculate rotation to point each monomer towards the center
                    rotation_angle = (angle * 180 / np.pi) + 90  # Convert to degrees and add offset

                    # Get base name without .lt extension for the instance
                    monomer_name = monomer_set[i][:-3] if monomer_set[i].endswith('.lt') else monomer_set[i]

                    # Create monomer with position and rotation
                    write_f.write(f"    monomer[{i}] = new {monomer_name}")
                    write_f.write(f".rot({rotation_angle},0,0,1)")  # Rotate around z-axis
                    write_f.write(f".move({x:.4f},{y:.4f},0)\n")

                # Add bonds between monomers including the ring closure
                write_f.write("\n    write('Data Bond List') {\n")
                for i in range(n_monomers):
                    next_i = (i + 1) % n_monomers  # Wrap around to 0 for last monomer
                    monomer_name_1 = monomer_set[i][:-3] if monomer_set[i].endswith('.lt') else monomer_set[i]
                    monomer_name_2 = monomer_set[next_i][:-3] if monomer_set[next_i].endswith('.lt') else monomer_set[next_i]
                    # Determine monomer bank path based on force field
                    monomer_bank = Path(self.poly.path_cwd)
                    merltfile_path_1 = monomer_bank / f"{monomer_name_1}.lt"
                    _, second_atom = read_lt_end_atoms(merltfile_path_1)
                    merltfile_path_2 = monomer_bank / f"{monomer_name_2}.lt"
                    first_atom, _ = read_lt_end_atoms(merltfile_path_2)
                    write_f.write(f"      $bond:b{i+1}  $atom:monomer[{i}]/{second_atom}  $atom:monomer[{next_i}]/{first_atom}\n")
                write_f.write("    }\n")

            else:
                # Original linear polymer code
                offset_cum = 0
                for indexii in range(len(monomer_set)):

                    monomer_name = monomer_set[indexii][:-3] if monomer_set[indexii].endswith('.lt') else monomer_set[indexii]

                    write_f.write(f"    monomer[{indexii}] = new {monomer_name}")
                    if indexii > 0:
                        write_f.write(f".rot({self.poly.rotate*(indexii%2)},1,0,0).move({offset_cum:.4f},0,0)")
                    write_f.write("\n")

                    # Call evaluate_offset through self.poly
                    self.poly.evaluate_offset(f"{monomer_name}.lt")
                    offset_cum += self.poly.offset

                # Add bonds for linear polymer
                write_f.write("\n    write('Data Bond List') {\n")
                for indexii in range(len(monomer_set)-1):

                    monomer_name_1 = monomer_set[indexii][:-3] if monomer_set[indexii].endswith('.lt') else monomer_set[indexii]
                    monomer_name_2 = monomer_set[indexii+1][:-3] if monomer_set[indexii+1].endswith('.lt') else monomer_set[indexii+1]
                    # Determine monomer bank path based on force field
                    monomer_bank = Path(self.poly.path_cwd)
                    merltfile_path_1 = monomer_bank / f"{monomer_name_1}.lt"
                    _, second_atom = read_lt_end_atoms(merltfile_path_1)
                    merltfile_path_2 = monomer_bank / f"{monomer_name_2}.lt"
                    first_atom, _ = read_lt_end_atoms(merltfile_path_2)
                    write_f.write(f"      $bond:b{indexii+1}  $atom:monomer[{indexii}]/{second_atom}  $atom:monomer[{indexii+1}]/{first_atom}\n")
                write_f.write("    }\n")

            write_f.write(f"\n}} # poly_{poly_index+1}\n")

    def make_system_lt_mc(self, placement_method: str = "mc_random") -> None:
        """
        Creates the system.lt file using Monte Carlo random placement.

        This method generates the main system.lt file with random positioning
        and orientation of polymers and molecules using collision detection.

        Args:
            placement_method: Placement method ("mc_random" for Monte Carlo)
        """
        from .mc import (
            CollisionDetector,
            MolecularPlacementMC,
            calculate_box_size,
        )

        output = Path(self.poly.path_cwd) / "system.lt"

        # Determine force field import
        ff_import = self._get_ff_import()

        # Calculate box size: use the larger of density-based and chain-length-based estimates
        total_monomers = self._count_total_monomers()
        monomer_density = getattr(self.poly, 'mc_monomer_density', 0.05)
        density_box = calculate_box_size(total_monomers, monomer_density)

        # Account for polymer extended length to ensure chains fit inside the box
        max_dop = max((m.dop for m in self.poly.model if hasattr(m, 'dop')), default=1)
        offset = getattr(self.poly, 'offset', 4.0)
        chain_length_box = max_dop**0.6 * offset * 3.0  # SAW end-to-end × safety factor

        box_size = max(density_box, chain_length_box)
        half_box = box_size / 2

        box_bounds = ((-half_box, half_box), (-half_box, half_box), (-half_box, half_box))

        # Initialize collision detector and placer
        cell_size = max(5.0, box_size / 20)
        collision_detector = CollisionDetector(box_bounds, cell_size)
        max_attempts = getattr(self.poly, 'mc_max_attempts', 10000)
        placer = MolecularPlacementMC(box_bounds, collision_detector, max_attempts)

        with open(output, "w") as write_f:
            write_f.write(ff_import)

            # Write imports
            polyindex = 0
            for modelii in self.poly.model:
                is_molecule = hasattr(modelii, '_is_molecule') and modelii._is_molecule

                if is_molecule:
                    write_f.write(f"import \"{modelii.sequenceSet[0][0]}\"\n")
                    write_f.write("\n")
                elif modelii.dop > 1:
                    n_poly = len(modelii.sequenceSet)
                    for indexi in range(n_poly):
                        write_f.write(f"import \"poly_{polyindex+1}.lt\"\n")
                        polyindex += 1
                    write_f.write("\n")
                else:
                    if len(modelii.merSet) > 1:
                        raise WorkflowError(
                            f"sequenceLen = {modelii.dop}, merSet should only have one mer type!"
                        )
                    unique_Sequence = [i[0] for i in modelii.sequenceSet]
                    for sequenceii in range(len(unique_Sequence)):
                        write_f.write("import \""+unique_Sequence[sequenceii]+"\"\n")
                    write_f.write("\n")

            # Place entities using MC
            polyindex = 0
            index = 0

            for modelii in self.poly.model:
                is_molecule = hasattr(modelii, '_is_molecule') and modelii._is_molecule

                if is_molecule:
                    # Estimate molecule radius (simple default)
                    mol_radius = 3.0
                    n_instances = modelii.Count

                    for moleii in range(n_instances):
                        placement = placer.place_molecule(
                            molecule_name=modelii.molecule_name,
                            radius=mol_radius,
                            instance_name=f"molecule_{index+1}"
                        )
                        if placement is None:
                            logger.warning(f"Failed to place molecule {index+1}, using fallback position")
                            # Fallback to grid position
                            pos_x = (index % 10) * 10.0
                            pos_y = (index // 10) * 10.0
                            pos_z = 0.0
                            write_f.write(f"molecule_{index+1} = new {modelii.molecule_name}")
                            write_f.write(f".move({pos_x:.4f},{pos_y:.4f},{pos_z:.4f})\n")
                        else:
                            cmd = placer.generate_molecule_lt_commands([placement])[0]
                            write_f.write(cmd + "\n")
                        index += 1

                elif modelii.dop > 1:
                    n_poly = len(modelii.sequenceSet)
                    # Estimate polymer radius based on DOP using random walk statistics.
                    # For a freely-jointed chain, R_g ≈ offset * sqrt(dop / 6).
                    # Use 2 * R_g as collision radius for inter-chain overlap avoidance.
                    is_ring = hasattr(modelii, 'topology') and modelii.topology == "ring"
                    if is_ring:
                        poly_radius = self.poly.offset * np.sqrt(modelii.dop / 12) * 2.0 + 2.0
                    else:
                        poly_radius = self.poly.offset * np.sqrt(modelii.dop / 6) * 2.0 + 2.0

                    for chain_idx in range(n_poly):
                        placement = placer.place_polymer(
                            poly_name=f"poly_{polyindex+1}",
                            radius=poly_radius
                        )
                        if placement is None:
                            logger.warning(f"Failed to place polymer {polyindex+1}, using fallback")
                            pos_x = (polyindex % 10) * poly_radius * 2.5
                            pos_y = (polyindex // 10) * poly_radius * 2.5
                            pos_z = 0.0
                            write_f.write(f"polymer_{polyindex+1} = new poly_{polyindex+1}")
                            write_f.write(f".move({pos_x:.4f},{pos_y:.4f},{pos_z:.4f})\n")
                        else:
                            cmd = placer.generate_polymer_lt_commands([placement])[0]
                            write_f.write(cmd + "\n")
                        polyindex += 1
                        index += 1

                else:
                    # Single monomers (DOP=1)
                    for chain_idx in range(len(modelii.sequenceSet)):
                        placement = placer.place_molecule(
                            molecule_name=modelii.merSet[0],
                            radius=3.0,
                            instance_name=f"molecule_{index+1}"
                        )
                        if placement is None:
                            pos_x = (index % 10) * 10.0
                            pos_y = (index // 10) * 10.0
                            pos_z = 0.0
                            write_f.write(f"molecule_{index+1} = new {modelii.merSet[0]}")
                            write_f.write(f".move({pos_x:.4f},{pos_y:.4f},{pos_z:.4f})\n")
                        else:
                            cmd = placer.generate_molecule_lt_commands([placement])[0]
                            write_f.write(cmd + "\n")
                        index += 1

                write_f.write("\n")

            # Write box boundaries
            write_f.write("write_once(\"Data Boundary\") {\n")
            write_f.write(f"   -{half_box:.4f}  {half_box:.4f}  xlo xhi\n")
            write_f.write(f"   -{half_box:.4f}  {half_box:.4f}  ylo yhi\n")
            write_f.write(f"   -{half_box:.4f}  {half_box:.4f}  zlo zhi\n")
            write_f.write("}\n")

        stats = placer.get_placement_stats()
        logger.info(f"MC placement complete: {stats['polymers']} polymers, "
                   f"{stats['molecules']} molecules placed")

    def make_poly_lt_mc(
        self,
        poly_index: int,
        monomer_set: list,
        model: object,
        collision_detector=None
    ) -> None:
        """
        Creates a poly.lt file using Monte Carlo chain growth.

        This method generates the polymer .lt file using self-avoiding random
        walk for monomer placement, providing more realistic chain conformations.

        Args:
            poly_index: The index of the polymer
            monomer_set: The list of monomers in the polymer
            model: The polymer model object
            collision_detector: Optional shared CollisionDetector
        """
        from .mc import CollisionDetector, ChainGrowthMC, calculate_box_size

        output = Path(self.poly.path_cwd) / f"poly_{poly_index+1}.lt"

        # Get force field settings
        ff_import, ff_inherits = self._get_ff_import_and_inherits()

        # Compute box bounds (needed for collision detector init and retries)
        n_monomers = len(monomer_set)
        estimated_chain_length = n_monomers**0.6 * 4.0
        box_size = max(estimated_chain_length * 3.0,
                      calculate_box_size(n_monomers, monomer_density=0.05))
        half_box = box_size / 2
        box_bounds = ((-half_box, half_box), (-half_box, half_box), (-half_box, half_box))

        # Initialize collision detector if not provided
        if collision_detector is None:
            collision_detector = CollisionDetector(box_bounds, cell_size=5.0)

        max_attempts = getattr(self.poly, 'mc_max_attempts', 1000)
        bond_angle_min = getattr(self.poly, 'mc_bond_angle_min', 50.0)
        bond_angle_max = getattr(self.poly, 'mc_bond_angle_max', 90.0)
        exclude_neighbors = getattr(self.poly, 'mc_intrachain_exclude_neighbors', 2)
        chain_mc = ChainGrowthMC(
            collision_detector,
            max_attempts,
            bond_angle_min=bond_angle_min,
            bond_angle_max=bond_angle_max,
            intrachain_exclude_neighbors=exclude_neighbors
        )

        # Build list of .lt file paths
        monomer_bank = Path(self.poly.path_cwd)
        lt_files = []
        for monomer in monomer_set:
            monomer_name = monomer[:-3] if monomer.endswith('.lt') else monomer
            lt_files.append(str(monomer_bank / f"{monomer_name}.lt"))

        is_ring = hasattr(model, 'topology') and model.topology == "ring"

        with open(output, "w") as write_f:
            write_f.write(ff_import)

            # Import unique monomers
            unique_monomers = list(dict.fromkeys(monomer_set))
            for monomer in unique_monomers:
                base_name = monomer[:-3] if monomer.endswith('.lt') else monomer
                write_f.write(f"import \"{base_name}.lt\"\n")

            write_f.write("\n")
            write_f.write(f"poly_{poly_index+1} inherits {ff_inherits} {{\n\n")
            write_f.write("    create_var {$mol}\n\n")

            if is_ring:
                # Ring topology: use circular placement
                self._write_ring_polymer_mc(write_f, monomer_set, chain_mc, lt_files)
            else:
                # Linear topology: use chain growth MC with retries.
                # Each retry resets the collision detector so stale monomer
                # registrations from a failed attempt don't block the next one.
                max_chain_retries = 5
                placed = False
                for retry in range(max_chain_retries):
                    if retry > 0:
                        # Fresh collision detector for each retry
                        collision_detector = CollisionDetector(box_bounds, cell_size=5.0)
                        chain_mc = ChainGrowthMC(
                            collision_detector,
                            max_attempts,
                            bond_angle_min=bond_angle_min,
                            bond_angle_max=bond_angle_max,
                            intrachain_exclude_neighbors=exclude_neighbors
                        )
                    try:
                        placements = chain_mc.grow_chain(lt_files, chain_id=poly_index)
                        commands = chain_mc.generate_lt_commands(placements)
                        for cmd in commands:
                            write_f.write(cmd + "\n")
                        placed = True
                        break
                    except RuntimeError as e:
                        logger.warning(
                            f"MC chain growth attempt {retry+1}/{max_chain_retries} "
                            f"failed for poly_{poly_index+1}: {e}"
                        )
                if not placed:
                    logger.warning("All MC retries exhausted, falling back to deterministic placement")
                    self._write_linear_polymer_deterministic(write_f, monomer_set)

            # Write bonds
            write_f.write("\n    write('Data Bond List') {\n")
            if is_ring:
                n_monomers = len(monomer_set)
                for i in range(n_monomers):
                    next_i = (i + 1) % n_monomers
                    monomer_name_1 = monomer_set[i][:-3] if monomer_set[i].endswith('.lt') else monomer_set[i]
                    monomer_name_2 = monomer_set[next_i][:-3] if monomer_set[next_i].endswith('.lt') else monomer_set[next_i]
                    merltfile_path_1 = monomer_bank / f"{monomer_name_1}.lt"
                    _, second_atom = read_lt_end_atoms(merltfile_path_1)
                    merltfile_path_2 = monomer_bank / f"{monomer_name_2}.lt"
                    first_atom, _ = read_lt_end_atoms(merltfile_path_2)
                    write_f.write(f"      $bond:b{i+1}  $atom:monomer[{i}]/{second_atom}  $atom:monomer[{next_i}]/{first_atom}\n")
            else:
                for indexii in range(len(monomer_set)-1):
                    monomer_name_1 = monomer_set[indexii][:-3] if monomer_set[indexii].endswith('.lt') else monomer_set[indexii]
                    monomer_name_2 = monomer_set[indexii+1][:-3] if monomer_set[indexii+1].endswith('.lt') else monomer_set[indexii+1]
                    merltfile_path_1 = monomer_bank / f"{monomer_name_1}.lt"
                    _, second_atom = read_lt_end_atoms(merltfile_path_1)
                    merltfile_path_2 = monomer_bank / f"{monomer_name_2}.lt"
                    first_atom, _ = read_lt_end_atoms(merltfile_path_2)
                    write_f.write(f"      $bond:b{indexii+1}  $atom:monomer[{indexii}]/{second_atom}  $atom:monomer[{indexii+1}]/{first_atom}\n")
            write_f.write("    }\n")

            write_f.write(f"\n}} # poly_{poly_index+1}\n")

    def _write_ring_polymer_mc(
        self,
        write_f,
        monomer_set: list,
        chain_mc,
        lt_files: list
    ) -> None:
        """Write ring polymer using MC with circular constraints."""
        n_monomers = len(monomer_set)
        radius = self.poly.offset * n_monomers / (2 * np.pi)

        for i in range(n_monomers):
            angle = 2 * np.pi * i / n_monomers
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            rotation_angle = (angle * 180 / np.pi) + 90

            monomer_name = monomer_set[i][:-3] if monomer_set[i].endswith('.lt') else monomer_set[i]

            write_f.write(f"    monomer[{i}] = new {monomer_name}")
            write_f.write(f".rot({rotation_angle:.4f},0,0,1)")
            write_f.write(f".move({x:.4f},{y:.4f},0)\n")

    def _write_linear_polymer_deterministic(self, write_f, monomer_set: list) -> None:
        """Fallback deterministic placement for linear polymers."""
        offset_cum = 0
        for indexii in range(len(monomer_set)):
            monomer_name = monomer_set[indexii][:-3] if monomer_set[indexii].endswith('.lt') else monomer_set[indexii]
            write_f.write(f"    monomer[{indexii}] = new {monomer_name}")
            if indexii > 0:
                write_f.write(f".rot({self.poly.rotate*(indexii%2)},1,0,0).move({offset_cum:.4f},0,0)")
            write_f.write("\n")
            self.poly.evaluate_offset(f"{monomer_name}.lt")
            offset_cum += self.poly.offset

    def _get_ff_import(self) -> str:
        """Get force field import statement."""
        ff_map = {
            "gaff": 'import "gaff.lt"\n\n',
            "gaff2": 'import "gaff2.lt"\n\n',
            "lopls": 'import "loplsaa.lt"\n\n',
            "dreiding": 'import "dreiding.lt"\n\n',
            "compass": 'import "compass_published.lt"\n\n',
        }
        return ff_map.get(self.poly.force_field, 'import "oplsaa.lt"\n\n')

    def _get_ff_import_and_inherits(self) -> Tuple[str, str]:
        """Get force field import and inheritance statements."""
        ff_map = {
            "gaff": ('import "gaff.lt"\n', "GAFF"),
            "gaff2": ('import "gaff2.lt"\n', "GAFF2"),
            "lopls": ('import "loplsaa.lt"\n', "OPLSAA"),
            "dreiding": ('import "dreiding.lt"\n', "DREIDING"),
            "compass": ('import "compass_published.lt"\n', "COMPASS"),
        }
        return ff_map.get(self.poly.force_field, ('import "oplsaa.lt"\n', "OPLSAA"))

    def _count_total_monomers(self) -> int:
        """Count total monomers across all models."""
        total = 0
        for modelii in self.poly.model:
            is_molecule = hasattr(modelii, '_is_molecule') and modelii._is_molecule
            if is_molecule:
                total += modelii.Count
            elif modelii.dop > 1:
                total += len(modelii.sequenceSet) * modelii.dop
            else:
                total += len(modelii.sequenceSet)
        return max(total, 1)
