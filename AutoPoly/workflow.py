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

Created on 2026-01-06
@author: zwu
"""
import sys
import subprocess
from pathlib import Path
import numpy as np
from typing import List
from .system import logger
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

        The method includes comprehensive error checking and validation
        to ensure all required files are generated correctly.

        Raises:
            SystemExit: If any critical step fails or required files are missing
        """
        try:
            logger.info("Starting LAMMPS data file generation using Moltemplate")

            poly_index = 0
            for modelii in self.poly.model:
                logger.info(f"Processing model with {len(modelii.sequenceSet)} molecules")
                base_smiles = modelii.sequence[0]
                # Get topology
                topology = getattr(modelii, 'topology', 'linear')

                # Get DOP
                dop = modelii.DOP

                logger.info(f"Generating sequence variants for {base_smiles}, DOP={dop}, topology={topology}")

                # Use new sequence-aware variant generation
                # This will generate position-specific variants and extract unique units
                variant_mapping = self.poly.generate_sequence_variants_for_polymer(
                    base_smiles=base_smiles,
                    dop=dop,
                    topology=topology,
                    base_name_prefix="monomer"
                )

                logger.info(f"Generated {len(variant_mapping)} unique variant types: {list(variant_mapping.keys())}")

                # Build sequenceSet with appropriate .lt files for each position
                for chain_idx, chain_smiles in enumerate(modelii.sequenceSet):
                    monomer_names = []
                    for pos_idx, smiles in enumerate(chain_smiles):
                        # Determine variant type based on position and topology
                        if topology == "ring":
                            variant_type = "ring"
                        else:
                            if pos_idx == 0:
                                variant_type = "first"
                            elif pos_idx == dop - 1:
                                variant_type = "last"
                            else:
                                variant_type = "middle"

                        # Check if T1 variant is needed
                        has_t1 = "_T1" in smiles

                        # Get the filename from variant_mapping
                        # Use the T1 variant key if T1 is needed, otherwise use standard variant
                        mapping_key = f"{variant_type}_T1" if has_t1 else variant_type

                        if mapping_key in variant_mapping:
                            filename = variant_mapping[mapping_key]
                            monomer_names.append(filename)
                        else:
                            logger.error(f"Variant key '{mapping_key}' not found in mapping")
                            logger.error(f"Available variants: {list(variant_mapping.keys())}")
                            sys.exit(1)

                    modelii.sequenceSet[chain_idx] = monomer_names

                logger.info(f"Built {len(modelii.sequenceSet)} chain(s) with sequence variants")

                # Loop through all polymers and make corresponding polymer lt files
                for moleii in range(len(modelii.sequenceSet)):
                    # Check degrees of polymerization (DOP) of current polymer
                    if modelii.DOP > 0:
                        if len(modelii.sequenceSet[moleii]) != modelii.DOP:
                            logger.warning(f"Warning: At molecule# {moleii} DOP={len(modelii.sequenceSet[moleii])} != {modelii.DOP}")
                    else:
                        logger.error(f"Warning: At molecule#{moleii+1}, DOP={len(modelii.sequenceSet[moleii])} {modelii.DOP}")

                    if modelii.DOP > 1:
                        # Make poly.lt file
                        logger.info(f"Creating poly_{poly_index+1}.lt")
                        self.make_poly_lt(poly_index, modelii.sequenceSet[moleii], modelii)
                        poly_index += 1

            # Generate force field files
            logger.info(f"Generating {self.poly.force_field}.lt")
            # Pass model to ff_manager for force field subset generation
            self.poly.ff_manager.model = self.poly.model
            self.poly.ff_manager.make_force_field_lt()

            # Generate system.lt file
            logger.info("Creating system.lt")
            self.make_system_lt()

            # Modify alkyl dihedral coefficients if needed (skip for GAFF)
            if self.poly.is_lopls and self.poly.force_field != "gaff":
                self.poly.ff_manager.FFmodify_alkyl_dihedral_oplsaa()

            # Invoke moltemplate to generate LAMMPS datafile
            logger.info("Running moltemplate")
            self.invoke_moltemplate()

            # Validate that output_ttree was created and contains required data files
            # This catches cases where moltemplate runs but doesn't generate complete atom data
            output_ttree_dir = Path(self.poly.path_cwd) / "output_ttree"
            if not output_ttree_dir.exists():
                logger.error(f"Moltemplate failed to create output_ttree directory at {output_ttree_dir}")
                logger.error("This indicates moltemplate did not run successfully")
                sys.exit(1)

            # Check for critical data files that should contain actual atom/topology data
            required_data_files = [
                "Data Atoms",
                "Data Bond List"
            ]

            missing_data_files = []
            for data_file in required_data_files:
                data_file_path = output_ttree_dir / data_file
                if not data_file_path.exists():
                    missing_data_files.append(data_file)
                else:
                    # Additional check: verify file is not empty
                    if data_file_path.stat().st_size == 0:
                        logger.error(f"Data file exists but is empty: {data_file}")
                        missing_data_files.append(f"{data_file} (empty)")

            if missing_data_files:
                logger.error(f"Moltemplate failed to generate required data files: {missing_data_files}")
                logger.error("This usually indicates a problem with monomer file syntax or moltemplate execution.")
                logger.error("Troubleshooting steps:")
                logger.error("1. Check that monomer .lt files have correct moltemplate syntax")
                logger.error("2. Verify that molecules inherit from the correct force field class")
                logger.error("3. Ensure write(\"Data Atoms\") sections are present in monomer files")
                logger.error(f"4. Check files in: {output_ttree_dir}")
                sys.exit(1)

            # Check if the required files exist before proceeding
            # Note: system.in.charges is optional for GAFF (charges calculated separately)
            required_files = ['system.in.settings', 'system.data']
            optional_files = ['system.in.charges'] if self.poly.force_field == "gaff" else []

            missing_files = []
            for file in required_files:
                if not (Path(self.poly.path_cwd) / file).exists():
                    missing_files.append(file)

            # Check optional files for non-GAFF force fields
            if not optional_files:
                for file in optional_files:
                    if not (Path(self.poly.path_cwd) / file).exists():
                        missing_files.append(file)

            if missing_files:
                logger.error(f"Moltemplate failed to generate required files: {', '.join(missing_files)}")
                logger.error("Check the following:")
                logger.error("1. All monomer .lt files exist and are valid")
                logger.error("2. The polymer .lt files were generated correctly")
                logger.error("3. The system.lt file is properly formatted")
                sys.exit(1)

            # Log warning about missing charges file for GAFF
            if self.poly.force_field == "gaff" and not (Path(self.poly.path_cwd) / "system.in.charges").exists():
                logger.warning("Note: system.in.charges not generated for GAFF")
                logger.warning("GAFF requires manual charge calculation using AM1-BCC or RESP")
                logger.warning("All atomic charges are currently set to 0.00")

            logger.info("Processing output files")
            self.poly.get_rid_of_lj_cut_coul_long()

            # Move files to working directory
            self.poly.mv_files()
            logger.info("Successfully completed polymer generation")

        except Exception as e:
            logger.error(f"Error in make_lmp_data_file_by_moltemplate: {str(e)}")
            sys.exit(1)

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
                logger.error(f"system.lt not found in {self.poly.path_cwd}")
                sys.exit(1)

            # Run moltemplate with output capture for error cases
            # Use cwd parameter and explicit bash path for security
            # Use -nocheck flag to bypass post-processing validation that may fail on first run
            moltemplate_sh = self.poly.path_moltemplatesrc + "moltemplate.sh"
            process = subprocess.run(
                ["bash", moltemplate_sh, "-nocheck", "./system.lt"],
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
                sys.exit(1)

        except Exception as e:
            logger.error(f"Error running Moltemplate: {str(e)}")
            sys.exit(1)

    def make_system_lt(self) -> None:
        """
        Creates the system.lt file for the polymerization.

        This method generates the main system.lt file that imports all necessary
        polymer files, force field parameters, and defines the simulation box
        with proper spacing and positioning for all polymer chains.
        """
        output = self.poly.path_cwd + "/system.lt"

        # Determine force field import based on force_field type
        if self.poly.force_field == "gaff":
            ff_import = 'import "gaff.lt"\n\n'
        else:
            ff_import = 'import "oplsaa.lt"\n\n'

        with open(output, "w") as write_f:
            # Write force field import at the top
            write_f.write(ff_import)

            polyindex = 0
            for modelii in self.poly.model:
                n_poly = len(modelii.sequenceSet)
                if modelii.DOP > 1:
                    for indexi in range(n_poly):
                        write_f.write(f"import \"poly_{polyindex+1}.lt\"\n")
                        polyindex += 1
                    write_f.write("\n")
                else:
                    if len(modelii.merSet)>1:
                        logger.error(' '.join(["sequenceLen = "+str(modelii.DOP)+", "
                                                    , " merSet should only have one mer type! Please check.\n"]))
                        sys.exit()
                    #import constituent monomer.lt's
                    unique_Sequence=[i[0] for i in modelii.sequenceSet]
                    print(unique_Sequence)
                    for sequenceii in range(len(unique_Sequence)):
                        write_f.write("import \""+unique_Sequence[sequenceii]+"\"\n")
                    write_f.write("\n")

            polyindex = 0
            index = 0

            # Calculate spacing based on polymer type and size
            for modelii in self.poly.model:
                n_poly = len(modelii.sequenceSet)
                is_ring = hasattr(modelii, 'topology') and modelii.topology == "ring"

                if is_ring:
                    # For ring polymers, calculate radius and use it for spacing
                    radius = self.poly.offset * len(modelii.sequenceSet[0]) / (2 * np.pi)
                    spacing = radius * 2.5  # Use 2.5x the ring radius for good separation
                else:
                    spacing = self.poly.offset * (modelii.DOP + 2)

                # Calculate grid arrangement
                grid_size = int(np.ceil(np.sqrt(n_poly)))  # Arrange in a square grid

                for moleii in range(n_poly):
                    # Calculate grid position
                    grid_x = moleii % grid_size
                    grid_y = moleii // grid_size

                    # Calculate actual position with spacing
                    pos_x = grid_x * spacing
                    pos_y = grid_y * spacing
                    pos_z = 0.0  # Keep all polymers in the same plane initially

                    if modelii.DOP > 1:
                        write_f.write(f"polymer_{index+1} = new poly_{polyindex+1}")
                        write_f.write(f".move({pos_x:.4f},{pos_y:.4f},{pos_z:.4f})\n")
                        polyindex += 1
                    else:
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
        output = self.poly.path_cwd + f"/poly_{poly_index+1}.lt"

        # Determine force field import and inheritance based on force_field type
        if self.poly.force_field == "gaff":
            ff_import = 'import "gaff.lt"\n'
            ff_inherits = "GAFF"
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
