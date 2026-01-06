#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Polymerization Module for AutoPoly Package

This module provides the core Polymerization class for generating polymer structures
using Moltemplate and preparing them for LAMMPS molecular dynamics simulations.

The Polymerization class handles:
- Polymer structure generation from monomer templates
- Moltemplate integration for LAMMPS data file creation
- Force field parameter management (OPLS-AA)
- Support for various polymer topologies and tacticity
- File organization and output management

Key Features:
- Atomistic polymer modeling with OPLS-AA force field
- Support for linear and ring polymer topologies
- Tacticity control (atactic, isotactic, syndiotactic)
- Automatic monomer bank management
- LAMMPS input file generation
- Comprehensive error handling and logging

Dependencies:
- Moltemplate: For generating LAMMPS data files
- OPLS-AA force field parameters
- Monomer bank with .lt template files

Created on Fri Dec 21 12:19:08 2018
@author: zwu
"""
import sys
import os
from pathlib import Path
import subprocess
import shutil
import re
import numpy as np
from typing import List, Optional, Dict, Any, Set, Tuple
from .system import logger
from .monomer_generator import MonomerGenerator

class Polymerization:
    """
    Core polymerization class for generating polymer structures using Moltemplate.
    
    This class manages the complete workflow for creating polymer structures
    from monomer templates and generating LAMMPS input files for molecular
    dynamics simulations. It handles file management, force field parameters,
    and integration with external tools like Moltemplate.
    
    Attributes:
        name (str): Name of the polymerization project
        system (object): System object containing path information
        path_cwd (str): Current working directory for the project
        path_master (str): Path to external dependencies
        path_moltemplatesrc (str): Path to Moltemplate source
        path_oplsaaprm (str): Path to OPLS-AA force field parameters
        is_lopls (bool): Whether to use LOPLS force field
        model (list): List of polymer models to generate
        rotate (float): Rotation angle for monomer placement
        offset_spacing (float): Spacing between polymer chains
        offset (float): Offset distance for monomer placement
        packingL_spacing (float): Packing length spacing
        moltemplate_box_size (float): Box size for Moltemplate
        FFmodify_alkylDihedral (np.array): Modified dihedral parameters
    """
    
    def __init__(self, name: str = None, system: object = None, model: list = None,
                 run: bool = True, is_lopls: bool = False,
                 force_field: str = "oplsaa") -> None:
        """
        Initialize the Polymerization class.

        Args:
            name (str, optional): Name of the polymerization project. Defaults to None.
            system (object, optional): System object containing folder path. Defaults to None.
            model (list, optional): List of models for polymerization. Defaults to None.
            run (bool, optional): Flag to run the process immediately. Defaults to True.
            is_lopls (bool, optional): Whether to use LOPLS force field. Defaults to False.
                                       DEPRECATED: Use force_field="lopls" instead.
            force_field (str, optional): Force field to use. Options: "oplsaa", "gaff", "lopls".
                                        Defaults to "oplsaa".

        Raises:
            SystemExit: If required directories or files are not found
        """
        # Deprecation warning for is_lopls
        if is_lopls:
            logger.warning("The 'is_lopls' parameter is deprecated. Use 'force_field=\"lopls\"' instead.")
            if force_field == "oplsaa":  # Only override if user hasn't explicitly set force_field
                force_field = "lopls"

        # Validate force_field parameter
        valid_force_fields = ["oplsaa", "gaff", "lopls"]
        if force_field not in valid_force_fields:
            logger.error(f"Invalid force_field '{force_field}'. Must be one of: {valid_force_fields}")
            sys.exit(1)

        self.name = name
        self.system = system
        self.path_cwd = f"{self.system.get_folder_path()}/{self.name}/moltemplate/"
        self.path_master = f"{Path(__file__).parent.resolve()}/extern/"
        self.path_moltemplatesrc = f"{self.path_master}moltemplate/src/"

        # SMILES to monomer name cache for dynamic generation
        self._generated_smiles = {}
        self._smiles_to_name_counter = 0

        self.is_lopls = is_lopls
        self.force_field = force_field

        # Set force field parameter path based on force_field type
        if force_field == "gaff":
            self.path_oplsaaprm = f"{self.path_master}moltemplate/common/gaff.lt"
        elif force_field == "lopls":
            self.path_oplsaaprm = f"{self.path_master}moltemplate/loplsaa.prm"
        else:  # oplsaa
            self.path_oplsaaprm = f"{self.path_master}moltemplate/oplsaa.prm"

        logger.info(f"\n'you are now using parameter set of {self.path_oplsaaprm}\n")
        self.model = model
        self.rotate = 90.0
        self.offset_spacing = 2.0
        self.offset = 4.0
        self.packingL_spacing = 5.0
        self.moltemplate_box_size = 400.0

        # Modified alkyl dihedral parameters (Kj/mol -> kcal/mol conversion)
        self.FFmodify_alkylDihedral = np.array([0.6446926386, -0.2143420172, 0.1782194073, 0.0])
        
        # Create working directory before proceeding
        self.create_working_directory()
        
        logger.info(f"\n'you are now using extern path of {self.path_master}\n")

        if run:
            self.make_lmp_data_file_by_moltemplate()

    def create_working_directory(self) -> None:
        """
        Create and manage the working directory structure for the polymerization.
        
        This method sets up the directory structure needed for the polymerization
        process. It creates the main project directory and the moltemplate
        subdirectory where all intermediate files will be stored.
        
        The method handles existing directories by prompting the user to either
        delete and recreate them or choose a different project name.
        
        Raises:
            SystemExit: If user chooses not to overwrite existing directory
        """
        base_path = Path(self.system.get_folder_path())
        polymer_path = base_path / self.name
        moltemplate_path = polymer_path / "moltemplate"
        
        # Check if base directory exists
        if polymer_path.exists():
            response = input(f"{polymer_path} folder exists, delete and make new?(y/n) ")
            if response.lower() == 'y':
                logger.info(f"removing {polymer_path}")
                shutil.rmtree(polymer_path)
            else:
                logger.error("Please remove the existing folder or choose a different name.")
                sys.exit(1)
        
        # Create directory structure
        moltemplate_path.mkdir(parents=True, exist_ok=True)

    def create_folder(self) -> None:
        """Creates the working directory for the polymerization."""
        path = Path(self.path_cwd)
        parent_path = path.parent
        
        if parent_path.exists():
            response = input(f"{parent_path} folder exists, delete and make new?(y/n) ")
            if response.lower() == 'y':
                logger.info(f"removing {parent_path}")
                import shutil
                shutil.rmtree(parent_path)
            else:
                logger.error("Please remove the existing folder or choose a different name.")
                sys.exit(1)
        
        # Create the directory structure
        path.mkdir(parents=True, exist_ok=True)

    def set_tacticity(self, tacticity: str) -> None:
        """Sets the tacticity of the polymer.

        Args:
            tacticity (str): The tacticity to set.
        """
        self.tacticity = tacticity

    def n_monomer_atoms(self, merltfile: str) -> int:
        """
        Count the number of monomer atoms in the specified .lt file.

        This method parses a Moltemplate monomer file (.lt) and counts the
        number of atoms defined in the "Data Atoms" block. This information
        is used for polymer structure generation and validation.

        Args:
            merltfile (str): The name of the monomer .lt file.

        Returns:
            int: The number of monomer atoms.

        Raises:
            SystemExit: If the monomer file cannot be opened or found.
        """
        n_monomer_atoms = 0
        merltfile_path = Path(self.path_cwd) / merltfile
        
        if merltfile_path.is_file():
            is_inside_block = False
            with open(merltfile_path) as f:
                while True: 
                    line = f.readline() 
                    if line.strip() == "write(\"Data Atoms\") {":
                        is_inside_block = True
                        line = f.readline() 
                    elif line.strip() == "}":
                        is_inside_block = False
                    
                    if is_inside_block:
                        n_monomer_atoms = n_monomer_atoms + 1

                    if not line: 
                        break
        else:
            logger.error(f"in MoltemplateLmpData::n_monomerAtoms(): {merltfile_path} file cannot open.")
            sys.exit(1)

        return n_monomer_atoms

    def extract_element_from_atom(self,atom_string):
        """
        Extract element name from atom identifier string.
        
        Args:
            atom_string (str): Atom identifier like "$atom:C1", "$atom:H16", etc.
        
        Returns:
            str: Element name (e.g., "C", "H", "Si", "Fe")
            None: If no match found
        """
        pattern = r'\$atom:([A-Z][a-z]?)\d*'
        match = re.search(pattern, atom_string)
        
        if match:
            return match.group(1)  # Return the captured element name
        else:
            return None

    def read_lt_end_atoms(self,lt_file):
        """Read the first and second atoms from a .lt file.
        
        Args:
            lt_file (str): Path to the .lt file
            
        Returns:
            tuple: (first_atom, second_atom) where each atom is a dict with:
                - element: atom element symbol
                - atom_type: OPLS atom type
                - x, y, z: coordinates
        """
        first_atom = None
        second_atom = None
        in_atoms_block = False
        
        with open(lt_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                if line == 'write("Data Atoms") {':
                    in_atoms_block = True
                    continue
                elif line == '}':
                    in_atoms_block = False
                    continue
                    
                if in_atoms_block and line:
                    # Parse atom line: $atom:atom_id $mol:... atom_type charge x y z # element
                    parts = line.split()
                    if len(parts) >= 7:  # Ensure we have all required fields
                        atom_type = parts[2]
                        x = float(parts[4])
                        y = float(parts[5])
                        z = float(parts[6])
                        #print(atom_type)
                        element = self.extract_element_from_atom(parts[0])  # Last field after #
                        #print(element)
                        atom_data = {
                            'element': element,
                            'atom_type': atom_type,
                            'x': x,
                            'y': y,
                            'z': z
                        }
                        
                        if first_atom is None:
                            first_atom = atom_data['element']
                        elif second_atom is None:
                            second_atom = atom_data['element']
                            break  # We have both atoms, no need to continue
        
        if first_atom is None or second_atom is None:
            raise ValueError(f"Could not find both end atoms in {lt_file}")
            
        return first_atom+"1", second_atom+"2"
    

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
            for modelii in self.model:
                logger.info(f"Processing model with {len(modelii.sequenceSet)} molecules")

                # Generate monomers from pSMILES/SMILES (once for all chains)
                logger.info("Generating monomers from pSMILES/SMILES...")

                # Get unique pSMILES/SMILES (preserving order)
                unique_psmiles = []
                seen = set()
                for psmiles in modelii.sequence:
                    if psmiles not in seen:
                        seen.add(psmiles)
                        unique_psmiles.append(psmiles)

                # Generate .lt files for each unique pSMILES/SMILES
                psmiles_to_base_name = {}
                for psmiles in unique_psmiles:
                    base_name = self.generate_monomer_from_psmiles(psmiles)
                    psmiles_to_base_name[psmiles] = base_name

                # Update sequence with generated monomer base names
                # set_Sequence() will add the le/re/i/_T1 suffixes
                modelii.sequence = [psmiles_to_base_name[psmiles] for psmiles in modelii.sequence]

                # Regenerate sequenceSet with new monomer names
                modelii.set_Sequence()

                logger.info(f"Generated {len(unique_psmiles)} unique monomer(s)")

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
            logger.info(f"Generating {self.force_field}.lt")
            self.make_force_field_lt()

            # Generate system.lt file
            logger.info("Creating system.lt")
            self.make_system_lt()

            # Modify alkyl dihedral coefficients if needed (skip for GAFF)
            if self.is_lopls and self.force_field != "gaff":
                self.FFmodify_alkyl_dihedral_oplsaa()

            # Invoke moltemplate to generate LAMMPS datafile
            logger.info("Running moltemplate")
            self.invoke_moltemplate()

            # Validate that output_ttree was created and contains required data files
            # This catches cases where moltemplate runs but doesn't generate complete atom data
            output_ttree_dir = Path(self.path_cwd) / "output_ttree"
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
            optional_files = ['system.in.charges'] if self.force_field == "gaff" else []

            missing_files = []
            for file in required_files:
                if not (Path(self.path_cwd) / file).exists():
                    missing_files.append(file)

            # Check optional files for non-GAFF force fields
            if not optional_files:
                for file in optional_files:
                    if not (Path(self.path_cwd) / file).exists():
                        missing_files.append(file)

            if missing_files:
                logger.error(f"Moltemplate failed to generate required files: {', '.join(missing_files)}")
                logger.error("Check the following:")
                logger.error("1. All monomer .lt files exist and are valid")
                logger.error("2. The polymer .lt files were generated correctly")
                logger.error("3. The system.lt file is properly formatted")
                sys.exit(1)

            # Log warning about missing charges file for GAFF
            if self.force_field == "gaff" and not (Path(self.path_cwd) / "system.in.charges").exists():
                logger.warning("Note: system.in.charges not generated for GAFF")
                logger.warning("GAFF requires manual charge calculation using AM1-BCC or RESP")
                logger.warning("All atomic charges are currently set to 0.00")

            logger.info("Processing output files")
            self.get_rid_of_lj_cut_coul_long()

            # Move files to working directory
            self.mv_files()
            logger.info("Successfully completed polymer generation")
            
        except Exception as e:
            logger.error(f"Error in make_lmp_data_file_by_moltemplate: {str(e)}")
            sys.exit(1)

    def get_rid_of_lj_cut_coul_long(self) -> None:
        """Removes lj/cut/coul/long from the settings file."""
        in_=self.path_cwd+"system.in.settings"
        out=self.path_cwd+"tmp.data"

        in_path=Path(in_)
        if not in_path.is_file():
            logger.error(' '.join(["system.in.setting does not exist plase check ",in_]))
            sys.exit()

        # Use context manager to ensure file is properly closed even if exception occurs
        with open(in_,'r') as read_f, open(out, "w") as write_f:
            while True:
                line = read_f.readline()
                if line.strip()=="":
                    write_f.write("\n")
                elif line.strip().split()[0]=="pair_coeff":
                    #write_f.write("    pair_coeff ")
                    space_i=0
                    for ii in line.split():
                        if ii=="lj/cut/coul/long":
                            continue
                        else:
                            if space_i==0:
                                write_f.write("    ")
                                space_i=space_i+1
                            write_f.write(ii+" ")
                    write_f.write("\n")
                else:
                    write_f.write(line)

                if not line:
                    break
        mv="rm "+in_+";mv "+out+" "+in_
        os.system(mv)

    def mv_files(self) -> None:
        """Moves generated files to the appropriate directories."""
        try:
            # Define paths
            moltemplate_dir = Path(self.path_cwd)
            parent_dir = moltemplate_dir.parent
            
            # Create output and input directories if they don't exist
            output_dir = parent_dir / "output"
            input_dir = parent_dir / "input"
            output_dir.mkdir(exist_ok=True)
            input_dir.mkdir(exist_ok=True)

            # Copy data files to parent directory
            for file in ["system.data", "system.in.charges", "system.in.settings", "system.in", "system.in.init"]:
                if (moltemplate_dir / file).exists():
                    shutil.copy2(moltemplate_dir / file, parent_dir)

            # Move files to output directory
            for pattern in ["system.in*", "system*data", "output_ttree"]:
                for file in moltemplate_dir.glob(pattern):
                    shutil.move(str(file), str(output_dir))

            # Move files to input directory
            for pattern in ["*.lt", "*.prm"]:
                for file in moltemplate_dir.glob(pattern):
                    shutil.move(str(file), str(input_dir))

        except Exception as e:
            logger.error(f"Error moving files: {str(e)}")
            sys.exit(1)

    def evaluate_box_len(self):
        in_=path_cwd+"system.data"
        dubVar=0
        lmin=0
        lmax=0

    def invoke_moltemplate(self) -> None:
        """Invokes Moltemplate to generate the LAMMPS data file."""
        try:
            # First check if system.lt exists
            system_lt = Path(self.path_cwd) / "system.lt"
            if not system_lt.exists():
                logger.error(f"system.lt not found in {self.path_cwd}")
                sys.exit(1)

            # Run moltemplate with output capture for error cases
            # Use cwd parameter and explicit bash path for security
            # Use -nocheck flag to bypass post-processing validation that may fail on first run
            moltemplate_sh = self.path_moltemplatesrc + "moltemplate.sh"
            process = subprocess.run(
                ["bash", moltemplate_sh, "-nocheck", "./system.lt"],
                cwd=self.path_cwd,
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
        """Creates the system.lt file for the polymerization."""
        output = self.path_cwd + "/system.lt"

        # Determine force field import based on force_field type
        if self.force_field == "gaff":
            ff_import = 'import "gaff.lt"\n\n'
        else:
            ff_import = 'import "oplsaa.lt"\n\n'

        with open(output, "w") as write_f:
            # Write force field import at the top
            write_f.write(ff_import)

            polyindex = 0
            for modelii in self.model:
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
            for modelii in self.model:
                n_poly = len(modelii.sequenceSet)
                is_ring = hasattr(modelii, 'topology') and modelii.topology == "ring"
                
                if is_ring:
                    # For ring polymers, calculate radius and use it for spacing
                    radius = self.offset * len(modelii.sequenceSet[0]) / (2 * np.pi)
                    spacing = radius * 2.5  # Use 2.5x the ring radius for good separation
                else:
                    spacing = self.offset * (modelii.DOP + 2)

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
                grid_size * spacing for modelii in self.model
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
        """Creates a poly.lt file for the specified polymer.

        Args:
            poly_index (int): The index of the polymer.
            monomer_set (list): The list of monomers in the polymer.
            model (object): The polymer model object containing topology information.
        """
        output = self.path_cwd + f"/poly_{poly_index+1}.lt"

        # Determine force field import and inheritance based on force_field type
        if self.force_field == "gaff":
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
                radius = self.offset * n_monomers / (2 * np.pi)
                
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
                    monomer_bank = Path(self.path_cwd)
                    merltfile_path_1 = monomer_bank / f"{monomer_name_1}.lt"
                    _, second_atom = self.read_lt_end_atoms(merltfile_path_1)
                    merltfile_path_2 = monomer_bank / f"{monomer_name_2}.lt"
                    first_atom, _ = self.read_lt_end_atoms(merltfile_path_2)
                    write_f.write(f"      $bond:b{i+1}  $atom:monomer[{i}]/{second_atom}  $atom:monomer[{next_i}]/{first_atom}\n")
                write_f.write("    }\n")

            else:
                # Original linear polymer code
                offset_cum = 0
                for indexii in range(len(monomer_set)):
                    
                    monomer_name = monomer_set[indexii][:-3] if monomer_set[indexii].endswith('.lt') else monomer_set[indexii]
                    
                    write_f.write(f"    monomer[{indexii}] = new {monomer_name}")
                    if indexii > 0:
                        write_f.write(f".rot({self.rotate*(indexii%2)},1,0,0).move({offset_cum:.4f},0,0)")
                    write_f.write("\n")

                    self.evaluate_offset(f"{monomer_name}.lt")
                    offset_cum += self.offset

                # Add bonds for linear polymer
                write_f.write("\n    write('Data Bond List') {\n")
                for indexii in range(len(monomer_set)-1):

                    monomer_name_1 = monomer_set[indexii][:-3] if monomer_set[indexii].endswith('.lt') else monomer_set[indexii]
                    monomer_name_2 = monomer_set[indexii+1][:-3] if monomer_set[indexii+1].endswith('.lt') else monomer_set[indexii+1]
                    # Determine monomer bank path based on force field
                    monomer_bank = Path(self.path_cwd)
                    merltfile_path_1 = monomer_bank / f"{monomer_name_1}.lt"
                    _, second_atom = self.read_lt_end_atoms(merltfile_path_1)
                    merltfile_path_2 = monomer_bank / f"{monomer_name_2}.lt"
                    first_atom, _ = self.read_lt_end_atoms(merltfile_path_2)
                    write_f.write(f"      $bond:b{indexii+1}  $atom:monomer[{indexii}]/{second_atom}  $atom:monomer[{indexii+1}]/{first_atom}\n")
                write_f.write("    }\n")

            write_f.write(f"\n}} # poly_{poly_index+1}\n")

    def evaluate_offset(self, merltfile: str) -> None:
        """Evaluates the offset distance based on the specified merlt file.

        Args:
            merltfile (str): The name of the merlt file.
        """
        # Determine monomer bank path based on force field
        MonomerBank = Path(self.path_cwd)
        merltfile_Path = MonomerBank / merltfile
        if merltfile_Path.is_file():
            C1=[]
            C2=[]
            dubVar=0
            with open(merltfile_Path) as f:
                while True: 
                    line = f.readline() 
                    if line.strip()=="write(\"Data Atoms\") {":
                        # C1 coordinates
                        line = f.readline() 
                        for i in range(3):
                            C1.append(float(line.split()[i+4]))
                        # C2 coordinates
                        line = f.readline() 
                        for i in range(3):
                            C2.append(float(line.split()[i+4]))
                        
                        # calculate C1-C2 distance
                        self.offset=np.linalg.norm(np.array(C1)-np.array(C2))+self.offset_spacing
                        
                        return
            
    def make_force_field_lt(self) -> None:
        """Creates the force field .lt file based on the force_field type."""
        try:
            if self.force_field == "gaff":
                self.make_gaff_subset()
            else:  # oplsaa or lopls
                self.make_oplsaa_subset()

                # Invoke oplsaa_moltemplate.py to make oplsaa.lt with suppressed output
                ff_name = "oplsaa.lt"
                ff_subset = self.path_cwd + "oplsaa_subset.prm"
                ff_py_script = self.path_moltemplatesrc + "oplsaa_moltemplate.py"

                # Redirect both stdout and stderr to devnull
                # Use cwd parameter instead of shell command to avoid shell injection
                with open(os.devnull, 'w') as devnull:
                    return_code = subprocess.call([ff_py_script, ff_subset],
                                                cwd=self.path_cwd,
                                                stdout=devnull,
                                                stderr=devnull)

                if return_code != 0:
                    logger.error(f"Failed to generate {ff_name} file. Check oplsaa_subset.prm for errors.")
                    sys.exit(1)

        except Exception as e:
            logger.error(f"Error in make_force_field_lt: {str(e)}")
            sys.exit(1)

    def _extract_gaff_atom_types_from_monomers(self) -> Set[str]:
        """Extract unique GAFF atom types used in monomer files.

        Parses monomer .lt files to find atom type references (e.g., @atom:c3, @atom:ce)
        and returns the set of unique type names.

        Returns:
            Set[str]: Set of atom type names (e.g., {'c3', 'ce', 'hc', 'o'})
        """
        atom_types = set()

        for modelii in self.model:
            for monomerii in range(len(modelii.sequenceSet)):
                monomerSet = modelii.sequenceSet[monomerii]

                for vecii in range(len(monomerSet)):
                    # Determine monomer bank path based on force field
                    MonomerBank = Path(self.path_cwd)

                    merltfile_Path = MonomerBank / monomerSet[vecii]

                    if merltfile_Path.is_file():
                        mono = str(merltfile_Path)
                        read_switch = False

                        try:
                            with open(mono) as f:
                                while True:
                                    line = f.readline()

                                    if line.strip() == 'write("Data Atoms") {':
                                        read_switch = True
                                        continue
                                    elif line.strip() == "}":
                                        read_switch = False
                                        break

                                    if read_switch:
                                        stringvector = line.split()

                                        if len(stringvector) >= 3:
                                            # Extract atom type from @atom:c3 format
                                            atom_type_full = stringvector[2]
                                            if atom_type_full.startswith("@atom:"):
                                                # Strip @atom: prefix to get type name
                                                atom_type = atom_type_full.split(":")[1]
                                                atom_types.add(atom_type)

                                    if not line:
                                        break
                        except Exception as e:
                            logger.warning(f"Error reading monomer file {mono}: {str(e)}")
                            continue
                    else:
                        logger.error(' '.join(["Monomer (" + monomerSet[vecii] + ") does NOT exist. \n",
                                               "Please check the following path to the file\n" + str(merltfile_Path) + "\n"]))
                        sys.exit()

        return atom_types

    def _parse_gaff_lt_sections(self, gaff_file: str) -> Dict[str, List[str]]:
        """Parse gaff.lt file into sections.

        Args:
            gaff_file: Path to gaff.lt file

        Returns:
            Dict with section names as keys and list of lines as values:
            {
                'header': [...],
                'masses': [],      # write_once("Data Masses")
                'pair_coeffs': [], # First write_once("In Settings")
                'bond_coeffs': [], # Second write_once("In Settings")
                'bond_definitions': [], # write_once("Data Bonds By Type")
                'angle_definitions': [], # write_once("Data Angles By Type")
                'angle_coeffs': [], # Third write_once("In Settings")
                'dihedral_definitions': [],
                'dihedral_coeffs': [],
                'improper_definitions': [],
                'improper_coeffs': [],
                'init': [],        # write_once("In Init")
            }
        """
        sections = {
            'header': [],
            'masses': [],
            'pair_coeffs': [],
            'bond_coeffs': [],
            'bond_definitions': [],
            'angle_definitions': [],
            'angle_coeffs': [],
            'dihedral_definitions': [],
            'dihedral_coeffs': [],
            'improper_definitions': [],
            'improper_coeffs': [],
            'init': []
        }

        current_section = 'header'
        in_settings_count = 0
        braces_depth = 0

        with open(gaff_file, 'r') as f:
            for line in f:
                stripped = line.strip()

                # Track braces to determine section boundaries
                if '{' in stripped:
                    braces_depth += 1
                if '}' in stripped:
                    braces_depth -= 1
                    if braces_depth == 0:
                        # End of section, go back to header
                        if current_section != 'header':
                            # Add closing brace to current section
                            sections[current_section].append(line)
                        current_section = 'header'
                        continue

                # Section detection
                if 'write_once("Data Masses")' in stripped:
                    current_section = 'masses'
                elif 'write_once("In Settings")' in stripped:
                    in_settings_count += 1
                    if in_settings_count == 1:
                        current_section = 'pair_coeffs'
                    elif in_settings_count == 2:
                        current_section = 'bond_coeffs'
                    elif in_settings_count == 3:
                        current_section = 'angle_coeffs'
                    elif in_settings_count == 4:
                        current_section = 'dihedral_coeffs'
                    elif in_settings_count == 5:
                        current_section = 'improper_coeffs'
                elif 'write_once("Data Bonds By Type")' in stripped:
                    current_section = 'bond_definitions'
                elif 'write_once("Data Angles By Type")' in stripped:
                    current_section = 'angle_definitions'
                elif 'write_once("Data Dihedrals By Type")' in stripped:
                    current_section = 'dihedral_definitions'
                elif 'write_once("Data Impropers By Type' in stripped:
                    current_section = 'improper_definitions'
                elif 'write_once("In Init")' in stripped:
                    current_section = 'init'

                # Add line to current section
                if current_section in sections:
                    sections[current_section].append(line)

        return sections

    def _filter_masses_section(self, masses_lines: List[str], atom_types: Set[str]) -> List[str]:
        """Filter Data Masses section to keep only used atom types.

        Args:
            masses_lines: Lines from Data Masses section
            atom_types: Set of atom type names (e.g., {'c3', 'ce'})

        Returns:
            Filtered lines
        """
        filtered = []
        for line in masses_lines:
            if '@atom:' in line:
                # Extract atom type
                parts = line.split()
                for part in parts:
                    if part.startswith('@atom:'):
                        atom_type = part.split(':')[1].split()[0]
                        if atom_type in atom_types:
                            filtered.append(line)
                            break
            else:
                # Keep header/footer lines, but NOT closing braces
                stripped = line.strip()
                if not (stripped.startswith('}') or (stripped.startswith('#') and 'end of' in stripped.lower())):
                    filtered.append(line)
        return filtered

    def _filter_pair_coeffs_section(self, pair_lines: List[str], atom_types: Set[str]) -> List[str]:
        """Filter In Settings (pair_coeff) section.

        Keeps pair_coeff for atom types in atom_types set.

        Args:
            pair_lines: Lines from pair_coeff In Settings section
            atom_types: Set of used atom types

        Returns:
            Filtered lines
        """
        filtered = []
        for line in pair_lines:
            if 'pair_coeff' in line and '@atom:' in line:
                # Extract atom type
                parts = line.split()
                for part in parts:
                    if part.startswith('@atom:'):
                        atom_type = part.split(':')[1]
                        if atom_type in atom_types:
                            filtered.append(line)
                            break
            else:
                # Keep non-pair_coeff lines but skip closing braces
                stripped = line.strip()
                if not (stripped.startswith('}') or (stripped.startswith('#') and 'end of' in stripped.lower())):
                    filtered.append(line)
        return filtered

    def _filter_bond_section(self,
                            bond_coeff_lines: List[str],
                            bond_def_lines: List[str],
                            used_bond_types: Set[str]) -> Tuple[List[str], List[str]]:
        """Filter bond coefficients and definitions.

        Keep bond if it's in the used_bond_types set.

        Args:
            bond_coeff_lines: Lines from bond_coeff In Settings section
            bond_def_lines: Lines from Data Bonds By Type section
            used_bond_types: Set of explicitly used bond type strings (e.g., "c2-c3", "ce-c2")

        Returns:
            Tuple of (filtered_coeff_lines, filtered_def_lines)
        """
        # First, identify which bonds to keep
        keep_bonds = set()

        for line in bond_def_lines:
            if '@bond:' in line and '@atom:' in line:
                # Extract atom types from definition
                parts = line.split()
                bond_types = []
                for part in parts:
                    if part.startswith('@atom:'):
                        atom_type = part.split(':')[1]
                        bond_types.append(atom_type)

                # Normalize bond type (alphabetically sort) and check if used
                if len(bond_types) == 2:
                    normalized_bond = '-'.join(sorted(bond_types))
                    if normalized_bond in used_bond_types:
                        # Extract bond name: @bond:c3-ce
                        bond_name = None
                        for part in parts:
                            if part.startswith('@bond:'):
                                bond_name = part.split(':')[1]
                                break
                        if bond_name:
                            keep_bonds.add(bond_name)

        # Filter coefficients
        filtered_coeffs = []
        for line in bond_coeff_lines:
            if '@bond:' in line:
                bond_name = None
                parts = line.split()
                for part in parts:
                    if part.startswith('@bond:'):
                        bond_name = part.split(':')[1]
                        break
                if bond_name and bond_name in keep_bonds:
                    filtered_coeffs.append(line)
            elif (not line.strip() or 'write_once' in line or
                  ('#' in line and not line.strip().startswith('#end of'))):
                # Keep empty lines, section headers, comments (but not closing braces)
                if not line.strip().startswith('}'):
                    filtered_coeffs.append(line)

        # Filter definitions
        filtered_defs = []
        for line in bond_def_lines:
            if '@bond:' in line:
                bond_name = None
                parts = line.split()
                for part in parts:
                    if part.startswith('@bond:'):
                        bond_name = part.split(':')[1]
                        break
                if bond_name and bond_name in keep_bonds:
                    filtered_defs.append(line)
            elif (not line.strip() or 'write_once' in line or
                  ('#' in line and not line.strip().startswith('#end of'))):
                # Keep empty lines, section headers, comments (but not closing braces)
                if not line.strip().startswith('}'):
                    filtered_defs.append(line)

        return filtered_coeffs, filtered_defs

    def _filter_angle_section(self,
                             angle_coeff_lines: List[str],
                             angle_def_lines: List[str],
                             used_angle_types: Set[str]) -> Tuple[List[str], List[str]]:
        """Filter angle coefficients and definitions.

        Keep angle if it's in the used_angle_types set.

        Args:
            angle_coeff_lines: Lines from angle_coeff In Settings section
            angle_def_lines: Lines from Data Angles By Type section
            used_angle_types: Set of explicitly used angle type strings (e.g., "c2-c3-c", "ce-c2-o")

        Returns:
            Tuple of (filtered_coeff_lines, filtered_def_lines)
        """
        # First, identify which angles to keep
        keep_angles = set()

        for line in angle_def_lines:
            if '@angle:' in line and '@atom:' in line:
                # Extract atom types from definition
                parts = line.split()
                angle_types = []
                for part in parts:
                    if part.startswith('@atom:'):
                        atom_type = part.split(':')[1]
                        angle_types.append(atom_type)

                # Normalize angle type (outer1-center-outer2 with outer atoms sorted)
                if len(angle_types) == 3:
                    # Middle atom is the center
                    center = angle_types[1]
                    outer1, outer2 = angle_types[0], angle_types[2]

                    # Sort outer atoms alphabetically
                    if outer1 > outer2:
                        outer1, outer2 = outer2, outer1

                    normalized_angle = f"{outer1}-{center}-{outer2}"
                    if normalized_angle in used_angle_types:
                        # Extract angle name: @angle:c3-ce-c3
                        angle_name = None
                        for part in parts:
                            if part.startswith('@angle:'):
                                angle_name = part.split(':')[1]
                                break
                        if angle_name:
                            keep_angles.add(angle_name)

        # Filter coefficients
        filtered_coeffs = []
        for line in angle_coeff_lines:
            if '@angle:' in line:
                angle_name = None
                parts = line.split()
                for part in parts:
                    if part.startswith('@angle:'):
                        angle_name = part.split(':')[1]
                        break
                if angle_name and angle_name in keep_angles:
                    filtered_coeffs.append(line)
            elif (not line.strip() or 'write_once' in line or
                  ('#' in line and not line.strip().startswith('#end of'))):
                # Keep empty lines, section headers, comments (but not closing braces)
                if not line.strip().startswith('}'):
                    filtered_coeffs.append(line)

        # Filter definitions
        filtered_defs = []
        for line in angle_def_lines:
            if '@angle:' in line:
                angle_name = None
                parts = line.split()
                for part in parts:
                    if part.startswith('@angle:'):
                        angle_name = part.split(':')[1]
                        break
                if angle_name and angle_name in keep_angles:
                    filtered_defs.append(line)
            elif (not line.strip() or 'write_once' in line or
                  ('#' in line and not line.strip().startswith('#end of'))):
                # Keep empty lines, section headers, comments (but not closing braces)
                if not line.strip().startswith('}'):
                    filtered_defs.append(line)

        return filtered_coeffs, filtered_defs

    def _filter_dihedral_section(self,
                                dihedral_coeff_lines: List[str],
                                dihedral_def_lines: List[str],
                                used_dihedral_types: Set[str]) -> Tuple[List[str], List[str]]:
        """Filter dihedral coefficients and definitions.

        Keep dihedral if it's in the used_dihedral_types set.

        Args:
            dihedral_coeff_lines: Lines from dihedral_coeff In Settings section
            dihedral_def_lines: Lines from Data Dihedrals By Type section
            used_dihedral_types: Set of explicitly used dihedral type strings (e.g., "hc-c3-c2-o")

        Returns:
            Tuple of (filtered_coeff_lines, filtered_def_lines)
        """
        # First, identify which dihedrals to keep
        keep_dihedrals = set()

        for line in dihedral_def_lines:
            if '@dihedral:' in line and '@atom:' in line:
                # Extract atom types from definition
                parts = line.split()
                dihedral_types = []
                for part in parts:
                    if part.startswith('@atom:'):
                        atom_type = part.split(':')[1]
                        dihedral_types.append(atom_type)

                # Check if this dihedral type is used (exact match, order matters)
                if len(dihedral_types) == 4:
                    dihedral_type_str = '-'.join(dihedral_types)
                    if dihedral_type_str in used_dihedral_types:
                        # Extract dihedral name: @dihedral:c3-ce-c3-h
                        dihedral_name = None
                        for part in parts:
                            if part.startswith('@dihedral:'):
                                dihedral_name = part.split(':')[1]
                                break
                        if dihedral_name:
                            keep_dihedrals.add(dihedral_name)

        # Filter coefficients
        filtered_coeffs = []
        for line in dihedral_coeff_lines:
            if '@dihedral:' in line:
                dihedral_name = None
                parts = line.split()
                for part in parts:
                    if part.startswith('@dihedral:'):
                        dihedral_name = part.split(':')[1]
                        break
                if dihedral_name and dihedral_name in keep_dihedrals:
                    filtered_coeffs.append(line)
            elif (not line.strip() or 'write_once' in line or
                  ('#' in line and not line.strip().startswith('#end of'))):
                # Keep empty lines, section headers, comments (but not closing braces)
                if not line.strip().startswith('}'):
                    filtered_coeffs.append(line)

        # Filter definitions
        filtered_defs = []
        for line in dihedral_def_lines:
            if '@dihedral:' in line:
                dihedral_name = None
                parts = line.split()
                for part in parts:
                    if part.startswith('@dihedral:'):
                        dihedral_name = part.split(':')[1]
                        break
                if dihedral_name and dihedral_name in keep_dihedrals:
                    filtered_defs.append(line)
            elif (not line.strip() or 'write_once' in line or
                  ('#' in line and not line.strip().startswith('#end of'))):
                # Keep empty lines, section headers, comments (but not closing braces)
                if not line.strip().startswith('}'):
                    filtered_defs.append(line)

        return filtered_coeffs, filtered_defs

    def _filter_improper_section(self,
                                improper_coeff_lines: List[str],
                                improper_def_lines: List[str],
                                used_improper_types: Set[str]) -> Tuple[List[str], List[str]]:
        """Filter improper coefficients and definitions.

        Keep improper if it's in the used_improper_types set.

        Args:
            improper_coeff_lines: Lines from improper_coeff In Settings section
            improper_def_lines: Lines from Data Impropers By Type section
            used_improper_types: Set of explicitly used improper type strings (e.g., "c-c2-c3-hc")

        Returns:
            Tuple of (filtered_coeff_lines, filtered_def_lines)
        """
        # First, identify which impropers to keep
        keep_impropers = set()

        for line in improper_def_lines:
            if '@improper:' in line and '@atom:' in line:
                # Extract atom types from definition
                parts = line.split()
                improper_types = []
                for part in parts:
                    if part.startswith('@atom:'):
                        atom_type = part.split(':')[1]
                        improper_types.append(atom_type)

                # Check if this improper type is used (exact match, order matters)
                if len(improper_types) == 4:
                    improper_type_str = '-'.join(improper_types)
                    if improper_type_str in used_improper_types:
                        # Extract improper name: @improper:X-c3-ce-c3
                        improper_name = None
                        for part in parts:
                            if part.startswith('@improper:'):
                                improper_name = part.split(':')[1]
                                break
                        if improper_name:
                            keep_impropers.add(improper_name)

        # Filter coefficients
        filtered_coeffs = []
        for line in improper_coeff_lines:
            if '@improper:' in line:
                improper_name = None
                parts = line.split()
                for part in parts:
                    if part.startswith('@improper:'):
                        improper_name = part.split(':')[1]
                        break
                if improper_name and improper_name in keep_impropers:
                    filtered_coeffs.append(line)
            elif (not line.strip() or 'write_once' in line or
                  ('#' in line and not line.strip().startswith('#end of'))):
                # Keep empty lines, section headers, comments (but not closing braces)
                if not line.strip().startswith('}'):
                    filtered_coeffs.append(line)

        # Filter definitions
        filtered_defs = []
        for line in improper_def_lines:
            if '@improper:' in line:
                improper_name = None
                parts = line.split()
                for part in parts:
                    if part.startswith('@improper:'):
                        improper_name = part.split(':')[1]
                        break
                if improper_name and improper_name in keep_impropers:
                    filtered_defs.append(line)
            elif (not line.strip() or 'write_once' in line or
                  ('#' in line and not line.strip().startswith('#end of'))):
                # Keep empty lines, section headers, comments (but not closing braces)
                if not line.strip().startswith('}'):
                    filtered_defs.append(line)

        return filtered_coeffs, filtered_defs

    def _extract_bonds_from_monomers(self, monomer_files: List[str]) -> Set[str]:
        """Extract bond type pairs from monomer .lt files.

        Parses the 'Data Bond List' sections to find which atom type pairs
        are actually bonded in the monomers.

        First builds an atom ID -> atom type mapping from 'Data Atoms' section,
        then extracts bonds and looks up the types.

        Args:
            monomer_files: List of paths to monomer .lt files

        Returns:
            Set of bond type strings like "c2-c3", "ce-c2", etc.
            Normalized alphabetically for consistency.
        """
        used_bonds = set()

        for filepath in monomer_files:
            try:
                # First, build atom ID -> atom type mapping
                atom_types_map = {}
                in_atoms_section = False
                with open(filepath, 'r') as f:
                    for line in f:
                        if 'write("Data Atoms")' in line or 'write(\'Data Atoms\')' in line:
                            in_atoms_section = True
                            continue
                        elif in_atoms_section and '}' in line and not line.strip().startswith('#'):
                            break

                        if in_atoms_section and '$atom:' in line and '@atom:' in line:
                            parts = line.split()
                            atom_id = None
                            atom_type = None
                            for part in parts:
                                if part.startswith('$atom:'):
                                    atom_id = part.split(':')[1]
                                elif part.startswith('@atom:'):
                                    atom_type = part.split(':')[1]
                            if atom_id and atom_type:
                                atom_types_map[atom_id] = atom_type

                # Now extract bonds using the atom type mapping
                in_bond_section = False
                with open(filepath, 'r') as f:
                    for line in f:
                        # Look for bond list section
                        if 'write(\'Data Bond List\')' in line or 'write("Data Bond List")' in line:
                            in_bond_section = True
                            continue

                        # End of section
                        if in_bond_section:
                            if '}' in line and not line.strip().startswith('#'):
                                break

                            # Extract bonded atom IDs
                            if '$bond:' in line and '$atom:' in line:
                                parts = line.split()
                                atom_ids = []
                                for part in parts:
                                    # Only extract $atom:, not $bond:
                                    if part.startswith('$atom:'):
                                        # Extract atom ID (after colon)
                                        atom_id = part.split(':')[1]
                                        atom_ids.append(atom_id)

                                # Look up atom types and create bond type
                                if len(atom_ids) >= 2:
                                    atom1_type = atom_types_map.get(atom_ids[0])
                                    atom2_type = atom_types_map.get(atom_ids[1])

                                    if atom1_type and atom2_type:
                                        # Normalize alphabetically for consistency
                                        bond_type = '-'.join(sorted([atom1_type, atom2_type]))
                                        used_bonds.add(bond_type)

            except Exception as e:
                logger.warning(f"  Warning: Could not parse bonds from {filepath}: {e}")
                continue

        return used_bonds

    def _build_bond_graph(self, filepath: str) -> Dict[str, Tuple[str, List[str]]]:
        """Build a bond connectivity graph from a monomer .lt file.

        Args:
            filepath: Path to monomer .lt file

        Returns:
            Dictionary mapping atom_id -> (atom_type, [neighbor_atom_ids])
        """
        graph = {}
        atom_types = {}  # atom_id -> atom_type

        try:
            with open(filepath, 'r') as f:
                in_atoms_section = False
                in_bond_section = False

                for line in f:
                    # Parse Data Atoms section
                    if 'write("Data Atoms")' in line or 'write(\'Data Atoms\')' in line:
                        in_atoms_section = True
                        continue
                    elif in_atoms_section and '}' in line and not line.strip().startswith('#'):
                        in_atoms_section = False

                    if in_atoms_section and '$atom:' in line and '@atom:' in line:
                        parts = line.split()
                        atom_id = None
                        atom_type = None
                        for part in parts:
                            if part.startswith('$atom:'):
                                atom_id = part.split(':')[1]
                            elif part.startswith('@atom:'):
                                atom_type = part.split(':')[1]
                        if atom_id and atom_type:
                            atom_types[atom_id] = atom_type
                            graph[atom_id] = (atom_type, [])

                    # Parse Data Bond List section
                    if 'write(\'Data Bond List\')' in line or 'write("Data Bond List")' in line:
                        in_bond_section = True
                        continue
                    elif in_bond_section and '}' in line and not line.strip().startswith('#'):
                        in_bond_section = False

                    if in_bond_section and '$bond:' in line:
                        parts = line.split()
                        bonded_atoms = []
                        for part in parts:
                            # Only extract $atom:, not $bond:
                            if part.startswith('$atom:'):
                                # Extract atom ID (after colon)
                                atom_id = part.split(':')[1]
                                bonded_atoms.append(atom_id)

                        # Add edges to graph (bidirectional)
                        if len(bonded_atoms) >= 2:
                            atom1, atom2 = bonded_atoms[0], bonded_atoms[1]
                            if atom1 in graph and atom2 in graph:
                                graph[atom1][1].append(atom2)
                                graph[atom2][1].append(atom1)

        except Exception as e:
            logger.warning(f"  Warning: Could not build bond graph from {filepath}: {e}")

        return graph

    def _infer_angles_from_graph(self, bond_graph: Dict) -> Set[str]:
        """Infer angle types from bond connectivity graph.

        For each atom with 2+ neighbors, generates all angle combinations
        with that atom as the center.

        Args:
            bond_graph: Dict mapping atom_id -> (atom_type, [neighbor_ids])

        Returns:
            Set of angle type strings like "c2-c3-c", "ce-c2-o", etc.
        """
        angles = set()

        for center_id, (center_type, neighbors) in bond_graph.items():
            if len(neighbors) >= 2:
                # Generate all angle combinations with this central atom
                for i in range(len(neighbors)):
                    for j in range(i+1, len(neighbors)):
                        atom1_type = bond_graph[neighbors[i]][0]
                        atom2_type = bond_graph[neighbors[j]][0]

                        # Create angle type string (atom1-center-atom2)
                        # Normalize: put smaller (alphabetically) atom type first
                        if atom1_type <= atom2_type:
                            angle_type = f"{atom1_type}-{center_type}-{atom2_type}"
                        else:
                            angle_type = f"{atom2_type}-{center_type}-{atom1_type}"
                        angles.add(angle_type)

        return angles

    def _infer_dihedrals_from_graph(self, bond_graph: Dict) -> Set[str]:
        """Infer dihedral types from bond connectivity graph.

        For each bond as the central bond, finds all atoms bonded to each end
        and generates dihedral combinations.

        Args:
            bond_graph: Dict mapping atom_id -> (atom_type, [neighbor_ids])

        Returns:
            Set of dihedral type strings like "hc-c3-c2-o", etc.
        """
        dihedrals = set()

        # For each bond as central bond
        for atom1_id, (atom1_type, neighbors1) in bond_graph.items():
            for atom2_id in neighbors1:
                atom2_type = bond_graph[atom2_id][0]

                # Find atoms bonded to atom1 (excluding atom2)
                outer_atoms_1 = [n for n in neighbors1 if n != atom2_id]

                # Find atoms bonded to atom2 (excluding atom1)
                neighbors2 = bond_graph[atom2_id][1]
                outer_atoms_2 = [n for n in neighbors2 if n != atom1_id]

                # Generate all dihedrals: outer1-atom1-atom2-outer2
                for outer1 in outer_atoms_1:
                    outer1_type = bond_graph[outer1][0]
                    for outer2 in outer_atoms_2:
                        outer2_type = bond_graph[outer2][0]
                        dihedral_type = f"{outer1_type}-{atom1_type}-{atom2_type}-{outer2_type}"
                        dihedrals.add(dihedral_type)

        return dihedrals

    def _infer_impropers_from_graph(self, bond_graph: Dict) -> Set[str]:
        """Infer improper types from bond connectivity graph.

        For each atom with 3+ neighbors (central atom of improper),
        generates all improper combinations.

        Args:
            bond_graph: Dict mapping atom_id -> (atom_type, [neighbor_ids])

        Returns:
            Set of improper type strings like "c-c2-c3-hc", etc.
        """
        impropers = set()

        # For each atom with 3+ neighbors as central atom
        for center_id, (center_type, neighbors) in bond_graph.items():
            if len(neighbors) >= 3:
                # Generate all improper combinations (3 neighbors at a time)
                for i in range(len(neighbors)):
                    for j in range(i+1, len(neighbors)):
                        for k in range(j+1, len(neighbors)):
                            # Get the three connected atom types
                            atom1_type = bond_graph[neighbors[i]][0]
                            atom2_type = bond_graph[neighbors[j]][0]
                            atom3_type = bond_graph[neighbors[k]][0]

                            # Create improper type: atom1-atom2-atom3-center
                            # (following GAFF convention where central atom is last)
                            improper_type = f"{atom1_type}-{atom2_type}-{atom3_type}-{center_type}"
                            impropers.add(improper_type)

        return impropers

    def make_gaff_subset(self) -> None:
        """Create a subset of GAFF parameters based on the models.

        This method analyzes the topology used in the monomer files and
        creates a filtered version of gaff.lt containing only the parameters
        relevant to those specific bond/angle/dihedral/improper types.

        Process:
            1. Extract monomer file paths from model
            2. Extract actual bond types from monomer files
            3. Build bond graphs and infer angle/dihedral/improper types
            4. Parse the full gaff.lt file into sections
            5. Filter each section to keep only used parameter types
            6. Write gaff_subset.lt to working directory and symlink to gaff.lt

        Falls back to full gaff.lt if no parameters found or errors occur.
        """
        try:
            gaff_src = f"{self.path_master}moltemplate/common/gaff.lt"
            gaff_dst = self.path_cwd + "gaff_subset.lt"

            # Check if source file exists
            if not Path(gaff_src).exists():
                logger.error(f"GAFF force field file not found: {gaff_src}")
                logger.error("Please ensure gaff.lt is installed in moltemplate/common/")
                sys.exit(1)

            logger.info("Creating GAFF parameter subset...")

            # Step 1: Collect monomer file paths
            logger.info("  Collecting monomer files...")
            monomer_files = []
            for modelii in self.model:
                for monomerii in range(len(modelii.sequenceSet)):
                    monomerSet = modelii.sequenceSet[monomerii]

                    for vecii in range(len(monomerSet)):
                        # Determine monomer bank path based on force field
                        MonomerBank = Path(self.path_cwd)
                        merltfile_Path = MonomerBank / monomerSet[vecii]

                        if merltfile_Path.is_file():
                            monomer_files.append(str(merltfile_Path))

            if not monomer_files:
                logger.warning("  No monomer files found!")
                logger.warning("  Falling back to full gaff.lt")
                shutil.copy(gaff_src, self.path_cwd + "gaff.lt")
                return

            logger.info(f"  Found {len(monomer_files)} monomer files")

            # Step 2: Extract atom types (still needed for masses and pair coeffs)
            logger.info("  Extracting atom types from monomers...")
            atom_types = self._extract_gaff_atom_types_from_monomers()
            logger.info(f"  Found {len(atom_types)} unique atom types: {sorted(atom_types)}")

            if not atom_types:
                logger.warning("  No atom types found in monomers!")
                logger.warning("  Falling back to full gaff.lt")
                shutil.copy(gaff_src, self.path_cwd + "gaff.lt")
                return

            # Step 3: Extract bond types from monomers and add inter-monomer possibilities
            logger.info("  Extracting bond types from monomers...")
            used_bond_types = self._extract_bonds_from_monomers(monomer_files)
            logger.info(f"  Found {len(used_bond_types)} unique intra-monomer bond types: {sorted(used_bond_types)}")

            # For bonds, include all possible atom type combinations to handle inter-monomer connections
            # This is conservative but ensures we don't miss bonds formed during polymerization
            all_atom_type_pairs = set()
            atom_list = sorted(atom_types)
            for i in range(len(atom_list)):
                for j in range(i, len(atom_list)):
                    all_atom_type_pairs.add(f"{atom_list[i]}-{atom_list[j]}")

            # Merge: explicit intra-monomer bonds + all possible combinations for inter-monomer
            used_bond_types.update(all_atom_type_pairs)
            logger.info(f"  Total bond types (including inter-monomer possibilities): {len(used_bond_types)}")

            # Step 4: Build bond graphs and infer angles/dihedrals/impropers
            logger.info("  Inferring angles, dihedrals, and impropers from bond connectivity...")
            used_angle_types = set()
            used_dihedral_types = set()
            used_improper_types = set()

            for monofile in monomer_files:
                try:
                    graph = self._build_bond_graph(monofile)
                    if graph:
                        # Infer angles
                        angles = self._infer_angles_from_graph(graph)
                        used_angle_types.update(angles)

                        # Infer dihedrals
                        dihedrals = self._infer_dihedrals_from_graph(graph)
                        used_dihedral_types.update(dihedrals)

                        # Infer impropers
                        impropers = self._infer_impropers_from_graph(graph)
                        used_improper_types.update(impropers)
                except Exception as e:
                    logger.warning(f"    Warning: Could not process {monofile}: {e}")
                    continue

            logger.info(f"  Found {len(used_angle_types)} unique angle types")
            logger.info(f"  Found {len(used_dihedral_types)} unique dihedral types")
            logger.info(f"  Found {len(used_improper_types)} unique improper types")

            # Step 5: Parse gaff.lt
            logger.info("  Parsing gaff.lt structure...")
            sections = self._parse_gaff_lt_sections(gaff_src)

            # Step 6: Filter each section
            logger.info("  Filtering parameter sections...")

            # Filter masses (still needs atom_types)
            filtered_masses = self._filter_masses_section(sections['masses'], atom_types)
            logger.info(f"    Masses: {len([l for l in filtered_masses if '@atom:' in l])} entries "
                       f"(from {len([l for l in sections['masses'] if '@atom:' in l])})")

            # Filter pair coefficients (still needs atom_types)
            filtered_pairs = self._filter_pair_coeffs_section(sections['pair_coeffs'], atom_types)
            logger.info(f"    Pair coeffs: {len([l for l in filtered_pairs if 'pair_coeff' in l])} entries")

            # Filter bonds (now uses explicit bond types)
            filtered_bond_coeffs, filtered_bond_defs = self._filter_bond_section(
                sections['bond_coeffs'], sections['bond_definitions'], used_bond_types
            )
            logger.info(f"    Bonds: {len([l for l in filtered_bond_coeffs if 'bond_coeff' in l])} coeffs, "
                       f"{len([l for l in filtered_bond_defs if '@bond:' in l])} defs")

            # Filter angles (now uses explicit angle types)
            filtered_angle_coeffs, filtered_angle_defs = self._filter_angle_section(
                sections['angle_coeffs'], sections['angle_definitions'], used_angle_types
            )
            logger.info(f"    Angles: {len([l for l in filtered_angle_coeffs if 'angle_coeff' in l])} coeffs, "
                       f"{len([l for l in filtered_angle_defs if '@angle:' in l])} defs")

            # Filter dihedrals (now uses explicit dihedral types)
            filtered_dihedral_coeffs, filtered_dihedral_defs = self._filter_dihedral_section(
                sections['dihedral_coeffs'], sections['dihedral_definitions'], used_dihedral_types
            )
            logger.info(f"    Dihedrals: {len([l for l in filtered_dihedral_coeffs if 'dihedral_coeff' in l])} coeffs, "
                       f"{len([l for l in filtered_dihedral_defs if '@dihedral:' in l])} defs")

            # Filter impropers (now uses explicit improper types)
            filtered_improper_coeffs, filtered_improper_defs = self._filter_improper_section(
                sections['improper_coeffs'], sections['improper_definitions'], used_improper_types
            )
            logger.info(f"    Impropers: {len([l for l in filtered_improper_coeffs if 'improper_coeff' in l])} coeffs, "
                       f"{len([l for l in filtered_improper_defs if '@improper:' in l])} defs")

            # Step 4: Write gaff_subset.lt
            logger.info(f"  Writing {gaff_dst}...")

            with open(gaff_dst, 'w') as f:
                # Header
                f.write("# GAFF Force Field Subset\n")
                f.write(f"# Generated by AutoPoly from {gaff_src}\n")
                f.write(f"# Atom types used: {', '.join(sorted(atom_types))}\n")
                f.write(f"# Total atom types: {len(atom_types)}\n")
                f.write("\n")
                f.write("GAFF {\n\n")

                # Masses
                f.write("  write_once(\"Data Masses\") {\n")
                for line in filtered_masses:
                    # Skip section headers, closing braces, and comments that are section markers
                    stripped = line.strip()
                    # Skip lines that are closing braces (with or without comments)
                    if stripped.startswith('}'):
                        continue
                    # Skip section headers
                    if stripped.startswith('write_once('):
                        continue
                    # Skip comment lines that reference "once" or "end"
                    if stripped.startswith('#') and ('once' in stripped or 'end of' in stripped.lower()):
                        continue
                    f.write(f"{line}")
                f.write("  } # (end of masses)\n\n")

                # Pair coeffs
                f.write("  write_once(\"In Settings\") {\n")
                for line in filtered_pairs:
                    # Skip section headers, closing braces, and comments that are section markers
                    stripped = line.strip()
                    # Skip lines that are closing braces (with or without comments)
                    if stripped.startswith('}'):
                        continue
                    # Skip section headers
                    if stripped.startswith('write_once('):
                        continue
                    # Skip comment lines that reference "once" or "end"
                    if stripped.startswith('#') and ('once' in stripped or 'end of' in stripped.lower()):
                        continue
                    f.write(f"{line}")
                f.write("  } # (end of pair_coeffs)\n\n")

                # Bond definitions and coeffs
                if filtered_bond_defs:
                    f.write("  write_once(\"Data Bonds By Type\") {\n")
                    for line in filtered_bond_defs:
                        # Skip section headers, closing braces, and comments that are section markers
                        stripped = line.strip()
                        if not (stripped.startswith('write_once(') or
                                stripped == '}' or
                                stripped.startswith('#end of') or
                                (stripped.startswith('#') and 'once' in stripped)):
                            f.write(f"{line}")
                    f.write("  } # (end of bonds by type)\n\n")

                if filtered_bond_coeffs:
                    f.write("  write_once(\"In Settings\") {\n")
                    for line in filtered_bond_coeffs:
                        # Skip section headers, closing braces, and comments that are section markers
                        stripped = line.strip()
                        if not (stripped.startswith('write_once(') or
                                stripped == '}' or
                                stripped.startswith('#end of') or
                                (stripped.startswith('#') and 'once' in stripped)):
                            f.write(f"{line}")
                    f.write("  } # (end of bond_coeffs)\n\n")

                # Angle definitions and coeffs
                if filtered_angle_defs:
                    f.write("  write_once(\"Data Angles By Type\") {\n")
                    for line in filtered_angle_defs:
                        # Skip section headers, closing braces, and comments that are section markers
                        stripped = line.strip()
                        if not (stripped.startswith('write_once(') or
                                stripped == '}' or
                                stripped.startswith('#end of') or
                                (stripped.startswith('#') and 'once' in stripped)):
                            f.write(f"{line}")
                    f.write("  } # (end of angles by type)\n\n")

                if filtered_angle_coeffs:
                    f.write("  write_once(\"In Settings\") {\n")
                    for line in filtered_angle_coeffs:
                        # Skip section headers, closing braces, and comments that are section markers
                        stripped = line.strip()
                        if not (stripped.startswith('write_once(') or
                                stripped == '}' or
                                stripped.startswith('#end of') or
                                (stripped.startswith('#') and 'once' in stripped)):
                            f.write(f"{line}")
                    f.write("  } # (end of angle_coeffs)\n\n")

                # Dihedral definitions and coeffs
                if filtered_dihedral_defs:
                    f.write("  write_once(\"Data Dihedrals By Type\") {\n")
                    for line in filtered_dihedral_defs:
                        # Skip section headers, closing braces, and comments that are section markers
                        stripped = line.strip()
                        if not (stripped.startswith('write_once(') or
                                stripped == '}' or
                                stripped.startswith('#end of') or
                                (stripped.startswith('#') and 'once' in stripped)):
                            f.write(f"{line}")
                    f.write("  } # (end of Dihedrals by type)\n\n")

                if filtered_dihedral_coeffs:
                    f.write("  write_once(\"In Settings\") {\n")
                    for line in filtered_dihedral_coeffs:
                        # Skip section headers, closing braces, and comments that are section markers
                        stripped = line.strip()
                        if not (stripped.startswith('write_once(') or
                                stripped == '}' or
                                stripped.startswith('#end of') or
                                (stripped.startswith('#') and 'once' in stripped)):
                            f.write(f"{line}")
                    f.write("  } # (end of dihedral_coeffs)\n\n")

                # Improper definitions and coeffs
                if filtered_improper_defs:
                    f.write("  write_once(\"Data Impropers By Type (gaff_imp.py)\") {\n")
                    for line in filtered_improper_defs:
                        # Skip section headers, closing braces, and comments that are section markers
                        stripped = line.strip()
                        if not (stripped.startswith('write_once(') or
                                stripped == '}' or
                                stripped.startswith('#end of') or
                                (stripped.startswith('#') and 'once' in stripped)):
                            f.write(f"{line}")
                    f.write("  } # (end of impropers by type)\n\n")

                if filtered_improper_coeffs:
                    f.write("  write_once(\"In Settings\") {\n")
                    for line in filtered_improper_coeffs:
                        # Skip section headers, closing braces, and comments that are section markers
                        stripped = line.strip()
                        if not (stripped.startswith('write_once(') or
                                stripped == '}' or
                                stripped.startswith('#end of') or
                                (stripped.startswith('#') and 'once' in stripped)):
                            f.write(f"{line}")
                    f.write("  } # (end of improp_coeffs)\n\n")

                # Init section (unchanged)
                if sections['init']:
                    f.write("  write_once(\"In Init\") {\n")
                    for line in sections['init']:
                        # Skip section headers, closing braces, and comments that are section markers
                        stripped = line.strip()
                        if not (stripped.startswith('write_once(') or
                                stripped == '}' or
                                stripped.startswith('#end of') or
                                (stripped.startswith('#') and 'once' in stripped)):
                            f.write(f"{line}")
                    f.write("  } #end of init parameters\n\n")

                f.write("} # GAFF\n")

            # Calculate file size reduction
            original_size = Path(gaff_src).stat().st_size
            subset_size = Path(gaff_dst).stat().st_size
            reduction = (1 - subset_size / original_size) * 100

            logger.info(f"  Successfully created gaff_subset.lt")
            logger.info(f"  File size: {original_size} -> {subset_size} bytes ({reduction:.1f}% reduction)")

            # Create symlink from gaff.lt to gaff_subset.lt for monomer compatibility
            gaff_link = self.path_cwd + "gaff.lt"
            if Path(gaff_link).exists():
                Path(gaff_link).unlink()
            Path(gaff_link).symlink_to("gaff_subset.lt")
            logger.info(f"  Created symlink: gaff.lt -> gaff_subset.lt")

        except Exception as e:
            logger.error(f"Error in make_gaff_subset: {str(e)}")
            import traceback
            traceback.print_exc()
            logger.warning("Falling back to full gaff.lt")
            shutil.copy(gaff_src, self.path_cwd + "gaff.lt")

    def make_gaff_lt(self) -> None:
        """Copy gaff.lt to working directory."""
        try:
            gaff_src = f"{self.path_master}moltemplate/common/gaff.lt"
            gaff_dst = self.path_cwd + "gaff.lt"

            # Check if source file exists
            if not Path(gaff_src).exists():
                logger.error(f"GAFF force field file not found: {gaff_src}")
                logger.error("Please ensure gaff.lt is installed in moltemplate/common/")
                sys.exit(1)

            shutil.copy(gaff_src, gaff_dst)
            logger.info(f"Copied gaff.lt to {gaff_dst}")

        except Exception as e:
            logger.error(f"Error in make_gaff_lt: {str(e)}")
            sys.exit(1)

    def FFmodify_alkyl_dihedral_oplsaa(self) -> None:
        """Modifies alkyl dihedral coefficients in the OPLSAA force field."""
        input_file = Path(self.path_cwd) / "oplsaa.lt"
        output_file = Path(self.path_cwd) / "oplsaa_tmp.lt"

        try:
            if not input_file.exists():
                logger.error("oplsaa.lt file cannot open.")
                sys.exit(1)
            logger.info(f"Start modifying alkyl dihedral coefficients")
            with open(input_file, 'r') as read_f, open(output_file, 'w') as write_f:
                is_inside_block = False
                
                for line in read_f:
                    line = line.strip()
                    
                    if not is_inside_block:
                        write_f.write(f"{line}\n")
                    
                    words = line.split()
                    if not words:
                        continue
                    
                    
                    if words[0] == 'write_once("In':
                        next_line = next(read_f).strip()
                        next_words = next_line.split()
                        
                        if next_words and next_words[0] == "dihedral_coeff":
                            is_inside_block = True
                            string_strip=next_line.strip()
                            if not string_strip.startswith('dihedral_coeff @dihedral:'):
                                write_f.write(f"{next_line}\n")
                                continue
                                
                            # Extract atom types from dihedral specification
                            dihedral_spec = next_line.split("@dihedral:")[1].split()[0]
                            atom_types = [int(x) for x in dihedral_spec.split('-')]

                            # Check if all atoms are CH3(80), CH2(81), or CH(82)
                            alkyl_atoms = {80, 81, 82}  # CH3, CH2, CH atoms
                            if all(atom in alkyl_atoms for atom in atom_types):
                                # Write modified dihedral coefficients
                                write_f.write("dihedral_coeff @dihedral:")
                                write_f.write('-'.join(str(x) for x in atom_types))
                                write_f.write(" opls")
                                
                                # Write new dihedral coefficients
                                # Note: FFmodify_alkylDihedral should be defined as a class attribute
                                if hasattr(self, 'FFmodify_alkylDihedral'):
                                    write_f.write(" " + " ".join(str(x) for x in self.FFmodify_alkylDihedral))
                                write_f.write("\n")
                            else:
                                write_f.write(f"{next_line}\n")
                            continue
                        else:
                            write_f.write(f"{next_line}\n")
                    
                    elif words[0] == "}":
                        if is_inside_block:
                            
                            write_f.write(f"{line}\n")
                        is_inside_block = False
                    
                    if is_inside_block:
                        # Parse dihedral specification
                        string_strip=line.strip()
                        if not string_strip.startswith('dihedral_coeff @dihedral:'):
                            
                            write_f.write(f"{line}\n")
                            continue
                            
                        # Extract atom types from dihedral specification
                        dihedral_spec = line.split("@dihedral:")[1].split()[0]
                        atom_types = [int(x) for x in dihedral_spec.split('-')]

                        # Check if all atoms are CH3(80), CH2(81), or CH(82)
                        alkyl_atoms = {80, 81, 82}  # CH3, CH2, CH atoms
                        if all(atom in alkyl_atoms for atom in atom_types):
                            # Write modified dihedral coefficients
                            write_f.write("dihedral_coeff @dihedral:")
                            write_f.write('-'.join(str(x) for x in atom_types))
                            write_f.write(" opls")
                            
                            # Write new dihedral coefficients
                            # Note: FFmodify_alkylDihedral should be defined as a class attribute
                            if hasattr(self, 'FFmodify_alkylDihedral'):
                                write_f.write(" " + " ".join(str(x) for x in self.FFmodify_alkylDihedral))
                            write_f.write("\n")
                        else:
                            write_f.write(f"{line}\n")

            # Replace original file with modified version
            shutil.move(str(output_file), str(input_file))

        except Exception as e:
            logger.error(f"Error modifying alkyl dihedral coefficients: {str(e)}")
            sys.exit(1)

    def make_oplsaa_subset(self) -> None:
        """Creates a subset of the oplsaa parameters based on the models."""
        # path to oplsaa_subset.prm file
        opls_subset_file = self.path_cwd+"oplsaa_subset.prm"

        atom_keys=[]
        for modelii in self.model:
            for monomerii in range(len(modelii.sequenceSet)):
                monomerSet=modelii.sequenceSet[monomerii]
                # vector to store all atom types including the repeats
                for vecii in range(len(monomerSet)):
                    # path to monomer.lt in monomer bank
                    # Determine monomer bank path based on force field
                    MonomerBank = Path(self.path_cwd)
                    merltfile_Path = MonomerBank / monomerSet[vecii]
                    if merltfile_Path.is_file():
                        mono = str(MonomerBank / monomerSet[vecii])
                        read_switch= False 
                        with open(mono) as f:
                            while True: 
                                line = f.readline() 
                                
                                if line.strip()=="write(\"Data Atoms\") {":
                                    read_switch=True
                                    continue
                                elif line.strip()=="}":
                                    read_switch=False
                                    break
                                
                                # Determine atom types, element names and the raw_charges as
                                # given in the opls table
                                
                                if read_switch:
                                    load_line=""
                                    stringvector=line.split()
                                    
                                    load_switch=False
                                    for readii in range(len(stringvector[2])):
                                        if stringvector[2][readii]==":":
                                            load_switch = True
                                            continue
                                        if load_switch:
                                            load_line += stringvector[2][readii]
                                    atom_keys.append(load_line)
                                
                                if not line: 
                                    break
                    else:
                        logger.error(' '.join(["Monomer ("+ monomerSet[vecii] + ") does NOT exist. \n",
                                                "Please check the following path to the file\n" + merltfile_Path + "\n"]))
                        sys.exit()
                        
        # Cleaning up the stored data. Remove duplicate atoms types
        atom_types=list(dict.fromkeys(atom_keys))
        #print(atom_types)
        # Convert the vectors string to vector int in order to sort the atom_types in ascending order
        atom_types=sorted([int(i) for i in atom_types])
        # Read the master opls file and store the ones that match the atom_types into new subset file
        write_f=open(opls_subset_file, "w")
        
        with open(self.path_oplsaaprm,'r') as read_f:
            path_oplsaaprm=Path(self.path_oplsaaprm)
            if path_oplsaaprm.is_file():
                check_switch=False
                while True:
                    prm_line = read_f.readline()
                    if len(prm_line.strip())!=0:
                       
                        if prm_line.strip() == "##  Atom Type Definitions  ##":
                            check_switch = True
                            write_f.write(prm_line+"\n")
                            prm_line = read_f.readline()
                            write_f.write(prm_line+"\n")
                            prm_line = read_f.readline()
                            write_f.write(prm_line+"\n")
                            continue
                        elif prm_line.strip()=="################################":
                            check_switch = False
                            write_f.write(prm_line+"\n")
                            continue
                        elif check_switch:
                            
                            stringvector=prm_line.split()
                            
                            
                            for checkii in range(len(atom_types)):
                                
                                if atom_types[checkii]==int(stringvector[1]):
                                    write_f.write(prm_line+"\n")
                                    break
                        else:
                            write_f.write(prm_line+"\n")
                    else:
                        write_f.write(prm_line+"\n")

                    if not prm_line:
                        break
        write_f.close()

    def generate_monomer_from_psmiles(self, psmiles: str) -> str:
        """
        Generate all 6 .lt files for a monomer from pSMILES or SMILES string.

        This method creates a complete set of monomer variant files (internal,
        left-end, right-end, with and without T1 chirality) from a pSMILES (for polymers)
        or SMILES (for molecules) string.

        Args:
            psmiles: pSMILES string (e.g., "[*]C=C[*]") or SMILES string (e.g., "CC(=O)C")

        Returns:
            str: Generated monomer base name (e.g., "monomer_0")

        Raises:
            SystemExit: If monomer generation fails
        """
        # Check cache (per-system caching)
        if psmiles in self._generated_smiles:
            logger.info(f"Using cached monomer for pSMILES/SMILES: {psmiles}")
            return self._generated_smiles[psmiles]

        # Generate unique monomer name
        monomer_name = f"monomer_{self._smiles_to_name_counter}"
        self._smiles_to_name_counter += 1

        try:
            # Preprocess pSMILES: strip wildcard atoms [*] to get core SMILES
            # pSMILES: [*]C=C[*] -> SMILES: C=C
            core_smiles = psmiles.replace("[*]", "").replace("*", "")

            # Create MonomerGenerator
            generator = MonomerGenerator(
                base_name=monomer_name,
                output_dir=self.path_cwd,
                is_gaff=(self.force_field == "gaff"),
                is_lopls=(self.force_field == "lopls"),
                verbose=False
            )

            # Generate variants from core SMILES
            # MonomerGenerator will handle capping the bonding sites with H atoms
            variants = generator.generate_variants(smiles=core_smiles)

            # Write .lt files to system folder
            generator.generate_lt_files(variants, generate_t1=True)

            # Cache the mapping
            self._generated_smiles[psmiles] = monomer_name

            logger.info(f"Generated monomer '{monomer_name}' from pSMILES/SMILES: {psmiles} (core: {core_smiles})")
            return monomer_name

        except Exception as e:
            logger.error(f"Failed to generate monomer from pSMILES/SMILES '{psmiles}': {e}")
            sys.exit(1)

    def create_ring_polymer_topology(self, poly_index: int, monomer_set: list) -> None:
        """
        Create a ring polymer topology and coordinates based on the provided monomer sequence.
        
        This method generates a Moltemplate .lt file for a ring polymer by:
        1. Importing required monomer templates and force field parameters
        2. Creating a circular arrangement of monomers
        3. Defining bonds between adjacent monomers to form a closed ring
        4. Calculating proper coordinates and rotations for each monomer
        
        The ring is created by placing monomers in a circle with appropriate
        spacing and rotation to ensure proper bonding geometry.
        
        Args:
            poly_index (int): The index of the polymer for file naming.
            monomer_set (list): The list of monomers in the ring polymer sequence.
        """
        output = self.path_cwd + f"/poly_{poly_index+1}.lt"

        with open(output, "w") as write_f:
            write_f.write("import \"oplsaa.lt\"\n")

            # Import unique monomers
            unique_monomers = list(dict.fromkeys(monomer_set))
            for monomer in unique_monomers:
                write_f.write(f"import \"{monomer}.lt\"\n")

            write_f.write("\n")

            # Define combined ring polymer
            write_f.write(f"poly_{poly_index+1} inherits OPLSAA {{\n\n")
            write_f.write("    create_var {$mol}\n\n")

            offset_cum = 0
            radius = self.offset * len(monomer_set) / (2 * np.pi)  # Calculate radius based on number of monomers
            
            # Place monomers in a circular arrangement
            for i in range(len(monomer_set)):
                angle = 2 * np.pi * i / len(monomer_set)  # Angle for current monomer
                x = radius * np.cos(angle)
                y = radius * np.sin(angle)
                
                # Calculate rotation to point each monomer towards the center
                rotation_angle = (angle * 180 / np.pi) + 90  # Convert to degrees and add 90° offset
                
                write_f.write(f"    monomer[{i}] = new {monomer_set[i]}")
                write_f.write(f".rot({rotation_angle},0,0,1)")  # Rotate around z-axis
                write_f.write(f".move({x:.4f},{y:.4f},0)")
                write_f.write("\n")

            # Add bonds between monomers to form the ring
            write_f.write("\n    write('Data Bond List') {\n")
            
            # Connect sequential monomers
            for i in range(len(monomer_set)-1):
                write_f.write(f"      $bond:b{i+1}  $atom:monomer[{i}]/C2  $atom:monomer[{i+1}]/C1\n")
            
            # Connect last monomer to first to close the ring
            write_f.write(f"      $bond:b{len(monomer_set)}  $atom:monomer[{len(monomer_set)-1}]/C2  $atom:monomer[0]/C1\n")
            
            write_f.write("    }\n")
            write_f.write(f"\n}} # poly_{poly_index+1}\n")