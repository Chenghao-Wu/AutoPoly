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
from typing import List, Optional, Dict, Any
from .system import logger

class Polymerization(object):
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
        path_monomer_bank (str): Path to monomer template bank
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
                 run: bool = True, path_monomer_bank: str = None, is_lopls: bool = False) -> None:
        """
        Initialize the Polymerization class.

        Args:
            name (str, optional): Name of the polymerization project. Defaults to None.
            system (object, optional): System object containing folder path. Defaults to None.
            model (list, optional): List of models for polymerization. Defaults to None.
            run (bool, optional): Flag to run the process immediately. Defaults to True.
            path_monomer_bank (str, optional): Path to the monomer bank. Defaults to None.
            is_lopls (bool, optional): Whether to use LOPLS force field. Defaults to False.
        
        Raises:
            SystemExit: If required directories or files are not found
        """
        self.name = name
        self.system = system
        self.path_cwd = f"{self.system.get_FolderPath}/{self.name}/moltemplate/"
        self.path_master = f"{Path(__file__).parent.resolve()}/extern/"
        self.path_moltemplatesrc = f"{self.path_master}moltemplate/src/"
        
        self.path_monomer_bank = (self.path_master + "Monomer_bank/"
                                  if path_monomer_bank is None else path_monomer_bank)
        self.is_lopls = is_lopls
        if is_lopls:
            self.path_oplsaaprm = f"{self.path_master}moltemplate/loplsaa.prm"
        else:
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
        base_path = Path(self.system.get_FolderPath)
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
        monomer_bank = Path(self.path_monomer_bank)
        merltfile_path = monomer_bank / merltfile
        
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
                
                # Loop through all polymers and make corresponding polymer lt files
                for moleii in range(len(modelii.sequenceSet)):
                    # Check degrees of polymerization (DOP) of current polymer
                    if modelii.DOP > 0:
                        if len(modelii.sequenceSet[moleii]) != modelii.DOP:
                            logger.warning(f"Warning: At molecule# {moleii} DOP={len(modelii.sequenceSet[moleii])} != {modelii.DOP}")
                    else:
                        logger.error(f"Warning: At molecule#{moleii+1}, DOP={len(modelii.sequenceSet[moleii])} {modelii.DOP}")

                    # Check monomer.lt's in the monomer bank
                    for indexii in range(len(modelii.sequenceSet[moleii])):
                        monomer_file = modelii.sequenceSet[moleii][indexii]
                        if not monomer_file.endswith('.lt'):
                            monomer_file += '.lt'
                        if self.check_monomer_bank(monomer_file):
                            source = Path(self.path_monomer_bank) / monomer_file
                            self.copy_to_cwd(source)
                        else:
                            logger.error(f"Error: Cannot find monomer file {monomer_file} in monomer bank at {self.path_monomer_bank}")
                            sys.exit(1)

                    if modelii.DOP > 1:
                        # Make poly.lt file
                        logger.info(f"Creating poly_{poly_index+1}.lt")
                        self.make_poly_lt(poly_index, modelii.sequenceSet[moleii], modelii)
                        poly_index += 1

            # Generate force field files
            logger.info("Generating oplsaa.lt")
            self.make_oplsaalt()

            # Generate system.lt file
            logger.info("Creating system.lt")
            self.make_system_lt()

            # Modify alkyl dihedral coefficients if needed
            if self.is_lopls:
                self.FFmodify_alkyl_dihedral_oplsaa()

            # Invoke moltemplate to generate LAMMPS datafile
            logger.info("Running moltemplate")
            self.invoke_moltemplate()

            # Check if the required files exist before proceeding
            required_files = ['system.in.settings', 'system.data', 'system.in.charges']
            missing_files = []
            for file in required_files:
                if not (Path(self.path_cwd) / file).exists():
                    missing_files.append(file)
            
            if missing_files:
                logger.error(f"Moltemplate failed to generate required files: {', '.join(missing_files)}")
                logger.error("Check the following:")
                logger.error("1. All monomer .lt files exist and are valid")
                logger.error("2. The polymer .lt files were generated correctly")
                logger.error("3. The system.lt file is properly formatted")
                sys.exit(1)

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
        write_f=open(out, "w")
        with open(in_,'r') as read_f:
            
            in_path=Path(in_)
            if in_path.is_file():
                
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
            else:
                logger.error(' '.join(["system.in.setting does not exist plase check ",in_]))
                sys.exit()
        write_f.close()
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
            process = subprocess.run(
                f"cd {self.path_cwd}; {self.path_moltemplatesrc}moltemplate.sh ./system.lt",
                shell=True,
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

        with open(output, "w") as write_f:
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

        with open(output, "w") as write_f:
            write_f.write("import \"oplsaa.lt\"\n")

            # Import unique monomers - ensure proper .lt extension
            unique_monomers = list(dict.fromkeys(monomer_set))
            for monomer in unique_monomers:
                # Remove .lt if it exists, then add it back
                base_name = monomer[:-3] if monomer.endswith('.lt') else monomer
                write_f.write(f"import \"{base_name}.lt\"\n")

            write_f.write("\n")

            # Define combined molecule (ex.polymer)
            write_f.write(f"poly_{poly_index+1} inherits OPLSAA {{\n\n")
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
                    monomer_bank = Path(self.path_monomer_bank)
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
                    monomer_bank = Path(self.path_monomer_bank)
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
        MonomerBank=Path(self.path_monomer_bank)
        merltfile_Path=MonomerBank/merltfile
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
            
    def make_oplsaalt(self) -> None:
        """Creates the oplsaa.lt file."""
        try:
            self.make_oplsaa_subset()

            # Invoke oplsaa_moltemplate.py to make oplsaa.lt with suppressed output
            oplsaa_subset = self.path_cwd + "oplsaa_subset.prm"
            oplsaa_py = self.path_moltemplatesrc + "oplsaa_moltemplate.py " + oplsaa_subset
            
            # Redirect both stdout and stderr to devnull
            with open(os.devnull, 'w') as devnull:
                return_code = subprocess.call(f"cd {self.path_cwd}; {oplsaa_py}", 
                                            shell=True,
                                            stdout=devnull,
                                            stderr=devnull)
            
            if return_code != 0:
                logger.error("Failed to generate oplsaa.lt file. Check oplsaa_subset.prm for errors.")
                sys.exit(1)
            
        except Exception as e:
            logger.error(f"Error in make_oplsaalt: {str(e)}")
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
                    MonomerBank=Path(self.path_monomer_bank)
                    merltfile_Path=MonomerBank/monomerSet[vecii]
                    if merltfile_Path.is_file():
                        mono = self.path_monomer_bank+monomerSet[vecii]
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

    def check_monomer_bank(self, monomer: str) -> bool:
        """
        Check if a monomer template exists in the monomer bank.
        
        This method verifies that a specific monomer .lt file exists in the
        monomer bank directory. It provides detailed logging when a monomer
        is not found to help with debugging.
        
        Args:
            monomer (str): Name of the monomer file to check (with or without .lt extension)
        
        Returns:
            bool: True if the monomer file exists, False otherwise
        """
        monomer_path = Path(self.path_monomer_bank) / monomer
        exists = monomer_path.is_file()
        if not exists:
            logger.warning(f"Could not find monomer file: {monomer_path}")
            logger.info(f"Looking in directory: {self.path_monomer_bank}")
            logger.info(f"Searching for file: {monomer}")
            logger.info(f"Available files in monomer bank: {sorted(list(Path(self.path_monomer_bank).glob('*.lt')))}")
            logger.info(f"Monomer bank path exists: {Path(self.path_monomer_bank).exists()}")
        return exists

    def copy_to_cwd(self, source: Path) -> None:
        """
        Copy the specified source file to the current working directory.
        
        This method copies monomer template files from the monomer bank
        to the current working directory for processing by Moltemplate.
        
        Args:
            source (Path): The path of the source file to copy.
        """
        bash = "cp "
        bash = bash + str(source) + " " + self.path_cwd
        os.system(bash)

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