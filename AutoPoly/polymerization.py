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
import os
from pathlib import Path
import subprocess
import shutil
import re
import numpy as np
from typing import List, Optional, Dict, Any, Set, Tuple
from .system import logger
from .exceptions import ValidationError
from .monomer_generator import MonomerGenerator
from .file_management import (
    create_working_directory,
    get_rid_of_lj_cut_coul_long,
    mv_files
)
from .gaff_analysis import GAFFAnalyzer
from . import monomer_processing
from .force_field import ForceFieldManager
from .workflow import WorkflowManager

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
        force_field (str): Force field to use ("oplsaa", "gaff", or "lopls")
        model (list): List of polymer models to generate
        rotate (float): Rotation angle for monomer placement
        offset_spacing (float): Spacing between polymer chains
        offset (float): Offset distance for monomer placement
        packingL_spacing (float): Packing length spacing
        moltemplate_box_size (float): Box size for Moltemplate
        FFmodify_alkylDihedral (np.array): Modified dihedral parameters
    """
    
    def __init__(self, name: str = None, system: object = None, model: list = None,
                 run: bool = True, force_field: str = "oplsaa") -> None:
        """
        Initialize the Polymerization class.

        Args:
            name (str, optional): Name of the polymerization project. Defaults to None.
            system (object, optional): System object containing folder path. Defaults to None.
            model (list, optional): List of models for polymerization. Defaults to None.
            run (bool, optional): Flag to run the process immediately. Defaults to True.
            force_field (str, optional): Force field to use. Options: "oplsaa", "gaff", "gaff2", "lopls".
                                        Defaults to "oplsaa".

        Raises:
            SystemExit: If required directories or files are not found
        """
        # Validate force_field parameter
        valid_force_fields = ["oplsaa", "gaff", "gaff2", "lopls"]
        if force_field not in valid_force_fields:
            raise ValidationError(
                f"Invalid force_field '{force_field}'. Must be one of: {valid_force_fields}"
            )

        self.name = name
        self.system = system
        self.path_cwd = str(Path(self.system.get_folder_path()) / self.name / "moltemplate/")
        self.path_master = str(Path(__file__).parent.resolve() / "extern/")
        self.path_moltemplatesrc = str(Path(self.path_master) / "moltemplate" / "scripts/")

        # SMILES to monomer name cache for dynamic generation
        self._generated_smiles = {}
        self._smiles_to_name_counter = 0

        self.force_field = force_field

        # Set force field parameter path based on force_field type
        # All force fields now use .lt files from moltemplate/force_fields/
        if force_field == "gaff":
            self.path_oplsaaprm = str(Path(self.path_master) / "moltemplate" / "force_fields" / "gaff.lt")
        elif force_field == "gaff2":
            self.path_oplsaaprm = str(Path(self.path_master) / "moltemplate" / "force_fields" / "gaff2.lt")
        elif force_field == "lopls":
            self.path_oplsaaprm = str(Path(self.path_master) / "moltemplate" / "force_fields" / "loplsaa.lt")
        else:  # oplsaa
            self.path_oplsaaprm = str(Path(self.path_master) / "moltemplate" / "force_fields" / "oplsaa.lt")

        logger.info(f"\n'you are now using parameter set of {self.path_oplsaaprm}\n")
        self.model = model
        self.rotate = 90.0
        self.offset_spacing = 2.0
        self.offset = 4.0
        self.packingL_spacing = 5.0
        self.moltemplate_box_size = 400.0

        # Modified alkyl dihedral parameters (Kj/mol -> kcal/mol conversion)
        self.FFmodify_alkylDihedral = np.array([0.6446926386, -0.2143420172, 0.1782194073, 0.0])

        # Initialize ForceFieldManager
        self.ff_manager = ForceFieldManager(
            path_cwd=self.path_cwd,
            path_master=self.path_master,
            path_moltemplatesrc=self.path_moltemplatesrc,
            force_field=self.force_field,
            ff_modify_dihedral=self.FFmodify_alkylDihedral
        )

        # Initialize GAFF analyzer if using GAFF or GAFF2 force field
        self.gaff_analyzer = None
        if force_field in ("gaff", "gaff2"):
            self.gaff_analyzer = GAFFAnalyzer(self.path_cwd, self.path_master, force_field)

        # Create working directory before proceeding
        self.create_working_directory()

        logger.info(f"\n'you are now using extern path of {self.path_master}\n")

        # Initialize WorkflowManager
        self.workflow = WorkflowManager(self)

        if run:
            self.workflow.make_lmp_data_file_by_moltemplate()

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
        self.path_cwd = str(create_working_directory(self.system, self.name))

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
        return monomer_processing.n_monomer_atoms(merltfile, self.path_cwd)

    def extract_element_from_atom(self, atom_string: str):
        """
        Extract element name from atom identifier string.

        Args:
            atom_string (str): Atom identifier like "$atom:C1", "$atom:H16", etc.

        Returns:
            str: Element name (e.g., "C", "H", "Si", "Fe")
            None: If no match found
        """
        return monomer_processing.extract_element_from_atom(atom_string)

    def read_lt_end_atoms(self, lt_file: str):
        """
        Read the first and second atoms from a .lt file.

        Args:
            lt_file (str): Path to the .lt file

        Returns:
            tuple: (first_atom, second_atom) where each atom is a string with
                   element name and position (e.g., "C1", "H2")
        """
        return monomer_processing.read_lt_end_atoms(lt_file)

    def generate_monomer_from_psmiles(self, psmiles: str):
        """
        Generate monomer from pSMILES/SMILES string.

        This method wraps the monomer_processing function with instance-specific
        parameters like force_field, generated cache, and counter.

        Args:
            psmiles (str): pSMILES or SMILES string

        Returns:
            tuple: (base_name, updated_counter) where:
                - base_name: Generated monomer base name
                - updated_counter: Incremented counter value
        """
        base_name, counter = monomer_processing.generate_monomer_from_psmiles(
            psmiles,
            self.path_cwd,
            self.force_field,
            self._generated_smiles,
            self._smiles_to_name_counter
        )
        # Update the internal counter state
        self._smiles_to_name_counter = counter
        return base_name, counter

    def generate_sequence_variants_for_polymer(
        self,
        smiles_list: list,
        topology: str = "linear",
        base_name_prefix: str = "monomer"
    ) -> dict:
        """
        Generate sequence variants for polymerization using complement SMILES.

        This method uses the complement SMILES approach where each position in the
        chain has its own SMILES with explicit terminal groups.

        Args:
            smiles_list (list): List of complement SMILES:
                - First: 1 wildcard (right connection) e.g., 'CC[*]'
                - Middle: 2 wildcards (left and right) e.g., '[*]CC[*]'
                - Last: 1 wildcard (left connection) e.g., '[*]CC'
            topology (str): Topology type ("linear" or "ring", default: "linear")
            base_name_prefix (str): Prefix for monomer names (default: "monomer")

        Returns:
            dict: Mapping from variant_type to .lt filename
            {
                'first': 'monomer_0_0le.lt',
                'middle': 'monomer_0_1i.lt',
                'last': 'monomer_0_2re.lt',
                ...
            }
        """
        variant_mapping, counter = monomer_processing.generate_sequence_variants_for_polymerization(
            smiles_list=smiles_list,
            topology=topology,
            path_cwd=self.path_cwd,
            force_field=self.force_field,
            generated_cache=self._generated_smiles,
            counter=self._smiles_to_name_counter,
            base_name_prefix=base_name_prefix
        )
        # Update the internal counter state
        self._smiles_to_name_counter = counter
        return variant_mapping

    def generate_molecule_from_smiles(self, smiles: str, molecule_name: str):
        """
        Generate molecule from SMILES string.

        This method wraps the monomer_processing function with instance-specific
        parameters like force_field, generated cache, and counter.

        Args:
            smiles (str): SMILES string WITHOUT wildcards (e.g., "O", "CCO", "c1ccccc1")
            molecule_name (str): Name for the molecule (e.g., "water", "ethanol", "benzene")

        Returns:
            tuple: (filename, updated_counter) where:
                - filename: Generated molecule .lt filename
                - updated_counter: Incremented counter value
        """
        filename, counter = monomer_processing.generate_molecule_from_smiles(
            smiles,
            molecule_name,
            self.path_cwd,
            self.force_field,
            self._generated_smiles,
            self._smiles_to_name_counter
        )
        # Update the internal counter state
        self._smiles_to_name_counter = counter
        return filename, counter

    def make_lmp_data_file_by_moltemplate(self) -> None:
        """
        Generate the LAMMPS data file using Moltemplate.

        This method delegates to the WorkflowManager to handle the complete
        polymer generation process.

        Raises:
            SystemExit: If any critical step fails or required files are missing
        """
        self.workflow.make_lmp_data_file_by_moltemplate()

    def get_rid_of_lj_cut_coul_long(self) -> None:
        """Removes lj/cut/coul/long from the settings file."""
        get_rid_of_lj_cut_coul_long(self.path_cwd)

    def mv_files(self) -> None:
        """Moves generated files to the appropriate directories."""
        mv_files(self.path_cwd)

    def evaluate_offset(self, merltfile: str) -> None:
        """
        Evaluates the offset distance based on the specified merlt file.

        Args:
            merltfile (str): The name of the merlt file.
        """
        self.offset = monomer_processing.evaluate_offset(
            merltfile,
            self.path_cwd,
            self.offset_spacing,
            self.offset
        )
