#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Force Field Management Module for AutoPoly Package

This module provides the ForceFieldManager class for managing force field
parameters and generating force field files for LAMMPS molecular dynamics simulations.

The ForceFieldManager class handles:
- OPLS-AA force field parameter management
- GAFF (General AMBER Force Field) support
- OPLS-AA subset generation based on monomer types
- Alkyl dihedral parameter modification
- Force field file generation for Moltemplate

Key Features:
- Support for OPLS-AA, LOPLS, and GAFF force fields
- Automatic parameter subset generation
- Dihedral coefficient customization
- Integration with Moltemplate force field system

Dependencies:
- Moltemplate: For generating LAMMPS data files
- OPLS-AA/GAFF force field parameters
- GAFFAnalyzer: For GAFF-specific analysis

Created on 2026-01-06
@author: zwu
"""
import sys
import os
from pathlib import Path
import subprocess
import shutil
import numpy as np
from typing import List, Optional
from .system import logger
from .gaff_analysis import GAFFAnalyzer


class ForceFieldManager:
    """
    Manages force field parameters and force field file generation.

    This class handles all force field-related operations including parameter
    subset generation, force field file creation, and parameter modification
    for LAMMPS molecular dynamics simulations.

    Attributes:
        path_cwd (str): Current working directory for the project
        path_master (str): Path to external dependencies
        path_moltemplatesrc (str): Path to Moltemplate source
        force_field (str): Force field type ("oplsaa", "gaff", or "lopls")
        ff_modify_dihedral (np.ndarray): Modified dihedral parameters
        path_oplsaaprm (str): Path to OPLS-AA force field parameters
        gaff_analyzer (GAFFAnalyzer): GAFF analyzer instance (for GAFF force field)
    """

    def __init__(self, path_cwd: str, path_master: str, path_moltemplatesrc: str,
                 force_field: str, ff_modify_dihedral: np.ndarray) -> None:
        """
        Initialize the ForceFieldManager.

        Args:
            path_cwd (str): Current working directory for the project
            path_master (str): Path to external dependencies
            path_moltemplatesrc (str): Path to Moltemplate source
            force_field (str): Force field type ("oplsaa", "gaff", or "lopls")
            ff_modify_dihedral (np.ndarray): Modified dihedral parameters
        """
        self.path_cwd = path_cwd
        self.path_master = path_master
        self.path_moltemplatesrc = path_moltemplatesrc
        self.force_field = force_field
        self.ff_modify_dihedral = ff_modify_dihedral

        # Set force field parameter path based on force_field type
        if force_field == "gaff":
            self.path_oplsaaprm = f"{self.path_master}moltemplate/common/gaff.lt"
        elif force_field == "lopls":
            self.path_oplsaaprm = f"{self.path_master}moltemplate/loplsaa.prm"
        else:  # oplsaa
            self.path_oplsaaprm = f"{self.path_master}moltemplate/oplsaa.prm"

        logger.info(f"\n'you are now using parameter set of {self.path_oplsaaprm}\n")

        # Initialize GAFF analyzer if using GAFF force field
        self.gaff_analyzer = None
        if force_field == "gaff":
            self.gaff_analyzer = GAFFAnalyzer(self.path_cwd, self.path_master)

    def make_force_field_lt(self) -> None:
        """
        Creates the force field .lt file based on the force_field type.

        This is the main entry point for force field file generation. It delegates
        to the appropriate method based on the force field type:
        - For GAFF: creates a filtered subset using GAFFAnalyzer
        - For OPLS-AA/LOPLS: creates a subset and invokes oplsaa_moltemplate.py

        Raises:
            SystemExit: If force field file generation fails
        """
        try:
            if self.force_field == "gaff":
                self.make_gaff_lt()
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

    def make_gaff_lt(self) -> None:
        """
        Create GAFF force field file using GAFFAnalyzer.

        This method delegates to the GAFFAnalyzer class to analyze the topology
        used in the monomer files and creates a filtered version of gaff.lt
        containing only the relevant parameters.

        The method:
        1. Validates that GAFF analyzer is initialized
        2. Calls GAFFAnalyzer.create_gaff_subset() to analyze monomers
        3. Generates gaff.lt with only necessary parameters

        Raises:
            SystemExit: If GAFF analyzer is not initialized or analysis fails

        Note:
            This requires force_field="gaff" to be set during initialization.
            The model attribute must be set before calling this method.
        """
        if self.gaff_analyzer is None:
            logger.error("GAFF analyzer not initialized. Set force_field='gaff' to use this method.")
            sys.exit(1)

        if not hasattr(self, 'model'):
            logger.error("Model attribute not set. Please set ff_manager.model before calling make_gaff_lt().")
            sys.exit(1)

        self.gaff_analyzer.create_gaff_subset(self.model)

    def make_oplsaa_subset(self) -> None:
        """
        Creates a subset of the OPLS-AA parameters based on the models.

        This method analyzes all monomers in the polymer models and extracts
        only the relevant OPLS-AA parameters from the master parameter file.
        This reduces the size of the force field file and improves performance.

        The process:
        1. Parses all monomer .lt files to extract atom types used
        2. Removes duplicate atom types and sorts them
        3. Reads the master OPLS-AA parameter file
        4. Writes only the relevant parameters to oplsaa_subset.prm

        The generated subset file includes:
        - Atom type definitions for all atoms in the system
        - All other parameter sections (bonds, angles, dihedrals, etc.)

        Raises:
            SystemExit: If a monomer file cannot be found or read
        """
        # path to oplsaa_subset.prm file
        opls_subset_file = self.path_cwd + "oplsaa_subset.prm"

        atom_keys = []
        for modelii in self.model:
            for monomerii in range(len(modelii.sequenceSet)):
                monomerSet = modelii.sequenceSet[monomerii]
                # vector to store all atom types including the repeats
                for vecii in range(len(monomerSet)):
                    # path to monomer.lt in monomer bank
                    MonomerBank = Path(self.path_cwd)
                    merltfile_Path = MonomerBank / monomerSet[vecii]
                    if merltfile_Path.is_file():
                        mono = str(MonomerBank / monomerSet[vecii])
                        read_switch = False
                        with open(mono) as f:
                            while True:
                                line = f.readline()

                                if line.strip() == "write(\"Data Atoms\") {":
                                    read_switch = True
                                    continue
                                elif line.strip() == "}":
                                    read_switch = False
                                    break

                                # Determine atom types, element names and the raw_charges as
                                # given in the opls table

                                if read_switch:
                                    load_line = ""
                                    stringvector = line.split()

                                    load_switch = False
                                    for readii in range(len(stringvector[2])):
                                        if stringvector[2][readii] == ":":
                                            load_switch = True
                                            continue
                                        if load_switch:
                                            load_line += stringvector[2][readii]
                                    atom_keys.append(load_line)

                                if not line:
                                    break
                    else:
                        logger.error(' '.join(["Monomer (" + monomerSet[vecii] + ") does NOT exist. \n",
                                                "Please check the following path to the file\n" + merltfile_Path + "\n"]))
                        sys.exit()

        # Cleaning up the stored data. Remove duplicate atoms types
        atom_types = list(dict.fromkeys(atom_keys))
        # Convert the vectors string to vector int in order to sort the atom_types in ascending order
        atom_types = sorted([int(i) for i in atom_types])

        # Read the master opls file and store the ones that match the atom_types into new subset file
        write_f = open(opls_subset_file, "w")

        with open(self.path_oplsaaprm, 'r') as read_f:
            path_oplsaaprm = Path(self.path_oplsaaprm)
            if path_oplsaaprm.is_file():
                check_switch = False
                while True:
                    prm_line = read_f.readline()
                    if len(prm_line.strip()) != 0:

                        if prm_line.strip() == "##  Atom Type Definitions  ##":
                            check_switch = True
                            write_f.write(prm_line + "\n")
                            prm_line = read_f.readline()
                            write_f.write(prm_line + "\n")
                            prm_line = read_f.readline()
                            write_f.write(prm_line + "\n")
                            continue
                        elif prm_line.strip() == "################################":
                            check_switch = False
                            write_f.write(prm_line + "\n")
                            continue
                        elif check_switch:

                            stringvector = prm_line.split()


                            for checkii in range(len(atom_types)):

                                if atom_types[checkii] == int(stringvector[1]):
                                    write_f.write(prm_line + "\n")
                                    break
                        else:
                            write_f.write(prm_line + "\n")
                    else:
                        write_f.write(prm_line + "\n")

                    if not prm_line:
                        break
        write_f.close()

    def FFmodify_alkyl_dihedral_oplsaa(self) -> None:
        """
        Modifies alkyl dihedral coefficients in the OPLS-AA force field.

        This method modifies the dihedral coefficients for alkyl chains (CH3-CH2-CH2-CH3
        and similar) in the generated oplsaa.lt file. The modification uses custom
        dihedral parameters stored in the ff_modify_dihedral attribute.

        The process:
        1. Reads the generated oplsaa.lt file
        2. Identifies dihedral coefficients for alkyl chains (atom types 80, 81, 82)
        3. Replaces them with modified parameters
        4. Writes the modified file back

        Alkyl atom types:
        - 80: CH3 (methyl)
        - 81: CH2 (methylene)
        - 82: CH (methine)

        Raises:
            SystemExit: If oplsaa.lt file cannot be opened or modification fails

        Note:
            The ff_modify_dihedral attribute should contain 4 dihedral coefficients
            in the format [K1, K2, K3, K4] in units of kcal/mol.
        """
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
                            string_strip = next_line.strip()
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
                                # Note: ff_modify_dihedral should be defined as a class attribute
                                if hasattr(self, 'ff_modify_dihedral'):
                                    write_f.write(" " + " ".join(str(x) for x in self.ff_modify_dihedral))
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
                        string_strip = line.strip()
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
                            # Note: ff_modify_dihedral should be defined as a class attribute
                            if hasattr(self, 'ff_modify_dihedral'):
                                write_f.write(" " + " ".join(str(x) for x in self.ff_modify_dihedral))
                            write_f.write("\n")
                        else:
                            write_f.write(f"{line}\n")

            # Replace original file with modified version
            shutil.move(str(output_file), str(input_file))

        except Exception as e:
            logger.error(f"Error modifying alkyl dihedral coefficients: {str(e)}")
            sys.exit(1)
