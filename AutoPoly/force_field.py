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
from typing import List, Optional, Set, Dict, Tuple
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
            self.path_oplsaaprm = str(Path(self.path_master) / "moltemplate" / "common" / "gaff.lt")
        elif force_field == "lopls":
            self.path_oplsaaprm = str(Path(self.path_master) / "moltemplate" / "loplsaa.prm")
        else:  # oplsaa
            self.path_oplsaaprm = str(Path(self.path_master) / "moltemplate" / "oplsaa.prm")

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
                ff_subset = str(Path(self.path_cwd) / "oplsaa_subset.prm")
                ff_py_script = str(Path(self.path_moltemplatesrc) / "oplsaa_moltemplate.py")

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
        2. Builds bond connectivity graphs from monomers
        3. Infers angles, dihedrals, and impropers from actual topology
        4. Maps atom type combinations to atom class combinations
        5. Filters parameters by actual topology, not just atom class presence

        The generated subset file includes:
        - Atom type definitions for atoms in the system
        - VDW parameters filtered by atom type
        - Bond parameters filtered by actual bonded atom class pairs
        - Angle parameters filtered by actual angle atom class triplets
        - Torsion parameters filtered by actual dihedral atom class quads (with wildcard support)
        - Improper torsion parameters filtered by actual improper atom class quads (with wildcard support)

        Raises:
            SystemExit: If a monomer file cannot be found or read
        """
        # path to oplsaa_subset.prm file
        opls_subset_file = str(Path(self.path_cwd) / "oplsaa_subset.prm")

        # Collect monomer file paths and atom types
        atom_keys = []
        monomer_files = []
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
                        monomer_files.append(mono)
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
                                                "Please check the following path to the file\n" + str(merltfile_Path) + "\n"]))
                        sys.exit()

        # Remove duplicates from monomer files
        monomer_files = list(dict.fromkeys(monomer_files))

        # Cleaning up the stored data. Remove duplicate atoms types
        atom_types = list(dict.fromkeys(atom_keys))
        # Convert the vectors string to vector int in order to sort the atom_types in ascending order
        atom_types_set = set(int(i) for i in atom_types)
        atom_types = sorted(atom_types_set)

        # FIRST PASS: Extract atom classes for the used atom types and build type->class mapping
        atom_classes = set()
        type_to_class: Dict[int, int] = {}
        path_oplsaaprm = Path(self.path_oplsaaprm)
        if not path_oplsaaprm.is_file():
            logger.error(f"OPLS-AA parameter file not found: {self.path_oplsaaprm}")
            sys.exit(1)

        with open(self.path_oplsaaprm, 'r') as read_f:
            for line in read_f:
                stripped = line.strip()
                if stripped.startswith('atom'):
                    parts = stripped.split()
                    if len(parts) >= 3:
                        try:
                            atom_type = int(parts[1])
                            atom_class = int(parts[2])
                            if atom_type in atom_types_set:
                                atom_classes.add(atom_class)
                                type_to_class[atom_type] = atom_class
                        except ValueError:
                            continue

        logger.info(f"OPLS-AA subset: {len(atom_types)} atom types, {len(atom_classes)} atom classes")

        # TOPOLOGY EXTRACTION: Build actual connectivity from monomers
        logger.info("  Extracting topology from monomers...")

        # Extract bonds from monomers
        used_bonds_types = self._extract_bonds_from_monomers_oplsaa(monomer_files)
        logger.info(f"    Found {len(used_bonds_types)} unique bond types from monomers")

        # Build bond graphs and infer angles/dihedrals/impropers
        used_angles_types: Set[Tuple[int, int, int]] = set()
        used_dihedrals_types: Set[Tuple[int, int, int, int]] = set()
        used_impropers_types: Set[Tuple[int, int, int, int]] = set()

        for monofile in monomer_files:
            try:
                graph = self._build_bond_graph_oplsaa(monofile)
                if graph:
                    used_angles_types.update(self._infer_angles_from_graph_oplsaa(graph))
                    used_dihedrals_types.update(self._infer_dihedrals_from_graph_oplsaa(graph))
                    used_impropers_types.update(self._infer_impropers_from_graph_oplsaa(graph))
            except Exception as e:
                logger.warning(f"    Could not process {monofile}: {e}")
                continue

        logger.info(f"    Found {len(used_angles_types)} unique angle types")
        logger.info(f"    Found {len(used_dihedrals_types)} unique dihedral types")
        logger.info(f"    Found {len(used_impropers_types)} unique improper types")

        # Convert atom type combinations to atom class combinations
        bond_class_pairs = self._convert_bonds_to_classes(used_bonds_types, type_to_class)
        angle_class_triplets = self._convert_angles_to_classes(used_angles_types, type_to_class)
        dihedral_class_quads = self._convert_dihedrals_to_classes(used_dihedrals_types, type_to_class)
        improper_class_quads = self._convert_impropers_to_classes(used_impropers_types, type_to_class)

        logger.info(f"    Mapped to {len(bond_class_pairs)} bond class pairs")
        logger.info(f"    Mapped to {len(angle_class_triplets)} angle class triplets")
        logger.info(f"    Mapped to {len(dihedral_class_quads)} dihedral class quads")
        logger.info(f"    Mapped to {len(improper_class_quads)} improper class quads")

        # SECOND PASS: Write filtered output
        # Section markers for OPLS-AA parameter file
        SECTION_MARKERS = {
            '##  Atom Type Definitions  ##': 'atom',
            '##  Van der Waals Parameters  ##': 'vdw',
            '##  Bond Stretching Parameters  ##': 'bond',
            '##  Angle Bending Parameters  ##': 'angle',
            '##  Torsional Parameters  ##': 'torsion',
            '##  Improper Torsional Parameters  ##': 'imptors',
            '##   Urey-Bradley Parameters  ##': 'ureybrad',
            '##  Atomic Partial Charge Parameters  ##': 'charge',
            '##  Biopolymer Atom Type Conversions  ##': 'biotype',
        }

        def should_include_parameter(line: str, section: str) -> bool:
            """Check if a parameter line should be included based on the current section."""
            parts = line.split()
            if len(parts) < 2:
                return True  # Empty or header lines

            keyword = parts[0].lower()

            # Filter based on section type
            if keyword == 'atom' and section == 'atom':
                # Atom type definitions: filter by atom type (column 1)
                if len(parts) >= 2:
                    try:
                        return int(parts[1]) in atom_types_set
                    except ValueError:
                        return True

            elif keyword == 'vdw' and section == 'vdw':
                # VDW parameters: filter by atom type (column 1)
                if len(parts) >= 2:
                    try:
                        return int(parts[1]) in atom_types_set
                    except ValueError:
                        return True

            elif keyword == 'bond' and section == 'bond':
                # Bond parameters: filter by actual bond class pairs
                if len(parts) >= 3:
                    try:
                        class1 = int(parts[1])
                        class2 = int(parts[2])
                        # Normalize to (min, max) for comparison
                        bond_pair = (min(class1, class2), max(class1, class2))
                        return bond_pair in bond_class_pairs
                    except ValueError:
                        return True

            elif keyword == 'angle' and section == 'angle':
                # Angle parameters: filter by actual angle class triplets
                if len(parts) >= 4:
                    try:
                        class1 = int(parts[1])
                        class2 = int(parts[2])
                        class3 = int(parts[3])
                        # Normalize: (min_outer, center, max_outer)
                        if class1 <= class3:
                            angle_triplet = (class1, class2, class3)
                        else:
                            angle_triplet = (class3, class2, class1)
                        return angle_triplet in angle_class_triplets
                    except ValueError:
                        return True

            elif keyword == 'torsion' and section == 'torsion':
                # Torsion parameters: filter by actual dihedral class quads
                # Class 0 is a wildcard that matches any class
                if len(parts) >= 5:
                    try:
                        classes = [int(parts[i]) for i in range(1, 5)]
                        # Check if this torsion matches any actual dihedral
                        for actual_dihedral in dihedral_class_quads:
                            match = True
                            for i in range(4):
                                if classes[i] != 0 and classes[i] != actual_dihedral[i]:
                                    match = False
                                    break
                            if match:
                                return True
                        # Also check reverse direction
                        for actual_dihedral in dihedral_class_quads:
                            rev_dihedral = (actual_dihedral[3], actual_dihedral[2],
                                          actual_dihedral[1], actual_dihedral[0])
                            match = True
                            for i in range(4):
                                if classes[i] != 0 and classes[i] != rev_dihedral[i]:
                                    match = False
                                    break
                            if match:
                                return True
                        return False
                    except ValueError:
                        return True

            elif keyword == 'imptors' and section == 'imptors':
                # Improper torsion parameters: filter by actual improper class quads
                # Class 0 is a wildcard that matches any class
                if len(parts) >= 5:
                    try:
                        classes = [int(parts[i]) for i in range(1, 5)]
                        # Check if this improper matches any actual improper
                        for actual_improper in improper_class_quads:
                            match = True
                            for i in range(4):
                                if classes[i] != 0 and classes[i] != actual_improper[i]:
                                    match = False
                                    break
                            if match:
                                return True
                        return False
                    except ValueError:
                        return True

            elif keyword == 'charge' and section == 'charge':
                # Charge parameters: filter by atom type (column 1)
                if len(parts) >= 2:
                    try:
                        return int(parts[1]) in atom_types_set
                    except ValueError:
                        return True

            # For other sections or non-parameter lines, include them
            return True

        current_section = None
        with open(self.path_oplsaaprm, 'r') as read_f, open(opls_subset_file, 'w') as write_f:
            for line in read_f:
                stripped = line.strip()

                # Check for section markers
                for marker, section_name in SECTION_MARKERS.items():
                    if marker in stripped:
                        current_section = section_name
                        break

                # Determine if this line should be written
                if current_section in ('atom', 'vdw', 'bond', 'angle', 'torsion', 'imptors', 'charge'):
                    # Check if it's a parameter line that needs filtering
                    if stripped and not stripped.startswith('#') and not stripped.startswith('##'):
                        if should_include_parameter(stripped, current_section):
                            write_f.write(line)
                    else:
                        # Header/comment lines - always write
                        write_f.write(line)
                else:
                    # For other sections (force field definition, literature, etc.), write as-is
                    write_f.write(line)

    # -------------------------------------------------------------------------
    # Topology extraction methods for OPLS-AA parameter filtering
    # -------------------------------------------------------------------------

    def _extract_bonds_from_monomers_oplsaa(self, monomer_files: List[str]) -> Set[Tuple[int, int]]:
        """Extract bond pairs (atom types) from OPLS-AA monomer .lt files.

        Parses the 'Data Bond List' sections to find which atom type pairs
        are actually bonded in the monomers.

        Args:
            monomer_files: List of paths to monomer .lt files

        Returns:
            Set of bond tuples (min_type, max_type) normalized for consistency.
        """
        used_bonds = set()

        for filepath in monomer_files:
            try:
                # First, build atom ID -> atom type mapping
                atom_types_map = {}
                in_atoms_section = False
                with open(filepath, 'r') as f:
                    for line in f:
                        if 'write("Data Atoms")' in line:
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
                                    # Extract integer atom type after colon
                                    type_str = part.split(':')[1]
                                    try:
                                        atom_type = int(type_str)
                                    except ValueError:
                                        continue
                            if atom_id and atom_type is not None:
                                atom_types_map[atom_id] = atom_type

                # Now extract bonds using the atom type mapping
                in_bond_section = False
                with open(filepath, 'r') as f:
                    for line in f:
                        if "write('Data Bond List')" in line or 'write("Data Bond List")' in line:
                            in_bond_section = True
                            continue

                        if in_bond_section:
                            if '}' in line and not line.strip().startswith('#'):
                                break

                            if '$bond:' in line and '$atom:' in line:
                                parts = line.split()
                                atom_ids = []
                                for part in parts:
                                    if part.startswith('$atom:'):
                                        atom_id = part.split(':')[1]
                                        atom_ids.append(atom_id)

                                if len(atom_ids) >= 2:
                                    atom1_type = atom_types_map.get(atom_ids[0])
                                    atom2_type = atom_types_map.get(atom_ids[1])

                                    if atom1_type is not None and atom2_type is not None:
                                        # Normalize: (min, max) for consistency
                                        bond = (min(atom1_type, atom2_type), max(atom1_type, atom2_type))
                                        used_bonds.add(bond)

            except Exception as e:
                logger.warning(f"Could not parse bonds from {filepath}: {e}")
                continue

        return used_bonds

    def _build_bond_graph_oplsaa(self, filepath: str) -> Dict[str, Tuple[int, List[str]]]:
        """Build a bond connectivity graph from an OPLS-AA monomer .lt file.

        Args:
            filepath: Path to monomer .lt file

        Returns:
            Dictionary mapping atom_id -> (atom_type, [neighbor_atom_ids])
        """
        graph = {}

        try:
            with open(filepath, 'r') as f:
                in_atoms_section = False
                in_bond_section = False

                for line in f:
                    # Parse Data Atoms section
                    if 'write("Data Atoms")' in line:
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
                                type_str = part.split(':')[1]
                                try:
                                    atom_type = int(type_str)
                                except ValueError:
                                    continue
                        if atom_id and atom_type is not None:
                            graph[atom_id] = (atom_type, [])

                    # Parse Data Bond List section
                    if "write('Data Bond List')" in line or 'write("Data Bond List")' in line:
                        in_bond_section = True
                        continue
                    elif in_bond_section and '}' in line and not line.strip().startswith('#'):
                        in_bond_section = False

                    if in_bond_section and '$bond:' in line:
                        parts = line.split()
                        bonded_atoms = []
                        for part in parts:
                            if part.startswith('$atom:'):
                                atom_id = part.split(':')[1]
                                bonded_atoms.append(atom_id)

                        # Add edges to graph (bidirectional)
                        if len(bonded_atoms) >= 2:
                            atom1, atom2 = bonded_atoms[0], bonded_atoms[1]
                            if atom1 in graph and atom2 in graph:
                                graph[atom1][1].append(atom2)
                                graph[atom2][1].append(atom1)

        except Exception as e:
            logger.warning(f"Could not build bond graph from {filepath}: {e}")

        return graph

    def _infer_angles_from_graph_oplsaa(self, bond_graph: Dict) -> Set[Tuple[int, int, int]]:
        """Infer angle types from bond connectivity graph.

        For each atom with 2+ neighbors, generates all angle combinations
        with that atom as the center.

        Args:
            bond_graph: Dict mapping atom_id -> (atom_type, [neighbor_ids])

        Returns:
            Set of angle tuples (type1, center_type, type2) normalized.
        """
        angles = set()

        for center_id, (center_type, neighbors) in bond_graph.items():
            if len(neighbors) >= 2:
                for i in range(len(neighbors)):
                    for j in range(i + 1, len(neighbors)):
                        atom1_type = bond_graph[neighbors[i]][0]
                        atom2_type = bond_graph[neighbors[j]][0]

                        # Normalize: (min_outer, center, max_outer)
                        if atom1_type <= atom2_type:
                            angle = (atom1_type, center_type, atom2_type)
                        else:
                            angle = (atom2_type, center_type, atom1_type)
                        angles.add(angle)

        return angles

    def _infer_dihedrals_from_graph_oplsaa(self, bond_graph: Dict) -> Set[Tuple[int, int, int, int]]:
        """Infer dihedral types from bond connectivity graph.

        For each bond as the central bond, finds all atoms bonded to each end
        and generates dihedral combinations.

        Args:
            bond_graph: Dict mapping atom_id -> (atom_type, [neighbor_ids])

        Returns:
            Set of dihedral tuples (type1, type2, type3, type4).
        """
        dihedrals = set()

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
                        dihedral = (outer1_type, atom1_type, atom2_type, outer2_type)
                        dihedrals.add(dihedral)

        return dihedrals

    def _infer_impropers_from_graph_oplsaa(self, bond_graph: Dict) -> Set[Tuple[int, int, int, int]]:
        """Infer improper types from bond connectivity graph.

        For each atom with 3+ neighbors (central atom of improper),
        generates all improper combinations.

        Args:
            bond_graph: Dict mapping atom_id -> (atom_type, [neighbor_ids])

        Returns:
            Set of improper tuples (type1, type2, type3, center_type).
        """
        impropers = set()

        for center_id, (center_type, neighbors) in bond_graph.items():
            if len(neighbors) >= 3:
                for i in range(len(neighbors)):
                    for j in range(i + 1, len(neighbors)):
                        for k in range(j + 1, len(neighbors)):
                            atom1_type = bond_graph[neighbors[i]][0]
                            atom2_type = bond_graph[neighbors[j]][0]
                            atom3_type = bond_graph[neighbors[k]][0]

                            # Create improper: (type1, type2, type3, center_type)
                            improper = (atom1_type, atom2_type, atom3_type, center_type)
                            impropers.add(improper)

        return impropers

    def _map_type_to_class(self, atom_type: int, type_to_class: Dict[int, int]) -> int:
        """Map atom type to atom class.

        Args:
            atom_type: Integer atom type
            type_to_class: Dict mapping atom types to atom classes

        Returns:
            Atom class (or atom type if not found in mapping)
        """
        return type_to_class.get(atom_type, atom_type)

    def _convert_bonds_to_classes(
        self, bonds: Set[Tuple[int, int]], type_to_class: Dict[int, int]
    ) -> Set[Tuple[int, int]]:
        """Convert bond type pairs to class pairs.

        Args:
            bonds: Set of (type1, type2) tuples
            type_to_class: Dict mapping atom types to atom classes

        Returns:
            Set of (class1, class2) tuples, normalized (min, max)
        """
        class_bonds = set()
        for t1, t2 in bonds:
            c1 = self._map_type_to_class(t1, type_to_class)
            c2 = self._map_type_to_class(t2, type_to_class)
            class_bonds.add((min(c1, c2), max(c1, c2)))
        return class_bonds

    def _convert_angles_to_classes(
        self, angles: Set[Tuple[int, int, int]], type_to_class: Dict[int, int]
    ) -> Set[Tuple[int, int, int]]:
        """Convert angle type triplets to class triplets.

        Args:
            angles: Set of (type1, center_type, type2) tuples
            type_to_class: Dict mapping atom types to atom classes

        Returns:
            Set of (class1, center_class, class2) tuples, normalized
        """
        class_angles = set()
        for t1, tc, t2 in angles:
            c1 = self._map_type_to_class(t1, type_to_class)
            cc = self._map_type_to_class(tc, type_to_class)
            c2 = self._map_type_to_class(t2, type_to_class)
            # Normalize: (min_outer, center, max_outer)
            if c1 <= c2:
                class_angles.add((c1, cc, c2))
            else:
                class_angles.add((c2, cc, c1))
        return class_angles

    def _convert_dihedrals_to_classes(
        self, dihedrals: Set[Tuple[int, int, int, int]], type_to_class: Dict[int, int]
    ) -> Set[Tuple[int, int, int, int]]:
        """Convert dihedral type quads to class quads.

        Args:
            dihedrals: Set of (type1, type2, type3, type4) tuples
            type_to_class: Dict mapping atom types to atom classes

        Returns:
            Set of (class1, class2, class3, class4) tuples
        """
        class_dihedrals = set()
        for t1, t2, t3, t4 in dihedrals:
            c1 = self._map_type_to_class(t1, type_to_class)
            c2 = self._map_type_to_class(t2, type_to_class)
            c3 = self._map_type_to_class(t3, type_to_class)
            c4 = self._map_type_to_class(t4, type_to_class)
            class_dihedrals.add((c1, c2, c3, c4))
        return class_dihedrals

    def _convert_impropers_to_classes(
        self, impropers: Set[Tuple[int, int, int, int]], type_to_class: Dict[int, int]
    ) -> Set[Tuple[int, int, int, int]]:
        """Convert improper type quads to class quads.

        Args:
            impropers: Set of (type1, type2, type3, center_type) tuples
            type_to_class: Dict mapping atom types to atom classes

        Returns:
            Set of (class1, class2, class3, center_class) tuples
        """
        class_impropers = set()
        for t1, t2, t3, tc in impropers:
            c1 = self._map_type_to_class(t1, type_to_class)
            c2 = self._map_type_to_class(t2, type_to_class)
            c3 = self._map_type_to_class(t3, type_to_class)
            cc = self._map_type_to_class(tc, type_to_class)
            class_impropers.add((c1, c2, c3, cc))
        return class_impropers

    # Alkyl atom types: CH3 (methyl), CH2 (methylene), CH (methine)
    ALKYL_ATOM_TYPES = {80, 81, 82}

    def _is_alkyl_dihedral(self, line: str) -> bool:
        """
        Check if a line contains an alkyl dihedral specification.

        Args:
            line: The line to check

        Returns:
            bool: True if the line specifies a dihedral with all alkyl atoms
        """
        if '@dihedral:' not in line:
            return False

        try:
            dihedral_spec = line.split("@dihedral:")[1].split()[0]
            atom_types = [int(x) for x in dihedral_spec.split('-')]
            return all(atom in self.ALKYL_ATOM_TYPES for atom in atom_types)
        except (IndexError, ValueError):
            return False

    def _write_modified_dihedral_line(self, write_f, line: str) -> None:
        """
        Write a modified dihedral coefficient line for alkyl chains.

        Args:
            write_f: File object to write to
            line: Original line containing dihedral specification
        """
        dihedral_spec = line.split("@dihedral:")[1].split()[0]
        atom_types = [int(x) for x in dihedral_spec.split('-')]

        write_f.write("dihedral_coeff @dihedral:")
        write_f.write('-'.join(str(x) for x in atom_types))
        write_f.write(" opls")

        if hasattr(self, 'ff_modify_dihedral'):
            write_f.write(" " + " ".join(str(x) for x in self.ff_modify_dihedral))
        write_f.write("\n")

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

        if not input_file.exists():
            logger.error("oplsaa.lt file cannot open.")
            sys.exit(1)

        logger.info("Start modifying alkyl dihedral coefficients")

        try:
            with open(input_file, 'r') as read_f, open(output_file, 'w') as write_f:
                for line in read_f:
                    stripped = line.strip()
                    words = stripped.split()

                    # Handle write_once block start
                    if words and words[0] == 'write_once("In':
                        write_f.write(stripped + "\n")
                        next_line = next(read_f).strip()
                        next_words = next_line.split()

                        if next_words and next_words[0] == "dihedral_coeff":
                            if self._is_alkyl_dihedral(next_line):
                                self._write_modified_dihedral_line(write_f, next_line)
                            else:
                                write_f.write(next_line + "\n")
                        else:
                            write_f.write(next_line + "\n")
                        continue

                    # Handle lines inside dihedral_coeff block
                    if words and words[0] == "dihedral_coeff":
                        if self._is_alkyl_dihedral(stripped):
                            self._write_modified_dihedral_line(write_f, stripped)
                        else:
                            write_f.write(stripped + "\n")
                    else:
                        write_f.write(stripped + "\n")

            # Replace original file with modified version
            shutil.move(str(output_file), str(input_file))

        except Exception as e:
            logger.error(f"Error modifying alkyl dihedral coefficients: {str(e)}")
            sys.exit(1)
