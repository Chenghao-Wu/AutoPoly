#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GAFF Force Field Analysis Module for AutoPoly Package

This module provides the GAFFAnalyzer class for analyzing and filtering
GAFF (General Amber Force Field) parameters based on molecular topology.
It extracts relevant force field parameters from monomer structures and
creates optimized subsets of the full GAFF parameter set.

The GAFFAnalyzer class handles:
- Extraction of atom types from monomer templates
- Parsing and filtering GAFF force field files
- Building bond connectivity graphs
- Inferring angles, dihedrals, and impropers from molecular topology
- Creating optimized GAFF parameter subsets

Key Features:
- Analyzes molecular topology from moltemplate .lt files
- Builds bond graphs to infer force field parameters
- Filters GAFF parameters to only those used in the system
- Significantly reduces parameter file size for improved performance
- Handles inter-monomer bond formation during polymerization

Dependencies:
- Moltemplate .lt files for monomer definitions
- GAFF force field parameters (gaff.lt)
- Standard library modules for file parsing

Created on 2025-01-06
@author: AutoPy Development Team
"""
import sys
import shutil
from pathlib import Path
from typing import Set, Dict, List, Tuple, Optional
from .system import logger


class GAFFAnalyzer:
    """
    Analyzer for GAFF force field parameters and molecular topology.

    This class provides methods to analyze molecular structures from
    monomer .lt files, extract relevant force field parameters, and
    create optimized subsets of the full GAFF parameter set.

    Attributes:
        path_cwd (str): Current working directory for the project
        path_master (str): Path to external dependencies
    """

    def __init__(self, path_cwd: str, path_master: str):
        """
        Initialize the GAFFAnalyzer class.

        Args:
            path_cwd (str): Current working directory for the project
            path_master (str): Path to external dependencies directory
        """
        self.path_cwd = path_cwd
        self.path_master = path_master

    def extract_atom_types(self, model) -> Set[str]:
        """Extract unique GAFF atom types used in monomer files.

        Parses monomer .lt files to find atom type references (e.g., @atom:c3, @atom:ce)
        and returns the set of unique type names.

        Args:
            model: Model object containing sequence information with monomer file names

        Returns:
            Set[str]: Set of atom type names (e.g., {'c3', 'ce', 'hc', 'o'})
        """
        atom_types = set()

        for modelii in model:
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

    def parse_gaff_lt_sections(self, gaff_file: str) -> Dict[str, List[str]]:
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

    # Helper methods for extracting names and filtering sections
    @staticmethod
    def _extract_reference_name(line: str, prefix: str) -> Optional[str]:
        """
        Extract a reference name from a line (e.g., @bond:c3-ce -> c3-ce).

        Args:
            line: The line to extract from
            prefix: The prefix to look for (e.g., '@bond:', '@angle:')

        Returns:
            The extracted name or None
        """
        for part in line.split():
            if part.startswith(prefix):
                return part.split(':')[1]
        return None

    @staticmethod
    def _extract_atom_types(line: str) -> List[str]:
        """
        Extract atom types from a line.

        Args:
            line: The line to parse

        Returns:
            List of atom type strings
        """
        atom_types = []
        for part in line.split():
            if part.startswith('@atom:'):
                atom_types.append(part.split(':')[1])
        return atom_types

    @staticmethod
    def _should_keep_header_line(line: str) -> bool:
        """
        Determine if a header/footer line should be kept.

        Args:
            line: The line to check

        Returns:
            True if the line should be kept
        """
        stripped = line.strip()
        if not stripped:
            return True
        if stripped.startswith('}'):
            return False
        if stripped.startswith('#') and 'end of' in stripped.lower():
            return False
        return True

    def _filter_section_by_atom_types(
        self,
        lines: List[str],
        atom_types: Set[str],
        keyword: str
    ) -> List[str]:
        """
        Filter a section keeping only lines with specified atom types.

        Args:
            lines: Lines to filter
            atom_types: Set of atom type names to keep
            keyword: Keyword that must be present (e.g., 'pair_coeff')

        Returns:
            Filtered lines
        """
        filtered = []
        for line in lines:
            if keyword in line and '@atom:' in line:
                atom_type = self._extract_reference_name(line, '@atom:')
                if atom_type in atom_types:
                    filtered.append(line)
            elif self._should_keep_header_line(line):
                filtered.append(line)
        return filtered

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
                atom_type = self._extract_reference_name(line, '@atom:')
                if atom_type in atom_types:
                    filtered.append(line)
            elif self._should_keep_header_line(line):
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
        return self._filter_section_by_atom_types(pair_lines, atom_types, 'pair_coeff')

    def _filter_section_by_topology(
        self,
        coeff_lines: List[str],
        def_lines: List[str],
        used_types: Set[str],
        reference_prefix: str,
        normalize_func: callable = None
    ) -> Tuple[List[str], List[str]]:
        """
        Generic filter for coefficients and definitions based on topology types.

        Args:
            coeff_lines: Lines from coefficient section
            def_lines: Lines from definition section
            used_types: Set of used topology type strings
            reference_prefix: Reference prefix (e.g., '@bond:', '@angle:')
            normalize_func: Optional function to normalize type strings for comparison

        Returns:
            Tuple of (filtered_coeff_lines, filtered_def_lines)
        """
        # Identify which references to keep
        keep_refs = set()

        for line in def_lines:
            if reference_prefix in line and '@atom:' in line:
                atom_types = self._extract_atom_types(line)
                type_str = normalize_func(atom_types) if normalize_func else '-'.join(atom_types)

                if type_str in used_types:
                    ref_name = self._extract_reference_name(line, reference_prefix)
                    if ref_name:
                        keep_refs.add(ref_name)

        # Filter using the common filtering logic
        filtered_coeffs = self._filter_lines_by_reference(coeff_lines, keep_refs, reference_prefix)
        filtered_defs = self._filter_lines_by_reference(def_lines, keep_refs, reference_prefix)

        return filtered_coeffs, filtered_defs

    def _normalize_bond_type(self, atom_types: List[str]) -> str:
        """Normalize bond type by sorting atom types alphabetically."""
        if len(atom_types) == 2:
            return '-'.join(sorted(atom_types))
        return '-'.join(atom_types)

    def _normalize_angle_type(self, atom_types: List[str]) -> str:
        """Normalize angle type with sorted outer atoms."""
        if len(atom_types) == 3:
            center = atom_types[1]
            outer1, outer2 = atom_types[0], atom_types[2]
            if outer1 > outer2:
                outer1, outer2 = outer2, outer1
            return f"{outer1}-{center}-{outer2}"
        return '-'.join(atom_types)

    def _normalize_dihedral_type(self, atom_types: List[str]) -> str:
        """Normalize dihedral type (order matters for dihedrals)."""
        return '-'.join(atom_types)

    def _filter_lines_by_reference(
        self,
        lines: List[str],
        keep_refs: Set[str],
        reference_prefix: str
    ) -> List[str]:
        """
        Filter lines keeping only those with references in keep_refs.

        Args:
            lines: Lines to filter
            keep_refs: Set of reference names to keep
            reference_prefix: Reference prefix to check (e.g., '@bond:')

        Returns:
            Filtered lines
        """
        filtered = []
        for line in lines:
            if reference_prefix in line:
                ref_name = self._extract_reference_name(line, reference_prefix)
                if ref_name and ref_name in keep_refs:
                    filtered.append(line)
            elif self._should_keep_header_line(line):
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
        return self._filter_section_by_topology(
            bond_coeff_lines, bond_def_lines, used_bond_types,
            '@bond:', self._normalize_bond_type
        )

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
        return self._filter_section_by_topology(
            angle_coeff_lines, angle_def_lines, used_angle_types,
            '@angle:', self._normalize_angle_type
        )

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
        return self._filter_section_by_topology(
            dihedral_coeff_lines, dihedral_def_lines, used_dihedral_types,
            '@dihedral:', self._normalize_dihedral_type
        )

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
        return self._filter_section_by_topology(
            improper_coeff_lines, improper_def_lines, used_improper_types,
            '@improper:', self._normalize_dihedral_type  # Same normalization as dihedral
        )

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

    def create_gaff_subset(self, model) -> None:
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

        Args:
            model: Model object containing sequence information with monomer file names
        """
        try:
            gaff_src = str(Path(self.path_master) / "moltemplate" / "common" / "gaff.lt")
            gaff_dst = str(Path(self.path_cwd) / "gaff_subset.lt")

            # Check if source file exists
            if not Path(gaff_src).exists():
                logger.error(f"GAFF force field file not found: {gaff_src}")
                logger.error("Please ensure gaff.lt is installed in moltemplate/common/")
                sys.exit(1)

            logger.info("Creating GAFF parameter subset...")

            # Step 1: Collect monomer file paths
            logger.info("  Collecting monomer files...")
            monomer_files = []
            for modelii in model:
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
                shutil.copy(gaff_src, str(Path(self.path_cwd) / "gaff.lt"))
                return

            logger.info(f"  Found {len(monomer_files)} monomer files")

            # Step 2: Extract atom types (still needed for masses and pair coeffs)
            logger.info("  Extracting atom types from monomers...")
            atom_types = self.extract_atom_types(model)
            logger.info(f"  Found {len(atom_types)} unique atom types: {sorted(atom_types)}")

            if not atom_types:
                logger.warning("  No atom types found in monomers!")
                logger.warning("  Falling back to full gaff.lt")
                shutil.copy(gaff_src, str(Path(self.path_cwd) / "gaff.lt"))
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
            sections = self.parse_gaff_lt_sections(gaff_src)

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
            gaff_link = str(Path(self.path_cwd) / "gaff.lt")
            if Path(gaff_link).exists():
                Path(gaff_link).unlink()
            Path(gaff_link).symlink_to("gaff_subset.lt")
            logger.info(f"  Created symlink: gaff.lt -> gaff_subset.lt")

        except Exception as e:
            logger.error(f"Error in create_gaff_subset: {str(e)}")
            import traceback
            traceback.print_exc()
            logger.warning("Falling back to full gaff.lt")
            shutil.copy(gaff_src, str(Path(self.path_cwd) / "gaff.lt"))
