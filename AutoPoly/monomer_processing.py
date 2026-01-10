#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Monomer Processing Module for AutoPoly Package

This module provides standalone functions for processing monomer data
including generating monomers from pSMILES/SMILES strings, reading monomer
files, counting atoms, and evaluating offset distances.

Created on 2026-01-06
@author: zwu
"""
import re
import numpy as np
from pathlib import Path
from typing import Tuple, List, Dict

from .system import logger
from .monomer_generator import MonomerGenerator
from .exceptions import GenerationError


def generate_monomer_from_psmiles(
    psmiles: str,
    path_cwd: str,
    force_field: str,
    generated_cache: dict,
    counter: int
) -> Tuple[str, int]:
    """
    Generate all 6 .lt files for a monomer from pSMILES or SMILES string.

    This function creates a complete set of monomer variant files (internal,
    left-end, right-end, with and without T1 chirality) from a pSMILES (for polymers)
    or SMILES (for molecules) string.

    Args:
        psmiles (str): pSMILES string (e.g., "[*]C=C[*]") or SMILES string (e.g., "CC(=O)C")
        path_cwd (str): Current working directory path
        force_field (str): Force field to use ("oplsaa", "gaff", or "lopls")
        generated_cache (dict): Cache mapping pSMILES to monomer names (modified in-place)
        counter (int): Current counter value for unique monomer naming

    Returns:
        Tuple[str, int]: (monomer_name, updated_counter) where:
            - monomer_name: Generated monomer base name (e.g., "monomer_0")
            - updated_counter: Incremented counter value

    Raises:
        SystemExit: If monomer generation fails
    """
    # Check cache (per-system caching)
    if psmiles in generated_cache:
        logger.info(f"Using cached monomer for pSMILES/SMILES: {psmiles}")
        return generated_cache[psmiles], counter

    # Generate unique monomer name
    monomer_name = f"monomer_{counter}"
    counter += 1

    try:
        # Create MonomerGenerator with correct API
        generator = MonomerGenerator(
            base_name=monomer_name,
            force_field=force_field,
            output_dir=path_cwd,
            verbose=False
        )

        # Generate variants from SMILES using from_smiles()
        # This builds a chain, assigns atom types, and splits into variants
        variants = generator.from_smiles(smiles=psmiles, n_monomers=3)

        # Write .lt files using write_lt_files()
        generator.write_lt_files(variants, generate_t1=True)

        # Cache the mapping
        generated_cache[psmiles] = monomer_name

        logger.info(f"Generated monomer '{monomer_name}' from pSMILES/SMILES: {psmiles}")
        return monomer_name, counter

    except Exception as e:
        raise GenerationError(
            f"Failed to generate monomer from pSMILES/SMILES '{psmiles}': {e}"
        ) from e


def generate_molecule_from_smiles(
    smiles: str,
    molecule_name: str,
    path_cwd: str,
    force_field: str,
    generated_cache: dict,
    counter: int
) -> Tuple[str, int]:
    """
    Generate .lt file for a single molecule from SMILES string.

    This function creates a single .lt file for a non-polymer molecule (e.g., water,
    benzene, ethanol) from a regular SMILES string WITHOUT wildcards. This is different
    from generate_monomer_from_psmiles() which is for polymer monomers with wildcards.

    Args:
        smiles (str): SMILES string WITHOUT wildcards (e.g., "O", "CCO", "c1ccccc1")
        molecule_name (str): Name for the molecule (e.g., "water", "ethanol", "benzene")
        path_cwd (str): Current working directory path
        force_field (str): Force field to use ("oplsaa", "gaff", or "lopls")
        generated_cache (dict): Cache mapping SMILES to molecule names (modified in-place)
        counter (int): Current counter value for unique molecule naming

    Returns:
        Tuple[str, int]: (molecule_filename, updated_counter) where:
            - molecule_filename: Generated molecule .lt filename (e.g., "water.lt")
            - updated_counter: Incremented counter value

    Raises:
        SystemExit: If molecule generation fails

    Example:
        >>> molecule_name, counter = generate_molecule_from_smiles(
        ...     smiles="O",
        ...     molecule_name="water",
        ...     path_cwd="./monomers",
        ...     force_field="gaff",
        ...     generated_cache={},
        ...     counter=0
        ... )
    """
    # Check cache (per-system caching)
    cache_key = f"molecule_{smiles}_{force_field}"
    if cache_key in generated_cache:
        logger.info(f"Using cached molecule for SMILES: {smiles}")
        return generated_cache[cache_key], counter

    # Use provided molecule_name for the .lt file
    base_name = molecule_name
    filename = f"{base_name}.lt"

    try:
        # Create MonomerGenerator with correct API
        generator = MonomerGenerator(
            base_name=base_name,
            force_field=force_field,
            output_dir=path_cwd,
            verbose=False
        )

        # Generate molecule variant from SMILES (no wildcards)
        # This uses from_single_molecule() which is designed for non-polymer molecules
        variant = generator.from_single_molecule(smiles=smiles, molecule_name=base_name)

        # Write .lt file using write_single_molecule()
        # For molecules, we typically don't need T1 variants (no tacticity)
        generator.write_single_molecule(variant, generate_t1=False)

        # Cache the mapping
        generated_cache[cache_key] = filename

        logger.info(f"Generated molecule '{filename}' from SMILES: {smiles}")
        return filename, counter

    except Exception as e:
        raise GenerationError(
            f"Failed to generate molecule from SMILES '{smiles}': {e}"
        ) from e


def generate_sequence_variants_for_polymerization(
    base_smiles: str,
    dop: int,
    topology: str,
    path_cwd: str,
    force_field: str,
    generated_cache: dict,
    counter: int,
    base_name_prefix: str = "monomer"
) -> Tuple[Dict[str, str], int]:
    """
    Generate monomer variants for polymerization using the MonomerGenerator API.

    This function creates position-aware variants for each monomer type in the
    polymer chain (first, middle, last for linear; all middle for ring).

    Args:
        base_smiles (str): Base monomer SMILES with wildcards (e.g., "[*]C=C[*]")
        dop (int): Degree of polymerization (number of monomers in chain)
        topology (str): Topology type ("linear" or "ring")
        path_cwd (str): Current working directory path
        force_field (str): Force field to use ("oplsaa", "gaff", or "lopls")
        generated_cache (dict): Cache mapping (base_smiles, dop, topology) to monomer names
        counter (int): Current counter value for unique monomer naming
        base_name_prefix (str): Prefix for monomer names (default: "monomer")

    Returns:
        Tuple[Dict[str, str], int]: (variant_name_mapping, updated_counter) where:
            - variant_name_mapping: Dict mapping variant keys to .lt filenames
              {
                  'first': 'monomer_0_0le.lt',
                  'middle': 'monomer_0_1i.lt',
                  'last': 'monomer_0_2re.lt',
                  'first_T1': 'monomer_0_0le_T1.lt',
                  ...
              }
            - updated_counter: Incremented counter value

    Raises:
        SystemExit: If monomer generation fails
    """
    # Check cache (per-system caching)
    cache_key = (base_smiles, dop, topology)
    if cache_key in generated_cache:
        logger.info(f"Using cached sequence variants for SMILES: {base_smiles}, DOP={dop}, topology={topology}")
        return generated_cache[cache_key], counter

    # Generate unique base name for this polymer
    base_name = f"{base_name_prefix}_{counter}"
    counter += 1

    try:
        # Create MonomerGenerator with correct API
        generator = MonomerGenerator(
            base_name=base_name,
            force_field=force_field,
            output_dir=path_cwd,
            verbose=False
        )

        # Generate variants using from_smiles()
        # Always use 3 monomers to get first/middle/last variants
        # This is a performance optimization - we only need these 3 positions
        n_monomers = 3
        logger.info(f"Generating 3 sequence variants (first, middle, last) for {base_name} ({topology} topology)")

        variants = generator.from_smiles(smiles=base_smiles, n_monomers=n_monomers)
        
        # For ring topology, we need all middle variants (both connections active)
        # For linear, we need first, middle, and last
        # The from_smiles already creates proper variant_types
        
        # Extract unique variants by variant_type
        # For linear: first (position 0), middle (position 1), last (position n-1)
        # For ring: all should be middle type (but from_smiles creates them as first/middle/last)
        unique_variants = {}
        for variant in variants:
            vtype = variant.variant_type
            # Keep first occurrence of each variant_type
            if vtype not in unique_variants:
                unique_variants[vtype] = variant
        
        # For ring topology, use middle variants for all positions
        if topology == "ring":
            if 'middle' in unique_variants:
                # Rename middle to 'ring' for clarity
                ring_variant = unique_variants['middle']
                unique_variants = {'ring': ring_variant}
            elif 'first' in unique_variants:
                # If no middle (DOP=2), use first as ring
                unique_variants = {'ring': unique_variants['first']}
        
        logger.info(f"Extracted {len(unique_variants)} unique variant types: {list(unique_variants.keys())}")

        # Generate .lt files for unique variants
        lt_files = generator.write_lt_files(list(unique_variants.values()), generate_t1=True)
        logger.info(f"Generated {len(lt_files)} .lt files for {base_name}")

        # Create mapping from variant_type to filename
        # The write_lt_files returns a list of file paths
        # File naming convention: {base_name}_{position}{suffix}.lt
        # where suffix is 'le' (first), 'i' (middle), 're' (last), 'single'
        variant_name_mapping = {}
        
        for variant in unique_variants.values():
            vtype = variant.variant_type
            pos = variant.position
            
            # Determine filename suffix based on variant_type
            if vtype == 'first':
                suffix = 'le'
            elif vtype == 'last':
                suffix = 're'
            elif vtype == 'single':
                suffix = 'single'
            elif vtype == 'ring':
                # Ring uses middle variant with 'i' suffix
                suffix = 'i'
            else:  # middle
                suffix = 'i'
            
            # Build filename (without _T1)
            filename = f"{base_name}_{pos}{suffix}.lt"
            filename_t1 = f"{base_name}_{pos}{suffix}_T1.lt"
            
            # Store mapping
            variant_name_mapping[vtype] = filename
            variant_name_mapping[f"{vtype}_standard"] = filename
            variant_name_mapping[f"{vtype}_T1"] = filename_t1

        # Cache the mapping
        generated_cache[cache_key] = variant_name_mapping

        logger.info(f"Generated sequence variants for '{base_name}' from SMILES: {base_smiles}")
        logger.info(f"Variant mapping keys: {list(variant_name_mapping.keys())}")
        return variant_name_mapping, counter

    except Exception as e:
        raise GenerationError(
            f"Failed to generate sequence variants from SMILES '{base_smiles}': {e}"
        ) from e


def n_monomer_atoms(merltfile: str, path_cwd: str) -> int:
    """
    Count the number of monomer atoms in the specified .lt file.

    This function parses a Moltemplate monomer file (.lt) and counts the
    number of atoms defined in the "Data Atoms" block. This information
    is used for polymer structure generation and validation.

    Args:
        merltfile (str): The name of the monomer .lt file.
        path_cwd (str): Current working directory path

    Returns:
        int: The number of monomer atoms.

    Raises:
        SystemExit: If the monomer file cannot be opened or found.
    """
    n_atoms = 0
    merltfile_path = Path(path_cwd) / merltfile

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
                    n_atoms = n_atoms + 1

                if not line:
                    break
    else:
        raise GenerationError(
            f"Monomer file not found: {merltfile_path}"
        )

    return n_atoms


def read_lt_end_atoms(lt_file: str) -> Tuple[str, str]:
    """
    Read the first and second atoms from a .lt file.

    Args:
        lt_file (str): Path to the .lt file

    Returns:
        tuple: (first_atom, second_atom) where each atom is a string with
               element name and position (e.g., "C1", "H2")

    Raises:
        ValueError: If both end atoms cannot be found in the file
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
                    element = extract_element_from_atom(parts[0])
                    atom_data = {
                        'element': element,
                        'atom_type': parts[2],
                        'x': float(parts[4]),
                        'y': float(parts[5]),
                        'z': float(parts[6])
                    }

                    if first_atom is None:
                        first_atom = atom_data['element']
                    elif second_atom is None:
                        second_atom = atom_data['element']
                        break  # We have both atoms, no need to continue

    if first_atom is None or second_atom is None:
        raise ValueError(f"Could not find both end atoms in {lt_file}")

    return first_atom + "1", second_atom + "2"


def extract_element_from_atom(atom_string: str) -> str:
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


def evaluate_offset(merltfile: str, path_cwd: str, offset_spacing: float, current_offset: float) -> float:
    """
    Evaluates the offset distance based on the specified merlt file.

    This function reads the first two atoms from the monomer file and calculates
    the distance between them, then adds the specified offset spacing to get
    the new offset value.

    Args:
        merltfile (str): The name of the merlt file.
        path_cwd (str): Current working directory path
        offset_spacing (float): Additional spacing to add to the calculated distance
        current_offset (float): Current offset value (not used in calculation but returned if file not found)

    Returns:
        float: The calculated offset distance (distance between first two atoms + offset_spacing)
    """
    monomer_bank = Path(path_cwd)
    merltfile_path = monomer_bank / merltfile

    if merltfile_path.is_file():
        C1 = []
        C2 = []

        with open(merltfile_path) as f:
            while True:
                line = f.readline()
                if line.strip() == "write(\"Data Atoms\") {":
                    # C1 coordinates
                    line = f.readline()
                    for i in range(3):
                        C1.append(float(line.split()[i + 4]))
                    # C2 coordinates
                    line = f.readline()
                    for i in range(3):
                        C2.append(float(line.split()[i + 4]))

                    # calculate C1-C2 distance
                    new_offset = np.linalg.norm(np.array(C1) - np.array(C2)) + offset_spacing

                    return new_offset

    return current_offset
