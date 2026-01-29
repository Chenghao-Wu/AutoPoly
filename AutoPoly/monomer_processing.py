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
from typing import Tuple, List, Dict, Optional

from .system import logger
from .monomer_generator import MonomerGenerator
from .exceptions import GenerationError


def _create_generator(base_name: str, force_field: str, output_dir: str) -> MonomerGenerator:
    """
    Create a MonomerGenerator instance with standard parameters.

    Args:
        base_name: Base name for generated files
        force_field: Force field to use ("oplsaa", "gaff", or "lopls")
        output_dir: Output directory path

    Returns:
        Configured MonomerGenerator instance
    """
    return MonomerGenerator(
        base_name=base_name,
        force_field=force_field,
        output_dir=output_dir,
        verbose=False
    )


def _check_cache(cache_key, cache: dict, counter: int) -> Optional[Tuple[str, int]]:
    """
    Check if a result is cached and return it if available.

    Args:
        cache_key: Key to look up in cache
        cache: Cache dictionary
        counter: Current counter value

    Returns:
        (cached_value, counter) if found, None otherwise
    """
    if cache_key in cache:
        logger.info(f"Using cached result for: {cache_key}")
        return cache[cache_key], counter
    return None


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
        psmiles: pSMILES string (e.g., "[*]C=C[*]") or SMILES string (e.g., "CC(=O)C")
        path_cwd: Current working directory path
        force_field: Force field to use ("oplsaa", "gaff", or "lopls")
        generated_cache: Cache mapping pSMILES to monomer names (modified in-place)
        counter: Current counter value for unique monomer naming

    Returns:
        Tuple[str, int]: (monomer_name, updated_counter)

    Raises:
        GenerationError: If monomer generation fails
    """
    cached = _check_cache(psmiles, generated_cache, counter)
    if cached:
        return cached

    monomer_name = f"monomer_{counter}"
    counter += 1

    try:
        generator = _create_generator(monomer_name, force_field, path_cwd)
        variants = generator.from_smiles(smiles=psmiles, n_monomers=3)
        generator.write_lt_files(variants, generate_t1=True)

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
    benzene, ethanol) from a regular SMILES string WITHOUT wildcards.

    Args:
        smiles: SMILES string WITHOUT wildcards (e.g., "O", "CCO", "c1ccccc1")
        molecule_name: Name for the molecule (e.g., "water", "ethanol", "benzene")
        path_cwd: Current working directory path
        force_field: Force field to use ("oplsaa", "gaff", or "lopls")
        generated_cache: Cache mapping SMILES to molecule names (modified in-place)
        counter: Current counter value for unique molecule naming

    Returns:
        Tuple[str, int]: (molecule_filename, updated_counter)

    Raises:
        GenerationError: If molecule generation fails

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
    cache_key = f"molecule_{smiles}_{force_field}"
    cached = _check_cache(cache_key, generated_cache, counter)
    if cached:
        return cached

    filename = f"{molecule_name}.lt"

    try:
        generator = _create_generator(molecule_name, force_field, path_cwd)
        variant = generator.from_single_molecule(smiles=smiles, molecule_name=molecule_name)
        generator.write_single_molecule(variant, generate_t1=False)

        generated_cache[cache_key] = filename
        logger.info(f"Generated molecule '{filename}' from SMILES: {smiles}")
        return filename, counter

    except Exception as e:
        raise GenerationError(
            f"Failed to generate molecule from SMILES '{smiles}': {e}"
        ) from e


def generate_sequence_variants_for_polymerization(
    smiles_list: List[str],
    topology: str,
    path_cwd: str,
    force_field: str,
    generated_cache: dict,
    counter: int,
    base_name_prefix: str = "monomer"
) -> Tuple[Dict[str, str], int]:
    """
    Generate monomer variants for polymerization using complement SMILES.

    This function creates position-aware variants for each unique SMILES in the
    complement SMILES list.

    Args:
        smiles_list: List of complement SMILES:
            - First: 1 wildcard (right connection) e.g., 'CC[*]'
            - Middle: 2 wildcards (left and right) e.g., '[*]CC[*]'
            - Last: 1 wildcard (left connection) e.g., '[*]CC'
        topology: Topology type ("linear" or "ring")
        path_cwd: Current working directory path
        force_field: Force field to use ("oplsaa", "gaff", or "lopls")
        generated_cache: Cache mapping (smiles_list, topology) to monomer names
        counter: Current counter value for unique monomer naming
        base_name_prefix: Prefix for monomer names (default: "monomer")

    Returns:
        Tuple[Dict[str, str], int]: (variant_name_mapping, updated_counter)

    Raises:
        GenerationError: If monomer generation fails
    """
    cache_key = (tuple(smiles_list), topology)
    cached = _check_cache(cache_key, generated_cache, counter)
    if cached:
        return cached

    base_name = f"{base_name_prefix}_{counter}"
    counter += 1

    try:
        generator = _create_generator(base_name, force_field, path_cwd)

        logger.info(f"Generating {len(smiles_list)} sequence variants for {base_name} ({topology} topology)")

        variants = generator.from_smiles(smiles_list)

        # Extract unique variants by variant_type
        unique_variants = {}
        for variant in variants:
            vtype = variant.variant_type
            if vtype not in unique_variants:
                unique_variants[vtype] = variant

        # For ring topology, use middle variants for all positions
        if topology == "ring":
            if 'middle' in unique_variants:
                unique_variants = {'ring': unique_variants['middle']}
            elif 'first' in unique_variants:
                unique_variants = {'ring': unique_variants['first']}

        logger.info(f"Extracted {len(unique_variants)} unique variant types: {list(unique_variants.keys())}")

        lt_files = generator.write_lt_files(list(unique_variants.values()), generate_t1=True)
        logger.info(f"Generated {len(lt_files)} .lt files for {base_name}")

        variant_name_mapping = _build_variant_mapping(unique_variants, base_name)
        generated_cache[cache_key] = variant_name_mapping

        logger.info(f"Generated sequence variants for '{base_name}' from {len(smiles_list)} complement SMILES")
        logger.info(f"Variant mapping keys: {list(variant_name_mapping.keys())}")
        return variant_name_mapping, counter

    except Exception as e:
        raise GenerationError(
            f"Failed to generate sequence variants from complement SMILES: {e}"
        ) from e


def _build_variant_mapping(unique_variants: dict, base_name: str) -> Dict[str, str]:
    """
    Build mapping from variant types to .lt filenames.

    Args:
        unique_variants: Dictionary of variant_type -> variant objects
        base_name: Base name for the files

    Returns:
        Dictionary mapping variant_type keys to filenames
    """
    variant_name_mapping = {}

    for variant in unique_variants.values():
        vtype = variant.variant_type
        pos = variant.position

        # Determine filename suffix based on variant_type
        suffix_map = {
            'first': 'le',
            'last': 're',
            'single': 'single',
            'ring': 'i',
            'middle': 'i'
        }
        suffix = suffix_map.get(vtype, 'i')

        filename = f"{base_name}_{pos}{suffix}.lt"
        filename_t1 = f"{base_name}_{pos}{suffix}_T1.lt"

        variant_name_mapping[vtype] = filename
        variant_name_mapping[f"{vtype}_standard"] = filename
        variant_name_mapping[f"{vtype}_T1"] = filename_t1

    return variant_name_mapping


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
