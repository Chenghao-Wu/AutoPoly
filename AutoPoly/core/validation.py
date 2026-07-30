#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validation Module for AutoPoly Package

This module provides validation functions for SMILES strings and other
user inputs to prevent injection attacks and ensure data integrity.

Created on 2026-01-10
@author: zwu
"""
from rdkit import Chem
from .exceptions import ValidationError


def validate_smiles(smiles: str, allow_wildcards: bool = True) -> bool:
    """
    Validate a SMILES string and check for wildcards if required.

    This function validates that a SMILES string is chemically valid and
    optionally checks for the presence of wildcard connection points ([*]).
    Wildcards are required for polymer monomers but not for standalone molecules.

    Args:
        smiles (str): SMILES string to validate
        allow_wildcards (bool): Whether wildcards are allowed. If True, requires
                               wildcards to be present. If False, requires wildcards
                               to be absent. Defaults to True.

    Returns:
        bool: True if validation passes

    Raises:
        ValidationError: If SMILES is invalid or wildcard check fails

    Examples:
        >>> # Validate polymer monomer SMILES (requires wildcards)
        >>> validate_smiles("[*]CC[*]")
        True

        >>> # Validate molecule SMILES (no wildcards)
        >>> validate_smiles("CCO", allow_wildcards=False)
        True

        >>> # This will raise ValidationError
        >>> validate_smiles("invalid_smiles")
        ValidationError: Invalid SMILES: invalid_smiles
    """
    if not smiles or not isinstance(smiles, str):
        raise ValidationError(
            f"SMILES must be a non-empty string, got: {type(smiles).__name__}"
        )

    # Try to parse the SMILES string
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValidationError(f"Invalid SMILES: {smiles}")
    except Exception as e:
        raise ValidationError(f"Failed to parse SMILES '{smiles}': {e}") from e

    # Check for wildcards if required
    has_wildcards = "[*]" in smiles

    if allow_wildcards and not has_wildcards:
        raise ValidationError(
            f"SMILES must contain [*] connection points for polymer monomers: {smiles}"
        )

    if not allow_wildcards and has_wildcards:
        raise ValidationError(
            f"SMILES for standalone molecules should not contain [*] wildcards: {smiles}"
        )

    return True


def validate_smiles_list(smiles_list: list, allow_wildcards: bool = True) -> bool:
    """
    Validate a list of SMILES strings.

    Args:
        smiles_list (list): List of SMILES strings to validate
        allow_wildcards (bool): Whether wildcards are allowed. Defaults to True.

    Returns:
        bool: True if all SMILES strings are valid

    Raises:
        ValidationError: If any SMILES string is invalid

    Examples:
        >>> validate_smiles_list(["[*]CC[*]", "[*]C=C[*]"])
        True
    """
    if not smiles_list or not isinstance(smiles_list, list):
        raise ValidationError(
            f"SMILES list must be a non-empty list, got: {type(smiles_list).__name__}"
        )

    for i, smiles in enumerate(smiles_list):
        try:
            validate_smiles(smiles, allow_wildcards=allow_wildcards)
        except ValidationError as e:
            raise ValidationError(
                f"SMILES at position {i} is invalid: {e}"
            ) from e

    return True


def validate_sequence_length(sequence_length: int, max_length: int) -> bool:
    """
    Validate that sequence length is within acceptable limits.

    Args:
        sequence_length (int): Length of the sequence to validate
        max_length (int): Maximum allowed sequence length

    Returns:
        bool: True if validation passes

    Raises:
        ValidationError: If sequence length exceeds maximum
    """
    if sequence_length > max_length:
        raise ValidationError(
            f"Sequence length ({sequence_length}) exceeds maximum {max_length}"
        )

    return True


def validate_unique_monomer_count(unique_count: int, max_count: int) -> bool:
    """
    Validate that the number of unique monomers is within acceptable limits.

    Args:
        unique_count (int): Number of unique monomers
        max_count (int): Maximum allowed unique monomers

    Returns:
        bool: True if validation passes

    Raises:
        ValidationError: If unique monomer count exceeds maximum
    """
    if unique_count > max_count:
        raise ValidationError(
            f"Number of unique monomers ({unique_count}) exceeds maximum {max_count}"
        )

    return True
