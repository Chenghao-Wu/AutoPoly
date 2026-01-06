"""
Atom Type Mappings Module for AutoPoly

Provides atom type change mappings for connection point modifications
in OPLS-AA and GAFF force fields.

When creating polymer end-caps, the atom type at the connection point
may need to change (e.g., CH3 → CH2 for carbon, OH → O for oxygen).

This module defines these mappings for common elements (C, O, N).

Author: AutoPoly Development Team
"""

from typing import Dict, Tuple, Optional


# OPLS-AA atom type mappings
# Format: (element, original_type_name, connection_type_name): (original_atom_type, connection_atom_type)
# If connection_atom_type is None, the atom type doesn't change

OPLS_CONNECTION_TYPE_MAP: Dict[Tuple[str, str, str], Tuple[Optional[str], Optional[str]]] = {
    # Carbon mappings
    ('C', 'CH3', 'CH2'): ('@atom:80', '@atom:82'),     # Methyl → Methylene
    ('C', 'CH2', 'CH'): ('@atom:82', '@atom:83'),       # Methylene → Methine
    ('C', 'CH', 'C'): ('@atom:83', '@atom:84'),         # Methine → Quaternary carbon
    # Aromatic carbons don't typically change at connection points
    ('C', 'Car', 'Car'): (None, None),                   # Aromatic carbon (no change)

    # Oxygen mappings
    ('O', 'OH', 'O'): ('@atom:96', '@atom:122'),        # Alcohol → Ether oxygen
    ('O', 'O', 'O'): ('@atom:122', '@atom:122'),        # Ether oxygen (no change)
    ('O', 'O=C', 'O=C'): (None, None),                   # Carbonyl oxygen (no change)
    ('O', 'O-C=O', 'O-C=O'): (None, None),               # Ester oxygen (no change)

    # Nitrogen mappings
    ('N', 'NH2', 'NH'): ('@atom:739', '@atom:740'),     # Primary amine → Secondary amine
    ('N', 'NH', 'N'): ('@atom:740', '@atom:741'),       # Secondary → Tertiary amine
    ('N', 'N', 'N'): (None, None),                       # Tertiary amine (no change)
    ('N', 'NH-C=O', 'N-C=O'): (None, None),              # Amide nitrogen (no change)
}


# GAFF atom type mappings
# GAFF uses different atom type names (ca, ha, oa, oh, etc.)
GAFF_CONNECTION_TYPE_MAP: Dict[Tuple[str, str, str], Tuple[Optional[str], Optional[str]]] = {
    # Carbon mappings
    ('C', 'ca', 'ca'): (None, None),                    # Aromatic carbon (no change)
    ('C', 'c3', 'c3'): (None, None),                    # sp3 carbon (no change for most cases)
    ('C', 'ct', 'ct'): (None, None),                    # sp3 carbon alkane (no change)

    # Oxygen mappings
    ('O', 'oh', 'os'): ('@atom:oh', '@atom:os'),        # Alcohol → Ether oxygen in GAFF
    ('O', 'os', 'os'): (None, None),                    # Ether oxygen (no change)
    ('O', 'o', 'o'): (None, None),                      # Carbonyl oxygen (no change)

    # Nitrogen mappings
    ('N', 'n', 'n'): (None, None),                      # Amine nitrogen (no change)
    ('N', 'na', 'na'): (None, None),                    # Aromatic nitrogen (no change)
    ('N', 'nh', 'n'): ('@atom:nh', '@atom:n'),          # Primary → Secondary amine in GAFF
}


def get_opls_connection_type(
    element: str,
    original_type: str,
    connection_type: str
) -> Tuple[Optional[str], Optional[str]]:
    """
    Get OPLS atom type change for a connection point modification.

    Args:
        element: Element symbol (e.g., 'C', 'O', 'N')
        original_type: Original atom type name (e.g., 'CH3', 'OH')
        connection_type: Connection atom type name (e.g., 'CH2', 'O')

    Returns:
        (original_atom_type, connection_atom_type) tuple
        Returns (None, None) if no change needed

    Example:
        >>> get_opls_connection_type('C', 'CH3', 'CH2')
        ('@atom:80', '@atom:82')
    """
    key = (element, original_type, connection_type)
    return OPLS_CONNECTION_TYPE_MAP.get(key, (None, None))


def get_gaff_connection_type(
    element: str,
    original_type: str,
    connection_type: str
) -> Tuple[Optional[str], Optional[str]]:
    """
    Get GAFF atom type change for a connection point modification.

    Args:
        element: Element symbol (e.g., 'C', 'O', 'N')
        original_type: Original atom type name (e.g., 'ct', 'oh')
        connection_type: Connection atom type name (e.g., 'c3', 'os')

    Returns:
        (original_atom_type, connection_atom_type) tuple
        Returns (None, None) if no change needed

    Example:
        >>> get_gaff_connection_type('O', 'oh', 'os')
        ('@atom:oh', '@atom:os')
    """
    key = (element, original_type, connection_type)
    return GAFF_CONNECTION_TYPE_MAP.get(key, (None, None))


def get_connection_type(
    element: str,
    original_type: str,
    connection_type: str,
    force_field: str = 'oplsaa'
) -> Tuple[Optional[str], Optional[str]]:
    """
    Get atom type change for a connection point modification.

    Generic function that dispatches to OPLS or GAFF mappings.

    Args:
        element: Element symbol (e.g., 'C', 'O', 'N')
        original_type: Original atom type name
        connection_type: Connection atom type name
        force_field: 'oplsaa', 'lopls', or 'gaff'

    Returns:
        (original_atom_type, connection_atom_type) tuple
        Returns (None, None) if no change needed

    Raises:
        ValueError: If force_field not supported

    Example:
        >>> get_connection_type('C', 'CH3', 'CH2', force_field='oplsaa')
        ('@atom:80', '@atom:82')
    """
    force_field_lower = force_field.lower()

    if 'gaff' in force_field_lower:
        return get_gaff_connection_type(element, original_type, connection_type)
    elif 'opls' in force_field_lower:
        return get_opls_connection_type(element, original_type, connection_type)
    else:
        raise ValueError(
            f"Unsupported force field: '{force_field}'. "
            f"Use 'oplsaa', 'lopls', or 'gaff'."
        )


def needs_type_change(
    element: str,
    original_type: str,
    connection_type: str,
    force_field: str = 'oplsaa'
) -> bool:
    """
    Check if an atom needs a type change at a connection point.

    Args:
        element: Element symbol
        original_type: Original atom type name
        connection_type: Connection atom type name
        force_field: Force field name

    Returns:
        True if atom type should change, False otherwise
    """
    original_type_str, connection_type_str = get_connection_type(
        element, original_type, connection_type, force_field
    )
    return original_type_str is not None and connection_type_str is not None
