"""
Polymerization Patterns Module for AutoPoly

Defines SMARTS patterns for different polymerization mechanisms.
This enables element-agnostic monomer generation for vinyl addition,
condensation, and step-growth polymerization.

Mechanisms:
- none: Non-polymerizable molecules (single molecules, solvents)
- vinyl_addition: C=C double bond opening (e.g., PE, PP, PS, PMMA)
- esterification: Carboxyl + alcohol → polyester (e.g., PLA, PET)
- amidation: Carboxyl + amine → polyamide (e.g., Nylon-6,6)
- etherification: Alcohol + alcohol → polyether (e.g., PEG)

Author: AutoPoly Development Team
"""

from typing import Dict, List, Optional


class PolymerizationPattern:
    """Represents a single polymerization mechanism pattern."""

    def __init__(
        self,
        name: str,
        smarts: Optional[str],
        connection_atoms: List[str],
        description: str,
        is_condensation: bool = False,
        requires_two_groups: bool = False
    ):
        """
        Initialize a polymerization pattern.

        Args:
            name: Mechanism name (e.g., 'vinyl_addition', 'esterification')
            smarts: SMARTS pattern to identify reactive sites (None for non-polymerizable)
            connection_atoms: Element symbols of atoms that form polymer bonds
            description: Human-readable description
            is_condensation: True if condensation polymerization (releases small molecule)
            requires_two_groups: True if monomer needs 2+ functional groups
        """
        self.name = name
        self.smarts = smarts
        self.connection_atoms = connection_atoms
        self.description = description
        self.is_condensation = is_condensation
        self.requires_two_groups = requires_two_groups

    def __repr__(self):
        return f"PolymerizationPattern({self.name})"


# Polymerization pattern definitions
POLYMERIZATION_PATTERNS: Dict[str, PolymerizationPattern] = {
    'none': PolymerizationPattern(
        name='none',
        smarts=None,
        connection_atoms=[],
        description='Non-polymerizable molecule (single molecules, solvents, non-reactive species)',
        is_condensation=False,
        requires_two_groups=False
    ),

    'vinyl_addition': PolymerizationPattern(
        name='vinyl_addition',
        smarts='[$([C]=[C])]',  # Carbon-carbon double bond
        connection_atoms=['C', 'C'],
        description='Vinyl addition polymerization (C=C double bond opening)',
        is_condensation=False,
        requires_two_groups=False
    ),

    'esterification': PolymerizationPattern(
        name='esterification',
        # Carboxyl group (C=O)OH + alcohol group
        smarts='[$([C](=[O])[OX2H1])][$([OX2H])]',  # Carboxyl carbon bonded to alcohol oxygen
        connection_atoms=['C', 'O'],
        description='Esterification polymerization (carboxyl + alcohol → polyester)',
        is_condensation=True,
        requires_two_groups=True
    ),

    'amidation': PolymerizationPattern(
        name='amidation',
        # Carboxyl group (C=O)OH + amine group (detected separately)
        smarts='[$([CX3](=[OX1])[OX2H1])]',  # Just match carboxyl group
        connection_atoms=['C', 'N'],
        description='Amidation polymerization (carboxyl + amine → polyamide)',
        is_condensation=True,
        requires_two_groups=True
    ),

    'etherification': PolymerizationPattern(
        name='etherification',
        # Two alcohol groups
        smarts='[$([OX2H])]',  # Alcohol oxygen
        connection_atoms=['O', 'O'],
        description='Etherification polymerization (alcohol + alcohol → polyether)',
        is_condensation=True,
        requires_two_groups=True
    ),
}


def get_pattern(mechanism: str) -> PolymerizationPattern:
    """
    Get a polymerization pattern by mechanism name.

    Args:
        mechanism: Mechanism name (e.g., 'vinyl_addition', 'esterification')

    Returns:
        PolymerizationPattern: Pattern definition

    Raises:
        ValueError: If mechanism not found
    """
    if mechanism not in POLYMERIZATION_PATTERNS:
        raise ValueError(
            f"Unknown mechanism: '{mechanism}'. "
            f"Available: {list(POLYMERIZATION_PATTERNS.keys())}"
        )
    return POLYMERIZATION_PATTERNS[mechanism]


def get_all_mechanisms() -> List[str]:
    """Get list of all available mechanism names."""
    return list(POLYMERIZATION_PATTERNS.keys())


def get_polymerization_mechanisms() -> List[str]:
    """Get list of polymer-capable mechanisms (excludes 'none')."""
    return [name for name, pattern in POLYMERIZATION_PATTERNS.items()
            if name != 'none']


def get_condensation_mechanisms() -> List[str]:
    """Get list of condensation polymerization mechanisms."""
    return [name for name, pattern in POLYMERIZATION_PATTERNS.items()
            if pattern.is_condensation]


def validate_mechanism(mechanism: str) -> bool:
    """Check if a mechanism name is valid."""
    return mechanism in POLYMERIZATION_PATTERNS


# Aliases for convenience
PATTERN_NONE = POLYMERIZATION_PATTERNS['none']
PATTERN_VINYL = POLYMERIZATION_PATTERNS['vinyl_addition']
PATTERN_ESTER = POLYMERIZATION_PATTERNS['esterification']
PATTERN_AMIDE = POLYMERIZATION_PATTERNS['amidation']
PATTERN_ETHER = POLYMERIZATION_PATTERNS['etherification']
