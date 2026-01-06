"""
Connection Point Modifier Module for AutoPoly

Provides generic connection point handling for polymer end-caps.
Handles atom type changes for different elements (C, O, N) in OPLS-AA and GAFF.

The ConnectionPointModifier class determines:
1. Which atom type should be used at a connection point
2. What the original atom type was
3. Element-agnostic type mapping

Author: AutoPoly Development Team
"""

from typing import Optional, Tuple
from rdkit import Chem

from .atom_type_mappings import get_connection_type
from .system import logger


class ConnectionPointModifier:
    """
    Handles connection point atom type modifications for polymer end-caps.

    This class provides element-agnostic methods to determine:
    - The original atom type at a connection site
    - The modified atom type for creating a polymer bond
    - Element-specific logic for C, O, N atoms

    Attributes:
        force_field: Force field name ('oplsaa', 'lopls', or 'gaff')
        verbose: Enable verbose logging

    Example:
        >>> modifier = ConnectionPointModifier(force_field='oplsaa')
        >>> original, connection = modifier.get_connection_atom_type(atom_carbon, mol)
        >>> print(f"{original} → {connection}")
        @atom:80 → @atom:82
    """

    # Atom type names for common patterns
    TYPE_NAMES = {
        'oplsaa': {
            'C': {
                'CH3': '@atom:80',   # Methyl
                'CH2': '@atom:82',   # Methylene
                'CH': '@atom:83',    # Methine
                'C': '@atom:84',     # Quaternary
                'Car': '@atom:145',  # Aromatic
            },
            'O': {
                'OH': '@atom:96',    # Alcohol
                'O': '@atom:122',    # Ether
                'O=C': '@atom:163',  # Carbonyl
                'O-C=O': '@atom:164' # Ester oxygen
            },
            'N': {
                'NH2': '@atom:739',  # Primary amine
                'NH': '@atom:740',   # Secondary amine
                'N': '@atom:741',    # Tertiary amine
            }
        },
        'gaff': {
            'C': {
                'ca': '@atom:ca',    # Aromatic
                'c3': '@atom:c3',    # sp3 carbon
                'ct': '@atom:ct',    # sp3 alkane
            },
            'O': {
                'oh': '@atom:oh',    # Alcohol
                'os': '@atom:os',    # Ether
                'o': '@atom:o',      # Carbonyl
            },
            'N': {
                'n': '@atom:n',      # Amine
                'na': '@atom:na',    # Aromatic
                'nh': '@atom:nh',    # Primary amine
            }
        }
    }

    def __init__(self, force_field: str = 'oplsaa', verbose: bool = False):
        """
        Initialize connection point modifier.

        Args:
            force_field: Force field name ('oplsaa', 'lopls', or 'gaff')
            verbose: Enable verbose logging
        """
        self.force_field = force_field.lower()
        self.verbose = verbose

        # Normalize force field name
        if 'gaff' in self.force_field:
            self.ff_key = 'gaff'
        else:
            self.ff_key = 'oplsaa'

    def _get_atom_type_name(self, atom: Chem.Atom) -> Optional[str]:
        """
        Determine the type name of an atom (e.g., 'CH3', 'OH').

        Analyzes the atom's environment to determine its type name.

        Args:
            atom: RDKit Atom object

        Returns:
            Type name string (e.g., 'CH3', 'OH', 'O')
            Returns None if type cannot be determined
        """
        element = atom.GetSymbol()
        degree = atom.GetDegree()
        num_hs = atom.GetTotalNumHs()
        total_valence = degree + num_hs

        # Carbon type determination
        if element == 'C':
            # Check for aromatic
            if atom.GetIsAromatic():
                return 'Car'
            # Check for carbonyl
            for bond in atom.GetBonds():
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    neighbor = bond.GetOtherAtomIdx(atom.GetIdx())
                    neighbor_elem = atom.GetOwningMol().GetAtomWithIdx(neighbor).GetSymbol()
                    if neighbor_elem == 'O':
                        return 'C=O'  # Carbonyl carbon
            # sp3 carbons
            if num_hs == 3:
                return 'CH3'
            elif num_hs == 2:
                return 'CH2'
            elif num_hs == 1:
                return 'CH'
            else:
                return 'C'  # Quaternary or other

        # Oxygen type determination
        elif element == 'O':
            # Check for carbonyl (double bond to carbon)
            for bond in atom.GetBonds():
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    return 'O=C'  # Carbonyl oxygen
            # Check if bonded to carbonyl carbon
            for bond in atom.GetBonds():
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    neighbor = bond.GetOtherAtomIdx(atom.GetIdx())
                    neighbor_atom = atom.GetOwningMol().GetAtomWithIdx(neighbor)
                    if neighbor_atom.GetSymbol() == 'C':
                        # Check if that carbon has double bond to O
                        for nbond in neighbor_atom.GetBonds():
                            if nbond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                nneighbor = nbond.GetOtherAtomIdx(neighbor)
                                if atom.GetOwningMol().GetAtomWithIdx(nneighbor).GetSymbol() == 'O':
                                    return 'O-C=O'  # Eester oxygen
            # Alcohol or ether
            if num_hs == 1:
                return 'OH'  # Alcohol
            else:
                return 'O'   # Ether

        # Nitrogen type determination
        elif element == 'N':
            # Check for amide (N bonded to carbonyl carbon)
            for bond in atom.GetBonds():
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    neighbor = bond.GetOtherAtomIdx(atom.GetIdx())
                    neighbor_atom = atom.GetOwningMol().GetAtomWithIdx(neighbor)
                    if neighbor_atom.GetSymbol() == 'C':
                        # Check if that carbon has double bond to O
                        for nbond in neighbor_atom.GetBonds():
                            if nbond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                nneighbor = nbond.GetOtherAtomIdx(neighbor)
                                if atom.GetOwningMol().GetAtomWithIdx(nneighbor).GetSymbol() == 'O':
                                    return 'N-C=O'  # Amide nitrogen
            # Amines
            if num_hs == 2:
                return 'NH2'  # Primary amine
            elif num_hs == 1:
                return 'NH'   # Secondary amine
            else:
                return 'N'    # Tertiary amine or other

        return None

    def get_connection_atom_type(
        self,
        atom: Chem.Atom,
        mol: Chem.Mol
    ) -> Tuple[Optional[str], Optional[str]]:
        """
        Get the atom type change for a connection point.

        Determines what atom type should be used when creating a polymer
        connection at this atom.

        Args:
            atom: RDKit Atom object at the connection point
            mol: RDKit Mol object (for context)

        Returns:
            (original_type, connection_type) tuple
            - original_type: Current atom type (e.g., '@atom:80')
            - connection_type: Atom type for connection point (e.g., '@atom:82')
            Returns (None, None) if no change needed

        Example:
            >>> carbon_atom = mol.GetAtomWithIdx(0)  # CH3 carbon
            >>> original, connection = modifier.get_connection_atom_type(carbon_atom, mol)
            >>> print(original, connection)
            @atom:80 @atom:82
        """
        element = atom.GetSymbol()

        # Get type name for this atom
        type_name = self._get_atom_type_name(atom)

        if type_name is None:
            if self.verbose:
                logger.warning(f"Could not determine type name for {element} atom {atom.GetIdx()}")
            return None, None

        # Determine connection type name
        if element == 'C':
            # CH3 → CH2, CH2 → CH, etc.
            if type_name == 'CH3':
                connection_type_name = 'CH2'
            elif type_name == 'CH2':
                connection_type_name = 'CH'
            elif type_name == 'CH':
                connection_type_name = 'C'
            elif type_name in ['Car', 'C=O', 'O-C=O']:
                # These don't change at connection points
                return None, None
            else:
                connection_type_name = type_name

        elif element == 'O':
            # OH → O (alcohol to ether)
            if type_name == 'OH':
                connection_type_name = 'O'
            elif type_name in ['O=C', 'O-C=O']:
                # Carbonyl and ester oxygens don't change
                return None, None
            else:
                connection_type_name = type_name

        elif element == 'N':
            # NH2 → NH, NH → N (primary to secondary to tertiary)
            if type_name == 'NH2':
                connection_type_name = 'NH'
            elif type_name == 'NH':
                connection_type_name = 'N'
            elif type_name == 'N-C=O':
                # Amide nitrogen doesn't change
                return None, None
            else:
                connection_type_name = type_name

        else:
            # Other elements typically don't change
            return None, None

        # Get actual atom types from mapping
        original_atom_type, connection_atom_type = get_connection_type(
            element, type_name, connection_type_name, self.force_field
        )

        # If no mapping found, try to get from TYPE_NAMES
        if original_atom_type is None:
            if self.ff_key in self.TYPE_NAMES and element in self.TYPE_NAMES[self.ff_key]:
                original_atom_type = self.TYPE_NAMES[self.ff_key][element].get(type_name)
                connection_atom_type = self.TYPE_NAMES[self.ff_key][element].get(connection_type_name)

        if self.verbose and original_atom_type and connection_atom_type:
            logger.info(
                f"Connection point: {element}({type_name}) "
                f"{original_atom_type} → {connection_atom_type}"
            )

        return original_atom_type, connection_atom_type

    def should_skip_hydrogen_removal(self, element: str, type_name: str) -> bool:
        """
        Determine if hydrogen should be skipped for this connection point.

        Some connection points (e.g., carbonyl carbon) don't have removable
        hydrogens and should be handled differently.

        Args:
            element: Element symbol
            type_name: Atom type name (e.g., 'CH3', 'O=C')

        Returns:
            True if hydrogen removal should be skipped, False otherwise
        """
        # Carbonyl carbons don't have removable H
        if element == 'C' and type_name in ['C=O', 'Car']:
            return True

        # Carbonyl oxygens and ester oxygens don't have removable H
        if element == 'O' and type_name in ['O=C', 'O-C=O']:
            return True

        # Amide nitrogen might keep its H
        if element == 'N' and type_name == 'N-C=O':
            return True

        return False
