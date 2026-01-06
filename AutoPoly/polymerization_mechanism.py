"""
Polymerization Mechanism Detector Module for AutoPoly

Automatically detects the polymerization mechanism from a molecule's structure.
Supports mechanism detection with DOP (Degree of Polymerization) awareness.

Key Features:
- Auto-detects vinyl addition, esterification, amidation, etherification
- DOP-aware: Returns 'none' for DOP=1 unless obvious vinyl pattern
- Preserves single molecule generation support

Author: AutoPoly Development Team
"""

from typing import Optional, List, Tuple
from rdkit import Chem

from .polymerization_patterns import (
    POLYMERIZATION_PATTERNS,
    PATTERN_NONE,
    PATTERN_VINYL,
    PATTERN_ESTER,
    PATTERN_AMIDE,
    PATTERN_ETHER
)
from .system import logger


class PolymerizationMechanism:
    """
    Detects polymerization mechanisms from molecular structure.

    This class analyzes a molecule and determines which polymerization
    mechanism it uses, with special handling for DOP=1 (single molecules).

    Attributes:
        verbose: Enable verbose logging

    Example:
        >>> detector = PolymerizationMechanism()
        >>> mechanism = detector.detect_mechanism(mol, dop=10)
        >>> print(mechanism)  # 'vinyl_addition', 'esterification', etc.
    """

    def __init__(self, verbose: bool = True):
        """
        Initialize mechanism detector.

        Args:
            verbose: Enable verbose logging
        """
        self.verbose = verbose

    def _has_pattern(self, mol: Chem.Mol, smarts: str) -> bool:
        """
        Check if molecule contains a SMARTS pattern.

        Args:
            mol: RDKit Mol object
            smarts: SMARTS pattern string

        Returns:
            True if pattern is found, False otherwise
        """
        try:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                return False
            return mol.HasSubstructMatch(pattern)
        except Exception as e:
            if self.verbose:
                logger.warning(f"Error matching SMARTS '{smarts}': {e}")
            return False

    def _count_functional_groups(self, mol: Chem.Mol, smarts: str) -> int:
        """
        Count the number of functional group matches in a molecule.

        Args:
            mol: RDKit Mol object
            smarts: SMARTS pattern string

        Returns:
            Number of matches
        """
        try:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                return 0
            matches = mol.GetSubstructMatches(pattern)
            return len(matches)
        except Exception:
            return 0

    def detect_mechanism(self, mol: Chem.Mol, dop: int = 1) -> str:
        """
        Auto-detect polymerization mechanism from molecule structure.

        Detection priority (when DOP>1):
        1. Vinyl addition (C=C detected)
        2. Condensation with both groups present (esterification, amidation)
        3. Diol/functional groups (etherification)
        4. Default: 'none' (non-polymerizable)

        Special handling for DOP=1:
        - Returns 'none' unless clear vinyl pattern is found
        - Allows single molecule generation for any compound

        Args:
            mol: RDKit molecule object
            dop: Degree of polymerization (1 = single molecule)

        Returns:
            str: Mechanism name ('none', 'vinyl_addition', 'esterification', etc.)

        Example:
            >>> detector = PolymerizationMechanism()
            >>> detector.detect_mechanism(mol, dop=1)  # Single molecule
            'none'
            >>> detector.detect_mechanism(mol, dop=10)  # Polymer
            'vinyl_addition'
        """
        # Add hydrogens for pattern matching
        mol_with_h = Chem.AddHs(mol)

        # For DOP=1, be conservative - only detect if obvious vinyl pattern
        if dop == 1:
            if self._has_pattern(mol_with_h, PATTERN_VINYL.smarts):
                if self.verbose:
                    logger.info("DOP=1 with vinyl pattern detected")
                return 'vinyl_addition'
            else:
                # For DOP=1 without clear pattern, treat as non-polymerizable
                # This preserves single molecule generation
                if self.verbose:
                    logger.info("DOP=1: treating as non-polymerizable (no clear pattern)")
                return 'none'

        # For DOP>1, check all patterns in priority order

        # 1. Check vinyl addition first (most common)
        if self._has_pattern(mol_with_h, PATTERN_VINYL.smarts):
            if self.verbose:
                logger.info("Detected: vinyl addition polymerization")
            return 'vinyl_addition'

        # 2. Check esterification (carboxyl + alcohol)
        # Need both carboxyl AND alcohol groups
        carboxyl_smarts = '[CX3](=[OX1])[OX2H1]'  # Carboxyl group
        alcohol_smarts = '[$([OX2H])]'  # Alcohol

        carboxyl_count = self._count_functional_groups(mol_with_h, carboxyl_smarts)
        alcohol_count = self._count_functional_groups(mol_with_h, alcohol_smarts)

        # For esterification, need carboxyl + separate alcohol (not the -OH from carboxyl)
        # Count alcohols that are NOT part of carboxyl groups
        # Carboxyl contributes 1 to alcohol_count, so need alcohol_count > carboxyl_count
        if carboxyl_count >= 1 and alcohol_count > carboxyl_count:
            if self.verbose:
                logger.info(f"Detected: esterification (carboxyl: {carboxyl_count}, alcohol: {alcohol_count})")
            return 'esterification'

        # 3. Check amidation (carboxyl + amine)
        amine_smarts = '[NX3]'  # Primary/secondary amine (any N with connectivity 3)

        amine_count = self._count_functional_groups(mol_with_h, amine_smarts)

        if carboxyl_count >= 1 and amine_count >= 1:
            if self.verbose:
                logger.info(f"Detected: amidation (carboxyl: {carboxyl_count}, amine: {amine_count})")
            return 'amidation'

        # 4. Check etherification (2+ alcohol groups)
        # Need at least 2 alcohol groups that are NOT from carboxyls
        # Carboxyl_count alcohols are from carboxyl -OH groups
        non_carboxyl_alcohols = alcohol_count - carboxyl_count
        if non_carboxyl_alcohols >= 2:
            if self.verbose:
                logger.info(f"Detected: etherification (alcohol groups: {alcohol_count}, non-carboxyl: {non_carboxyl_alcohols})")
            return 'etherification'

        # 5. Default: non-polymerizable
        if self.verbose:
            logger.info("No clear polymerization pattern detected, treating as non-polymerizable")

        return 'none'

    def get_connection_atoms(self, mol: Chem.Mol, mechanism: str) -> List[int]:
        """
        Identify the atoms that form polymer connections.

        Uses SMARTS patterns to find the specific atoms involved in
        polymerization bonds.

        Args:
            mol: RDKit Mol object
            mechanism: Mechanism name (e.g., 'vinyl_addition', 'esterification')

        Returns:
            List of atom indices that are connection points
            Empty list for 'none' mechanism

        Example:
            >>> atoms = detector.get_connection_atoms(mol, 'vinyl_addition')
            >>> print(atoms)  # [0, 1] - the two carbons in C=C
        """
        if mechanism == 'none':
            return []

        # Special handling for amidation: find carboxyl carbon and amine nitrogen separately
        if mechanism == 'amidation':
            # Find carboxyl carbon
            carboxyl_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]')
            carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)

            # Find amine nitrogen
            amine_pattern = Chem.MolFromSmarts('[NX3]')
            amine_matches = mol.GetSubstructMatches(amine_pattern)

            if carboxyl_matches and amine_matches:
                carboxyl_c = carboxyl_matches[0][0]  # First atom in carboxyl match
                amine_n = amine_matches[0][0]        # First atom in amine match
                return [carboxyl_c, amine_n]
            return []

        pattern = POLYMERIZATION_PATTERNS.get(mechanism)
        if pattern is None or pattern.smarts is None:
            return []

        try:
            smarts_pattern = Chem.MolFromSmarts(pattern.smarts)
            if smarts_pattern is None:
                return []

            matches = mol.GetSubstructMatches(smarts_pattern)

            # Extract atom indices from matches
            # For vinyl (C=C), returns both carbons
            # For condensation, returns the connecting atoms
            connection_atoms = []
            for match in matches:
                connection_atoms.extend(list(match))

            # Remove duplicates and return first 2
            connection_atoms = list(set(connection_atoms))[:2]

            return connection_atoms

        except Exception as e:
            if self.verbose:
                logger.warning(f"Error finding connection atoms: {e}")
            return []

    def get_backbone_atoms(self, mol: Chem.Mol, mechanism: str) -> List[int]:
        """
        Identify backbone atoms for the polymerization mechanism.

        Backbone atoms are the main atoms that form the polymer chain.

        Args:
            mol: RDKit Mol object
            mechanism: Mechanism name

        Returns:
            List of backbone atom indices (typically 2 atoms)

        Example:
            >>> atoms = detector.get_backbone_atoms(mol, 'vinyl_addition')
            >>> print(atoms)  # [0, 1] - the two backbone carbons
        """
        return self.get_connection_atoms(mol, mechanism)

    def requires_variant_generation(self, mechanism: str) -> bool:
        """
        Check if a mechanism requires variant generation (le/re/i).

        Args:
            mechanism: Mechanism name

        Returns:
            True if variants needed, False otherwise (for 'none')
        """
        return mechanism != 'none'

    def is_condensation_polymerization(self, mechanism: str) -> bool:
        """
        Check if mechanism is condensation polymerization.

        Args:
            mechanism: Mechanism name

        Returns:
            True if condensation, False if addition
        """
        pattern = POLYMERIZATION_PATTERNS.get(mechanism)
        return pattern.is_condensation if pattern else False

    def get_mechanism_info(self, mechanism: str) -> dict:
        """
        Get detailed information about a mechanism.

        Args:
            mechanism: Mechanism name

        Returns:
            dict with keys: name, smarts, connection_atoms, description, is_condensation
            Returns None if mechanism not found
        """
        pattern = POLYMERIZATION_PATTERNS.get(mechanism)
        if pattern is None:
            return None

        return {
            'name': pattern.name,
            'smarts': pattern.smarts,
            'connection_atoms': pattern.connection_atoms,
            'description': pattern.description,
            'is_condensation': pattern.is_condensation,
            'requires_two_groups': pattern.requires_two_groups
        }


def detect_mechanism(mol: Chem.Mol, dop: int = 1, verbose: bool = True) -> str:
    """
    Convenience function to detect polymerization mechanism.

    Args:
        mol: RDKit Mol object
        dop: Degree of polymerization (1 = single molecule)
        verbose: Enable verbose logging

    Returns:
        Mechanism name string

    Example:
        >>> from rdkit import Chem
        >>> mol = Chem.MolFromSmiles('C=C')
        >>> detect_mechanism(mol, dop=10)
        'vinyl_addition'
    """
    detector = PolymerizationMechanism(verbose=verbose)
    return detector.detect_mechanism(mol, dop)
