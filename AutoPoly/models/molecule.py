#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Molecule Definition Module

This module provides the Molecule class for defining small molecule structures
for molecular dynamics simulations.

The Molecule class handles:
- Small molecule definitions (water, benzene, ethanol, etc.)
- Molecule count and concentration control
- Integration with Polymer class for mixed systems
- SMILES validation for non-polymer molecules

Created on 2025
@author: zwu
"""
from typing import List, Dict, Union

from ..core.system import logger


class Molecule:
    """
    Molecule class for defining small molecule structures.

    This class manages small molecule definitions for molecular dynamics simulations.
    It follows the same API pattern as the Polymer class but is designed for
    non-polymer molecules using regular SMILES (without wildcards).

    Attributes:
        Count (int): Number of molecules to generate
        Smiles (str): SMILES string WITHOUT wildcards (e.g., "O", "CCO", "c1ccccc1")
        molecule_name (str): Name for the molecule
        sequenceSet (list): List of molecule identifiers for each molecule instance
        sequenceName (list): List of molecule names for each molecule instance
        merSet (list): Single-element list containing the molecule name
        DOP (int): Always 1 for molecules (for pipeline compatibility)
        _is_molecule (bool): Type detection flag (always True for Molecule class)

    Example:
        >>> water = Molecule(Count=100, Smiles="O", Name="water")
        >>> benzene = Molecule(Count=50, Smiles="c1ccccc1", Name="benzene")
    """

    def __init__(self, Count: int = None, Smiles: str = None, Name: str = None) -> None:
        """
        Initialize the Molecule class.

        Args:
            Count (int): Number of molecules to generate
            Smiles (str): SMILES string WITHOUT wildcards (e.g., "O" for water,
                         "CCO" for ethanol, "c1ccccc1" for benzene)
            Name (str, optional): Name for the molecule. Defaults to None (auto-generate
                                  from SMILES).

        Raises:
            ValueError: If Count is None or zero, Smiles is None/empty, or if Smiles
                       contains wildcards (* or [*])
        """
        # Validate required parameters
        if Count is None or Count <= 0:
            raise ValueError("Count must be a positive integer")

        if Smiles is None or len(Smiles.strip()) == 0:
            raise ValueError("Smiles cannot be None or empty")

        self.Count = Count
        self.Smiles = Smiles.strip()

        # Set molecule name
        if Name is not None:
            self.molecule_name = Name
        else:
            # Auto-generate name from SMILES
            # Create a safe name by replacing special characters
            safe_name = self.Smiles.replace('(', '').replace(')', '') \
                                 .replace('[', '').replace(']', '') \
                                 .replace('=', '').replace('#', '') \
                                 .replace('/', '').replace('\\', '')
            self.molecule_name = f"molecule_{safe_name}"

        # SMILES validation - reject wildcards (used for polymers)
        if '*' in self.Smiles or '[*]' in self.Smiles:
            raise ValueError(
                "Molecule SMILES should not contain wildcards (* or [*]). "
                "Use Polymer class for pSMILES with wildcards."
            )

        # Initialize attributes for pipeline compatibility
        self.DOP = 1  # Always 1 for molecules
        self._is_molecule = True  # Type detection flag

        # Initialize sequence-related attributes
        self.sequenceSet = []
        self.sequenceName = []
        self.merSet = []

        # Set up the molecule structure
        self._set_molecule_structure()

        logger.info(f"Created Molecule: {self.molecule_name} (Count={self.Count}, Smiles={self.Smiles})")

    def _set_molecule_structure(self) -> None:
        """
        Set up the molecule structure for compatibility with the pipeline.

        This method creates sequenceSet and sequenceName lists that mirror the
        Polymer class structure, allowing the GeometryBuilder to handle molecules
        and polymers seamlessly.
        """
        # Set merSet (unique set of molecules)
        self.merSet = [self.molecule_name]

        # Create sequenceSet - one entry per molecule instance
        # Each entry is a list containing the molecule name
        # This mirrors Polymer's sequenceSet structure where:
        # - Polymer: sequenceSet[chain_idx][position_idx]
        # - Molecule: sequenceSet[molecule_idx][0] (always position 0)
        for _ in range(self.Count):
            self.sequenceSet.append([self.molecule_name])
            self.sequenceName.append([self.molecule_name])

        logger.debug(
            f"Set up molecule structure: {len(self.sequenceSet)} instances, "
            f"merSet={self.merSet}"
        )

    def get_count(self) -> int:
        """
        Get the number of molecules.

        Returns:
            int: Number of molecules
        """
        return self.Count

    def get_smiles(self) -> str:
        """
        Get the SMILES string.

        Returns:
            str: SMILES string without wildcards
        """
        return self.Smiles

    def get_name(self) -> str:
        """
        Get the molecule name.

        Returns:
            str: Molecule name
        """
        return self.molecule_name

    def get_sequence_set(self) -> List[List[str]]:
        """
        Get the sequence set for all molecule instances.

        This method mirrors Polymer.get_sequence_set() for compatibility.

        Returns:
            List[List[str]]: List of molecule identifiers for each instance
        """
        return self.sequenceSet

    def get_sequence_names(self) -> List[List[str]]:
        """
        Get the sequence names for all molecule instances.

        This method mirrors Polymer.get_sequence_names() for compatibility.

        Returns:
            List[List[str]]: List of molecule names for each instance
        """
        return self.sequenceName

    def get_mer_set(self) -> List[str]:
        """
        Get the unique set of molecules (always single-element list).

        This method mirrors Polymer.get_mer_set() for compatibility.

        Returns:
            List[str]: Single-element list containing the molecule name
        """
        return self.merSet

    def get_molecule_info(self) -> Dict[str, Union[int, str, List]]:
        """
        Get comprehensive information about the molecule.

        This method mirrors Polymer.get_chain_info() for compatibility.

        Returns:
            dict: Dictionary containing molecule properties including:
                  - count: Number of molecules
                  - smiles: SMILES string
                  - name: Molecule name
                  - dop: Degree of polymerization (always 1)
                  - mer_set: List of molecule names
                  - sequence_set: List of molecule identifiers
                  - sequence_names: List of molecule names
        """
        return {
            'count': self.Count,
            'smiles': self.Smiles,
            'name': self.molecule_name,
            'dop': self.DOP,
            'mer_set': self.merSet,
            'sequence_set': self.sequenceSet,
            'sequence_names': self.sequenceName
        }

    def __repr__(self) -> str:
        """String representation of the Molecule object."""
        return (
            f"Molecule(Count={self.Count}, Smiles='{self.Smiles}', "
            f"Name='{self.molecule_name}')"
        )
