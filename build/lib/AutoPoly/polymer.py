#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Polymer Definition Module

This module provides the Polymer class for defining polymer properties,
sequences, and structural characteristics for molecular dynamics simulations.

The Polymer class handles:
- Polymer chain definitions with various topologies (linear, ring)
- Tacticity control (atactic, isotactic, syndiotactic)
- Monomer sequence management
- Chain generation with proper bonding patterns

Created on Fri Dec 21 12:19:08 2018
@author: zwu
"""
import random
from typing import List, Optional, Union

from .system import logger
from .exceptions import ValidationError
from .conf import MAX_DOP, MAX_SEQUENCE_LENGTH, MAX_UNIQUE_MONOMERS
from .validation import validate_smiles


class Polymer:
    """
    Polymer class for defining polymer structures and properties.
    
    This class manages polymer chain definitions including topology, tacticity,
    monomer sequences, and chain generation for molecular dynamics simulations.
    
    Attributes:
        chain_num (int): Number of polymer chains to generate
        sequence (list): Explicit monomer sequence (pSMILES strings)
        dop (int): Degree of polymerization (derived from sequence length)
        topology (str): Polymer topology ('linear' or 'ring')
        tacticity (str): Polymer tacticity ('atactic', 'isotactic', 'syndiotactic')
        sequence_set (list): List of monomer identifiers for each chain (includes _T1 markers)
        sequence_name (list): List of monomer names for each chain
        tacticity_set (list): List of tacticity choices (bool) for each position in each chain
        mer_set (list): Unique set of monomers used
    """

    def __init__(self, chain_num: int = None, sequence: list = None,
                 topology: str = "linear", tacticity: str = 'atactic') -> None:
        """
        Initialize the Polymer with explicit monomer sequence.

        Args:
            chain_num (int): Number of polymer chains to generate
            sequence (list): List of pSMILES strings specifying exact monomer at each position.
                           Example: ["[*]CC[*]", "[*]C=C[*]", "[*]CC[*]"] for ABA copolymer
            topology (str): Polymer topology, either "linear" (default) or "ring"
            tacticity (str): Polymer tacticity ('atactic', 'isotactic', or 'syndiotactic')

        Note:
            DOP is automatically derived from len(sequence). No separate DOP parameter needed.
            This is a BREAKING CHANGE from the old API.

        Raises:
            ValueError: If topology or tacticity is invalid
            ValidationError: If sequence is empty or exceeds limits

        Examples:
            >>> # Block copolymer (ABA triblock)
            >>> poly = Polymer(
            ...     chain_num=5,
            ...     sequence=["[*]CC[*]", "[*]CC[*]", "[*]C=C[*]", "[*]C=C[*]", "[*]C=C[*]", "[*]CC[*]", "[*]CC[*]"],
            ...     topology="linear",
            ...     tacticity="atactic"
            ... )
        """
        self.chain_num = chain_num
        # Handle None, empty sequences, and nested lists properly
        if sequence and isinstance(sequence[0], list):
            self.sequence = sequence[0]
        elif sequence:
            self.sequence = sequence
        else:
            raise ValidationError("sequence cannot be empty")
        self.topology = topology
        self.tacticity = tacticity

        # Initialize empty lists with Pythonic naming
        self.sequence_set = []
        self.sequence_name = []
        self.tacticity_set = []  # Store tacticity choices separately
        self.mer_set = []

        # Validate topology
        if topology not in ["linear", "ring"]:
            raise ValueError("topology must be either 'linear' or 'ring'")

        # Validate tacticity
        if tacticity not in ["isotactic", "syndiotactic", "atactic"]:
            raise ValueError("tacticity must be 'isotactic', 'syndiotactic', or 'atactic'")

        # Set up unique monomers
        self.set_merSet(self.sequence)

        # Validate sequence length
        sequence_length = len(self.sequence)
        if sequence_length > MAX_SEQUENCE_LENGTH:
            raise ValidationError(
                f"Sequence length ({sequence_length}) exceeds maximum {MAX_SEQUENCE_LENGTH}. "
                f"This limit prevents resource exhaustion."
            )

        # Validate unique monomer count
        if len(self.mer_set) > MAX_UNIQUE_MONOMERS:
            raise ValidationError(
                f"Number of unique monomers ({len(self.mer_set)}) exceeds maximum {MAX_UNIQUE_MONOMERS}. "
                f"This limit prevents resource exhaustion."
            )

        # Validate all unique SMILES in the sequence (only if they look like SMILES)
        # Skip validation for monomer names (like "PE", "PS") that don't contain wildcards
        for smiles in self.mer_set:
            try:
                # Remove .lt extension if present before validation
                smiles_clean = smiles.replace('.lt', '')
                # Only validate if it looks like a SMILES string (contains wildcards or brackets)
                # Monomer names like "PE", "PS" will skip validation
                if "[" in smiles_clean:  # Looks like a SMILES string
                    validate_smiles(smiles_clean, allow_wildcards=True)
            except ValidationError as e:
                raise ValidationError(
                    f"Invalid monomer SMILES in sequence: {e}"
                ) from e

        # DOP is derived from sequence length (single source of truth)
        self.dop = sequence_length

        # Validate DOP
        if self.dop > MAX_DOP:
            raise ValidationError(
                f"DOP ({self.dop}) exceeds maximum {MAX_DOP}. "
                f"This limit prevents resource exhaustion."
            )

        self.set_Sequence()

    def set_merSet(self, merSet: Union[List[str], str]) -> None:
        """
        Set the unique set of monomers used in the polymer.

        Args:
            merSet (Union[List[str], str]): List of monomers or single monomer
        """
        if isinstance(merSet, list):
            self.mer_set = list(dict.fromkeys(merSet))  # Remove duplicates
        else:
            self.mer_set = [merSet]

    def set_dop(self, dop: int) -> None:
        """
        Set the degree of polymerization.

        Note: This method is kept for backward compatibility but should not be used
        in new code. DOP is automatically derived from sequence length.

        Args:
            dop (int): Degree of polymerization
        """
        self.dop = dop

    def set_Sequence(self) -> None:
        """
        Generate monomer identifier sequences using explicit sequence (no cycling).

        This method generates identifier sequences for each chain. Each identifier
        consists of the base SMILES and a tacticity marker (_T1) if applicable.

        The identifiers are NOT pure SMILES - they are strings used by the workflow
        to determine which variant files to use:
        - Base identifier: exact SMILES from the sequence at each position
        - Tacticity marker: "_T1" suffix indicates use of T1 chirality variant

        Key changes from old API:
        - NO CYCLING: Each position in the sequence is used exactly once
        - DOP is len(sequence), not independently specified
        - Use explicit sequences for block copolymers and arbitrary patterns

        Handles:
        - Linear vs ring topology
        - Atactic, isotactic, and syndiotactic tacticity
        - Block copolymers and arbitrary sequences

        Raises:
            ValidationError: If chain_num is 0 (no chains specified)
        """
        # Clear existing sequences before regenerating
        self.sequence_set = []
        self.sequence_name = []
        self.tacticity_set = []

        sequence = self.sequence
        self.set_merSet(sequence)

        if self.chain_num == 0:
            raise ValidationError("chain_num must be greater than 0")

        # For isotactic polymers, make the chirality choice once per polymer instance
        if self.tacticity == 'isotactic' and not hasattr(self, '_isotactic_use_t1'):
            self._isotactic_use_t1 = random.choice([True, False])

        for chain_idx in range(self.chain_num):
            identifier_sequence = []
            tacticity_choices = []

            for position in range(self.dop):
                # Get monomer at position (explicit sequence - NO CYCLING)
                base_smiles = sequence[position]

                # Remove any existing .lt extension if present
                base_smiles = base_smiles.replace('.lt', '')

                # Determine tacticity for this position
                use_t1 = self._get_tacticity_choice(position)
                tacticity_choices.append(use_t1)

                # Create identifier with tacticity marker
                # Note: This is an identifier string, not a pure SMILES
                identifier = base_smiles + ("_T1" if use_t1 else "")
                identifier_sequence.append(identifier)

            self.sequence_set.append(identifier_sequence)
            self.sequence_name.append(identifier_sequence)
            self.tacticity_set.append(tacticity_choices)

        logger.debug(f"Generated {len(self.sequence_set)} chains with DOP={self.dop}")

    def _get_tacticity_choice(self, position: int) -> bool:
        """
        Determine tacticity (T1 variant) choice for a given position.

        Args:
            position: Position in the polymer chain (0-based)

        Returns:
            bool: True if T1 variant should be used, False otherwise

        Tacticity rules:
        - Isotactic: All monomers have same chirality (decided once per instance)
        - Syndiotactic: Alternating chirality (regular, T1, regular, T1, ...)
        - Atactic: Random chirality assignment
        """
        if self.tacticity == 'isotactic':
            return self._isotactic_use_t1

        elif self.tacticity == 'syndiotactic':
            return position % 2 == 1

        else:  # atactic
            return random.choice([True, False])

    def get_tacticity_for_chain(self, chain_idx: int) -> List[bool]:
        """
        Get tacticity choices for a specific chain.
        
        Args:
            chain_idx: Index of the chain (0-based)
            
        Returns:
            List[bool]: List of T1 choices for each position in the chain
        """
        if 0 <= chain_idx < len(self.tacticitySet):
            return self.tacticitySet[chain_idx]
        return []

    def get_sequence_set(self) -> List[List[str]]:
        """
        Get the sequence set for all chains.
        
        Returns:
            List[List[str]]: List of monomer identifiers for each chain
        """
        return self.sequenceSet
    
    def get_sequence_names(self) -> List[List[str]]:
        """
        Get the sequence names for all chains.
        
        Returns:
            List[List[str]]: List of monomer names for each chain
        """
        return self.sequenceName
    
    def get_mer_set(self) -> List[str]:
        """
        Get the unique set of monomers used.

        Returns:
            List[str]: Unique list of monomers
        """
        return self.mer_set
    
    def get_chain_info(self) -> dict:
        """
        Get comprehensive information about the polymer.

        Returns:
            dict: Dictionary containing polymer properties with Pythonic naming
        """
        return {
            'chain_num': self.chain_num,
            'sequence': self.sequence,
            'dop': self.dop,
            'topology': self.topology,
            'tacticity': self.tacticity,
            'mer_set': self.mer_set,
            'sequence_set': self.sequence_set,
            'sequence_names': self.sequence_name,
            'tacticity_set': self.tacticity_set
        }
