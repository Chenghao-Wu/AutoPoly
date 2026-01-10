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
        ChainNum (int): Number of polymer chains to generate
        sequence (list): Original monomer sequence (base SMILES)
        DOP (int): Degree of polymerization
        topology (str): Polymer topology ('linear' or 'ring')
        tacticity (str): Polymer tacticity ('atactic', 'isotactic', 'syndiotactic')
        sequenceSet (list): List of monomer identifiers for each chain (includes _T1 markers)
        sequenceName (list): List of monomer names for each chain
        tacticitySet (list): List of tacticity choices (bool) for each position in each chain
        merSet (list): Unique set of monomers used
        SequenceLen (int): Length of the monomer sequence
    """
    
    def __init__(self, ChainNum: int = None, Sequence: list = None, DOP: int = 0,
                 topology: str = "linear", tacticity: str = 'atactic') -> None:
        """
        Initialize the Polymer class.

        Args:
            ChainNum (int, optional): Number of polymer chains. Defaults to None.
            Sequence (list, optional): List of monomer sequences. Defaults to None.
            DOP (int, optional): Degree of polymerization. If 0, uses sequence length.
                               Defaults to 0.
            topology (str, optional): Polymer topology, either "linear" (default)
                                    or "ring". Defaults to "linear".
            tacticity (str, optional): Polymer tacticity ('atactic', 'isotactic',
                                     or 'syndiotactic'). Defaults to 'atactic'.

        Raises:
            ValueError: If topology is not 'linear' or 'ring'
        """
        self.ChainNum = ChainNum
        # Handle None, empty sequences, and nested lists properly
        if Sequence and isinstance(Sequence[0], list):
            self.sequence = Sequence[0]
        elif Sequence:
            self.sequence = Sequence
        else:
            raise ValueError("Sequence cannot be None or empty")
        self.topology = topology
        self.tacticity = tacticity

        # Initialize empty lists
        self.sequenceSet = []
        self.sequenceName = []
        self.tacticitySet = []  # Store tacticity choices separately
        self.merSet = []

        # Validate topology
        if topology not in ["linear", "ring"]:
            raise ValueError("Topology must be either 'linear' or 'ring'")

        # Set up the sequence
        self.SequenceLen = len(self.sequence)
        self.set_merSet(self.sequence)

        # Validate sequence length
        if self.SequenceLen > MAX_SEQUENCE_LENGTH:
            raise ValidationError(
                f"Sequence length ({self.SequenceLen}) exceeds maximum {MAX_SEQUENCE_LENGTH}. "
                f"This limit prevents resource exhaustion."
            )

        # Validate unique monomer count
        if len(self.merSet) > MAX_UNIQUE_MONOMERS:
            raise ValidationError(
                f"Number of unique monomers ({len(self.merSet)}) exceeds maximum {MAX_UNIQUE_MONOMERS}. "
                f"This limit prevents resource exhaustion."
            )

        # Validate all unique SMILES in the sequence (only if they look like SMILES)
        # Skip validation for monomer names (like "PE", "PS") that don't contain wildcards
        for smiles in self.merSet:
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

        # Set DOP before calling set_Sequence() so it's not overwritten
        # DOP represents chain length, not number of unique monomer types
        self.DOP = DOP if DOP > 0 else len(self.sequence)

        # Validate DOP
        if self.DOP > MAX_DOP:
            raise ValidationError(
                f"DOP ({self.DOP}) exceeds maximum {MAX_DOP}. "
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
            self.merSet = list(dict.fromkeys(merSet))  # Remove duplicates
        else:
            self.merSet = [merSet]

    def set_dop(self, dop: int) -> None:
        """
        Set the degree of polymerization.
        
        Args:
            dop (int): Degree of polymerization
        """
        self.DOP = dop

    def set_Sequence(self) -> None:
        """
        Generate monomer identifier sequences for polymer chains.

        This method generates identifier sequences for each chain. Each identifier
        consists of the base SMILES and a tacticity marker (_T1) if applicable.
        
        The identifiers are NOT pure SMILES - they are strings used by the workflow
        to determine which variant files to use:
        - Base identifier: position-based connection pattern (for workflow's variant selection)
        - Tacticity marker: "_T1" suffix indicates use of T1 chirality variant

        Handles:
        - Linear vs ring topology
        - Atactic, isotactic, and syndiotactic tacticity
        - Copolymers (mixed sequences)

        Raises:
            SystemExit: If ChainNum is 0 (no chains specified)
        """
        # Clear existing sequences before regenerating
        self.sequenceSet = []
        self.sequenceName = []
        self.tacticitySet = []

        sequence = self.sequence
        self.SequenceLen = len(sequence)
        self.set_merSet(sequence)

        if self.ChainNum == 0:
            raise ValidationError("ChainNum must be greater than 0")

        # For isotactic polymers, make the chirality choice once per polymer instance
        if self.tacticity == 'isotactic' and not hasattr(self, '_isotactic_use_t1'):
            self._isotactic_use_t1 = random.choice([True, False])

        for chainii in range(self.ChainNum):
            identifier_sequence = []
            tacticity_choices = []

            for i in range(self.DOP):
                # Cycle through base sequence (for copolymers)
                seq_idx = i % self.SequenceLen
                base_smiles = sequence[seq_idx]

                # Remove any existing .lt extension if present
                base_smiles = base_smiles.replace('.lt', '')

                # Determine tacticity for this position
                use_t1 = self._get_tacticity_choice(i)
                tacticity_choices.append(use_t1)

                # Create identifier with tacticity marker
                # Note: This is an identifier string, not a pure SMILES
                identifier = base_smiles + ("_T1" if use_t1 else "")
                identifier_sequence.append(identifier)

            self.sequenceSet.append(identifier_sequence)
            self.sequenceName.append(identifier_sequence)
            self.tacticitySet.append(tacticity_choices)

        logger.debug(f"Generated {len(self.sequenceSet)} chains with DOP={self.DOP}")

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
        return self.merSet
    
    def get_chain_info(self) -> dict:
        """
        Get comprehensive information about the polymer.
        
        Returns:
            dict: Dictionary containing polymer properties
        """
        return {
            'chain_num': self.ChainNum,
            'sequence': self.sequence,
            'dop': self.DOP,
            'topology': self.topology,
            'tacticity': self.tacticity,
            'sequence_length': self.SequenceLen,
            'mer_set': self.merSet,
            'sequence_set': self.sequenceSet,
            'sequence_names': self.sequenceName,
            'tacticity_set': self.tacticitySet
        }
