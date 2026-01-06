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
import sys
import random
from typing import List, Optional, Union

from .system import logger

class Polymer:
    """
    Polymer class for defining polymer structures and properties.
    
    This class manages polymer chain definitions including topology, tacticity,
    monomer sequences, and chain generation for molecular dynamics simulations.
    
    Attributes:
        ChainNum (int): Number of polymer chains to generate
        sequence (list): Original monomer sequence
        DOP (int): Degree of polymerization
        topology (str): Polymer topology ('linear' or 'ring')
        tacticity (str): Polymer tacticity ('atactic', 'isotactic', 'syndiotactic')
        sequenceSet (list): List of monomer file names for each chain
        sequenceName (list): List of monomer names for each chain
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
        self.merSet = []

        # Validate topology
        if topology not in ["linear", "ring"]:
            raise ValueError("Topology must be either 'linear' or 'ring'")

        # Set up the sequence
        self.SequenceLen = len(self.sequence)
        self.set_merSet(self.sequence)

        # Set DOP before calling set_Sequence() so it's not overwritten
        # DOP represents chain length, not number of unique monomer types
        self.DOP = DOP if DOP > 0 else len(self.sequence)
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
        Set up the polymer sequence based on tacticity and chain number.

        This method generates the monomer file names and names for each chain
        based on the specified topology and tacticity. It handles:
        - Linear vs ring topology
        - Atactic, isotactic, and syndiotactic tacticity
        - Proper file naming conventions for different monomer positions

        The sequence defines the monomer pattern that repeats DOP times.
        For example, sequence=['A', 'B'] with DOP=6 creates A-B-A-B-A-B.

        Raises:
            SystemExit: If ChainNum is 0 (no chains specified)
        """
        # Clear existing sequences before regenerating
        self.sequenceSet = []
        self.sequenceName = []

        sequence = self.sequence
        self.SequenceLen = len(sequence)
        self.set_merSet(sequence)
        # Don't overwrite DOP - it represents chain length, not sequence length

        if self.ChainNum == 0:
            logger.error("Error : Please set number of chains ")
            sys.exit()

        for chainii in range(self.ChainNum):
            merSet = []
            merSet_ = []

            if self.topology == "ring":
                # For ring polymers, all monomers are equivalent
                for merii in range(self.DOP):
                    # Cycle through sequence elements
                    seq_idx = merii % self.SequenceLen
                    # Remove any existing .lt extension and add it cleanly
                    base_name = sequence[seq_idx].replace('.lt', '')
                    merSet.append(f"{base_name}i.lt")  # Add internal monomer suffix
                    merSet_.append(f"{base_name}i")    # Name without extension
                self.sequenceSet.append(merSet)
                self.sequenceName.append(merSet_)
            else:
                # Original logic for linear polymers
                if self.tacticity == 'atactic':
                    if self.DOP > 1:
                        for merii in range(self.DOP):
                            # Cycle through sequence elements
                            seq_idx = merii % self.SequenceLen
                            if merii == 0:
                                if bool(random.choice([True, False])):
                                    merSet.append(sequence[seq_idx]+"le_T1.lt")
                                    merSet_.append(sequence[seq_idx]+"le_T1")
                                else:
                                    merSet.append(sequence[seq_idx]+"le.lt")
                                    merSet_.append(sequence[seq_idx]+"le")
                            elif merii == self.DOP-1:
                                if bool(random.choice([True, False])):
                                    merSet.append(sequence[seq_idx]+"re_T1.lt")
                                    merSet_.append(sequence[seq_idx]+"re_T1")
                                else:
                                    merSet.append(sequence[seq_idx]+"re.lt")
                                    merSet_.append(sequence[seq_idx]+"re")
                            else:
                                if bool(random.choice([True, False])):
                                    merSet.append(sequence[seq_idx]+"i_T1.lt")
                                    merSet_.append(sequence[seq_idx]+"i_T1")
                                else:
                                    merSet.append(sequence[seq_idx]+"i.lt")
                                    merSet_.append(sequence[seq_idx]+"i")
                        self.sequenceSet.append(merSet)
                        self.sequenceName.append(merSet_)
                    elif self.DOP == 1:
                        for merii in range(self.DOP):
                            merSet.append(sequence[merii]+".lt")
                            merSet_.append(sequence[merii])
                        self.sequenceSet.append(merSet)
                        self.sequenceName.append(merSet_)
                elif self.tacticity == 'isotactic':
                    chosenTac =".lt"
                    chosenTac_name=''
                    if self.DOP>1:
                        for merii in range(self.DOP):
                            # Cycle through sequence elements
                            seq_idx = merii % self.SequenceLen
                            if merii==0:
                                merSet.append(sequence[seq_idx]+"le"+chosenTac)
                                merSet_.append(sequence[seq_idx]+"le"+chosenTac_name)
                            elif merii == self.DOP-1:
                                merSet.append(sequence[seq_idx]+"re"+chosenTac)
                                merSet_.append(sequence[seq_idx]+"re"+chosenTac_name)
                            else:
                                merSet.append(sequence[seq_idx]+"i"+chosenTac)
                                merSet_.append(sequence[seq_idx]+"i"+chosenTac_name)
                        self.sequenceSet.append(merSet)
                        self.sequenceName.append(merSet_)
                    elif self.DOP==1:
                        for merii in range(self.DOP):
                            merSet.append(sequence[merii]+".lt")
                            merSet_.append(sequence[merii])
                        self.sequenceSet.append(merSet)
                        self.sequenceName.append(merSet_)

                elif self.tacticity == 'syndiotactic':

                    randbool = bool(random.choice([True, False]))
                    if randbool:
                        startTac="_T1.lt"
                        nextTac =".lt"
                        startTac_name="_T1"
                        nextTac_name =""
                    else:
                        startTac=".lt"
                        nextTac ="_T1.lt"
                        startTac_name=""
                        nextTac_name ="_T1"

                    if self.DOP>1:
                        for merii in range(self.DOP):
                            # Cycle through sequence elements
                            seq_idx = merii % self.SequenceLen
                            if merii%2==0:
                                currentTac=startTac
                                currentTac_name = startTac_name
                            else:
                                currentTac=nextTac
                                currentTac_name=nextTac_name

                            if merii==0:
                                merSet.append(sequence[seq_idx]+"le"+currentTac)
                                merSet_.append(sequence[seq_idx]+"le"+currentTac_name)
                            elif merii == self.DOP-1:
                                merSet.append(sequence[seq_idx]+"re"+currentTac)
                                merSet_.append(sequence[seq_idx]+"re"+currentTac_name)
                            else:
                                merSet.append(sequence[seq_idx]+"i"+currentTac)
                                merSet_.append(sequence[seq_idx]+"i"+currentTac_name)
                        self.sequenceSet.append(merSet)
                        self.sequenceName.append(merSet_)
                    elif self.DOP==1:
                        for merii in range(self.DOP):
                            merSet.append(sequence[merii]+".lt")
                            merSet_.append(sequence[merii])
                        self.sequenceSet.append(merSet)
                        self.sequenceName.append(merSet_)
        print(self.sequenceSet)
        print(self.sequenceName)
    
    def get_sequence_set(self) -> List[List[str]]:
        """
        Get the sequence set for all chains.
        
        Returns:
            List[List[str]]: List of monomer file names for each chain
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
            'sequence_names': self.sequenceName
        }