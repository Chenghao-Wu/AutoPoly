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
import os
import sys
import random
from typing import List, Optional, Union

from .system import logger

class Polymer(object):
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
        SequnceLen (int): Length of the monomer sequence
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
        self.sequence = Sequence[0] if isinstance(Sequence[0], list) else Sequence  # Store original sequence
        self.DOP = DOP if DOP > 0 else len(self.sequence)
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
        self.SequnceLen = len(self.sequence)
        self.set_merSet(self.sequence)
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
        
        Raises:
            SystemExit: If ChainNum is 0 (no chains specified)
        """
        sequence = self.sequence
        self.SequnceLen = len(sequence)
        self.set_merSet(sequence)
        self.set_dop(self.SequnceLen)

        if self.ChainNum == 0:
            logger.error("Error : Please set number of chains ")
            sys.exit()

        for chainii in range(self.ChainNum):
            merSet = []
            merSet_ = []

            if self.topology == "ring":
                # For ring polymers, all monomers are equivalent
                for merii in range(self.SequnceLen):
                    # Remove any existing .lt extension and add it cleanly
                    base_name = sequence[merii].replace('.lt', '')
                    merSet.append(f"{base_name}i.lt")  # Add internal monomer suffix
                    merSet_.append(f"{base_name}i")    # Name without extension
                self.sequenceSet.append(merSet)
                self.sequenceName.append(merSet_)
            else:
                # Original logic for linear polymers
                if self.tacticity == 'atactic':
                    if self.SequnceLen > 1:
                        for merii in range(self.SequnceLen):
                            if merii == 0:
                                if bool(random.choice([True, False])):
                                    merSet.append(sequence[merii]+"le_T1.lt")
                                    merSet_.append(sequence[merii]+"le_T1")
                                else:
                                    merSet.append(sequence[merii]+"le.lt")
                                    merSet_.append(sequence[merii]+"le")
                            elif merii == self.SequnceLen-1:
                                if bool(random.choice([True, False])):
                                    merSet.append(sequence[merii]+"re_T1.lt")
                                    merSet_.append(sequence[merii]+"re_T1")
                                else:
                                    merSet.append(sequence[merii]+"re.lt")
                                    merSet_.append(sequence[merii]+"re")
                            else:
                                if bool(random.choice([True, False])):
                                    merSet.append(sequence[merii]+"i_T1.lt")
                                    merSet_.append(sequence[merii]+"i_T1")
                                else:
                                    merSet.append(sequence[merii]+"i.lt")
                                    merSet_.append(sequence[merii]+"i")
                        self.sequenceSet.append(merSet)
                        self.sequenceName.append(merSet_)
                    elif self.SequnceLen == 1:
                        for merii in range(self.SequnceLen):
                            merSet.append(sequence[merii]+".lt")
                            merSet_.append(sequence[merii])
                        self.sequenceSet.append(merSet)
                        self.sequenceName.append(merSet_)
                elif self.tacticity == 'isotactic':
                    chosenTac =".lt"
                    chosenTac_name=''
                    if self.SequnceLen>1:
                        for merii in range(self.SequnceLen):
                            if merii==0:
                                merSet.append(sequence[merii]+"le"+chosenTac)
                                merSet_.append(sequence[merii]+"le"+chosenTac_name)
                            elif merii == self.SequnceLen-1:
                                merSet.append(sequence[merii]+"re"+chosenTac)
                                merSet_.append(sequence[merii]+"re"+chosenTac_name)
                            else:
                                merSet.append(sequence[merii]+"i"+chosenTac)
                                merSet_.append(sequence[merii]+"i"+chosenTac_name)
                        self.sequenceSet.append(merSet)
                        self.sequenceName.append(merSet_)
                    elif self.SequnceLen==1:
                        for merii in range(self.SequnceLen):
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

                    if self.SequnceLen>1:
                        for merii in range(self.SequnceLen):
                            if merii%2==0:
                                currentTac=startTac
                                currentTac_name = startTac_name
                            else:
                                currentTac=nextTac
                                currentTac_name=nextTac_name
                                
                            if merii==0:
                                merSet.append(sequence[merii]+"le"+currentTac)
                                merSet_.append(sequence[merii]+"le"+currentTac_name)
                            elif merii == self.SequnceLen-1:
                                merSet.append(sequence[merii]+"re"+currentTac)
                                merSet_.append(sequence[merii]+"re"+currentTac_name)
                            else:
                                merSet.append(sequence[merii]+"i"+currentTac)
                                merSet_.append(sequence[merii]+"i"+currentTac_name)
                        self.sequenceSet.append(merSet)
                        self.sequenceName.append(merSet_)
                    elif self.SequnceLen==1:
                        for merii in range(self.SequnceLen):
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
            'sequence_length': self.SequnceLen,
            'mer_set': self.merSet,
            'sequence_set': self.sequenceSet,
            'sequence_names': self.sequenceName
        }