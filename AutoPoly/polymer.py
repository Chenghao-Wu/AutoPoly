#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 12:19:08 2018

@author: zwu
"""
import os
import sys

from .system import logger

class Polymer:
    def __init__(self,ChainNum=None,Sequence=None):

        self.ChainNum=ChainNum
        self.sequence=Sequence

        self.DOP=1
        self.sequenceSet=[]
        self.sequenceName=[]
        self.merSet=[]
        self.set_Sequence()
        
    
    def set_merSet(self,merSet):
        self.merSet.append(merSet)
        if isinstance(merSet,list):
            self.merSet=merSet

    def set_dop(self,dop):
        self.DOP=dop

    def set_Sequence(self):
        sequence=self.sequence
        self.SequnceLen=len(sequence)
        self.set_merSet(sequence)
        self.set_dop(self.SequnceLen)

        if self.ChainNum==0:
            logger.error("Error : Please set number of chains ")
            sys.exit()

        for chainii in range(self.ChainNum):
            merSet=[]
            merSet_=[]
            if self.SequnceLen>1:
                for merii in range(self.SequnceLen):
                    if merii==0:
                        merSet.append(sequence[merii]+"le.lt")
                        merSet_.append(sequence[merii]+"le")
                    elif merii == self.SequnceLen-1:
                        merSet.append(sequence[merii]+"re.lt")
                        merSet_.append(sequence[merii]+"re")
                    else:
                        merSet.append(sequence[merii]+"i.lt")
                        merSet_.append(sequence[merii]+"i")
                self.sequenceSet.append(merSet)
                self.sequenceName.append(merSet_)
            elif self.SequnceLen==1:
                for merii in range(self.SequnceLen):
                    merSet.append(sequence[merii]+".lt")
                    merSet_.append(sequence[merii])
                self.sequenceSet.append(merSet)
                self.sequenceName.append(merSet_)
