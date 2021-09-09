#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 12:19:08 2018

@author: zwu
"""
import sys
sys.path.append('/home/zwu/AutoPoly')
import AutoPoly

# Define the system
# out is the folder name for the output
system=AutoPoly.System(out="test")

# create polymers
# Just an example of polypropylene and polyethylene with all-atom and united-atom resolutions
PE1=AutoPoly.Polymer(ChainNum=2,Sequence=["PPUA"]*20)
PE2=AutoPoly.Polymer(ChainNum=2,Sequence=["PEUA"]*10)
PE3=AutoPoly.Polymer(ChainNum=3,Sequence=["PEAA"]*15)

# polymerization
# Name is the output folder for this polymer
poly=AutoPoly.Polymerization(Name="Polymer",System=system,Model=[PE1,PE2,PE3],run=True)
