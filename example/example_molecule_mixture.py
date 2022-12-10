#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 12:19:08 2018

@author: zwu
"""
import sys
sys.path.append('/home/zwq2834/development/AutoPoly')
import AutoPoly

# Define the system
# out is the folder name for the output
system=AutoPoly.System(out="test_mole_mix")

# create polymers
# Just an example of polypropylene and polyethylene with all-atom and united-atom resolutions
naphthalene=AutoPoly.Polymer(ChainNum=50,Sequence=["Naphthalene"])
benzene=AutoPoly.Polymer(ChainNum=50,Sequence=["Benzene"])

# polymerization
# Name is the output folder for this polymer
poly=AutoPoly.Polymerization(Name="Molecular_Mixture",System=system,Model=[naphthalene,benzene],run=True)
