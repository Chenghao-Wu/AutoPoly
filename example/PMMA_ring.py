#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 12:19:08 2018

@author: zwu
"""

import AutoPoly

# Define the system
# out is the folder name for the output
system = AutoPoly.System(out="test")
system.made_folder = True

# Create a ring polymer with 10 PMMA monomers
# Note: Using just "PMMA" as the base name - the code will append the necessary suffixes
ring_polymer = AutoPoly.Polymer(ChainNum=10, Sequence=["PMMA"]*50, topology="ring")

# Polymerization with ring topology
poly = AutoPoly.Polymerization(name="RingPolymer", system=system, model=[ring_polymer], run=True)
