#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 12:19:08 2018

@author: zwu
"""

import AutoPoly

# Define the system
# out is the folder name for the output
system = AutoPoly.System(out="LinearPLA")
system.made_folder = True

# Create a ring polymer with 10 PMMA monomers
# Note: Using just "PMMA" as the base name - the code will append the necessary suffixes
linear_polymer = AutoPoly.Polymer(ChainNum=2, Sequence=["PLA"]*3, topology="linear")

# Polymerization with ring topology
poly = AutoPoly.Polymerization(name="LinearPLA", system=system, model=[linear_polymer], run=True)
