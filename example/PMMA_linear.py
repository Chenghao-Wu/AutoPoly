#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 12:19:08 2018

@author: zwu
"""

import AutoPoly

# Define the system
# out is the folder name for the output
system = AutoPoly.System(out="test2")
system.made_folder = True

# Create a ring polymer with 10 PMMA monomers
# Note: Using just "PMMA" as the base name - the code will append the necessary suffixes
linear_polymer = AutoPoly.Polymer(ChainNum=10, Sequence=["PMMA"]*50, topology="linear")

# Polymerization with ring topology
poly = AutoPoly.Polymerization(name="LinearPolymer", system=system, model=[linear_polymer], run=True)
