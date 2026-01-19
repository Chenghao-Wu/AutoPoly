#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example: Creating an ABA triblock copolymer with explicit monomer sequences.

This example demonstrates how to create block copolymers with precise
monomer placement at each position in the polymer chain.

The new explicit sequence API allows you to:
- Create block copolymers (e.g., ABA triblock)
- Specify exact monomer at each position
- Create arbitrary sequences

Created on 2026-01-13
@author: zwu
"""

from AutoPoly import System, Polymer, Polymerization

# Create system
system = System(out="aba_triblock")

# Define explicit monomer sequence
# ABA triblock: 2 PE-like, 3 PS-like, 2 PE-like
# This creates a block copolymer structure
sequence = [
    "[*]CC[*]",      # Position 0: Ethylene (Block A)
    "[*]CC[*]",      # Position 1: Ethylene (Block A)
    "[*]C=C[*]",     # Position 2: Styrene (Block B)
    "[*]C=C[*]",     # Position 3: Styrene (Block B)
    "[*]C=C[*]",     # Position 4: Styrene (Block B)
    "[*]CC[*]",      # Position 5: Ethylene (Block A)
    "[*]CC[*]"       # Position 6: Ethylene (Block A)
]

# Create polymer with explicit sequence
poly = Polymer(
    chain_num=5,
    sequence=sequence,  # DOP is automatically 7
    topology="linear",
    tacticity="atactic"
)

# Print polymer information
print(f"Created ABA triblock copolymer:")
print(f"  Number of chains: {poly.chain_num}")
print(f"  DOP (degree of polymerization): {poly.dop}")
print(f"  Unique monomer types: {poly.mer_set}")
print(f"  Topology: {poly.topology}")
print(f"  Tacticity: {poly.tacticity}")

# You can also use get_chain_info() for complete information
info = poly.get_chain_info()
print(f"\nComplete chain info:")
for key, value in info.items():
    print(f"  {key}: {value}")

# Note: Polymerization would be done here if needed
# polyz = Polymerization(
#     name="aba_triblock",
#     system=system,
#     model=[poly],
#     force_field="oplsaa",
#     run=True
# )
