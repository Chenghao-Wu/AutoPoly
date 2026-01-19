#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example: Creating uniform polymers using a helper function.

For long uniform chains, you can use a simple helper function to
create the sequence list. This demonstrates the explicit sequence
API where you specify each monomer position.

Created on 2026-01-13
@author: zwu
"""

from AutoPoly import System, Polymer, Polymerization


def create_uniform_sequence(smiles: str, length: int) -> list:
    """
    Create a uniform sequence of given length.

    Args:
        smiles: pSMILES string with [*] connection points
        length: Number of monomer units

    Returns:
        List of SMILES strings

    Example:
        >>> sequence = create_uniform_sequence("[*]CC[*]", 50)
        >>> len(sequence)
        50
    """
    return [smiles] * length


# Create system
system = System(out="uniform_pe")

# Create polyethylene with 50 units using helper function
sequence = create_uniform_sequence("[*]CC[*]", 50)

poly = Polymer(
    chain_num=10,
    sequence=sequence,  # DOP is automatically 50
    topology="linear",
    tacticity="isotactic"
)

# Print polymer information
print(f"Created uniform polyethylene:")
print(f"  Number of chains: {poly.chain_num}")
print(f"  DOP (degree of polymerization): {poly.dop}")
print(f"  Unique monomer types: {len(poly.mer_set)}")
print(f"  Topology: {poly.topology}")
print(f"  Tacticity: {poly.tacticity}")

# You can also use get_chain_info() for complete information
info = poly.get_chain_info()
print(f"\nComplete chain info:")
for key, value in info.items():
    if key not in ['sequence_set', 'sequence_names', 'tacticity_set']:
        print(f"  {key}: {value}")

# Note: Polymerization would be done here if needed
# polyz = Polymerization(
#     name="polyethylene",
#     system=system,
#     model=[poly],
#     force_field="oplsaa",
#     run=True
# )
