#!/usr/bin/env python3
"""
Example: Create a system of 100 benzene molecules

Demonstrates the Molecule class for building a pure molecular liquid
for LAMMPS simulations. Benzene (C6H6) is defined with regular SMILES
notation (no [*] wildcards — those are only for polymers).

GAFF is recommended for small organic molecules with aromatic rings.

Requires: pip install -e .  (from the AutoPoly repo root)
"""

from AutoPoly import System, Molecule, Polymerization

system = System(out="benzene_system_100")

benzene = Molecule(
    Count=100,               # Number of benzene molecules
    Smiles="c1ccccc1",       # Aromatic ring
    Name="benzene",
)

Polymerization(
    name="benzene_100",
    system=system,
    model=[benzene],
    force_field="gaff",
)

print("Benzene system created successfully!")
print("Output directory: benzene_system_100/benzene_100/")
print("Total atoms: 100 x 12 = 1200 (C6H6)")
