#!/usr/bin/env python3
"""
Example: PEO Chains in Explicit Water — a Polymer Solution
===========================================================

Demonstrates mixing the Polymer and Molecule classes in a single system:
poly(ethylene oxide) chains solvated in explicit water. PEO/water is a
classic biocompatible polymer solution (hydrogels, drug delivery,
battery electrolytes).

  - Polymer: complement SMILES with [*] wildcards, chain_num + sequence
  - Molecule: regular SMILES, Count
  - Both are passed together as generate(system, name, [polymer, solvent])

Force field: GAFF covers both components; Gasteiger charges are assigned
automatically. For production runs consider AM1-BCC/RESP charges and a
water-specific model if quantitative aqueous properties matter.

Requires: pip install -e .  (from the AutoPoly repo root)
"""

from AutoPoly import System, Molecule, Polymer, generate

# Complement SMILES for PEO (methyl-terminated, matching the other examples)
PEO_FIRST = "CCO[*]"
PEO_MIDDLE = "[*]CCO[*]"
PEO_LAST = "[*]CCO"

DOP = 20
sequence = [PEO_FIRST] + [PEO_MIDDLE] * (DOP - 2) + [PEO_LAST]

system = System(out="peo_solution")

peo = Polymer(
    chain_num=5,
    sequence=sequence,
    topology="linear",
    tacticity="atactic",
)

water = Molecule(
    Count=200,
    Smiles="O",
    Name="water",
)

generate(
    system,
    "peo_water",
    [peo, water],   # polymer + solvent in one box
    force_field="gaff",
)

print()
print("Done! LAMMPS input files written to peo_solution/peo_water/")
print(f"System: 5 PEO chains (DOP={DOP}) + 200 water molecules")
print()
print("Note: the box is built at a low density to allow overlap-free")
print("placement. Equilibrate with NPT compression to reach the target")
print("density before production runs.")
