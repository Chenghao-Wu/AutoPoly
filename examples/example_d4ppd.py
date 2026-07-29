#!/usr/bin/env python3
"""
Example: Generate a single D4PPD molecule with GAFF2

D4PPD is a diaryl-p-phenylenediamine antioxidant used in rubber and
lubricant formulations. This example demonstrates the Molecule class
with a larger, drug-like organic molecule and the GAFF2 force field
(extended GAFF with additional atom types such as "nq").

Requires: pip install -e .  (from the AutoPoly repo root)
"""

from AutoPoly import System, Molecule, Polymerization

# D4PPD: a diaryl-p-phenylenediamine antioxidant
D4PPD_SMILES = "CC(C)Nc1c(C)cc(Nc2ccc(C)cc2)cc1"

system = System(out="d4ppd")

d4ppd = Molecule(
    Count=1,
    Smiles=D4PPD_SMILES,
    Name="d4ppd",
)

Polymerization(
    name="d4ppd",
    system=system,
    model=[d4ppd],
    force_field="gaff2",  # GAFF2 (extended atom types, e.g. nq)
)

print("D4PPD system created successfully!")
print("Output directory: d4ppd/d4ppd/")
print(f"SMILES: {D4PPD_SMILES}")
