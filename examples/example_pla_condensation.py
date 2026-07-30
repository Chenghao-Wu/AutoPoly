#!/usr/bin/env python3
"""
Example: PLA — a Condensation Polymer with GAFF
================================================

Polylactic acid (PLA) is a biodegradable polyester formed by step-growth
condensation, unlike the vinyl (chain-growth) polymers in most other
examples.

  Vinyl addition (PE, PP, PS, PMMA):
      all-carbon backbone, no byproducts
  Condensation (PLA, Nylon, PET):
      heteroatom backbone (C-O ester links), small molecule byproducts

  PLA:  n HO-CH(CH3)-COOH -> [-O-CH(CH3)-C(=O)-]n + n H2O

Complement SMILES for the PLA repeat unit:
  first:  "OC(C)C(=O)[*]"      alcohol end, connects via the carbonyl C
  middle: "[*]OC(C)C(=O)[*]"   ester repeat unit
  last:   "[*]OC(C)C(=O)O"     carboxylic acid end

Force field: GAFF is recommended for polyesters — its ester parameters
are well validated. Gasteiger charges are assigned automatically.

Requires: pip install -e .  (from the AutoPoly repo root)
"""

from AutoPoly import System, Polymer, generate

PLA_FIRST = "OC(C)C(=O)[*]"
PLA_MIDDLE = "[*]OC(C)C(=O)[*]"
PLA_LAST = "[*]OC(C)C(=O)O"

DOP = 20
sequence = [PLA_FIRST] + [PLA_MIDDLE] * (DOP - 2) + [PLA_LAST]

system = System(out="pla_condensation")

pla = Polymer(
    chain_num=5,
    sequence=sequence,
    topology="linear",
    tacticity="atactic",
)

generate(
    system,
    "pla",
    [pla],
    force_field="gaff",  # well-validated ester parameters
)

print()
print("Done! LAMMPS input files written to pla_condensation/pla/")
print()
print("PLA reference values for validation:")
print("  Density: ~1.24-1.26 g/cm^3")
print("  Tg:      ~330 K (57 degC)")
