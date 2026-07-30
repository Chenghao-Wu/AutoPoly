#!/usr/bin/env python3
"""
Beginner Tutorial: Linear PMMA from pSMILES
============================================

This is the recommended starting point for new AutoPoly users. It walks
through the complete workflow for generating linear poly(methyl
methacrylate) — PMMA — using complement SMILES (pSMILES) notation.

Key concepts:
  - Complement SMILES: monomers are defined with [*] wildcards marking
    the atoms that connect to neighbors during polymerization.
      first  variant: 1 wildcard (right)  -> chain start
      middle variant: 2 wildcards         -> chain interior
      last   variant: 1 wildcard (left)   -> chain end
  - From these, AutoPoly's MonomerGenerator automatically creates the
    moltemplate .lt files (including mirror-tacticity _T1 variants).
  - OPLS-AA is the recommended force field for vinyl polymers.

Workflow:
  System -> Polymer -> generate -> moltemplate -> LAMMPS files

Output (pmma_tutorial/pmma/):
  system.data          topology + coordinates (LAMMPS read_data)
  system.in.init       units / atom / bond styles
  system.in.settings   force-field parameters
  system.in.charges    atomic charges (verify before production runs)

Requires: pip install -e .  (from the AutoPoly repo root)
"""

from AutoPoly import System, Polymer, generate

# ---------------------------------------------------------------------------
# Step 1: Complement SMILES for methyl methacrylate (vinyl addition, C-C backbone)
#
#   Monomer:  CH2=C(CH3)COOCH3
#   Repeat:   [-CH2-C(CH3)(COOCH3)-]
#
# The two [*] wildcards mark the backbone carbons that link repeat units.
# ---------------------------------------------------------------------------
PMMA_FIRST = "CC(C)(C(=O)OC)[*]"      # chain start (connects on the right)
PMMA_MIDDLE = "[*]CC([*])(C)C(=O)OC"  # interior repeat unit
PMMA_LAST = "[*]CC(C)(C(=O)OC)"       # chain end (connects on the left)

# ---------------------------------------------------------------------------
# Step 2: Build the explicit monomer sequence for one chain.
# DOP (degree of polymerization) is simply the sequence length.
# Start small while testing: 4 chains x 10 monomers.
# ---------------------------------------------------------------------------
DOP = 10
sequence = [PMMA_FIRST] + [PMMA_MIDDLE] * (DOP - 2) + [PMMA_LAST]

# ---------------------------------------------------------------------------
# Step 3: Create the output System and the Polymer model.
# ---------------------------------------------------------------------------
system = System(out="pmma_tutorial")

pmma = Polymer(
    chain_num=4,            # 4 chains in the simulation box
    sequence=sequence,      # DOP=10 derived from sequence length
    topology="linear",      # distinct chain start/end (use "ring" for cyclic)
    tacticity="atactic",    # random stereochemistry (commercial PMMA)
)

info = pmma.get_chain_info()
print("Polymer defined:")
for key, value in info.items():
    print(f"  {key}: {value}")

# ---------------------------------------------------------------------------
# Step 4: Run the generation pipeline.
# This generates monomer .lt templates, grows the chains, invokes
# moltemplate, and writes the LAMMPS input files.
# ---------------------------------------------------------------------------
generate(
    system,
    "pmma",
    [pmma],
    force_field="oplsaa",   # recommended for vinyl polymers
)

print()
print("Done! LAMMPS input files written to pmma_tutorial/pmma/")
print()
print("Before running LAMMPS, remember to:")
print("  1. Check system.in.charges — Gasteiger charges are assigned")
print("     automatically, but AM1-BCC or RESP charges are more accurate.")
print("  2. Equilibrate: minimize -> NVT heating -> NPT compression.")
print("  3. Validate density (~1.18 g/cm^3 for amorphous PMMA).")
