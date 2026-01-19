#!/usr/bin/env python3
"""
Example: Automatic Gasteiger Charges for PMMA

Demonstrates that GAFF force field now works automatically
with Gasteiger charges - no manual configuration needed.

This example creates PMMA (Poly(methyl methacrylate)) using the
GAFF force field. Gasteiger charges are calculated automatically
during the LT file generation process.

Before this implementation:
- GAFF charges would default to 0.0
- System would not run correctly

After this implementation:
- Gasteiger charges are calculated automatically
- System is ready for simulation
"""

from AutoPoly import System, Polymer, Polymerization

# Create system
system = System(out="pmma_gasteiger")

# Define PMMA polymer
# SMILES: [*]CC([*])(C)C(=O)OC
# - [*] marks connection points for polymerization
# - Atactic (random stereochemistry)
polymer = Polymer(
    chain_num=2,  # Number of polymer chains
    sequence=["[*]CC([*])(C)C(=O)OC"] * 10,  # 10 monomers per chain
    topology="linear",  # Linear polymer
    tacticity="atactic"  # Random stereochemistry
)

# Run polymerization with GAFF
# Gasteiger charges are calculated automatically!
polymerization = Polymerization(
    name="pmma_gasteiger",
    system=system,
    model=[polymer],
    force_field="gaff"  # GAFF force field - no special configuration needed
)

print("✓ PMMA polymerization with GAFF complete!")
print("✓ Gasteiger charges were calculated automatically")
print(f"✓ Output files in: {system.out}")
print()
print("To verify non-zero charges:")
print(f"  grep -A 30 'Data Atoms' {system.out}/output_ttree/system.in.charges | head -35")
