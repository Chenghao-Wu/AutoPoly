#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reactor example: step-growth polymerization of a monomer melt.

Builds a melt of ethylene glycol + adipic acid and prepares LAMMPS
``fix bond/react`` inputs so the two monomers polymerize into a polyester
during MD (with water split off as a byproduct).

The reactor:
  1. detects the diol + di-carboxylic-acid polyesterification reaction,
  2. builds the pre/post reaction molecule templates + map file (typed with
     the same force field as the system),
  3. writes supplementary force-field parameters for the types the reaction
     creates (ester linkage, rehybridized carbonyl, water),
  4. writes a ready-to-run ``in.bond_react`` input script.

Run the generated script from the project directory:

    cd reactor_out/melt
    lmp -in in.bond_react
"""
from AutoPoly import System, Molecule, generate, Reactor

# 1. Build a monomer melt (no pre-built chains; the monomers react in MD)
system = System(out="reactor_out")
eg = Molecule(Count=20, Smiles="OCCO", Name="eg")
adipic = Molecule(Count=20, Smiles="O=C(O)CCCCC(=O)O", Name="adipic")

generate(system, "melt", [eg, adipic], force_field="gaff2")

# 2. Detect the polymerization reactions between the packed monomers
project_dir = system.get_folder_path() + "/melt"
reactor = Reactor(project_dir, monomers=[eg, adipic])
instances = reactor.detect_reactions()
print("Detected reactions:")
for inst in instances:
    print("  -", inst.description)

# 3. Build the fix bond/react templates + input script
result = reactor.build(temperature=300.0, nevery=100, rmin=0.0, rmax=3.5)

print("\nBuilt reaction template sets:")
for name in result.reactions:
    print("  -", name)
print("\nFiles written:")
for f in result.template_files:
    print("  ", f)
print(f"\nRun with:  cd {project_dir} && lmp -in in.bond_react")
