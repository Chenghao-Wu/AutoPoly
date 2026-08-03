#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example: Polymer film on a physical substrate, with optional carve subtract.

Builds a PE film on top of an ordered small-molecule slab (fully periodic,
two-interface slab model), and shows the subtract post-pass (whole-instance
removal — no covalent bonds are cut).

Run:
    python example_film_on_substrate.py
"""
from AutoPoly import System, Polymer, Molecule, generate
from AutoPoly.packing import SubstrateSpec, CutAbove, Cylinder
from AutoPoly.pipeline.geometry import GeometryConfig

PE_SEQUENCE = ["CC[*]"] + ["[*]CC[*]"] * 8 + ["[*]CC"]  # DOP 10

system = System(out="film_on_substrate_out")

film = Polymer(chain_num=6, sequence=PE_SEQUENCE, tacticity="atactic")

# Model-built substrate: 200 ethanol instances grid-packed into a 10 A slab
# at the bottom of the box. (For a crystalline surface, pass an external
# slab instead: SubstrateSpec(lt_file="au111.lt", class_name="Au111", ...).)
substrate = SubstrateSpec(
    model=Molecule(Count=64, Smiles="CCO", Name="etoh_sub"),  # or count="auto" w/ density
    thickness=10.0,   # slab z-extent in Angstrom
    packing="grid",   # ordered slab; "mc" gives an amorphous one
    gap=3.0,          # empty space between slab top and film
)

generate(
    system, "pe_film", [film],
    force_field="gaff",
    substrate=substrate,               # auto-selects strategy="on_substrate"
    box_dims=(50.0, 50.0, 50.0),       # explicit lz (or None: slab+gap+film at monomer_density)
    geometry_config=GeometryConfig(use_mc_chain_growth=False),
    rng_seed=42,
)

# Variant: pattern the film after placement (whole-instance subtract).
# A cylindrical hole through the film; substrate untouched (apply_to="film"
# is the default). Use CutAbove(z=...) to trim the film to a target
# thickness, or apply_to="substrate" to carve the slab.
generate(
    system, "pe_film_patterned", [film],
    force_field="gaff",
    substrate=substrate,
    box_dims=(50.0, 50.0, 50.0),
    subtract=[Cylinder(axis="z", center=(0.0, 0.0), radius=8.0)],
    geometry_config=GeometryConfig(use_mc_chain_growth=False),
    rng_seed=42,
)

print("Done — see film_on_substrate_out/pe_film*/ for system.data etc.")
