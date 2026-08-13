#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example: Polymer film on a built-in beta-cristobalite(111) silica substrate.

beta-cristobalite(111) carries isolated Q3 silanols at ~4.5/nm^2 — the
crystalline face whose silanol density matches hydroxylated amorphous
silica (Zhuravlev 2000). The lateral box is snapped to integer surface
cells (10.1258 x 17.5383 A) so the slab is seamlessly periodic.

The slab .lt is self-typed (LJ + charges, no bonded terms): freeze or
`fix rigid` the slab atoms in MD.

Run:
    python example_film_on_cristobalite.py
"""
from AutoPoly import System, Polymer, generate
from AutoPoly.packing import SubstrateSpec

PE_SEQUENCE = ["CC[*]"] + ["[*]CC[*]"] * 8 + ["[*]CC"]  # DOP 10

system = System(out="film_on_cristobalite_out")

film = Polymer(chain_num=3, sequence=PE_SEQUENCE, tacticity="atactic")

substrate = SubstrateSpec(
    builder="beta_cristobalite",
    thickness=13.0,        # slab envelope (A), incl. hydroxyl coatings
    gap=3.0,
    slab_ff="interface",   # INTERFACE FF v1.5; "clayff" also built in
    oh_density=4.6,        # ~ the intrinsic Q3 density of this face
)

generate(
    system, "pe_on_cristobalite", [film],
    force_field="gaff",
    substrate=substrate,
    box_dims=(55.7, 52.7, 57.0),   # lateral snapped to 6x3 surface cells;
                               # must fit the chains' bounding spheres
                               # (r = 12.3 A for DOP-10 PE)
    rng_seed=42,
)

print("Done — see film_on_cristobalite_out/pe_on_cristobalite*/ for "
      "system.data etc.")
