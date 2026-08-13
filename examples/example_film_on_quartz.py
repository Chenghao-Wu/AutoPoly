#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example: Polymer film on a built-in alpha-quartz(0001) silica substrate.

Builds a PE film on a hydroxylated crystalline quartz slab (generated at
pack time, INTERFACE FF types) — the standard fully periodic, two-interface
slab model. The lateral box is snapped to integer quartz surface cells
(4.9019 x 8.4903 A) so the slab is seamlessly periodic.

The slab .lt is self-typed (LJ + charges, no bonded terms): freeze or
`fix rigid` the slab atoms in MD. Both faces are hydroxylated (geminal
Q2 silanols, ~9.6/nm^2 — the intrinsic alpha-quartz(0001) termination).

Run:
    python example_film_on_quartz.py
"""
from AutoPoly import System, Polymer, generate
from AutoPoly.packing import SubstrateSpec

PE_SEQUENCE = ["CC[*]"] + ["[*]CC[*]"] * 8 + ["[*]CC"]  # DOP 10

system = System(out="film_on_quartz_out")

film = Polymer(chain_num=3, sequence=PE_SEQUENCE, tacticity="atactic")

substrate = SubstrateSpec(
    builder="alpha_quartz",
    thickness=12.0,        # slab envelope (A), incl. hydroxyl coatings
    gap=3.0,
    slab_ff="interface",   # INTERFACE FF v1.5; "clayff" also built in
    # oh_density=4.6,      # target silanol density (note: quartz(0001) is
                           # intrinsically Q2, ~9.6/nm^2; see docs)
    # hydroxylate_bottom=False,  # leave the bottom face bare
)

generate(
    system, "pe_on_quartz", [film],
    force_field="gaff",
    substrate=substrate,
    box_dims=(54.0, 51.0, 57.0),   # lateral snapped to 11x6 quartz cells;
                               # must fit the chains' bounding spheres
                               # (r = 12.3 A for DOP-10 PE)
    rng_seed=42,
)

print("Done — see film_on_quartz_out/pe_on_quartz*/ for system.data etc.")
