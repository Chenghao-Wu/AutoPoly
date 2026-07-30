#!/usr/bin/env python3
"""
Three-Stage Pipeline Example: Geometry -> Typing -> Packing

Demonstrates the stage-level API introduced with the three-stage pipeline:

    Stage 1  GeometryBuilder  models -> geometry/geometry.json (FF-agnostic)
    Stage 2  UnitTyper        geometry + force field -> build/<ff>/ + units.json
    Stage 3  BoxPacker        units -> moltemplate/ -> system.data

Key capability: ONE geometry can be typed under MULTIPLE force fields.
Typing is the cheap stage, so comparing force fields no longer means
rebuilding every chain (6x less conformer work for a 6-FF comparison).

Requires: pip install -e .  (from the AutoPoly repo root)
"""

from AutoPoly import (
    System,
    Polymer,
    Molecule,
    GeometryBuilder,
    GeometryConfig,
    UnitTyper,
    BoxPacker,
)

FORCE_FIELDS = ["oplsaa", "gaff2"]  # extend as needed


def main():
    system = System(out="peo_staged")

    # PEO, 5 chains x 10 monomers
    sequence = ["CCO[*]"] + ["[*]CCO[*]"] * 8 + ["[*]CCO"]
    polymer = Polymer(
        chain_num=5,
        sequence=sequence,
        topology="linear",
        tacticity="atactic",
    )

    # ------------------------------------------------------------------
    # Stage 1: build the geometry ONCE (conformers + chain placements).
    # No force field is involved at this stage.
    # ------------------------------------------------------------------
    geometry = GeometryBuilder(
        system,
        name="peo",
        config=GeometryConfig(use_mc_chain_growth=True, rng_seed=42),
    ).build([polymer])
    print(f"Geometry written to {geometry.dir}")

    for ff in FORCE_FIELDS:
        # --------------------------------------------------------------
        # Stage 2: type the SAME geometry under each force field.
        # Produces build/<ff>/*.lt + units.json.
        # --------------------------------------------------------------
        units = UnitTyper(geometry.dir, ff).type()
        print(f"Typed {len(units.units)} units under {ff} -> {units.source_dir}")

        # --------------------------------------------------------------
        # Stage 3: pack into a box and run moltemplate.
        # (Re-running pack overwrites <name>/moltemplate/, so in practice
        #  pick one force field per packed system, or copy the output.)
        # --------------------------------------------------------------
        BoxPacker(system, "peo", strategy="mc_random", rng_seed=42).pack(units)
        print(f"Packed with {ff}: peo_staged/peo/system.data")


if __name__ == "__main__":
    main()
