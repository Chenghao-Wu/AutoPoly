# Substrates & Films

AutoPoly can build a polymer or small-molecule **film on top of a physical
substrate slab** — the standard fully periodic, two-interface slab model used
for surface, coating, and interfacial simulations:

```
zhi ┌───────────────────────────┐
    │    vacuum (optional)      │  spec.vacuum, 0 by default
    ├── ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ┤
    │    film region            │  your models, MC-placed with
    │                           │  collision detection
    ├── ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ┤  film_zmin = slab_top + gap
    │    gap                    │
    ├───────────────────────────┤  slab_top = zlo + thickness
    │    substrate slab         │  model-built (grid or mc) or
    │                           │  one external .lt instance
zlo └───────────────────────────┘
```

Passing `substrate=` to `generate()` auto-selects the `on_substrate` packing
strategy — no other changes to your workflow.

## Quick start

```python
from AutoPoly import System, Polymer, Molecule, generate
from AutoPoly.packing import SubstrateSpec

film = Polymer(chain_num=10, sequence=["CC[*]"] + ["[*]CC[*]"] * 18 + ["[*]CC"])

substrate = SubstrateSpec(
    model=Molecule(Count=64, Smiles="CCO", Name="etoh_sub"),
    thickness=10.0,   # slab z-extent (Å), must be > 0
    packing="grid",   # "grid" = ordered slab, "mc" = amorphous
    gap=3.0,          # empty space between slab top and film (Å)
)

generate(
    system, "pe_film", [film],
    force_field="gaff",
    substrate=substrate,
    box_dims=(50.0, 50.0, 50.0),
)
```

## Choosing the substrate source

### In-pipeline model: `SubstrateSpec(model=...)`

Any `Molecule` or `Polymer` works. It flows through the same three pipeline
stages as the film (same force field), and its units are tagged
`role="substrate"` in `units.json`. The slab is packed two ways:

| `packing` | Result |
|---|---|
| `"grid"` (default) | Ordered layers — instances on a 3D grid spanning the slab. Spacing is `1.8 ×` the unit's bounding radius; the run errors if your count exceeds grid capacity (grow `box_dims`/`thickness` or lower the count). |
| `"mc"` | Amorphous — MC placement with collision detection inside the slab region. Errors if the slab is too dense to place (grow the slab or use `"grid"`). |

The instance count comes from the model (`Count` / `chain_num`), an explicit
`count=N` override, or `count="auto"`:

```python
SubstrateSpec(
    model=Molecule(Count=1, Smiles="CCO", Name="etoh_sub"),  # Count ignored
    count="auto", density=0.01,   # particles/Å³
    thickness=10.0,
)
# -> count = density × lx × ly × thickness  (requires explicit box_dims lateral)
```

!!! note
    `density` here is *particles* per Å³ (whole molecules), not the
    `monomer_density` used for film sizing. For a small-molecule liquid,
    ~0.01/Å³ is a typical magnitude — 0.085 would be atomic density.

### External pre-built surface: `SubstrateSpec(lt_file=..., class_name=...)`

For crystalline surfaces (Au(111), SiO2, ...) built outside AutoPoly:

```python
SubstrateSpec(
    lt_file="surfaces/au111.lt", class_name="Au111",
    thickness=12.0, gap=3.0,
)
```

The file is copied into the moltemplate directory, imported in `system.lt`,
and instantiated **once**, centered laterally at the slab mid-plane. The
slab's atoms are not collision-checked individually — instead the film is
kept out geometrically (every film bounding-sphere center stays above
`slab_top + gap`). Make sure your external slab's coordinates are centered
on its own origin so the `.move()` placement lands it correctly, and that it
imports whatever force-field file it needs (it is *not* part of AutoPoly's
force-field subsetting).

## Box geometry: `box_dims`

Substrate systems use rectangular boxes. `box_dims=(lx, ly, lz)` takes
precedence per axis over the cubic `box_size`; any element may be `None`:

- **Lateral (`lx`, `ly`)** — auto: the standard melt box estimate.
- **`lz`** — auto: `thickness + gap + film_thickness + vacuum`, where
  `film_thickness = film_particles / (monomer_density × lx × ly)`.

Validation catches impossible stacks early: `lz` smaller than
`thickness + gap`, or a film region thinner than the largest film chain's
bounding-sphere diameter, both raise `ValidationError` with the numbers.

With the default `vacuum=0` the film fills the rest of the box, giving the
usual fully periodic slab model (film touches the slab's periodic image —
that's the second interface). Set `vacuum > 0` if you plan to run with
`boundary p p f`.

## Subtract: carving after placement

`subtract=` removes **whole instances** (chains/molecules) whose centers fall
inside a carve region. Because moltemplate instantiates whole `.lt` classes,
removal never cuts a covalent bond — the carved boundary is fuzzy at the
scale of one bounding radius (tighten it with `conservative=True`, which
also removes instances merely *touching* the region):

```python
from AutoPoly.packing import CutAbove, CutBelow, Cylinder, BoxRegion

generate(
    system, "pe_film_patterned", [film],
    substrate=substrate,
    box_dims=(50.0, 50.0, 50.0),
    subtract=[
        Cylinder(axis="z", center=(0, 0), radius=8.0),       # hole
        BoxRegion(x=(-25, -5), y=None, z=None),              # trench
        CutAbove(z=25.0),                                    # thickness trim
    ],
)
```

| Region | Carves |
|---|---|
| `CutAbove(z)` / `CutBelow(z)` | everything above/below a plane |
| `Cylinder(axis, center, radius)` | infinite cylinder (`axis` = "x"/"y"/"z") |
| `BoxRegion(x=, y=, z=)` | axis-aligned box; `None` = full span on that axis |

Every region takes `apply_to` — `"film"` (default), `"substrate"`, or
`"all"` — so a film pattern never touches the slab unless you ask. Regions
apply in order; every removal is logged.

Subtract is strategy-agnostic: it also works with plain `mc_random` melts
(`generate(..., subtract=[Cylinder(...)])` → nanopore). The `grid` strategy
rejects it.

## What subtract does *not* do

Whole-instance granularity means no dangling bonds, no charge drift, and a
topology that always matches the `.lt` definitions — at the price of a
carved surface that is smooth only at the chain-radius scale. If your model
needs atomically sharp interfaces (severed bonds, capped fragments), that is
a different, post-processing operation on `system.data`; it is not part of
this feature.

## See also

- [example_film_on_substrate.py](https://github.com/WuGroup-XJTLU/AutoPoly/blob/v2.0/examples/example_film_on_substrate.py) — runnable film + patterned variant
- [MC Placement & Chain Growth](mc-placement.md) — the placement engine the film uses
- [generate API](../reference/generate.md) — `substrate`, `subtract`, `box_dims` parameters
