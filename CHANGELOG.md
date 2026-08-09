# Changelog

All notable changes to AutoPoly are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Bead-spring architectures (graph core)** — bead-spring chains are now
  built on an explicit graph representation (`BeadArchitecture`: nodes =
  beads, edges = bonds) instead of an implicit linear path, so arbitrary
  topologies are supported uniformly:
  - New module `AutoPoly.models.architectures` with factories `linear`,
    `ring`, `star` (incl. miktoarm/asymmetric arms), `comb`, `graft`
    (explicit graft points / side groups), `tadpole` (ring + tail),
    `dendrimer`, and `custom` (explicit bead + bond lists).
  - `MonomerTemplate` + `BeadArchitecture.from_monomers(...)` for
    explicitly defined multi-bead monomers (backbone bead, side-group
    bead) with named connection points (`head`/`tail`/`side`).
  - Copolymer sequence generators: `block_sequence`,
    `alternating_sequence`, `random_sequence` (seeded, weighted),
    `gradient_sequence`; usable in any sequence slot.
  - `BeadSpringPolymer(..., architecture=...)` accepts any architecture;
    legacy `sequence=`/`topology=` remain fully backward compatible.
  - Branch-aware angles: triplets centered on branch points are included
    by default with their own canonical triplet types (configurable via
    `AngleType`); disable with `include_branch_angles=False`.
  - Graph-based SAW generation (`saw_grow_graph`, `saw_generate_graphs`):
    DFS spanning-tree growth with cycle-closing constraints.
  - Branched MC equilibration: tree-pivot moves across bridge edges
    (`mc_tree_pivot_move`) and segment crankshaft
    (`mc_segment_crankshaft_move`); reptation is linear-only.
- **Bead-spring mixtures** — new `BeadSpringSystem` packs multiple
  species (`(architecture, n_chains)` pairs — e.g. rings + linear +
  combs) into one box and one LAMMPS data file with a shared bead-type
  table.
- **Moltemplate is now the default backend for bead-spring models** —
  new unified `generate()` entry point on `BeadSpringPolymer` and
  `BeadSpringSystem` dispatches to the moltemplate backend by default
  (constructor `backend="moltemplate"`; override per call with
  `generate(backend="direct")` for the lightweight direct writer, e.g.
  very large melts). `generate_moltemplate()` and `generate_data_file()`
  remain available for explicit control.
- **SAW whole-system retries** — `SAWConfig.system_retries` (default 3)
  retries the whole multi-chain SAW generation when a single crowded
  chain fails, instead of immediately falling back to geometric
  placement.
- **Example: bead-spring with side groups** —
  `examples/example_bead_spring_side_groups.py` demonstrates a comb with
  single-bead side groups + branch angles, a graft copolymer with
  oligomeric side chains, a comb built from an explicit
  `MonomerTemplate`, and the moltemplate backend.
- **Moltemplate backend for bead-spring models** — `generate_moltemplate()`
  on both `BeadSpringPolymer` and `BeadSpringSystem` emits
  `bead_spring.lt` (CG force field), `bead_<Type>.lt` monomer objects,
  `chains.lt` (one object per chain at generated coordinates, explicit
  bond list, typed angle list), and `system.lt`, then optionally runs the
  bundled moltemplate to produce `system.data` + `system.in.init/settings`.

- **Physical substrates** — build a polymer/molecule film on top of a
  substrate slab via `generate(..., substrate=SubstrateSpec(...))`:
  - `SubstrateSpec(model=...)` packs a `Molecule`/`Polymer` into a slab of
    given `thickness` at the bottom of the box (`packing="grid"` for an
    ordered slab, `"mc"` for amorphous), with a `gap` to the film region.
    The substrate is typed in-pipeline alongside the film models;
    `UnitSpec.role` ("film"/"substrate") records the partition.
  - `SubstrateSpec(lt_file=..., class_name=...)` instantiates a pre-built
    external surface (e.g. an Au(111) or SiO2 slab) once, centered
    laterally at the slab mid-plane.
  - `count="auto"` derives a `Molecule` substrate's instance count from
    slab volume × `density` (requires explicit lateral `box_dims`).
- **New packing strategy `on_substrate`** (auto-selected when `substrate`
  is passed): z-layered assembly — slab at the bottom, film MC-placed
  above `slab_top + gap`, shared collision detector so film and slab
  never interpenetrate. Box sides are set per axis with
  `box_dims=(lx, ly, lz)` (any element `None` = auto;
  `lz = thickness + gap + film thickness at monomer_density + vacuum`).
- **Subtract (carve regions)** — whole-instance removal after placement,
  no covalent bonds cut: `CutAbove`, `CutBelow`, `Cylinder`, `BoxRegion`
  from `AutoPoly.packing`, passed as `generate(..., subtract=[...])`.
  Each region has `apply_to` ("film"/"substrate"/"all") and an optional
  `conservative=True` sphere-inflated boundary. Supported by the
  `mc_random` and `on_substrate` strategies (e.g. carve a `Cylinder`
  through a bulk melt to build a nanopore); removals are logged.
- `BoxSpec.box_dims` for rectangular (non-cubic) boxes.

## [2.0.0] - 2026-07-30

AutoPoly 2.0 is a ground-up rearchitecture of the generation pipeline. It is a
**breaking release**: the `Polymerization` class, the agent/CLI surface, and
loose monomer inputs are gone. See the migration guide below.

### Breaking changes

| 1.x | 2.0 | Notes |
|---|---|---|
| `Polymerization(name=..., system=..., model=..., force_field=...)` | `generate(system, name, models, force_field=...)` | One-shot function composing the three pipeline stages; returns a `PlacementResult` |
| `Polymerization(..., run=False)` then manual method calls | `GeometryBuilder` → `UnitTyper` → `BoxPacker` | The stages are now public, independently callable APIs |
| `Polymerization(placement_method="mc_random"/"grid")` | `generate(..., strategy="mc_random"/"grid")` | Strategies are pluggable via `register_strategy()` |
| `Polymerization(mc_monomer_density=...)` | `generate(..., monomer_density=...)` | Renamed |
| `Polymerization(mc_bond_angle_min/max=..., mc_intrachain_exclude_neighbors=...)` | `generate(..., geometry_config=GeometryConfig(...))` | Chain-growth parameters live in `GeometryConfig` |
| `Polymer(chain_num=..., sequence=["PE", "PS"], ...)` with monomer **names** | `Polymer(chain_num=..., sequence=["CC[*]", "[*]C=C"], ...)` with **pSMILES** | Sequences must be complement SMILES with explicit wildcards |
| Implicit DOP / sequence cycling | `DOP = len(sequence)`, no cycling | The sequence you pass is the chain you get |
| `AutoPoly.agent`, `AutoPoly.cli`, `AutoPoly.tools` (LangChain agent & CLI) | removed | Library-only package; drive it from Python |
| `from AutoPoly.polymer import Polymer` (flat layout) | `from AutoPoly import Polymer` (unchanged) or `AutoPoly.models.polymer` | Package reorganized into `core/`, `models/`, `monomers/`, `forcefields/`, `pipeline/`, `packing/`, `mc/` subpackages |

### Stricter pSMILES validation

`Polymer` now validates wildcard counts per position at construction time:

- **First monomer:** exactly 1 wildcard (right connection), e.g. `"CC[*]"`
- **Middle monomers:** exactly 2 wildcards, e.g. `"[*]CC[*]"`
- **Last monomer:** exactly 1 wildcard (left connection), e.g. `"[*]CC"`
- **Single-monomer sequence (DOP=1):** 0 wildcards, e.g. `"CC"`

Monomer names such as `"PE"` or `"PS"` are no longer accepted and raise
`ValidationError` with a message identifying the offending position.

### Added

- **Three-stage pipeline**, each stage usable standalone:
  - `GeometryBuilder` (stage 1): force-field-agnostic chain growth → `geometry/geometry.json`
  - `UnitTyper` (stage 2): atom typing → `build/<ff>/*.lt` + `units.json`; type one geometry for multiple force fields without regrowing chains
  - `BoxPacker` (stage 3): packing + moltemplate → `system.data`, `system.in.*`
- `generate()` one-shot convenience function (`AutoPoly.pipeline.workflow`)
- `UnitLibrary` / `UnitSpec` manifest contract between stages 2 and 3
- Pluggable packing strategies: `PlacementStrategy`, `register_strategy()`, `get_strategy()`; built-ins `mc_random` and `grid`
- `GeometryConfig` for MC chain-growth parameters
- GAFF `gaff.lt` now covers the `nu`/`nv` amine-to-aromatic nitrogen types
  (parameters aliased from `nh`), closing a runtime gap where rdlt could
  assign types that moltemplate could not parameterize

### Removed

- `AutoPoly.agent` (LangChain agent), `AutoPoly.cli`, `AutoPoly.tools`
- `Polymerization` class (superseded by `generate()` and the stage APIs)
- Implicit DOP parameter and monomer cycling in `Polymer`

### Migration

**1.x:**

```python
from AutoPoly import System, Polymer, Polymerization

system = System(out="my_polymer")
poly = Polymer(chain_num=4, sequence=["PE"] * 50)
Polymerization(name="pe", system=system, model=[poly], force_field="oplsaa")
```

**2.0:**

```python
from AutoPoly import System, Polymer, generate

system = System(out="my_polymer")
sequence = ["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"]  # DOP = 50
poly = Polymer(chain_num=4, sequence=sequence)
generate(system, "pe", [poly], force_field="oplsaa")
```

For multi-force-field work, build the geometry once and type it per force field:

```python
from AutoPoly import GeometryBuilder, UnitTyper, BoxPacker

GeometryBuilder(system, "pe").build([poly])  # writes <out>/pe/geometry/
for ff in ("oplsaa", "gaff", "gaff2"):
    units = UnitTyper("<out>/pe/geometry", force_field=ff).type()
    BoxPacker(system, "pe").pack(units)
```

See the [quickstart](docs/getting-started/quickstart.md) and
[workflow guide](docs/guides/workflow.md) for full examples.

## [1.0.0]

Initial stable release: `Polymerization` workflow, OPLS-AA/GAFF/GAFF2/l-OPLS
typing, MC chain growth and placement, bead-spring models, LangChain agent
and CLI.
