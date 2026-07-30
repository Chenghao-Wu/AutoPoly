# Changelog

All notable changes to AutoPoly are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
