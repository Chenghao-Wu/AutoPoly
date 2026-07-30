# API Reference

AutoPoly's public API, generated from the source docstrings.

## Core classes

| Class | Module | Purpose |
|---|---|---|
| [`System`](system.md) | `AutoPoly.system` | Output directory and path management |
| [`Polymer`](polymer.md) | `AutoPoly.polymer` | Polymer chains from complement SMILES sequences |
| [`Molecule`](molecule.md) | `AutoPoly.molecule` | Small molecules (solvents, additives) |
| [`Polymerization`](polymerization.md) | `AutoPoly.polymerization` | Full pipeline: templates → moltemplate → LAMMPS files |

## Coarse-grained

| Class | Module | Purpose |
|---|---|---|
| [`BeadSpringPolymer`](bead-spring.md) | `AutoPoly.bead_spring` | Direct LAMMPS data files for bead-spring models |
| [`BeadType` / `AngleType`](bead-spring.md) | `AutoPoly.bead_spring` | Bead species and angle parameters |
| [`MCConfig` / `SAWConfig`](bead-spring.md) | `AutoPoly.bead_spring` | Generation and equilibration tuning |

## Building blocks

| Module | Purpose |
|---|---|
| [`AutoPoly.monomer_generator`](monomer-generator.md) | SMILES → moltemplate `.lt` monomer templates |
| [`AutoPoly.mc`](mc.md) | Monte Carlo placement: collision detection, SAW chain growth |
| [`AutoPoly.exceptions`](exceptions.md) | Exception hierarchy |

## Agent interface

| Module | Purpose |
|---|---|
| [`AutoPoly.agent`](agent.md) | Config-driven JSON API: `info`, `validate`, `generate`, `describe_smiles`, `suggest_force_field` |
| [`AutoPoly.tools`](tools.md) | LangChain tool wrappers |
| [`AutoPoly.cli`](cli.md) | The `autopoly` command-line interface |

!!! note "RDKit-dependent modules"
    `Polymer`, `Molecule`, `Polymerization`, `MonomerGenerator`, and `mc` require RDKit. `System`, `BeadSpringPolymer`, and `agent` work without it.
