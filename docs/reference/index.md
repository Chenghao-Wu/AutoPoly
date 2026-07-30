# API Reference

AutoPoly's public API, generated from the source docstrings.

## Core classes

| Class | Module | Purpose |
|---|---|---|
| [`System`](system.md) | `AutoPoly.core.system` | Output directory and path management |
| [`Polymer`](polymer.md) | `AutoPoly.models.polymer` | Polymer chains from complement SMILES sequences |
| [`Molecule`](molecule.md) | `AutoPoly.models.molecule` | Small molecules (solvents, additives) |
| [`generate`](generate.md) | `AutoPoly.pipeline.workflow` | One-shot pipeline: geometry → typing → packing → LAMMPS files |

## Coarse-grained

| Class | Module | Purpose |
|---|---|---|
| [`BeadSpringPolymer`](bead-spring.md) | `AutoPoly.models.bead_spring` | Direct LAMMPS data files for bead-spring models |
| [`BeadType` / `AngleType`](bead-spring.md) | `AutoPoly.models.bead_spring` | Bead species and angle parameters |
| [`MCConfig` / `SAWConfig`](bead-spring.md) | `AutoPoly.models.bead_spring` | Generation and equilibration tuning |

## Building blocks

| Module | Purpose |
|---|---|
| [`AutoPoly.monomers.monomer_generator`](monomer-generator.md) | SMILES → moltemplate `.lt` monomer templates |
| [`AutoPoly.mc`](mc.md) | Monte Carlo placement: collision detection, SAW chain growth |
| [`AutoPoly.core.exceptions`](exceptions.md) | Exception hierarchy |

!!! note "RDKit-dependent modules"
    `Polymer`, `Molecule`, `generate`, `MonomerGenerator`, and `mc` require RDKit. `System` and `BeadSpringPolymer` work without it.
