# The 3-Step Workflow

AutoPoly turns SMILES strings into LAMMPS input files through a small, fixed set of objects. Every script — from a single solvent box to a multi-block copolymer melt — follows the same three steps.

```
System  →  Polymer / Molecule  →  generate  →  LAMMPS files
 (where)        (what)             (how)         (output)
```

## Step 1 — System: where the output goes

```python
from AutoPoly import System

system = System(out="my_simulation")
```

`System` is a lightweight container for the output directory path. Everything AutoPoly writes lands under this folder. Create one per project.

## Step 2 — Polymer / Molecule: what to build

Define each chemical species in the box. There are two classes, and you can mix them freely:

- **`Polymer`** — polymer chains described by an explicit **complement SMILES** sequence (one entry per monomer). Controls chain count, topology (linear/ring), and tacticity.
- **`Molecule`** — small molecules (solvents, additives) described by a regular SMILES string and a count.

```python
from AutoPoly import Polymer, Molecule

polymer = Polymer(
    chain_num=10,
    sequence=["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"],
    topology="linear",
    tacticity="atactic",
)

water = Molecule(Count=100, Smiles="O", Name="water")
```

The degree of polymerization is never specified directly — it is always `len(sequence)`.

## Step 3 — generate: how to build it

```python
from AutoPoly import generate

generate(
    system,
    "polyethylene",
    [polymer, water],
    force_field="oplsaa",
)
```

`generate` takes the list of models and runs the full three-stage pipeline (Geometry → Typing → Packing):

1. **Monomer generation** — each unique complement SMILES becomes a moltemplate `.lt` template (six variants per monomer: first/middle/last × two chiralities), with 3D geometry from RDKit.
2. **Atom typing & parameters** — the force field manager assigns atom types from SMARTS patterns and collects bond/angle/dihedral/LJ parameters for the six supported force fields.
3. **Chain placement** — chains are grown monomer-by-monomer with the Monte Carlo self-avoiding-walk algorithm (or placed on a deterministic grid) inside a sized simulation box.
4. **Moltemplate** — the assembled `system.lt` is compiled into a LAMMPS data file.
5. **Output organization** — `system.data`, `system.in.init`, `system.in.settings`, and `system.in.charges` are written to `<System out>/<name>/`.

!!! tip "Need more than one force field?"
    `generate` is the one-shot path. For stage-level control — e.g. building
    the geometry once and typing it under several force fields — use the
    `GeometryBuilder` / `UnitTyper` / `BoxPacker` classes from
    `AutoPoly.pipeline` directly.

See [Output Files](output-files.md) for what each file contains.

## The coarse-grained exception

`BeadSpringPolymer` skips this pipeline entirely: it needs no SMILES, no atom typing, and no moltemplate. It writes a LAMMPS data file directly. See the [Bead-Spring guide](bead-spring.md).

## Where to go next

- [Complement SMILES](complement-smiles.md) — the monomer sequence format
- [Polymers vs Molecules](polymers-vs-molecules.md) — choosing the right model class
- [Force Fields](force-fields.md) — the six supported parameter sets
- [MC Placement & Chain Growth](mc-placement.md) — how chains fill the box
