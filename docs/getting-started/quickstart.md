# Quickstart

This page walks through the complete AutoPoly workflow on the simplest possible system: a box of polyethylene (PE) chains with the OPLS-AA force field.

## The 3-step pattern

Every AutoPoly script follows the same shape:

```
System  →  Polymer / Molecule  →  Polymerization  →  LAMMPS files
```

1. **System** — decides *where* output goes
2. **Polymer / Molecule** — defines *what* to build
3. **Polymerization** — generates the files with a chosen force field

## Step by step

### 1. Create the system

```python
from AutoPoly import System

system = System(out="my_polymer")
```

`out` names the output directory that will hold everything AutoPoly generates.

### 2. Define the polymer

```python
from AutoPoly import Polymer

polymer = Polymer(
    chain_num=10,                                         # 10 chains in the box
    sequence=["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"],   # PE, DOP = 50
    topology="linear",                                    # or "ring"
    tacticity="atactic"                                   # or "isotactic" / "syndiotactic"
)
```

The `sequence` is a list of **complement SMILES** — one entry per monomer, in chain order. Wildcard atoms (`[*]`) mark the connection points:

- `"CC[*]"` — first monomer (1 wildcard, right side)
- `"[*]CC[*]"` — middle monomers (2 wildcards)
- `"[*]CC"` — last monomer (1 wildcard, left side)

The degree of polymerization is implicit: `DOP = len(sequence)` (here 1 + 48 + 1 = 50). See the [Complement SMILES guide](../guides/complement-smiles.md) for the full format.

### 3. Generate the LAMMPS files

```python
from AutoPoly import Polymerization

Polymerization(
    name="polyethylene",
    system=system,
    model=[polymer],
    force_field="oplsaa"
)
```

AutoPoly now generates monomer `.lt` templates, performs atom typing for OPLS-AA, grows chains with the Monte Carlo self-avoiding-walk placer, runs moltemplate, and writes the LAMMPS inputs.

## What you get

```
my_polymer/polyethylene/
├── moltemplate/           # Intermediate files (.lt inputs, moltemplate output)
├── system.data            # LAMMPS data file (topology & coordinates)
├── system.in.init         # Units, atom/bond/angle styles
├── system.in.settings     # Force field parameters
└── system.in.charges      # Atomic charges
```

Run it with LAMMPS by including the pieces from your own input script:

```bash
lmp -in your_run.in   # with: include system.in.init / system.in.settings
```

See [Output Files](../guides/output-files.md) for a description of each file.

## Next steps

- Build a **block copolymer** by mixing monomers in the sequence — [tutorial](../tutorials/block-copolymer.md)
- Add a **solvent** with the `Molecule` class — [Polymers vs Molecules](../guides/polymers-vs-molecules.md)
- Try a different **force field** — [selection guide](../guides/force-fields.md)
- Control **chain placement** — [MC Placement & Chain Growth](../guides/mc-placement.md)
