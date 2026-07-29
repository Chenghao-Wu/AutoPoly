# AutoPoly: Automated Polymer Generation for Molecular Simulation

![Python](https://img.shields.io/badge/python-3.7+-blue.svg)
![LAMMPS](https://img.shields.io/badge/LAMMPS-compatible-orange.svg)

AutoPoly generates polymer structures and prepares them for molecular dynamics simulations with LAMMPS. Build everything from simple homopolymers to complex block copolymers with explicit sequence control.

**Key Features:**
- **6 Force Fields** - OPLS-AA, LOPLS, GAFF, GAFF2, DREIDING, COMPASS
- **Block Copolymers** - Explicit sequence control for any block arrangement
- **Complement SMILES** - Unique format for precise positional control
- **Small Molecules** - Built-in support for solvents and additives
- **Ring & Linear** - Both topologies supported
- **SAW Placement** - Monte Carlo self-avoiding walk for realistic initial configurations
- **Automatic Setup** - Generates complete LAMMPS input files

Note: The current atom typing system (for all force fields) relies on manually built SMARTS patterns. A data-driven approach such as BESMARTS could generate more robust SMARTS patterns, and we are exploring this for a better atom typing system.

## Quick Start

### Installation

```bash
git clone https://github.com/WuGroup-XJTLU/AutoPoly.git
cd AutoPoly
pip install -e .
```

### Your First Polymer (3 Steps)

```python
from AutoPoly import System, Polymer, Polymerization

# Step 1: Create system
system = System(out="my_polymer")

# Step 2: Define polymer
polymer = Polymer(
    chain_num=10,                    # 10 chains
    sequence=["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"],  # PE, DOP=50
    topology="linear",
    tacticity="atactic"
)

# Step 3: Generate LAMMPS files
Polymerization(
    name="polyethylene",
    system=system,
    model=[polymer],
    force_field="oplsaa"
)
```

**Output:** Ready-to-run LAMMPS files in `my_polymer/` directory!

## Core Concepts

### The 3-Step Workflow

```
System → Polymer/Molecule → Polymerization → LAMMPS Files
```

1. **System** - Defines output directory
2. **Polymer/Molecule** - Defines what to build
3. **Polymerization** - Generates files with chosen force field

### Complement SMILES (Unique to AutoPoly)

AutoPoly uses **complement SMILES** to control each monomer's position:

- **First position:** `"CC[*]"` (1 wildcard, right)
- **Middle positions:** `"[*]CC[*]"` (2 wildcards, both sides)
- **Last position:** `"[*]CC"` (1 wildcard, left)

This enables precise block copolymer design:

```python
# ABA triblock: PE(2)-PS(3)-PE(2)
sequence = [
    "CC[*]",                    # PE first
    "[*]CC[*]",                 # PE middle
    "[*]CC([*])c1ccccc1",       # PS middle
    "[*]CC([*])c1ccccc1",       # PS middle
    "[*]CC([*])c1ccccc1",       # PS middle
    "[*]CC[*]",                 # PE middle
    "[*]CC"                     # PE last
]
```

**Why it matters:** Standard SMILES can't distinguish first/middle/last positions. Complement SMILES gives you control over every monomer.

[Deep dive: Complement SMILES Guide →](docs/COMPLEMENT_SMILES.md)

### Polymer vs Molecule

| Feature | Polymer | Molecule |
|---------|---------|----------|
| **Use for** | Polymers, chains | Solvents, small molecules |
| **SMILES** | Complement SMILES with `[*]` | Regular SMILES, no `[*]` |
| **Parameters** | `chain_num`, `sequence`, `topology` | `Count`, `Smiles`, `Name` |
| **Example** | `Polymer(chain_num=10, sequence=["CC[*]"]+["[*]CC[*]"]*48+["[*]CC"])` | `Molecule(Count=100, Smiles="O", Name="water")` |

## Examples

### Block Copolymers (Complement SMILES)

```python
# PE-PS-PE triblock (10-20-10)
sequence = (
    ["CC[*]"] + ["[*]CC[*]"] * 9 +                  # PE block (10)
    ["[*]CC([*])c1ccccc1"] * 20 +                   # PS block (20)
    ["[*]CC[*]"] * 9 + ["[*]CC"]                    # PE block (10)
)

polymer = Polymer(
    chain_num=5,
    sequence=sequence,  # DOP = 40
    topology="linear",
    tacticity="atactic"
)
```

[Full example: examples/example_block_copolymer.py →](examples/example_block_copolymer.py)

### Small Molecules (Solvent)

```python
from AutoPoly import Molecule

# Water molecules
water = Molecule(
    Count=100,
    Smiles="O",     # Regular SMILES (no wildcards)
    Name="water"
)

# Ethanol molecules
ethanol = Molecule(
    Count=20,
    Smiles="CCO",
    Name="ethanol"
)

# Generate system
Polymerization(
    name="solvent_mixture",
    system=system,
    model=[water, ethanol],
    force_field="gaff"
)
```

[Full example: examples/example_molecules.py →](examples/example_molecules.py)

### Mixed System (Polymer + Solvent)

```python
# Polymer in water
polymer = Polymer(
    chain_num=5,
    sequence=["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"]
)
water = Molecule(Count=100, Smiles="O", Name="water")

Polymerization(
    name="polymer_solution",
    system=system,
    model=[polymer, water],
    force_field="gaff"
)
```

### Ring Polymer

```python
# Ring: use only middle variants (all 2 wildcards)
polymer = Polymer(
    chain_num=5,
    sequence=["[*]CC[*]"] * 30,  # All middle
    topology="ring"              # Specify ring
)
```

### Monte Carlo Placement with Self-Avoiding Walk

AutoPoly uses a Monte Carlo (MC) self-avoiding walk (SAW) algorithm to generate realistic initial polymer configurations. Instead of placing chains on a grid, the SAW method grows each chain monomer-by-monomer with collision detection, producing coiled conformations that better approximate equilibrium structures.

```python
Polymerization(
    name="polymer_mc",
    system=system,
    model=[polymer],
    force_field="oplsaa",
    placement_method="mc_random",       # Monte Carlo placement (default)
    use_mc_chain_growth=True,           # SAW chain growth (default)
    mc_max_attempts=10000,              # Max placement attempts
    mc_monomer_density=0.085,           # Target density (monomers/Å³)
    mc_bond_angle_min=50.0,             # Min deflection angle (degrees)
    mc_bond_angle_max=90.0              # Max deflection angle (degrees)
)
```

**Placement methods:**
- `"mc_random"` (default) — Monte Carlo with SAW chain growth
- `"grid"` — Deterministic grid placement

Box sizing uses SAW scaling (`N^0.6 × bond_length`) rather than fully-extended chain length, producing compact, realistic simulation boxes.

### More Examples

The [examples directory](examples/) contains 12 runnable scripts covering:

- **Beginner tutorial** — PMMA step by step ([example_pmma_linear.py](examples/example_pmma_linear.py))
- **Condensation polymers** — PLA with GAFF ([example_pla_condensation.py](examples/example_pla_condensation.py))
- **Polymer solutions** — PEO in explicit water ([example_peo_solution.py](examples/example_peo_solution.py))
- **Batch generation** — 10 commodity polymers ([example_commodity_polymers_10.py](examples/example_commodity_polymers_10.py))
- **Force field comparison** — all 6 force fields on PEO ([example_peo_all_forcefields.py](examples/example_peo_all_forcefields.py))
- **Placement methods** — grid vs MC random vs MC chain growth ([example_peo_mc_placement.py](examples/example_peo_mc_placement.py))
- **Bead-spring models** — coarse-grained homo/block/ring polymers ([example_bead_spring.py](examples/example_bead_spring.py))

[See the full list with descriptions →](examples/README.md)

## API Quick Reference

### System

```python
System(out="folder_name")
```

Creates output directory for simulation files.

### Polymer

```python
Polymer(
    chain_num=10,               # Number of chains
    sequence=["[*]CC[*]"] * 50, # Explicit sequence (DOP=50)
    topology="linear",          # or "ring"
    tacticity="atactic"         # or "isotactic", "syndiotactic"
)
```

**Key points:**
- DOP is automatic from `len(sequence)`
- Use complement SMILES: first=`"CC[*]"`, middle=`"[*]CC[*]"`, last=`"[*]CC"`
- [Complete API →](docs/API.md#polymer)

### Molecule

```python
Molecule(
    Count=100,       # Number of molecules
    Smiles="O",      # Regular SMILES (no wildcards)
    Name="water"     # Identifier
)
```

**Key point:** Use regular SMILES without `[*]` wildcards.

[Complete API →](docs/API.md#molecule)

### Polymerization

```python
Polymerization(
    name="project",
    system=system,
    model=[polymer1, polymer2, molecule1],  # Mix polymers and molecules
    force_field="oplsaa"  # See force fields below
)
```

**Supported force fields:**

| Force Field | Value | Best For |
|------------|-------|----------|
| OPLS-AA | `"oplsaa"` | General organic polymers |
| LOPLS | `"lopls"` | Liquid-phase, better densities |
| GAFF | `"gaff"` | Small molecules, drug-like |
| GAFF2 | `"gaff2"` | Updated GAFF |
| DREIDING | `"dreiding"` | Generic, metals, inorganics |
| COMPASS | `"compass"` | Commercial polymers |

[Force field selection guide →](docs/FORCE_FIELDS.md)

[Complete API →](docs/API.md)

## Force Field Selection

Quick guide:

- **Organic polymers** → OPLS-AA or LOPLS
- **Small molecules/solvents** → GAFF or GAFF2
- **Accurate densities** → LOPLS or COMPASS
- **Exploratory/generic** → DREIDING
- **Commercial polymers** → COMPASS

**Note:** Gasteiger charges are assigned automatically for GAFF/GAFF2. For production runs, replace them with AM1-BCC or RESP charges in `system.in.charges`.

[Detailed comparison →](docs/FORCE_FIELDS.md)

## Troubleshooting

### Common Issues

**API errors:**
- `TypeError: 'ChainNum'` → Use `chain_num` (v1.0 uses snake_case)
- `TypeError: 'DOP'` → Removed, DOP = `len(sequence)` automatically

**Sequence errors:**
- `ValidationError: sequence cannot be empty` → Provide at least one monomer
- `ValidationError: Sequence length exceeds maximum` → Max DOP is 10000

**SMILES errors:**
- `ValidationError: Invalid SMILES` → Check wildcard count (first=1, middle=2, last=1)
- Ring polymers → Use only middle variants: `["[*]CC[*]"] * 50`

**Moltemplate:**
- `Moltemplate not found` → Install: `pip install moltemplate`

[Full troubleshooting guide →](docs/TROUBLESHOOTING.md)

## Output Structure

Each run writes into `<System out>/<Polymerization name>/`:

```
my_polymer/polyethylene/
├── moltemplate/           # Intermediate files (.lt inputs, moltemplate output)
├── system.data            # LAMMPS data file (topology & coordinates)
├── system.in.init         # Units, atom/bond/angle styles
├── system.in.settings     # Force field parameters
└── system.in.charges      # Atomic charges
```

Run with LAMMPS by including the pieces from your own input script:
```bash
lmp -in your_run.in   # with: include system.in.init / system.in.settings
```

## More Information

**Documentation:**
- 📖 [Complete API Reference](docs/API.md) - All classes and methods
- 🧬 [Complement SMILES Guide](docs/COMPLEMENT_SMILES.md) - Deep dive on SMILES format
- ⚙️ [Force Field Guide](docs/FORCE_FIELDS.md) - Detailed comparison of all 6 force fields
- 🐛 [Troubleshooting](docs/TROUBLESHOOTING.md) - Solutions to common issues
- 📝 [Examples Directory](examples/) - 12 working examples

**Quick links:**
- [Installation details](docs/TROUBLESHOOTING.md#installation-issues)
- [Ring vs linear polymers](docs/COMPLEMENT_SMILES.md#ring-vs-linear-polymers)
- [Common monomer SMILES](docs/COMPLEMENT_SMILES.md#common-monomers-reference)

## Citation

If you use AutoPoly in your research, please cite:

```bibtex
@software{autopoly,
  title={AutoPoly: Automated Polymer Generation for Molecular Simulation},
  author={Wu, Zhenghao},
  url={https://github.com/WuGroup-XJTLU/AutoPoly}
}
```

## License

MIT License - see [license.md](license.md) for details.
