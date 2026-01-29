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
- **Automatic Setup** - Generates complete LAMMPS input files

Note: The current atom typing systems (for all force fields) relies on SMARTS pattern which is built mannuly. The correct SMARTS pattern can be built via a data-driven method as BESMARTS. We are currently exploring this for better atom typing system.

## Quick Start

### Installation

```bash
git clone <repository-url>
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

### Other Examples

- **Bead-spring models** - Coarse-grained simulations
- **Multiple polymers** - Blends and mixtures
- **Custom tacticity** - Isotactic, syndiotactic, atactic

[See all examples →](examples/)

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

**Note:** GAFF/GAFF2 require explicit charge calculation (not automatic).

[Detailed comparison →](docs/FORCE_FIELDS.md)

## Troubleshooting

### Common Issues

**API errors:**
- `TypeError: 'ChainNum'` → Use `chain_num` (v1.0+ uses snake_case)
- `TypeError: 'DOP'` → Removed, DOP = `len(sequence)` automatically
- [Migration guide →](MIGRATION.md)

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

```
project_name/
├── moltemplate/           # Intermediate files
├── system.data            # LAMMPS data file (topology & coordinates)
├── system.in              # LAMMPS input script
└── system.in.settings     # Force field parameters
```

Run with LAMMPS:
```bash
lmp -in system.in
```

## More Information

**Documentation:**
- 📖 [Complete API Reference](docs/API.md) - All classes and methods
- 🔄 [v1.0 Migration Guide](MIGRATION.md) - Upgrading from v0.x
- 🧬 [Complement SMILES Guide](docs/COMPLEMENT_SMILES.md) - Deep dive on SMILES format
- ⚙️ [Force Field Guide](docs/FORCE_FIELDS.md) - Detailed comparison of all 6 force fields
- 🐛 [Troubleshooting](docs/TROUBLESHOOTING.md) - Solutions to common issues
- 📝 [Examples Directory](examples/) - Working examples

**Quick links:**
- [Installation details](docs/TROUBLESHOOTING.md#installation-issues)
- [Ring vs linear polymers](docs/COMPLEMENT_SMILES.md#ring-vs-linear-polymers)
- [Common monomer SMILES](docs/COMPLEMENT_SMILES.md#common-monomers-reference)

## Citation

If you use AutoPoly in your research, please cite:

```bibtex
@software{autopoly2024,
  title={AutoPoly: Automated Polymer Generation for Molecular Simulation},
  author={Wu, Zhenghao},
  year={2024},
  url={https://github.com/your-repo/autopoly}
}
```

## License

MIT License - see [license.md](license.md) for details.

---

**Note for v0.x users:** Version 1.0 has breaking API changes (snake_case, explicit sequences). See [MIGRATION.md](MIGRATION.md) for upgrade guide.
