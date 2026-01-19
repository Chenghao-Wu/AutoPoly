# AutoPoly: Automated Polymer Generation for Molecular Simulation

![Python](https://img.shields.io/badge/python-3.7+-blue.svg)
![LAMMPS](https://img.shields.io/badge/LAMMPS-compatible-orange.svg)

AutoPoly is a Python package for generating polymer structures and preparing them for molecular simulations.

## ⚠️ BREAKING CHANGE NOTICE (v1.0)

**Version 1.0 is a major breaking release** with complete API modernization:

- **Pythonic naming**: `chain_num` instead of `ChainNum`, `sequence` instead of `Sequence`
- **Explicit sequences**: Specify exact monomer at each position (no cycling mode)
- **DOP derived from sequence**: No separate `DOP` parameter needed
- **Block copolymers**: Now supported with explicit monomer placement

**See [Migration Guide](#migration-guide-to-v10) below for details.**

## Quick Start

```bash
# Install AutoPoly
git clone <repository-url>
cd AutoPoly
pip install -e .

# For development
pip install -e ".[dev]"
```

**Note:** The `psmiles` package (for pSMILES canonicalization) is automatically installed from GitHub during installation.

```python
from AutoPoly import System, Polymer, Polymerization

# Create a linear polyethylene system
system = System(out="my_simulation")

# Define polymer using explicit pSMILES sequence
# Use a helper function for uniform polymers
def create_uniform_sequence(smiles: str, length: int) -> list:
    """Create a uniform sequence of given length."""
    return [smiles] * length

polymer = Polymer(
    chain_num=10,
    sequence=create_uniform_sequence("[*]CC[*]", 100),
    topology="linear",
    tacticity="atactic"
)

# Generate LAMMPS files
polymerization = Polymerization(
    name="polyethylene",
    system=system,
    model=[polymer]
)
```

## Key Features

- **Explicit Monomer Sequences** - Create block copolymers with precise monomer placement
- **Multiple Topologies** - Linear and ring polymer structures
- **Tacticity Control** - Atactic, isotactic, and syndiotactic configurations
- **Force Fields** - OPLS-AA, LOPLS, and GAFF support
- **Bead-Spring Models** - Coarse-grained simulations
- **LAMMPS Integration** - Complete input file generation
- **Automatic pSMILES Canonicalization** - Ensures consistent polymer SMILES representation
- **Resource Protection** - Built-in validation prevents resource exhaustion

## Basic Usage

### Linear Polymer (Uniform)

```python
from AutoPoly import System, Polymer, Polymerization

system = System(out="output_dir")

# Helper function for uniform polymers
def create_uniform_sequence(smiles: str, length: int) -> list:
    return [smiles] * length

# Use pSMILES with [*] connection points
polymer = Polymer(
    chain_num=10,
    sequence=create_uniform_sequence("[*]CC[*]", 50),
    topology="linear",
    tacticity="isotactic"
)
polymerization = Polymerization(name="polyethylene", system=system, model=[polymer])
```

### Linear Polymer (Block Copolymer)

**NEW in v1.0**: Create block copolymers with explicit monomer placement:

```python
# ABA triblock copolymer: 2 PE, 3 PS, 2 PE
sequence = [
    "[*]CC[*]",      # Position 0: Ethylene (Block A)
    "[*]CC[*]",      # Position 1: Ethylene (Block A)
    "[*]C=C[*]",     # Position 2: Styrene (Block B)
    "[*]C=C[*]",     # Position 3: Styrene (Block B)
    "[*]C=C[*]",     # Position 4: Styrene (Block B)
    "[*]CC[*]",      # Position 5: Ethylene (Block A)
    "[*]CC[*]"       # Position 6: Ethylene (Block A)
]

polymer = Polymer(
    chain_num=5,
    sequence=sequence,  # DOP is automatically 7
    topology="linear",
    tacticity="atactic"
)
```

### Ring Polymer

```python
# Ring polymer with explicit sequence
polymer = Polymer(
    chain_num=5,
    sequence=["[*]CC[*]"] * 30,  # 30 monomer units
    topology="ring",
    tacticity="atactic"
)
polymerization = Polymerization(name="ring_polyethylene", system=system, model=[polymer])
```

### Bead-Spring Model

```python
from AutoPoly import BeadSpringPolymer

bead_polymer = BeadSpringPolymer(
    name="coarse_grained",
    system=system,
    n_chains=5,
    n_beads=20,
    topology="linear"
)
bead_polymer.generate_data_file()
```

## pSMILES Canonicalization

AutoPoly automatically canonicalizes all pSMILES input using the [psmiles package](https://github.com/kuennethgroup/psmiles), ensuring consistent representation of polymer structures.

**Benefits:**
- Equivalent pSMILES always produce identical results
- Automatic validation catches invalid input
- Follows polymer SMILES best practices

**Example:**
```python
# All of these are equivalent after canonicalization:
generator.generate_variants("[*]CC([*])c1ccccc1")
generator.generate_variants("C(c1ccccc1)(C[*])[*]")
generator.generate_variants("*CC(C*)c1ccccc1")  # all produce same result
```

## API Reference

### System
Manages file paths and output directories.

```python
system = System(out="simulation_name")
output_path = system.get_folder_path()
```

### Polymer
Defines polymer structure and properties with explicit monomer sequences.

**Parameters:**
- `chain_num` (int): Number of polymer chains
- `sequence` (list): Explicit monomer sequence as pSMILES strings (e.g., `["[*]CC[*]", "[*]C=C[*]"]`)
- `topology` (str): "linear" or "ring"
- `tacticity` (str): "atactic", "isotactic", or "syndiotactic"

**Important Notes:**
- DOP is automatically derived from `len(sequence)` - no separate DOP parameter
- All monomers must use **pSMILES notation** with `[*]` wildcards for connection points
- Sequence is explicit - each position corresponds to one monomer unit
- For uniform polymers, use the helper function: `[smiles] * length`

**Common Monomer pSMILES:**
- Ethylene: `"[*]CC[*]"`
- Propylene: `"[*]CC(C)[*]"`
- Styrene: `"[*]Cc1ccccc1[*]"`
- Methyl methacrylate: `"[*]CC([*])(C)C(=O)OC"`

For a complete SMILES reference guide with more monomers, see [docs/SMILES_GUIDE.md](docs/SMILES_GUIDE.md)

**Examples:**

```python
# Uniform polymer
poly = Polymer(
    chain_num=10,
    sequence=["[*]CC[*]"] * 50,  # DOP = 50
    topology="linear",
    tacticity="isotactic"
)

# Block copolymer (NEW)
poly = Polymer(
    chain_num=5,
    sequence=[
        "[*]CC[*]",    # Position 0
        "[*]CC[*]",    # Position 1
        "[*]C=C[*]",   # Position 2
        "[*]CC[*]"     # Position 3
    ],  # DOP = 4
    topology="linear",
    tacticity="atactic"
)
```

### Polymerization
Generates polymer structures using Moltemplate.

**Parameters:**
- `name` (str): Project name
- `system` (System): System object
- `model` (list): List of Polymer objects
- `force_field` (str): "oplsaa", "lopls", or "gaff"

### BeadSpringPolymer
Creates coarse-grained bead-spring models.

**Parameters:**
- `n_chains` (int): Number of chains
- `n_beads` (int): Beads per chain
- `topology` (str): "linear" or "ring"
- `bond_length` (float): Equilibrium bond length
- `mass` (float): Bead mass

## Migration Guide (to v1.0)

### Breaking Changes Summary

**This is a major breaking release.** All existing code must be updated.

| Old API (Pre-1.0) | New API (v1.0) | Notes |
|------------------|----------------|-------|
| `ChainNum=10` | `chain_num=10` | Pythonic naming |
| `Sequence=["PE"]` | `sequence=["[*]CC[*]"] * 50` | Explicit SMILES |
| `DOP=50` | *removed* | Derived from `len(sequence)` |
| Cycling behavior | *removed* | No automatic sequence repetition |

### Migration Examples

#### Example 1: Uniform Polymer

**Before (Old API - No Longer Supported):**
```python
from AutoPoly import Polymer

poly = Polymer(
    ChainNum=10,      # OLD: CamelCase
    Sequence=["PE"],  # OLD: Monomer names
    DOP=50,           # OLD: Separate DOP parameter
    topology="linear",
    tacticity="atactic"
)
# OLD BEHAVIOR: Cycled ["PE"] 50 times
```

**After (New API):**
```python
from AutoPoly import Polymer

# Option 1: Direct list multiplication
poly = Polymer(
    chain_num=10,  # NEW: Pythonic snake_case
    sequence=["[*]CC[*]"] * 50,  # NEW: Explicit pSMILES
    topology="linear",
    tacticity="atactic"
)

# Option 2: Using helper function (recommended for clarity)
def create_uniform_sequence(smiles: str, length: int) -> list:
    return [smiles] * length

poly = Polymer(
    chain_num=10,
    sequence=create_uniform_sequence("[*]CC[*]", 50),
    topology="linear",
    tacticity="atactic"
)
```

#### Example 2: Block Copolymers (Now Possible!)

**Before (Old API - Not Possible):**
```python
# Could not create block copolymers
# Only uniform polymers supported
```

**After (New API):**
```python
# ABA triblock copolymer: 2 PE, 3 PS, 2 PE
sequence = [
    "[*]CC[*]",      # Position 0: Ethylene
    "[*]CC[*]",      # Position 1: Ethylene
    "[*]C=C[*]",     # Position 2: Styrene
    "[*]C=C[*]",     # Position 3: Styrene
    "[*]C=C[*]",     # Position 4: Styrene
    "[*]CC[*]",      # Position 5: Ethylene
    "[*]CC[*]"       # Position 6: Ethylene
]

poly = Polymer(
    chain_num=5,
    sequence=sequence,  # DOP is automatically 7
    topology="linear",
    tacticity="atactic"
)
```

#### Example 3: Multiple Unique Monomers

**Before (Old API - Cycling Mode):**
```python
# Would cycle through monomers
poly = Polymer(
    ChainNum=5,
    Sequence=["PE", "PS"],  # Would cycle: PE, PS, PE, PS, ...
    DOP=8,
    topology="linear"
)
# Result: PE, PS, PE, PS, PE, PS, PE, PS (4 cycles)
```

**After (New API - Explicit Sequence):**
```python
# Specify exact sequence (no cycling)
sequence = [
    "[*]CC[*]",      # PE
    "[*]C=C[*]",     # PS
    "[*]CC[*]",      # PE
    "[*]C=C[*]",     # PS
    "[*]CC[*]",      # PE
    "[*]C=C[*]",     # PS
    "[*]CC[*]",      # PE
    "[*]C=C[*]"      # PS
]

poly = Polymer(
    chain_num=5,
    sequence=sequence,  # Explicit sequence, DOP = 8
    topology="linear",
    tacticity="isotactic"
)
```

### Key Behavioral Changes

1. **No More Cycling Mode**
   - Old: `Sequence=["A", "B"]` with `DOP=6` → `A, B, A, B, A, B`
   - New: Must specify full sequence explicitly

2. **Explicit SMILES Required**
   - Old: Monomer names like `"PE"`, `"PS"`
   - New: pSMILES strings like `"[*]CC[*]"`, `"[*]C=C[*]"`

3. **DOP is Derived**
   - Old: `DOP` parameter controls sequence repetition
   - New: `DOP = len(sequence)` automatically

4. **Pythonic Naming**
   - Old: `ChainNum`, `Sequence`, `DOP`
   - New: `chain_num`, `sequence`, `dop` (read-only)

### Validation and Resource Limits

**New in v1.0**: Built-in validation protects against resource exhaustion:

```python
from AutoPoly.exceptions import ValidationError

# Sequence too long
try:
    poly = Polymer(
        chain_num=1,
        sequence=["[*]CC[*]"] * 15000  # Exceeds MAX_SEQUENCE_LENGTH
    )
except ValidationError as e:
    print(f"Protected: {e}")
    # "Sequence length (15000) exceeds maximum 10000"

# Too many unique monomers
try:
    poly = Polymer(
        chain_num=1,
        sequence=[f"[*]C{i}C[*]" for i in range(150)]  # Exceeds MAX_UNIQUE_MONOMERS
    )
except ValidationError as e:
    print(f"Protected: {e}")
    # "Number of unique monomers (150) exceeds maximum 100"
```

**Resource Limits:**
- `MAX_DOP = 10000` - Maximum degree of polymerization
- `MAX_SEQUENCE_LENGTH = 10000` - Maximum sequence length
- `MAX_UNIQUE_MONOMERS = 100` - Maximum unique monomer types

## Advanced Documentation

- **Full Examples** → [example/README.md](example/README.md)
- **API Documentation** → [docs/API.md](docs/API.md)
- **Migration Guide** → [MIGRATION.md](MIGRATION.md)
- **Troubleshooting** → See common issues below

## Troubleshooting

**API-related errors:**
- `TypeError: __init__() got an unexpected keyword argument 'ChainNum'`: Use `chain_num` instead (Pythonic naming)
- `TypeError: __init__() got an unexpected keyword argument 'DOP'`: DOP is derived from sequence length
- `TypeError: __init__() got an unexpected keyword argument 'Sequence'`: Use `sequence` instead

**Sequence-related errors:**
- `ValidationError: sequence cannot be empty`: Provide at least one monomer in sequence
- `ValidationError: Sequence length exceeds maximum`: Reduce sequence length below 10000
- `ValidationError: Invalid SMILES`: Ensure all monomers use valid pSMILES notation with `[*]` wildcards

**psmiles package errors:**
- The psmiles package is required and should be installed automatically
- If installation fails: `pip install git+https://github.com/kuennethgroup/psmiles.git`
- Ensure you have Git installed if using direct Git dependency

**Invalid pSMILES errors:**
- All monomers must use pSMILES notation with exactly 2 wildcard atoms ([*])
- Example valid: `"[*]CC[*]"`
- Example invalid: `"CC"` (missing wildcards)
- Check the [SMILES_GUIDE.md](docs/SMILES_GUIDE.md) for more details

**Moltemplate errors:**
- Ensure Moltemplate is installed and in PATH
- Check monomer .lt file syntax

**File permission errors:**
- Verify write permissions for output directory

## Output Structure

```
project_name/
├── moltemplate/     # Intermediate files
├── system.data      # LAMMPS data file
├── system.in        # LAMMPS input script
└── system.in.settings # Force field parameters
```

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