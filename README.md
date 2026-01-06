# AutoPoly: Automated Polymer Generation for Molecular Simulation

![Python](https://img.shields.io/badge/python-3.7+-blue.svg)
![LAMMPS](https://img.shields.io/badge/LAMMPS-compatible-orange.svg)

AutoPoly is a Python package for generating polymer structures and preparing them for molecular simulations.

## Quick Start

```bash
# Install AutoPoly
git clone <repository-url>
cd AutoPoly
pip install -e .

# For development
pip install -e ".[dev]"
```

```python
from AutoPoly import System, Polymer, Polymerization

# Create a linear polyethylene system
system = System(out="my_simulation")

# Define polymer using SMILES notation
# "C=C" is the SMILES string for ethylene (polyethylene monomer)
polymer = Polymer(
    ChainNum=10,
    Sequence=["C=C"],
    DOP=100,
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

- **Multiple Topologies** - Linear and ring polymer structures
- **Tacticity Control** - Atactic, isotactic, and syndiotactic configurations
- **Force Fields** - OPLS-AA, LOPLS, and GAFF support
- **Bead-Spring Models** - Coarse-grained simulations
- **LAMMPS Integration** - Complete input file generation

## Basic Usage

### Linear Polymer

```python
from AutoPoly import System, Polymer, Polymerization

system = System(out="output_dir")
# Use SMILES strings for monomers (C=C = ethylene)
polymer = Polymer(
    ChainNum=10,
    Sequence=["C=C"],
    DOP=50,
    topology="linear",
    tacticity="isotactic"
)
polymerization = Polymerization(name="polyethylene", system=system, model=[polymer])
```

### Ring Polymer

```python
# For ring polymers, use SMILES with [*] wildcards (pSMILES notation)
polymer = Polymer(
    ChainNum=5,
    Sequence=["[*]C=C[*]"],  # pSMILES for ring closure
    DOP=30,
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

## API Reference

### System
Manages file paths and output directories.

```python
system = System(out="simulation_name")
output_path = system.get_folder_path()
```

### Polymer
Defines polymer structure and properties.

**Parameters:**
- `ChainNum` (int): Number of chains
- `Sequence` (list): Monomer sequence as SMILES or pSMILES strings
- `DOP` (int): Degree of polymerization
- `topology` (str): "linear" or "ring"
- `tacticity` (str): "atactic", "isotactic", or "syndiotactic"

**Note**: All monomers must be specified using SMILES notation. Common examples:
- Ethylene: `"C=C"` or `"[*]C=C[*]"` (for linear/ring polymers)
- Propylene: `"C=C(C)"` or `"[*]C=C(C)[*]"`
- Styrene: `"C=C(C1=CC=CC=C1)"` or `"[*]C=C(C1=CC=CC=C1)[*]"`
- Methyl methacrylate: `"C=C(C)C(=O)OC"` or `"[*]C=C(C)C(=O)OC[*]"`

For a complete SMILES reference guide with more monomers, see [docs/SMILES_GUIDE.md](docs/SMILES_GUIDE.md)

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

## Advanced Documentation

- **Full Examples** → [example/README.md](example/README.md)
- **API Documentation** → [docs/API.md](docs/API.md)
- **Migration Guide** → [MIGRATION.md](MIGRATION.md)
- **Troubleshooting** → See common issues below

## Troubleshooting

**Monomer not found:**
- Check `extern/Monomer_bank/` for available monomers
- Verify monomer names match your sequence

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