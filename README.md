# AutoPoly: Automated Polymer Generation and Simulation Package

AutoPoly is a comprehensive Python package for generating polymer structures and preparing them for molecular dynamics simulations using LAMMPS and Moltemplate.

## Overview

AutoPoly provides tools for:
- **Polymer Structure Generation**: Create atomistic polymer models with various topologies and tacticity
- **Force Field Integration**: Seamless integration with OPLS-AA and LOPLS force fields
- **LAMMPS Preparation**: Generate complete LAMMPS input files and data structures
- **Bead-Spring Models**: Simplified coarse-grained polymer models for quick simulations
- **File Management**: Automated organization of simulation files and outputs

## Key Features

### 🧬 Polymer Topologies
- **Linear Polymers**: Standard linear chain structures
- **Ring Polymers**: Circular polymer topologies
- **Custom Sequences**: User-defined monomer sequences

### 🎯 Tacticity Control
- **Atactic**: Random stereochemistry
- **Isotactic**: All monomers with same stereochemistry
- **Syndiotactic**: Alternating stereochemistry

### 🔬 Force Fields
- **OPLS-AA**: All-atom force field for accurate simulations
- **LOPLS**: Lipid-optimized force field variant
- **Custom Parameters**: Support for modified force field parameters

### 📁 File Management
- **Monomer Bank**: Centralized monomer template library
- **Output Organization**: Structured file organization
- **Error Handling**: Comprehensive validation and error reporting

## Installation

### Prerequisites
- Python 3.7+
- LAMMPS (for simulations)
- Moltemplate (for structure generation)

### Install AutoPoly
```bash
# Clone the repository
git clone <repository-url>
cd AutoPoly

# Install the package
pip install -e .
```

## Quick Start

### 1. Basic Polymer Generation

```python
from AutoPoly import System, Polymer, Polymerization

# Create system
system = System(out="my_polymer")

# Define polymer
polymer = Polymer(
    ChainNum=10,
    Sequence=["ethylene"],
    DOP=100,
    topology="linear",
    tacticity="atactic"
)

# Generate structure
polymerization = Polymerization(
    name="ethylene_polymer",
    system=system,
    model=[polymer]
)
```

### 2. Bead-Spring Model

```python
from AutoPoly import BeadSpringPolymer, System

# Create system
system = System(out="bead_spring")

# Generate bead-spring polymer
bead_polymer = BeadSpringPolymer(
    name="test_polymer",
    system=system,
    n_chains=5,
    n_beads=20,
    topology="linear"
)

# Generate LAMMPS files
bead_polymer.generate_data_file()
```

### 3. Ring Polymer

```python
# Create ring polymer
ring_polymer = Polymer(
    ChainNum=5,
    Sequence=["styrene"],
    topology="ring",
    tacticity="atactic"
)

# Generate with polymerization
polymerization = Polymerization(
    name="ring_styrene",
    system=system,
    model=[ring_polymer]
)
```

## Package Structure

```
AutoPoly/
├── AutoPoly/
│   ├── __init__.py          # Package initialization
│   ├── system.py            # System management
│   ├── polymer.py           # Polymer definition
│   ├── polymerization.py    # Core polymerization logic
│   ├── bead_spring.py       # Bead-spring models
│   ├── conf.py              # Configuration settings
│   ├── logger.py            # Logging utilities
│   └── extern/              # External dependencies
│       ├── Monomer_bank/    # Monomer templates
│       ├── moltemplate/     # Moltemplate files
│       └── rdlt.py          # LAMMPS data utilities
├── examples/                # Usage examples
├── tests/                   # Test suite
└── docs/                    # Documentation
```

## Core Classes

### System
Manages file paths and directory operations for polymer simulations.

```python
system = System(out="output_directory")
```

### Polymer
Defines polymer properties including topology, tacticity, and monomer sequences.

```python
polymer = Polymer(
    ChainNum=10,           # Number of chains
    Sequence=["monomer"],  # Monomer sequence
    DOP=100,              # Degree of polymerization
    topology="linear",     # "linear" or "ring"
    tacticity="atactic"    # "atactic", "isotactic", "syndiotactic"
)
```

### Polymerization
Core class for generating polymer structures using Moltemplate.

```python
polymerization = Polymerization(
    name="project_name",
    system=system,
    model=[polymer],
    is_lopls=False  # Use LOPLS force field
)
```

### BeadSpringPolymer
Simplified bead-spring model generator for coarse-grained simulations.

```python
bead_polymer = BeadSpringPolymer(
    name="model_name",
    system=system,
    n_chains=5,
    n_beads=20,
    topology="linear"
)
```

## Configuration

### Logging
Configure logging levels in `conf.py`:

```python
LOG = {
    'ROOT_LEVEL': logging.INFO,
    'CONSOLE_LEVEL': logging.INFO,
    'FILE_LEVEL': logging.INFO,
    'TO_FILE': False
}
```

### Output Paths
Set output directory in `conf.py`:

```python
OUT_PATH = os.path.join(os.path.expanduser("~"))
```

## Examples

### Linear Polyethylene
```python
from AutoPoly import System, Polymer, Polymerization

system = System(out="polyethylene")
polymer = Polymer(ChainNum=10, Sequence=["ethylene"], DOP=50)
polymerization = Polymerization(name="PE", system=system, model=[polymer])
```

### Ring Polystyrene
```python
polymer = Polymer(
    ChainNum=5,
    Sequence=["styrene"],
    topology="ring",
    tacticity="atactic"
)
polymerization = Polymerization(name="ring_PS", system=system, model=[polymer])
```

### Copolymer
```python
polymer = Polymer(
    ChainNum=10,
    Sequence=["ethylene", "propylene"],
    DOP=100,
    tacticity="syndiotactic"
)
```

## Output Files

AutoPoly generates the following files:

### LAMMPS Files
- `system.data`: Atom coordinates and connectivity
- `system.in.settings`: Force field parameters
- `system.in.charges`: Atomic charges
- `system.in`: LAMMPS input script
- `system.in.init`: Initialization script

### Organization
```
project_name/
├── moltemplate/     # Intermediate files
├── input/          # Input templates
├── output/         # LAMMPS files
└── logs/           # Log files
```

## Troubleshooting

### Common Issues

1. **Monomer Not Found**
   - Check monomer bank path in `extern/Monomer_bank/`
   - Verify monomer file names match sequence

2. **Moltemplate Errors**
   - Ensure Moltemplate is properly installed
   - Check monomer .lt file syntax
   - Verify force field parameter files

3. **File Permission Errors**
   - Check write permissions for output directory
   - Ensure sufficient disk space

### Debug Mode
Enable detailed logging:

```python
from AutoPoly.conf import LOG
LOG['ROOT_LEVEL'] = logging.DEBUG
LOG['TO_FILE'] = True
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests
5. Submit a pull request

## License

This project is licensed under the BSD License - see the `license.md` file for details.

## Citation

If you use AutoPoly in your research, please cite:

```bibtex
@software{autopoly2024,
  title={AutoPoly: Automated Polymer Generation and Simulation Package},
  author={Wu, Zhenghao},
  year={2024},
  url={https://github.com/your-repo/autopoly}
}
```

## Support

For questions and support:
- Check the documentation
- Review example files
- Open an issue on GitHub
- Contact the maintainers

## Acknowledgments

- Moltemplate developers for the structure generation framework
- LAMMPS developers for the molecular dynamics engine
- OPLS-AA force field developers
