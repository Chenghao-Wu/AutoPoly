# AutoPoly.py: Automatic Generation of Data File for Polymers

AutoPoly is a Python package for automatically generating LAMMPS data files for polymer systems. It supports both atomistic (OPLS force field) and coarse-grained (bead-spring) polymer models.

## Why AutoPoly?

- Create polymers with any chain length
- Support for linear and ring polymer topologies
- Easy integration with RDKit for custom monomer creation
- Automated handling of tacticity
- Built-in support for common polymer types
- Coarse-grained bead-spring polymer models

## Installation

Install AutoPoly via pip:
```bash
git clone https://github.com/Chenghao-Wu/AutoPoly.git
cd AutoPoly
pip install .
```

## Supported Models

### Atomistic Polymers (OPLS)
1. PMMA (Tacticity supported)
2. PS (Tacticity supported)
3. PE (Tacticity supported)
4. PI (cis)
5. PP (Tacticity supported)
6. PVA (Tacticity supported)

### Coarse-Grained Models
- Bead-spring polymer chains
  - Linear topology
  - Ring topology
  - Configurable parameters (mass, bond length, LJ parameters)
  - 3D spatial distribution for non-overlapping configurations

## Examples

### 1. Linear Polymer (PMMA)
```python
import AutoPoly

# Define the system
system = AutoPoly.System(out="pmma_linear")

# Create a linear PMMA polymer with 10 chains, each with 50 monomers
linear_polymer = AutoPoly.Polymer(ChainNum=10, Sequence=["PMMA"]*50, topology="linear")

# Generate the polymer system
poly = AutoPoly.Polymerization(name="LinearPolymer", 
                              system=system, 
                              model=[linear_polymer], 
                              run=True)
```

### 2. Bead-Spring Polymer System
```python
import AutoPoly
from AutoPoly.bead_spring import BeadSpringPolymer

# Define the system
system = AutoPoly.System(out="bead_spring_test")

# Create a linear bead-spring polymer
linear_polymer = BeadSpringPolymer(
    name="linear_polymer",
    system=system,
    n_chains=10,
    n_beads=50,
    topology="linear",
    bond_length=1.0,
    mass=1.0,
    epsilon=1.0,
    sigma=1.0
)
linear_polymer.generate_data_file()

# Create a ring bead-spring polymer
ring_polymer = BeadSpringPolymer(
    name="ring_polymer",
    system=system,
    n_chains=60,
    n_beads=50,
    topology="ring",
    bond_length=1.0,
    mass=1.0,
    epsilon=1.0,
    sigma=1.0
)
ring_polymer.generate_data_file()
```

The bead-spring model features:
- Customizable number of chains and beads per chain
- Linear or ring topology
- Lennard-Jones non-bonded interactions
- Harmonic bond potentials
- Non-overlapping 3D configurations for ring polymers
- Automatic LAMMPS input script generation

### 3. Mixed Polymer System (PP and PE)
```python
# Define the system
system = AutoPoly.System(out="mixed_polymers")

# Create different polymer types
pp_ua = AutoPoly.Polymer(ChainNum=2, Sequence=["PPUA"]*20)  # United-atom PP
pe_ua = AutoPoly.Polymer(ChainNum=2, Sequence=["PEUA"]*10)  # United-atom PE
pe_aa = AutoPoly.Polymer(ChainNum=3, Sequence=["PEAA"]*15)  # All-atom PE

# Generate the mixed polymer system
poly = AutoPoly.Polymerization(name="MixedSystem",
                              system=system,
                              model=[pp_ua, pe_ua, pe_aa],
                              run=True)
```

### 4. Using RDKit to Create Custom Monomers
```python
# Create a custom monomer using SMILES
rdlt = AutoPoly.RDlt(smiles='c1c2ccccc2ccc1')  # Naphthalene
rdlt.run(to_file='Naphthalene.lt')
rdlt.store_bank()  # Store in monomer bank for future use

# Create a system using the custom monomer
system = AutoPoly.System(out="custom_system")
custom_polymer = AutoPoly.Polymer(ChainNum=50, Sequence=["Naphthalene"])

poly = AutoPoly.Polymerization(name="CustomPolymer",
                              system=system,
                              model=[custom_polymer],
                              run=True)
```

### 5. Small Molecule System
```python
# Create a system of 50 benzene molecules
system = AutoPoly.System(out="benzene_system")
benzene = AutoPoly.Polymer(ChainNum=50, Sequence=["Benzene"])

poly = AutoPoly.Polymerization(name="BenzeneSystem",
                              system=system,
                              model=[benzene],
                              run=True)
```

## Output Files

AutoPoly generates:
- LAMMPS data files (.data)
- Moltemplate files (.lt) for atomistic models
- LAMMPS input scripts
- Force field parameters
- System configuration files

## Features

### Atomistic Models
- OPLS force field support
- Tacticity control
- Custom monomer creation via RDKit
- Moltemplate integration

### Coarse-Grained Models
- Bead-spring representation
- Configurable force field parameters
- 3D spatial distribution
- Multiple chain topologies
- Ready-to-run LAMMPS configurations

## Monomer Bank

The package includes a built-in monomer bank with common polymer units. Custom monomers can be created using:
1. Avogadro (https://avogadro.cc/)
2. RDKit integration (for custom molecular structures)

## Acknowledgments

- Moltemplate: "A Tool for Coarse-Grained Modeling of Complex Biological Matter and Soft Condensed Matter Physics", J. Mol. Biol., 2021, 433(11):166841
- rdlt: Script for automating OPLS atom type assignment and .lt file generation

## Future Development

- [ ] Enhanced coarse-grained model support
- [ ] Additional polymer architectures
- [ ] Pre-relaxation using LAMMPS
- [ ] Integration with LigParGen
- [ ] Extended force field support
- [ ] More polymer topologies

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
