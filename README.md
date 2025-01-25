# AutoPoly.py: Automatic Generation of Data File for Polymers

AutoPoly is a Python package for automatically generating LAMMPS data files for polymer systems. It supports the OPLS (Optimized Potentials for Liquid Simulations) force field developed by Prof. William L. Jorgensen.

## Why AutoPoly?

- Create polymers with any chain length
- Support for linear and ring polymer topologies
- Easy integration with RDKit for custom monomer creation
- Automated handling of tacticity
- Built-in support for common polymer types

## Installation

Install AutoPoly via pip:
```bash
pip install AutoPoly
```

## Supported Polymers

1. PMMA (Tacticity supported)
2. PS (Tacticity supported)
3. PE (Tacticity supported)
4. PI (cis)
5. PP (Tacticity supported)
6. PVA (Tacticity supported)

## Examples

### 1. Linear Polymer (PMMA)
```python
import AutoPoly

# Define the system
system = AutoPoly.System(out="pmma_linear")
system.made_folder = True

# Create a linear PMMA polymer with 10 chains, each with 50 monomers
linear_polymer = AutoPoly.Polymer(ChainNum=10, Sequence=["PMMA"]*50, topology="linear")

# Generate the polymer system
poly = AutoPoly.Polymerization(name="LinearPolymer", 
                              system=system, 
                              model=[linear_polymer], 
                              run=True)
```

### 2. Mixed Polymer System (PP and PE)
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

### 3. Using RDKit to Create Custom Monomers
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

### 4. Small Molecule System
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
- Moltemplate files (.lt)
- Force field parameters
- System configuration files

## Monomer Bank

The package includes a built-in monomer bank with common polymer units. Custom monomers can be created using:
1. Avogadro (https://avogadro.cc/)
2. RDKit integration (for custom molecular structures)

## Acknowledgments

- Moltemplate: "A Tool for Coarse-Grained Modeling of Complex Biological Matter and Soft Condensed Matter Physics", J. Mol. Biol., 2021, 433(11):166841
- rdlt: Script for automating OPLS atom type assignment and .lt file generation

## Future Development

- [ ] Support for coarse-grained polymers
- [ ] Enhanced support for branched polymers
- [ ] Pre-relaxation using LAMMPS
- [ ] Integration with LigParGen for small molecules
- [ ] Additional polymer architectures
- [ ] Extended force field support

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
