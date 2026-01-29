# AutoPoly API Reference

Complete API documentation for all AutoPoly classes and functions (v1.0+).

## Table of Contents

- [System](#system)
- [Polymer](#polymer)
- [Molecule](#molecule)
- [Polymerization](#polymerization)
- [BeadSpringPolymer](#beadspringpolymer)
- [MonomerGenerator](#monomergenerator)

---

## System

Manages file paths and output directories for simulation setup.

### Class Definition

```python
class System:
    """Utility class for managing file paths and system operations."""
```

### Constructor

```python
System(out: str)
```

**Parameters:**
- `out` (str): Name of the output directory where simulation files will be generated

**Example:**
```python
from AutoPoly import System

system = System(out="my_simulation")
output_path = system.get_folder_path()
print(f"Files will be saved to: {output_path}")
```

### Methods

#### `get_folder_path()`

Returns the full path to the output directory.

**Returns:**
- `str`: Absolute path to the output folder

**Example:**
```python
system = System(out="pe_simulation")
path = system.get_folder_path()
```

---

## Polymer

Defines polymer structure and properties with explicit monomer sequences.

### Class Definition

```python
class Polymer:
    """
    Represents a polymer with explicit monomer sequence.

    Uses complement SMILES format where:
    - First position: 1 wildcard on right (e.g., "CC[*]")
    - Middle positions: 2 wildcards (e.g., "[*]CC[*]")
    - Last position: 1 wildcard on left (e.g., "[*]CC")
    """
```

### Constructor

```python
Polymer(
    chain_num: int,
    sequence: list,
    topology: str = "linear",
    tacticity: str = "atactic"
)
```

**Parameters:**

- `chain_num` (int): Number of polymer chains to generate
  - Must be positive integer
  - Example: `chain_num=10` creates 10 independent chains

- `sequence` (list): Explicit monomer sequence as complement SMILES strings
  - Each element is a pSMILES string representing one monomer unit
  - DOP is automatically `len(sequence)`
  - Use complement SMILES format with proper wildcards
  - Example: `["CC[*]", "[*]CC[*]", "[*]CC[*]", "[*]CC"]`

- `topology` (str, optional): Chain topology. Default: `"linear"`
  - `"linear"`: Linear polymer chain with two end groups
  - `"ring"`: Cyclic polymer (no end groups)

- `tacticity` (str, optional): Stereochemistry configuration. Default: `"atactic"`
  - `"atactic"`: Random stereochemistry
  - `"isotactic"`: Same stereochemistry at all stereocenters
  - `"syndiotactic"`: Alternating stereochemistry

**Examples:**

```python
from AutoPoly import Polymer

# Uniform linear polymer
poly = Polymer(
    chain_num=10,
    sequence=["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"],  # DOP = 50
    topology="linear",
    tacticity="isotactic"
)

# Block copolymer with complement SMILES
sequence = [
    "CC[*]",         # First: ethylene
    "[*]CC[*]",      # Middle: ethylene
    "[*]C=C[*]",     # Middle: styrene
    "[*]C=C[*]",     # Middle: styrene
    "[*]CC"          # Last: ethylene
]
poly = Polymer(
    chain_num=5,
    sequence=sequence,  # DOP = 5
    topology="linear",
    tacticity="atactic"
)

# Ring polymer
ring_poly = Polymer(
    chain_num=5,
    sequence=["[*]CC[*]"] * 30,  # All middle positions for ring
    topology="ring",
    tacticity="atactic"
)
```

### Properties

#### `dop` (read-only)

Degree of polymerization, automatically derived from sequence length.

```python
poly = Polymer(
    chain_num=10,
    sequence=["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"]
)
print(poly.dop)  # Output: 50
```

#### `mer_set` (read-only)

Set of unique monomer types in the sequence.

```python
sequence = ["[*]CC[*]", "[*]CC[*]", "[*]C=C[*]"]
poly = Polymer(chain_num=5, sequence=sequence)
print(poly.mer_set)  # Output: {"[*]CC[*]", "[*]C=C[*]"}
```

### Methods

#### `get_chain_info()`

Returns complete information about the polymer chain.

**Returns:**
- `dict`: Dictionary containing chain properties

**Example:**
```python
poly = Polymer(
    chain_num=10,
    sequence=["[*]CC[*]"] * 50,
    topology="linear",
    tacticity="isotactic"
)

info = poly.get_chain_info()
print(info)
# {
#     'chain_num': 10,
#     'dop': 50,
#     'topology': 'linear',
#     'tacticity': 'isotactic',
#     'unique_monomers': 1,
#     'sequence_length': 50
# }
```

### Validation

The Polymer class enforces resource limits to prevent excessive memory usage:

```python
from AutoPoly.exceptions import ValidationError

# Exceeding MAX_SEQUENCE_LENGTH (10000)
try:
    poly = Polymer(chain_num=1, sequence=["[*]CC[*]"] * 15000)
except ValidationError as e:
    print(e)  # "Sequence length (15000) exceeds maximum 10000"

# Exceeding MAX_UNIQUE_MONOMERS (100)
try:
    sequence = [f"[*]C{i}C[*]" for i in range(150)]
    poly = Polymer(chain_num=1, sequence=sequence)
except ValidationError as e:
    print(e)  # "Number of unique monomers (150) exceeds maximum 100"
```

**Resource Limits:**
- `MAX_SEQUENCE_LENGTH = 10000`
- `MAX_UNIQUE_MONOMERS = 100`

---

## Molecule

Defines small molecule structures (solvents, additives) using regular SMILES notation.

### Class Definition

```python
class Molecule:
    """
    Represents small molecules like solvents or additives.

    Unlike Polymer class which uses complement SMILES with wildcards,
    Molecule uses regular SMILES notation without [*] wildcards.
    """
```

### Constructor

```python
Molecule(
    Count: int,
    Smiles: str,
    Name: str
)
```

**Parameters:**

- `Count` (int): Number of molecules to generate
  - Must be positive integer
  - Example: `Count=100` creates 100 water molecules

- `Smiles` (str): Regular SMILES string (no wildcards)
  - Must be valid SMILES notation
  - Do NOT use `[*]` wildcards
  - Example: `"O"` for water, `"CCO"` for ethanol

- `Name` (str): Identifier for the molecule type
  - Used for file naming and identification
  - Example: `"water"`, `"benzene"`, `"ethanol"`

**Examples:**

```python
from AutoPoly import Molecule

# Water molecules
water = Molecule(
    Count=100,
    Smiles="O",
    Name="water"
)

# Ethanol molecules
ethanol = Molecule(
    Count=20,
    Smiles="CCO",
    Name="ethanol"
)

# Benzene molecules
benzene = Molecule(
    Count=50,
    Smiles="c1ccccc1",
    Name="benzene"
)
```

### Common SMILES Strings

**Solvents:**
- Water: `"O"`
- Methanol: `"CO"`
- Ethanol: `"CCO"`
- Acetone: `"CC(=O)C"`
- DMSO: `"CS(C)=O"`
- Benzene: `"c1ccccc1"`
- Toluene: `"Cc1ccccc1"`

**Hydrocarbons:**
- Methane: `"C"`
- Ethane: `"CC"`
- Propane: `"CCC"`
- Cyclohexane: `"C1CCCCC1"`

**Acids/Bases:**
- Acetic acid: `"CC(=O)O"`
- Ammonia: `"N"`

### Properties

#### `molecule_name` (read-only)

Returns the name identifier for the molecule.

```python
water = Molecule(Count=100, Smiles="O", Name="water")
print(water.molecule_name)  # Output: "water"
```

---

## Polymerization

Generates polymer structures and LAMMPS input files using Moltemplate.

### Class Definition

```python
class Polymerization:
    """
    Core class for generating polymer structures.

    Coordinates the workflow:
    1. Monomer generation from SMILES
    2. Moltemplate file creation
    3. LAMMPS data file generation
    4. Force field parameter assignment
    """
```

### Constructor

```python
Polymerization(
    name: str,
    system: System,
    model: list,
    run: bool = True,
    force_field: str = "oplsaa"
)
```

**Parameters:**

- `name` (str): Project name for the simulation
  - Used for file naming
  - Example: `"polyethylene_simulation"`

- `system` (System): System object defining output directory
  - Created with `System(out="folder_name")`

- `model` (list): List of Polymer and/or Molecule objects
  - Can contain only Polymer objects
  - Can contain only Molecule objects
  - Can mix Polymer and Molecule objects
  - Example: `[polymer1, polymer2, water_molecule]`

- `run` (bool, optional): Execute generation immediately. Default: `True`
  - `True`: Generate files immediately upon object creation
  - `False`: Defer generation (call manually later)

- `force_field` (str, optional): Force field to use. Default: `"oplsaa"`
  - Options: `"oplsaa"`, `"lopls"`, `"gaff"`, `"gaff2"`, `"dreiding"`, `"compass"`
  - See [Force Field Guide](FORCE_FIELDS.md) for selection guidance

**Examples:**

```python
from AutoPoly import System, Polymer, Molecule, Polymerization

# Example 1: Single polymer
system = System(out="pe_simulation")
polymer = Polymer(
    chain_num=10,
    sequence=["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"]
)

polymerization = Polymerization(
    name="polyethylene",
    system=system,
    model=[polymer],
    force_field="oplsaa"
)

# Example 2: Polymer mixture
pe = Polymer(
    chain_num=5,
    sequence=["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"]
)
ps = Polymer(
    chain_num=5,
    sequence=["Cc1ccccc1[*]"] + ["[*]Cc1ccccc1[*]"] * 48 + ["[*]Cc1ccccc1"]
)

polymerization = Polymerization(
    name="pe_ps_blend",
    system=system,
    model=[pe, ps],
    force_field="oplsaa"
)

# Example 3: Polymer + solvent
pe = Polymer(
    chain_num=10,
    sequence=["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"]
)
water = Molecule(Count=100, Smiles="O", Name="water")

polymerization = Polymerization(
    name="pe_in_water",
    system=system,
    model=[pe, water],
    force_field="gaff"
)

# Example 4: Molecule-only system
water = Molecule(Count=100, Smiles="O", Name="water")
ethanol = Molecule(Count=20, Smiles="CCO", Name="ethanol")

polymerization = Polymerization(
    name="water_ethanol",
    system=system,
    model=[water, ethanol],
    force_field="gaff"
)
```

### Force Fields

Six force fields are supported:

| Force Field | String Value | Best For |
|------------|--------------|----------|
| OPLS-AA | `"oplsaa"` | General organic polymers |
| LOPLS | `"lopls"` | Liquid-phase simulations |
| GAFF | `"gaff"` | Small molecules, drug-like compounds |
| GAFF2 | `"gaff2"` | Updated GAFF with improved parameters |
| DREIDING | `"dreiding"` | Generic force field, metals/inorganics |
| COMPASS | `"compass"` | Commercial polymers, condensed phases |

See [docs/FORCE_FIELDS.md](FORCE_FIELDS.md) for detailed comparison and selection guide.

### Output Files

After successful execution, the output directory contains:

```
project_name/
├── moltemplate/           # Intermediate Moltemplate files
│   ├── monomer1.lt
│   ├── monomer2.lt
│   └── system.lt
├── system.data            # LAMMPS data file (coordinates, topology)
├── system.in              # LAMMPS input script
└── system.in.settings     # Force field parameters
```

### Exceptions

```python
from AutoPoly.exceptions import ValidationError

# Invalid force field
try:
    polymerization = Polymerization(
        name="test",
        system=system,
        model=[polymer],
        force_field="invalid"
    )
except ValidationError as e:
    print(e)  # "Invalid force_field 'invalid'. Must be one of: ..."
```

---

## BeadSpringPolymer

Creates coarse-grained bead-spring polymer models for simplified simulations.

### Class Definition

```python
class BeadSpringPolymer:
    """
    Coarse-grained bead-spring polymer model.

    Uses simple harmonic bonds instead of detailed atomistic interactions.
    Useful for:
    - Rapid prototyping
    - Large-scale simulations
    - Studying polymer physics
    """
```

### Constructor

```python
BeadSpringPolymer(
    name: str,
    system: System,
    n_chains: int,
    n_beads: int,
    topology: str = "linear",
    bond_length: float = 1.0,
    mass: float = 1.0
)
```

**Parameters:**

- `name` (str): Name for the bead-spring system

- `system` (System): System object defining output directory

- `n_chains` (int): Number of polymer chains

- `n_beads` (int): Number of beads per chain
  - Equivalent to degree of polymerization in atomistic models

- `topology` (str, optional): Chain topology. Default: `"linear"`
  - `"linear"`: Linear chain
  - `"ring"`: Ring polymer

- `bond_length` (float, optional): Equilibrium bond length. Default: `1.0`
  - In reduced units (LJ units typically)

- `mass` (float, optional): Bead mass. Default: `1.0`
  - In reduced units

**Example:**

```python
from AutoPoly import System, BeadSpringPolymer

system = System(out="bead_spring_simulation")

# Linear bead-spring polymer
bead_polymer = BeadSpringPolymer(
    name="coarse_grained",
    system=system,
    n_chains=5,
    n_beads=20,
    topology="linear",
    bond_length=1.0,
    mass=1.0
)

# Generate data file
bead_polymer.generate_data_file()

# Ring bead-spring polymer
ring_bead_polymer = BeadSpringPolymer(
    name="ring_cg",
    system=system,
    n_chains=10,
    n_beads=50,
    topology="ring",
    bond_length=0.97,
    mass=1.0
)
ring_bead_polymer.generate_data_file()
```

### Methods

#### `generate_data_file()`

Generates LAMMPS data file for the bead-spring system.

**Returns:**
- None (writes file to disk)

**Output:**
Creates `system.data` file in the output directory with bead-spring topology.

---

## MonomerGenerator

Advanced class for generating monomer variants with tacticity and complement SMILES.

### Class Definition

```python
class MonomerGenerator:
    """
    Generates monomer structure variants.

    Handles:
    - Tacticity (isotactic, syndiotactic, atactic)
    - Complement SMILES (first/middle/last positions)
    - 3D structure generation
    - Force field typing
    """
```

### Usage

This is an internal class typically used by the `Polymerization` workflow. Advanced users can use it directly for custom monomer generation.

**Example:**

```python
from AutoPoly import MonomerGenerator

generator = MonomerGenerator()

# Generate isotactic styrene variants
variants = generator.generate_variants(
    "[*]Cc1ccccc1[*]",
    tacticity="isotactic"
)

# variants contains first/middle/last SMILES strings with proper stereochemistry
```

---

## Exceptions

AutoPoly defines custom exceptions for error handling:

### ValidationError

Raised when input validation fails.

```python
from AutoPoly.exceptions import ValidationError

try:
    poly = Polymer(chain_num=0, sequence=["[*]CC[*]"])
except ValidationError as e:
    print(f"Validation failed: {e}")
```

**Common causes:**
- Empty sequence
- Invalid SMILES
- Sequence length exceeds limit
- Too many unique monomers
- Invalid force field
- Zero or negative chain_num

---

## Helper Functions

### Creating Uniform Sequences

For uniform polymers, use list multiplication:

```python
# Create 100 units of the same monomer
sequence = ["[*]CC[*]"] * 100

# Or define a helper function
def create_uniform_sequence(smiles: str, length: int) -> list:
    """Create uniform sequence of given length."""
    return [smiles] * length

sequence = create_uniform_sequence("[*]CC[*]", 100)
```

### Creating Block Copolymers

Build sequences by concatenating blocks:

```python
# ABA triblock
block_a = ["[*]CC[*]"] * 10
block_b = ["[*]C=C[*]"] * 20
sequence = block_a + block_b + block_a

# Or more explicitly
sequence = (
    ["[*]CC[*]"] * 10 +      # Block A (ethylene)
    ["[*]C=C[*]"] * 20 +     # Block B (styrene)
    ["[*]CC[*]"] * 10        # Block A (ethylene)
)
```

---

## Complete Workflow Example

Here's a complete example showing all major components:

```python
from AutoPoly import System, Polymer, Molecule, Polymerization

# Step 1: Create system
system = System(out="complete_example")

# Step 2: Define polymers
# Block copolymer using complement SMILES
pe_block = ["CC[*]"] + ["[*]CC[*]"] * 8 + ["[*]CC"]
ps_block = ["[*]C=C[*]"] * 5
aba_sequence = pe_block + ps_block + pe_block

aba_copolymer = Polymer(
    chain_num=5,
    sequence=aba_sequence,
    topology="linear",
    tacticity="atactic"
)

# Ring polymer
ring_polymer = Polymer(
    chain_num=3,
    sequence=["[*]CC[*]"] * 30,
    topology="ring",
    tacticity="isotactic"
)

# Step 3: Define solvent
water = Molecule(Count=100, Smiles="O", Name="water")

# Step 4: Generate system
polymerization = Polymerization(
    name="complex_system",
    system=system,
    model=[aba_copolymer, ring_polymer, water],
    force_field="gaff"
)

print("System generated successfully!")
print(f"Output: {system.get_folder_path()}")
```

---

## Additional Resources

- [Complement SMILES Guide](COMPLEMENT_SMILES.md) - Deep dive on SMILES format
- [Force Field Selection](FORCE_FIELDS.md) - Detailed force field comparison
- [Troubleshooting](TROUBLESHOOTING.md) - Common issues and solutions
- [Migration Guide](../MIGRATION.md) - Upgrading from v0.x to v1.0
- [Examples Directory](../examples/) - Full working examples
