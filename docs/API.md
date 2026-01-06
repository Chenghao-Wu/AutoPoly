# AutoPoly API Documentation

Complete API reference for AutoPoly package.

## Table of Contents

- [System](#system)
- [Polymer](#polymer)
- [Polymerization](#polymerization)
- [BeadSpringPolymer](#beadspringpolymer)

---

## System

Manages file paths and directory operations for polymer simulations.

### Class: `System`

```python
from AutoPoly import System

system = System(out="output_directory")
```

### Constructor

```python
System(out: str = None) -> None
```

**Parameters:**
- `out` (str, optional): Output directory name. If None, uses current working directory.

### Methods

#### `get_folder_path()`

Get the full path to the output directory.

```python
path = system.get_folder_path()
# Returns: "/path/to/output_directory"
```

**Returns:** `str` - Full path to the output directory

---

#### `get_output_path()`

Alias for `get_folder_path()`.

**Returns:** `str` - Full path to the output directory

---

#### `change_output_directory(new_out: str)`

Change the output directory to a new location.

```python
system.change_output_directory("new_output")
```

**Parameters:**
- `new_out` (str): New output directory name

---

#### `cleanup_output_directory()`

Remove the output directory and all its contents.

```python
system.cleanup_output_directory()
```

**Warning:** This permanently deletes all files in the output directory.

---

#### `get_FolderPath` (Deprecated)

**Deprecated:** Use `get_folder_path()` instead.

---

## Polymer

Defines polymer structures and properties including topology, tacticity, and monomer sequences.

### Class: `Polymer`

```python
from AutoPoly import Polymer

polymer = Polymer(
    ChainNum=10,
    Sequence=["PE", "PE"],
    DOP=100,
    topology="linear",
    tacticity="atactic"
)
```

### Constructor

```python
Polymer(
    ChainNum: int = None,
    Sequence: list = None,
    DOP: int = 0,
    topology: str = "linear",
    tacticity: str = "atactic"
) -> None
```

**Parameters:**
- `ChainNum` (int, optional): Number of polymer chains. Defaults to None.
- `Sequence` (list, optional): List of monomer names. Defaults to None.
- `DOP` (int, optional): Degree of polymerization. If 0, uses sequence length. Defaults to 0.
- `topology` (str, optional): Polymer topology, either "linear" (default) or "ring". Defaults to "linear".
- `tacticity` (str, optional): Polymer tacticity - "atactic" (default), "isotactic", or "syndiotactic". Defaults to "atactic".

**Raises:**
- `ValueError`: If topology is not 'linear' or 'ring', or if Sequence is None/empty.

### Methods

#### `set_merSet(merSet: Union[List[str], str])`

Set the unique set of monomers used in the polymer.

```python
polymer.set_merSet(["PE", "PS"])
```

**Parameters:**
- `merSet` (Union[List[str], str]): List of monomers or single monomer

---

#### `set_dop(dop: int)`

Set the degree of polymerization.

```python
polymer.set_dop(200)
```

**Parameters:**
- `dop` (int): Degree of polymerization

---

#### `set_Sequence()`

Set up the polymer sequence based on tacticity and chain number.

This method generates the monomer file names and names for each chain based on the specified topology and tacticity.

```python
polymer.set_Sequence()
```

---

#### `get_sequence_set()`

Get the sequence set for all chains.

```python
sequences = polymer.get_sequence_set()
# Returns: [["PEle.lt", "PEre.lt"], ["PEi.lt", ...]]
```

**Returns:** `List[List[str]]` - List of monomer file names for each chain

---

#### `get_sequence_names()`

Get the sequence names for all chains.

```python
names = polymer.get_sequence_names()
# Returns: [["PEle", "PEre"], ["PEi", ...]]
```

**Returns:** `List[List[str]]` - List of monomer names for each chain

---

#### `get_mer_set()`

Get the unique set of monomers used.

```python
monomers = polymer.get_mer_set()
# Returns: ["PE", "PS"]
```

**Returns:** `List[str]` - Unique list of monomers

---

#### `get_chain_info()`

Get comprehensive information about the polymer.

```python
info = polymer.get_chain_info()
# Returns: {
#   'chain_num': 10,
#   'sequence': ['PE', 'PE'],
#   'dop': 100,
#   'topology': 'linear',
#   'tacticity': 'atactic',
#   'sequence_length': 2,
#   'mer_set': ['PE'],
#   'sequence_set': [...],
#   'sequence_names': [...]
# }
```

**Returns:** `dict` - Dictionary containing polymer properties

---

## Polymerization

Core class for generating polymer structures using Moltemplate and preparing them for LAMMPS simulations.

### Class: `Polymerization`

```python
from AutoPoly import Polymerization

polymerization = Polymerization(
    name="project_name",
    system=system,
    model=[polymer],
    force_field="oplsaa"
)
```

### Constructor

```python
Polymerization(
    name: str = None,
    system: object = None,
    model: list = None,
    run: bool = True,
    path_monomer_bank: str = None,
    is_lopls: bool = False,
    force_field: str = "oplsaa"
) -> None
```

**Parameters:**
- `name` (str, optional): Name of the polymerization project. Defaults to None.
- `system` (object, optional): System object containing folder path. Defaults to None.
- `model` (list, optional): List of Polymer objects for polymerization. Defaults to None.
- `run` (bool, optional): Flag to run the process immediately. Defaults to True.
- `path_monomer_bank` (str, optional): Path to the monomer bank. Defaults to None.
- `is_lopls` (bool, optional): **Deprecated:** Use `force_field="lopls"` instead. Whether to use LOPLS force field. Defaults to False.
- `force_field` (str, optional): Force field to use - "oplsaa" (default), "gaff", or "lopls". Defaults to "oplsaa".

**Raises:**
- `SystemExit`: If required directories or files are not found, or if invalid force_field is specified.

### Methods

#### `create_working_directory()`

Create and manage the working directory structure for the polymerization.

```python
polymerization.create_working_directory()
```

---

#### `set_tacticity(tacticity: str)`

Set the tacticity of the polymer.

```python
polymerization.set_tacticity("isotactic")
```

**Parameters:**
- `tacticity` (str): The tacticity to set

---

## BeadSpringPolymer

Simplified bead-spring polymer model generator for coarse-grained LAMMPS simulations.

### Class: `BeadSpringPolymer`

```python
from AutoPoly import BeadSpringPolymer

bead_polymer = BeadSpringPolymer(
    name="coarse_grained",
    system=system,
    n_chains=5,
    n_beads=20,
    topology="linear"
)
```

### Constructor

```python
BeadSpringPolymer(
    name: str = None,
    system: object = None,
    n_chains: int = 1,
    n_beads: int = 10,
    topology: str = "linear",
    bond_length: float = 1.0,
    mass: float = 1.0,
    epsilon: float = 1.0,
    sigma: float = 1.0
) -> None
```

**Parameters:**
- `name` (str, optional): Name for the output files. Defaults to None.
- `system` (object, optional): System object containing path information. Defaults to None.
- `n_chains` (int, optional): Number of polymer chains. Defaults to 1.
- `n_beads` (int, optional): Number of beads per chain. Defaults to 10.
- `topology` (str, optional): "linear" (default) or "ring". Defaults to "linear".
- `bond_length` (float, optional): Equilibrium bond length. Defaults to 1.0.
- `mass` (float, optional): Mass of each bead. Defaults to 1.0.
- `epsilon` (float, optional): LJ energy parameter. Defaults to 1.0.
- `sigma` (float, optional): LJ distance parameter. Defaults to 1.0.

**Raises:**
- `ValueError`: If topology is not 'linear' or 'ring'

### Methods

#### `generate_data_file()`

Generate LAMMPS data file for bead-spring polymer.

```python
bead_polymer.generate_data_file()
```

This method creates:
- `polymer.data` - LAMMPS data file with atoms and bonds
- `in.polymer` - LAMMPS input script

---

#### `get_system_info()`

Get comprehensive information about the bead-spring polymer system.

```python
info = bead_polymer.get_system_info()
# Returns: {
#   'name': 'coarse_grained',
#   'n_chains': 5,
#   'n_beads_per_chain': 20,
#   'topology': 'linear',
#   'total_beads': 100,
#   'total_bonds': 95,
#   'bond_length': 1.0,
#   'mass': 1.0,
#   'epsilon': 1.0,
#   'sigma': 1.0,
#   'output_path': '/path/to/output'
# }
```

**Returns:** `dict` - Dictionary containing system properties

---

#### `modify_parameters(**kwargs)`

Modify polymer parameters after initialization.

```python
bead_polymer.modify_parameters(
    n_chains=10,
    n_beads=30,
    bond_length=1.5,
    mass=2.0
)
```

**Parameters:**
- `**kwargs`: Keyword arguments for parameters to modify
  - Valid params: `n_chains`, `n_beads`, `topology`, `bond_length`, `mass`, `epsilon`, `sigma`

**Raises:**
- `ValueError`: If topology is modified to invalid value

---

## Common Workflows

### Basic Atomistic Polymer

```python
from AutoPoly import System, Polymer, Polymerization

# 1. Create system
system = System(out="pe_simulation")

# 2. Define polymer
polymer = Polymer(
    ChainNum=10,
    Sequence=["PE"],
    DOP=100,
    topology="linear",
    tacticity="atactic"
)

# 3. Generate structure
polymerization = Polymerization(
    name="polyethylene",
    system=system,
    model=[polymer],
    force_field="oplsaa"
)
```

### Coarse-Grained Simulation

```python
from AutoPoly import System, BeadSpringPolymer

system = System(out="cg_simulation")
polymer = BeadSpringPolymer(
    name="cg_polymer",
    system=system,
    n_chains=5,
    n_beads=50,
    topology="ring"
)
polymer.generate_data_file()
```

### Multiple Polymers

```python
# Linear PE
pe = Polymer(ChainNum=5, Sequence=["PE"], DOP=50, topology="linear")

# Ring PS
ps = Polymer(ChainNum=3, Sequence=["PS"], DOP=30, topology="ring")

# Generate both
polymerization = Polymerization(
    name="mixed_system",
    system=system,
    model=[pe, ps]
)
```

## Type Hints

AutoPoly uses Python type hints for better IDE support:

```python
from typing import List, Optional, Union

def get_folder_path(self) -> str: ...
def set_merSet(self, merSet: Union[List[str], str]) -> None: ...
```

## Error Handling

```python
# Invalid topology
try:
    polymer = Polymer(topology="invalid")
except ValueError as e:
    print(f"Error: {e}")

# Missing monomer files
try:
    polymerization = Polymerization(name="test", system=system, model=[polymer])
except SystemExit as e:
    print("Check monomer bank path")
```
