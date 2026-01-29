# Troubleshooting Guide

Common issues and solutions for AutoPoly.

## Table of Contents

- [API Errors](#api-errors)
- [Sequence Errors](#sequence-errors)
- [SMILES Errors](#smiles-errors)
- [Force Field Errors](#force-field-errors)
- [Moltemplate Errors](#moltemplate-errors)
- [Installation Issues](#installation-issues)
- [Performance Issues](#performance-issues)

---

## API Errors

### TypeError: Unexpected keyword argument 'ChainNum'

**Error:**
```python
TypeError: __init__() got an unexpected keyword argument 'ChainNum'
```

**Cause:** Using old v0.x API with v1.0+ AutoPoly

**Solution:** Update to snake_case naming:
```python
# OLD (v0.x)
poly = Polymer(ChainNum=10, Sequence=["PE"], DOP=50)

# NEW (v1.0+)
poly = Polymer(
    chain_num=10,
    sequence=["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"]
)
```

**See:** [Migration Guide](../MIGRATION.md)

---

### TypeError: Unexpected keyword argument 'DOP'

**Error:**
```python
TypeError: __init__() got an unexpected keyword argument 'DOP'
```

**Cause:** DOP parameter removed in v1.0

**Solution:** DOP is automatically `len(sequence)`:
```python
# OLD (v0.x)
poly = Polymer(ChainNum=10, Sequence=["PE"], DOP=50)

# NEW (v1.0+)
# DOP = 50 automatically (first + 48 middle + last)
sequence = ["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"]
poly = Polymer(chain_num=10, sequence=sequence)
print(poly.dop)  # Output: 50
```

---

### TypeError: Unexpected keyword argument 'Sequence'

**Error:**
```python
TypeError: __init__() got an unexpected keyword argument 'Sequence'
```

**Cause:** Using CamelCase instead of snake_case

**Solution:** Use lowercase `sequence`:
```python
# OLD
poly = Polymer(ChainNum=10, Sequence=["PE"])

# NEW
poly = Polymer(chain_num=10, sequence=["[*]CC[*]"] * 50)
```

---

## Sequence Errors

### ValidationError: sequence cannot be empty

**Error:**
```python
ValidationError: sequence cannot be empty
```

**Cause:** Empty sequence list provided

**Solution:** Provide at least one monomer:
```python
# BAD
sequence = []

# GOOD
sequence = ["[*]CC[*]"] * 10
```

---

### ValidationError: Sequence length exceeds maximum

**Error:**
```python
ValidationError: Sequence length (15000) exceeds maximum 10000
```

**Cause:** Sequence longer than MAX_SEQUENCE_LENGTH (10000)

**Solution:** Reduce sequence length:
```python
# BAD
sequence = ["[*]CC[*]"] * 15000  # Too long

# GOOD
sequence = ["[*]CC[*]"] * 5000  # Within limit
```

**Note:** 10000 is a safety limit to prevent memory exhaustion.

---

### ValidationError: Number of unique monomers exceeds maximum

**Error:**
```python
ValidationError: Number of unique monomers (150) exceeds maximum 100
```

**Cause:** More than 100 unique monomer types in sequence

**Solution:** Reduce variety of monomers:
```python
# BAD - 150 different monomers
sequence = [f"[*]C{i}C[*]" for i in range(150)]

# GOOD - Reuse monomers
pe = "[*]CC[*]"
ps = "[*]Cc1ccccc1[*]"
sequence = [pe] * 50 + [ps] * 50  # Only 2 unique types
```

---

## SMILES Errors

### ValidationError: Invalid SMILES

**Error:**
```python
ValidationError: Invalid SMILES in sequence at position 5
```

**Causes & Solutions:**

#### 1. Wrong wildcard count
```python
# BAD - First position with 2 wildcards
sequence = ["[*]CC[*]", ...]  # Should be "CC[*]"

# GOOD - Correct wildcard count
sequence = ["CC[*]", "[*]CC[*]", ..., "[*]CC"]
```

#### 2. Missing brackets around wildcard
```python
# BAD
sequence = ["*CC*"]  # Wrong syntax

# GOOD
sequence = ["[*]CC[*]"]  # Correct syntax
```

#### 3. SMILES syntax error
```python
# BAD - Invalid SMILES
sequence = ["[*]C(C[*]"]  # Missing closing parenthesis

# GOOD
sequence = ["[*]CC([*])C"]  # Valid SMILES
```

**Validation tip:**
```python
from rdkit import Chem

def validate_smiles(smiles: str) -> bool:
    """Validate SMILES syntax."""
    # Replace wildcards with carbons for validation
    test_smiles = smiles.replace("[*]", "C")
    mol = Chem.MolFromSmiles(test_smiles)
    return mol is not None

# Check before using
if not validate_smiles("[*]CC[*]"):
    print("Invalid SMILES!")
```

---

## Force Field Errors

### ValidationError: Invalid force_field

**Error:**
```python
ValidationError: Invalid force_field 'invalid'. Must be one of: ['oplsaa', 'lopls', 'gaff', 'gaff2', 'dreiding', 'compass']
```

**Cause:** Typo or unsupported force field name

**Solution:** Use valid force field name:
```python
# BAD
force_field="opls"    # Wrong name

# GOOD
force_field="oplsaa"  # Correct name
```

**Valid options:**
- `"oplsaa"` - OPLS-AA
- `"lopls"` - Liquid-OPLS
- `"gaff"` - GAFF
- `"gaff2"` - GAFF2
- `"dreiding"` - DREIDING
- `"compass"` - COMPASS

**See:** [Force Field Guide](FORCE_FIELDS.md)

---

### GAFF Charge Calculation Error

**Problem:** GAFF simulations produce incorrect energies

**Cause:** GAFF requires explicit charge calculation, which AutoPoly doesn't do automatically

**Solution:**

1. Generate initial files with AutoPoly:
```python
polymerization = Polymerization(
    name="system",
    system=system,
    model=[molecules],
    force_field="gaff"
)
```

2. Calculate charges using Antechamber:
```bash
cd output_folder
antechamber -i molecule.mol2 -fi mol2 -o charged.mol2 -fo mol2 -c bcc
```

3. Update `system.in.charges` file with calculated charges

**See:** [Force Field Guide - GAFF section](FORCE_FIELDS.md#gaff-general-amber-force-field)

---

## Moltemplate Errors

### Moltemplate not found

**Error:**
```
SystemExit: Moltemplate not found in PATH
```

**Cause:** Moltemplate not installed or not in PATH

**Solution:**

**Option 1: Install via pip**
```bash
pip install moltemplate
```

**Option 2: Install from source**
```bash
git clone https://github.com/jewettaij/moltemplate
cd moltemplate
pip install .
```

**Verify installation:**
```bash
which moltemplate.sh
# Should print path to moltemplate
```

---

### Monomer .lt file generation failed

**Error:**
```
Error: Failed to generate monomer .lt file for [*]CC[*]
```

**Causes & Solutions:**

#### 1. Invalid SMILES
Check SMILES syntax using RDKit or online SMILES validator

#### 2. Missing dependencies
Ensure OpenBabel is installed:
```bash
# Ubuntu/Debian
sudo apt-get install openbabel

# macOS
brew install open-babel

# Conda
conda install -c conda-forge openbabel
```

#### 3. Unsupported functional groups
Some exotic groups may not be supported. Try:
- Simplifying the structure
- Using a different force field
- Using DREIDING (more general)

---

## Installation Issues

### psmiles package installation failed

**Error:**
```
ERROR: Could not install psmiles from GitHub
```

**Cause:** Git not installed or GitHub access issues

**Solution:**

**Option 1: Install Git**
```bash
# Ubuntu/Debian
sudo apt-get install git

# macOS (install Xcode Command Line Tools)
xcode-select --install

# Windows
# Download from https://git-scm.com/
```

**Option 2: Manual install**
```bash
git clone https://github.com/kuennethgroup/psmiles.git
cd psmiles
pip install .
```

---

### Import Error: No module named 'AutoPoly'

**Error:**
```python
ImportError: No module named 'AutoPoly'
```

**Cause:** AutoPoly not installed or not in Python path

**Solution:**

**Option 1: Install in development mode**
```bash
cd /path/to/AutoPoly
pip install -e .
```

**Option 2: Install from source**
```bash
cd /path/to/AutoPoly
pip install .
```

**Verify installation:**
```python
import AutoPoly
print(AutoPoly.__version__)  # Should print version number
```

---

### Dependency conflicts

**Error:**
```
ERROR: Cannot install AutoPoly due to conflicting dependencies
```

**Solution:** Create fresh virtual environment:
```bash
# Create new environment
python -m venv autopoly_env

# Activate
source autopoly_env/bin/activate  # Linux/macOS
autopoly_env\Scripts\activate     # Windows

# Install AutoPoly
pip install -e /path/to/AutoPoly
```

---

## Performance Issues

### Generation takes too long

**Problem:** Polymerization takes hours for large systems

**Solutions:**

#### 1. Reduce system size for testing
```python
# Start small
poly = Polymer(chain_num=5, sequence=["[*]CC[*]"] * 20)

# Then scale up
poly = Polymer(chain_num=100, sequence=["[*]CC[*]"] * 1000)
```

#### 2. Use bead-spring models for large systems
```python
# Atomistic (slow)
poly = Polymer(chain_num=100, sequence=["[*]CC[*]"] * 1000)

# Bead-spring (fast)
bead_poly = BeadSpringPolymer(
    name="cg",
    system=system,
    n_chains=100,
    n_beads=1000
)
```

#### 3. Parallelize multiple systems
Run independent systems in parallel rather than one large system.

---

### Memory errors

**Error:**
```
MemoryError: Unable to allocate array
```

**Cause:** System too large for available RAM

**Solutions:**

#### 1. Reduce chain_num or DOP
```python
# Too large
poly = Polymer(chain_num=1000, sequence=["[*]CC[*]"] * 10000)

# Reasonable
poly = Polymer(chain_num=100, sequence=["[*]CC[*]"] * 1000)
```

#### 2. Check sequence length
```python
# Respect limits
MAX_SEQUENCE_LENGTH = 10000
sequence = ["[*]CC[*]"] * min(desired_dop, MAX_SEQUENCE_LENGTH)
```

#### 3. Use 64-bit Python
Ensure you're using 64-bit Python for large systems:
```bash
python -c "import sys; print(sys.maxsize > 2**32)"
# Should print: True
```

---

## File Permission Errors

### Permission denied writing to output directory

**Error:**
```
PermissionError: [Errno 13] Permission denied: 'output_folder/system.data'
```

**Solutions:**

#### 1. Check directory permissions
```bash
ls -ld output_folder
# Should show write permissions
```

#### 2. Create directory first
```python
import os
os.makedirs("output_folder", exist_ok=True)

system = System(out="output_folder")
```

#### 3. Use different output location
```python
# Use home directory
import os
home = os.path.expanduser("~")
system = System(out=os.path.join(home, "autopoly_output"))
```

---

## Ring Polymer Issues

### Ring polymer generation fails

**Problem:** Ring polymers produce errors or incorrect structures

**Solutions:**

#### 1. Use only middle variants
```python
# BAD - Using first/last for ring
sequence = ["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"]

# GOOD - All middle variants
sequence = ["[*]CC[*]"] * 50
```

#### 2. Specify topology correctly
```python
poly = Polymer(
    chain_num=10,
    sequence=["[*]CC[*]"] * 50,
    topology="ring"  # Must specify!
)
```

---

## Common Workflow Errors

### Molecule class vs Polymer class confusion

**Problem:** Using Polymer for small molecules or Molecule for polymers

**Solution:**

**For polymers** (use Polymer):
```python
# Correct
poly = Polymer(
    chain_num=10,
    sequence=["[*]CC[*]"] * 50,  # Complement SMILES with wildcards
    topology="linear"
)
```

**For small molecules** (use Molecule):
```python
# Correct
water = Molecule(
    Count=100,
    Smiles="O",  # Regular SMILES, NO wildcards
    Name="water"
)
```

**See:** [API Documentation - Molecule vs Polymer](API.md#molecule)

---

## Getting Help

If your issue isn't covered here:

1. Check the [API Documentation](API.md)
2. Review [examples/](../examples/) for working code
3. Read the [Complement SMILES Guide](COMPLEMENT_SMILES.md)
4. Check the [Force Field Guide](FORCE_FIELDS.md)
5. Review the [Migration Guide](../MIGRATION.md) if upgrading
6. Open an issue on GitHub with:
   - Error message (full traceback)
   - Minimal code to reproduce
   - AutoPoly version (`AutoPoly.__version__`)
   - Python version (`python --version`)
   - Operating system

---

## Quick Diagnostic Commands

```python
# Check AutoPoly version
import AutoPoly
print(f"AutoPoly version: {AutoPoly.__version__}")

# Check Python version
import sys
print(f"Python version: {sys.version}")

# Check if Moltemplate is available
import shutil
print(f"Moltemplate found: {shutil.which('moltemplate.sh') is not None}")

# Check dependencies
try:
    import rdkit
    print("RDKit: OK")
except ImportError:
    print("RDKit: NOT FOUND")

try:
    import psmiles
    print("pSMILES: OK")
except ImportError:
    print("pSMILES: NOT FOUND")
```

---

**Still stuck?** Open an issue on GitHub with diagnostic output and we'll help!
