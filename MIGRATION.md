# Migration Guide: Upgrading to AutoPoly v1.0

## ⚠️ Breaking Changes Overview

**Version 1.0 is a major breaking release** with complete API modernization. All existing code must be updated.

### Summary of Changes

- **Pythonic naming**: `chain_num` instead of `ChainNum`, `sequence` instead of `Sequence`
- **Explicit sequences**: Specify exact monomer at each position (no cycling mode)
- **DOP derived from sequence**: No separate `DOP` parameter needed
- **Block copolymers**: Now supported with explicit monomer placement

## Quick Reference

| Old API (Pre-1.0) | New API (v1.0) | Notes |
|------------------|----------------|-------|
| `ChainNum=10` | `chain_num=10` | Pythonic naming |
| `Sequence=["PE"]` | `sequence=["[*]CC[*]"] * 50` | Explicit SMILES |
| `DOP=50` | *removed* | Derived from `len(sequence)` |
| Cycling behavior | *removed* | No automatic sequence repetition |

## Migration Examples

### Example 1: Uniform Polymer

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

# Option 1: Direct list construction with first/middle/last
poly = Polymer(
    chain_num=10,  # NEW: Pythonic snake_case
    sequence=["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"],  # Complement SMILES
    topology="linear",
    tacticity="atactic"
)

# Option 2: Using helper function (recommended for clarity)
def create_uniform_sequence(first: str, middle: str, last: str, dop: int) -> list:
    """Create uniform sequence with proper first/middle/last."""
    if dop == 1:
        return [middle]
    elif dop == 2:
        return [first, last]
    return [first] + [middle] * (dop - 2) + [last]

poly = Polymer(
    chain_num=10,
    sequence=create_uniform_sequence("CC[*]", "[*]CC[*]", "[*]CC", 50),
    topology="linear",
    tacticity="atactic"
)
```

### Example 2: Block Copolymers (Now Possible!)

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

### Example 3: Multiple Unique Monomers

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

## Key Behavioral Changes

### 1. No More Cycling Mode
- **Old**: `Sequence=["A", "B"]` with `DOP=6` → `A, B, A, B, A, B`
- **New**: Must specify full sequence explicitly

### 2. Explicit SMILES Required
- **Old**: Monomer names like `"PE"`, `"PS"`
- **New**: pSMILES strings like `"[*]CC[*]"`, `"[*]C=C[*]"`

### 3. DOP is Derived
- **Old**: `DOP` parameter controls sequence repetition
- **New**: `DOP = len(sequence)` automatically

### 4. Pythonic Naming
- **Old**: `ChainNum`, `Sequence`, `DOP`
- **New**: `chain_num`, `sequence`, `dop` (read-only)

## Validation and Resource Limits

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

## Common Migration Errors

### API-related errors

**TypeError: `__init__() got an unexpected keyword argument 'ChainNum'`**
- **Fix**: Use `chain_num` instead (Pythonic naming)

**TypeError: `__init__() got an unexpected keyword argument 'DOP'`**
- **Fix**: DOP is derived from sequence length, remove this parameter

**TypeError: `__init__() got an unexpected keyword argument 'Sequence'`**
- **Fix**: Use `sequence` instead (lowercase)

### Sequence-related errors

**ValidationError: `sequence cannot be empty`**
- **Fix**: Provide at least one monomer in sequence

**ValidationError: `Sequence length exceeds maximum`**
- **Fix**: Reduce sequence length below 10000

**ValidationError: `Invalid SMILES`**
- **Fix**: Ensure all monomers use valid pSMILES notation with `[*]` wildcards

### psmiles package errors

The psmiles package is required and should be installed automatically during AutoPoly installation.

**If installation fails:**
```bash
pip install git+https://github.com/kuennethgroup/psmiles.git
```

Ensure you have Git installed if using direct Git dependency.

## Common Monomer pSMILES

Here are common monomers in the new pSMILES format:

- **Ethylene**: `"[*]CC[*]"`
- **Propylene**: `"[*]CC(C)[*]"`
- **Styrene**: `"[*]Cc1ccccc1[*]"`
- **Methyl methacrylate**: `"[*]CC([*])(C)C(=O)OC"`
- **Vinyl chloride**: `"[*]CCCl[*]"`
- **Acrylonitrile**: `"[*]CC#N[*]"`

For more monomer examples, see [docs/COMPLEMENT_SMILES.md](docs/COMPLEMENT_SMILES.md).

## Additional Resources

- [Complete API Reference](docs/API.md)
- [Complement SMILES Guide](docs/COMPLEMENT_SMILES.md)
- [Force Field Selection](docs/FORCE_FIELDS.md)
- [Troubleshooting](docs/TROUBLESHOOTING.md)
- [Examples Directory](examples/)

## Need Help?

If you encounter issues during migration:

1. Check the [Troubleshooting Guide](docs/TROUBLESHOOTING.md)
2. Review the [API Documentation](docs/API.md)
3. Examine example files in the `examples/` directory
4. Open an issue on GitHub with your migration question
