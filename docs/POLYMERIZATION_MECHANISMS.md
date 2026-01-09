# Polymerization Mechanisms in AutoPoly

This document explains the mechanism-based polymerization system in AutoPoly, which enables element-agnostic monomer generation for both vinyl addition and condensation polymers.

## Table of Contents
1. [Overview](#overview)
2. [Available Mechanisms](#available-mechanisms)
3. [Mechanism Detection](#mechanism-detection)
4. [DOP-Aware Detection](#dop-aware-detection)
5. [Using Mechanisms](#using-mechanisms)
6. [SMARTS Patterns](#smarts-patterns)
7. [Connection Points](#connection-points)

---

## Overview

AutoPoly uses a **SMARTS-based** polymerization system that:
- **Detects** the polymerization mechanism from molecular structure
- **Identifies** connection points for any element (C, O, N, etc.)
- **Modifies** atom types appropriately for end-caps
- **Preserves** single molecule generation (DOP=1)

This replaces the previous carbon-only hard-coded approach with a flexible, extensible system.

### Key Benefits

1. **Element-Agnostic**: Works with C, O, N backbones and beyond
2. **Extensible**: Easy to add new mechanisms
3. **SMARTS-Based**: Uses chemical pattern matching
4. **DOP-Aware**: Respects Degree of Polymerization
5. **Backward Compatible**: All existing vinyl polymers still work

---

## Available Mechanisms

### 1. None (`none`)
**Type:** Non-polymerizable
**Connection Atoms:** None
**Description:** Single molecules, solvents, non-reactive species

**Examples:**
- Water (O)
- Ethanol (CCO)
- Methane (C)

**Use Case:**
```python
generator = MonomerGenerator(
    base_name="ethanol",
    output_dir="./monomers",
    mechanism='none',
    is_gaff=False
)
variants = generator.generate_variants(smiles="CCO")
```

### 2. Vinyl Addition (`vinyl_addition`)
**Type:** Addition polymerization
**Connection Atoms:** C, C
**SMARTS:** `[$([C]=[C])]`
**Description:** C=C double bond opening (most common)

**Examples:**
- Polyethylene (C=C)
- Polypropylene (C=C(C))
- Polystyrene (C=C(C1=CC=CC=C1))
- PMMA (C=C(C)C(=O)OC)

**Use Case:**
```python
generator = MonomerGenerator(
    base_name="PE",
    output_dir="./monomers",
    mechanism='vinyl_addition',
    is_gaff=False
)
variants = generator.generate_variants(smiles="C=C")
```

### 3. Esterification (`esterification`)
**Type:** Condensation polymerization
**Connection Atoms:** C, O
**SMARTS:** `[$([C](=[O])[OX2H0])][$([OX2H])]`
**Description:** Carboxyl + alcohol → polyester (releases H₂O)

**Examples:**
- PLA from lactic acid
- PET from ethylene glycol + terephthalic acid
- PCL from caprolactone

**Use Case:**
```python
generator = MonomerGenerator(
    base_name="PLA",
    output_dir="./monomers",
    mechanism='esterification',
    is_gaff=True  # GAFF recommended for polyesters
)
variants = generator.generate_variants(smiles="CC(C(=O)O)O")
```

### 4. Amidation (`amidation`)
**Type:** Condensation polymerization
**Connection Atoms:** C, N
**SMARTS:** `[$([C](=[O])[OX2H0])][$([NX3H])]`
**Description:** Carboxyl + amine → polyamide (releases H₂O)

**Examples:**
- Nylon-6 from caprolactam
- Nylon-6,6 from hexamethylenediamine + adipic acid
- Kevlar precursors

**Use Case:**
```python
generator = MonomerGenerator(
    base_name="Nylon66",
    output_dir="./monomers",
    mechanism='amidation',
    is_gaff=True  # GAFF recommended for polyamides
)
variants = generator.generate_variants(smiles="NCC(=O)O")
```

### 5. Etherification (`etherification`)
**Type:** Condensation polymerization
**Connection Atoms:** O, O
**SMARTS:** `[$([OX2H])]`
**Description:** Alcohol + alcohol → polyether (releases H₂O)

**Examples:**
- PEG from ethylene glycol
- PPG from propylene glycol
- PTMG from tetrahydrofuran

**Use Case:**
```python
generator = MonomerGenerator(
    base_name="PEG",
    output_dir="./monomers",
    mechanism='etherification',
    is_gaff=False  # OPLS-AA works for polyethers
)
variants = generator.generate_variants(smiles="OCCO")
```

---

## Mechanism Detection

### Auto-Detection

AutoPoly can automatically detect the polymerization mechanism from SMILES:

```python
from AutoPy.polymerization_mechanism import detect_mechanism
from rdkit import Chem

# Detect mechanism for lactic acid
mol = Chem.MolFromSmiles("CC(C(=O)O)O")
mechanism = detect_mechanism(mol, dop=10)  # dop>1 for polymerization
print(mechanism)  # 'esterification'

# Detect mechanism for ethylene
mol = Chem.MolFromSmiles("C=C")
mechanism = detect_mechanism(mol, dop=10)
print(mechanism)  # 'vinyl_addition'
```

### Detection Priority Order

When `dop>1`, mechanisms are checked in this order:

1. **Vinyl addition** - If C=C double bond present
2. **Esterification** - If both carboxyl and alcohol groups present
3. **Amidation** - If both carboxyl and amine groups present
4. **Etherification** - If 2+ alcohol groups present
5. **None** - Default (no clear pattern)

### Explicit Mechanism Specification

You can also explicitly specify the mechanism:

```python
generator = MonomerGenerator(
    base_name="PLA",
    output_dir="./monomers",
    mechanism='esterification',  # Explicit mechanism
    is_gaff=True
)
```

---

## DOP-Aware Detection

The mechanism detector is **DOP-aware** (Degree of Polymerization):

### DOP = 1 (Single Molecule)

```python
mol = Chem.MolFromSmiles("CC(C(=O)O)O")  # Lactic acid

# DOP=1: Returns 'none' (don't polymerize)
mechanism = detect_mechanism(mol, dop=1)
print(mechanism)  # 'none'

# Generate single molecule
generator = MonomerGenerator(
    base_name="lactic_acid",
    output_dir="./monomers",
    mechanism='none',
    is_gaff=True
)
variants = generator.generate_variants(smiles="CC(C(=O)O)O")
# Generates: plain molecule with all functional groups intact
```

### DOP > 1 (Polymer)

```python
mol = Chem.MolFromSmiles("CC(C(=O)O)O")  # Lactic acid

# DOP>1: Returns 'esterification'
mechanism = detect_mechanism(mol, dop=10)
print(mechanism)  # 'esterification'

# Generate polymer
generator = MonomerGenerator(
    base_name="PLA",
    output_dir="./monomers",
    mechanism='esterification',
    is_gaff=True
)
variants = generator.generate_variants(smiles="CC(C(=O)O)O")
# Generates: internal, left_end, right_end variants
```

### Special Case: Vinyl with DOP=1

Vinyl compounds are still detected as vinyl even with DOP=1:

```python
mol = Chem.MolFromSmiles("C=C")

# Still detects vinyl_addition even with DOP=1
mechanism = detect_mechanism(mol, dop=1)
print(mechanism)  # 'vinyl_addition'
```

---

## Using Mechanisms

### Basic Usage

```python
from AutoPy.monomer_generator import MonomerGenerator

# Auto-detect mechanism (default)
generator = MonomerGenerator(
    base_name="PLA",
    output_dir="./monomers",
    is_gaff=True
)
# Mechanism will be auto-detected based on SMILES and DOP

# Explicit mechanism
generator = MonomerGenerator(
    base_name="PLA",
    output_dir="./monomers",
    mechanism='esterification',
    is_gaff=True
)
# Forces esterification mechanism
```

### With Polymer Class

```python
from AutoPy import System, Polymer, Polymerization

system = System(out="pla_output")

polymer = Polymer(
    ChainNum=10,
    Sequence=["CC(C(=O)O)O"] * 50,  # PLA SMILES
    DOP=50,                         # Degree of polymerization
    topology="linear"
)

poly = Polymerization(
    name="pla",
    system=system,
    model=[polymer],
    force_field="gaff"  # GAFF for heteroatom polymers
)
```

---

## SMARTS Patterns

### What are SMARTS?

SMARTS (SMILES Arbitrary Target Specification) is a language for specifying molecular patterns:

- `C` - Carbon atom
- `[C]=[C]` - Carbon double-bonded to carbon
- `[$([C](=[O])]` - Carbon with double bond to oxygen (carboxyl carbon)
- `[$([OX2H])]` - Oxygen with single bonds and hydrogen (alcohol oxygen)

### Vinyl Addition SMARTS

```python
smarts = '[$([C]=[C])]'
# Matches: C=C double bonds
```

**Matches:**
- Ethylene: C=C ✓
- Propylene: C=CC ✓
- Butadiene: C=CC=C ✓

**Does not match:**
- Ethane: CC ✗
- Ethanol: CCO ✗

### Esterification SMARTS

```python
smarts = '[$([C](=[O])[OX2H0])][$([OX2H])]'
# Matches: Carboxyl carbon + Alcohol oxygen
```

**Matches:**
- Lactic acid: CC(C(=O)O)O ✓
- Acetic acid + ethanol: CC(=O)O + CCO ✓

**Requires:**
- Carboxyl group: C(=O)OH
- Alcohol group: O-H

### Amidation SMARTS

```python
smarts = '[$([C](=[O])[OX2H0])][$([NX3H])]'
# Matches: Carboxyl carbon + Amine nitrogen
```

**Matches:**
- Glycine: NCC(=O)O ✓
- Amino acids ✓

**Requires:**
- Carboxyl group: C(=O)OH
- Amine group: N-H

### Etherification SMARTS

```python
smarts = '[$([OX2H])]'
# Matches: Alcohol oxygen
```

**Matches:**
- Ethylene glycol: OCCO ✓
- Glycerol: OCC(C)O ✓

**Requires:**
- 2+ alcohol groups

---

## Connection Points

### What are Connection Points?

Connection points are the atoms that form bonds between monomers during polymerization:

| Mechanism | Connection Atoms | Bond Formed |
|-----------|-----------------|-------------|
| vinyl_addition | C, C | C-C single bond |
| esterification | C, O | C-O ester bond |
| amidation | C, N | C-N amide bond |
| etherification | O, O | C-O-C ether bond |
| none | None | No bond |

### Connection Point Modification

When creating end-cap variants (le/re), AutoPoly:

1. **Identifies** the connection atom
2. **Removes** one hydrogen (if applicable)
3. **Changes** the atom type for polymer bonding

**Example (Vinyl - Carbon):**
```
CH3 (methyl) → CH2 (methylene)
@atom:80     → @atom:82
```

**Example (Ester - Oxygen):**
```
OH (alcohol) → O (ether)
@atom:96      → @atom:122
```

**Example (Amide - Nitrogen):**
```
NH2 (amine) → NH (secondary amine)
@atom:739    → @atom:740
```

### Element-Agnostic Handling

The `ConnectionPointModifier` class handles all elements:

```python
from AutoPy.connection_point import ConnectionPointModifier

modifier = ConnectionPointModifier(force_field='oplsaa')

# Works with carbon
carbon_atom = mol.GetAtomWithIdx(0)
original, connection = modifier.get_connection_atom_type(carbon_atom, mol)

# Works with oxygen
oxygen_atom = mol.GetAtomWithIdx(5)
original, connection = modifier.get_connection_atom_type(oxygen_atom, mol)

# Works with nitrogen
nitrogen_atom = mol.GetAtomWithIdx(3)
original, connection = modifier.get_connection_atom_type(nitrogen_atom, mol)
```

---

## API Reference

### PolymerizationMechanism Class

```python
from AutoPy.polymerization_mechanism import PolymerizationMechanism

detector = PolymerizationMechanism(verbose=True)

# Detect mechanism
mechanism = detector.detect_mechanism(mol, dop=10)

# Get connection atoms
atoms = detector.get_connection_atoms(mol, mechanism)

# Check if condensation
is_cond = detector.is_condensation_polymerization(mechanism)

# Get mechanism info
info = detector.get_mechanism_info(mechanism)
```

### Convenience Function

```python
from AutoPy.polymerization_mechanism import detect_mechanism

mechanism = detect_mechanism(mol, dop=10, verbose=True)
```

### ConnectionPointModifier Class

```python
from AutoPy.connection_point import ConnectionPointModifier

modifier = ConnectionPointModifier(force_field='oplsaa', verbose=True)

# Get connection atom type
original, connection = modifier.get_connection_atom_type(atom, mol)

# Check if should skip H removal
skip = modifier.should_skip_hydrogen_removal(element, type_name)
```

---

## Best Practices

### 1. Force Field Selection

**For Vinyl Polymers:**
- Use OPLS-AA (most common)
- PE, PP, PS work well with OPLS-AA

**For Heteroatom Polymers:**
- Use GAFF (recommended)
- PLA, Nylon, PEG have better coverage with GAFF

**Example:**
```python
# Vinyl: OPLS-AA
generator = MonomerGenerator(
    base_name="PE",
    output_dir="./monomers",
    mechanism='vinyl_addition',
    is_gaff=False  # OPLS-AA
)

# Polyester: GAFF
generator = MonomerGenerator(
    base_name="PLA",
    output_dir="./monomers",
    mechanism='esterification',
    is_gaff=True  # GAFF
)
```

### 2. DOP Parameter

**For Single Molecules:**
- Set `DOP=1`
- Mechanism will auto-detect as 'none' (unless vinyl)

**For Polymers:**
- Set `DOP>1` (e.g., DOP=10, DOP=100)
- Mechanism will be properly detected

### 3. Mechanism Specification

**Auto-Detect (Recommended):**
```python
generator = MonomerGenerator(
    base_name="PLA",
    output_dir="./monomers",
    # mechanism omitted - will auto-detect
    is_gaff=True
)
```

**Explicit (When Needed):**
```python
generator = MonomerGenerator(
    base_name="PLA",
    output_dir="./monomers",
    mechanism='esterification',  # Force specific mechanism
    is_gaff=True
)
```

---

## Troubleshooting

### Problem: "Mechanism not detected"

**Possible Causes:**
1. SMILES doesn't match any pattern
2. Functional groups not present
3. DOP=1 for non-vinyl molecule

**Solutions:**
```python
# Check detected mechanism
from AutoPy.polymerization_mechanism import detect_mechanism
mechanism = detect_mechanism(mol, dop=10)
print(f"Detected: {mechanism}")

# Force mechanism if needed
generator = MonomerGenerator(
    base_name="custom",
    output_dir="./monomers",
    mechanism='esterification',  # Force specific mechanism
    is_gaff=True
)
```

### Problem: "Atom typing failed"

**Possible Causes:**
1. Unusual atom types not in force field
2. Missing GAFF parameters

**Solutions:**
```python
# Try different force field
generator = MonomerGenerator(
    base_name="PLA",
    output_dir="./monomers",
    mechanism='esterification',
    is_gaff=True  # Try GAFF instead of OPLS
)
```

### Problem: "Wrong connection atoms"

**Possible Causes:**
1. SMARTS pattern matching multiple sites
2. Symmetric molecules

**Solutions:**
```python
# Explicitly specify mechanism
generator = MonomerGenerator(
    base_name="custom",
    output_dir="./monomers",
    mechanism='esterification',
    is_gaff=True
)
```

---

## Additional Resources

- **SMILES Guide:** See `SMILES_GUIDE.md`
- **API Documentation:** See `API.md`
- **Test Files:** See `tests/test_integration/` for examples

---

## Version History

- **v1.0** - Initial implementation with vinyl addition only
- **v2.0** - Added heteroatom support (PLA, Nylon, PEG)
  - SMARTS-based mechanism detection
  - Element-agnostic connection point handling
  - DOP-aware single molecule support

---

For questions or issues, please refer to the main documentation or open an issue on GitHub.
