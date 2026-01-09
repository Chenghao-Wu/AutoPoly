# SMILES Guide for AutoPoly

This guide explains how to use SMILES (Simplified Molecular-Input Line-Entry System) notation to define polymer monomers in AutoPoly.

## Table of Contents
1. [What is SMILES?](#what-is-smiles)
2. [Basic SMILES Syntax](#basic-smiles-syntax)
3. [Common Polymer Monomers](#common-polymer-monomers)
4. [Heteroatom Polymers (NEW)](#heteroatom-polymers)
5. [Single Molecule Generation (NEW)](#single-molecule-generation)
6. [pSMILES for Ring Polymers](#psmiles-for-ring-polymers)
7. [Finding SMILES for Custom Monomers](#finding-smiles-for-custom-monomers)
8. [Examples](#examples)

---

## What is SMILES?

SMILES is a string notation for describing molecular structures. For example:
- `C=C` represents ethylene (the monomer for polyethylene)
- `CC` represents ethane
- `C=C(C)C(=O)OC` represents methyl methacrylate (the monomer for PMMA)

AutoPoly automatically generates 3D structures and all necessary monomer variants from SMILES strings.

---

## Basic SMILES Syntax

### Atoms
- `C` - Carbon (aliphatic)
- `c` - Carbon (aromatic)
- `O` - Oxygen
- `N` - Nitrogen
- `F`, `Cl`, `Br`, `I` - Halogens

### Bonds
- `-` or omitted: Single bond (e.g., `CC` or `C-C`)
- `=`: Double bond (e.g., `C=C`)
- `#`: Triple bond (e.g., `C#N`)

### Branches
Use parentheses to indicate branching:
- `C=C(C)` means a carbon double-bonded to another carbon, which has a methyl group attached
- `C(C)(C)C` represents a carbon with three methyl groups attached (tert-butyl)

### Rings
- Ring numbers indicate closures: `C1CCCCC1` is cyclohexane
- Aromatic rings: `c1ccccc1` is benzene

---

## Common Polymer Monomers

### Polyethylene (PE)
**SMILES:** `C=C`
**Structure:** CH₂=CH₂
**Notes:** Simplest vinyl monomer

```python
polymer = Polymer(ChainNum=10, Sequence=["C=C"]*100, topology="linear")
```

### Polypropylene (PP)
**SMILES:** `C=C(C)`
**Structure:** CH₂=CH(CH₃)
**Notes:** Has a chiral center (use tacticity parameter)

```python
polymer = Polymer(ChainNum=10, Sequence=["C=C(C)"]*100, topology="linear", tacticity="atactic")
```

### Polystyrene (PS)
**SMILES:** `C=C(C1=CC=CC=C1)`
**Structure:** CH₂=CH(C₆H₅)
**Notes:** Contains phenyl ring

```python
polymer = Polymer(ChainNum=10, Sequence=["C=C(C1=CC=CC=C1)"]*50, topology="linear")
```

### Poly(methyl methacrylate) (PMMA)
**SMILES:** `C=C(C)C(=O)OC`
**Structure:** CH₂=C(CH₃)C(=O)OCH₃
**Notes:** Contains ester group; requires GAFF force field

```python
polymer = Polymer(ChainNum=10, Sequence=["C=C(C)C(=O)OC"]*50,
                  topology="linear", tacticity="atactic",
                  force_field="gaff")
```

### Poly(vinyl acetate) (PVAc)
**SMILES:** `C=C(OC(=O)C)`
**Structure:** CH₂=CH(OC(=O)CH₃)
**Notes:** Contains acetate group

```python
polymer = Polymer(ChainNum=10, Sequence=["C=C(OC(=O)C)"]*50,
                  topology="linear", force_field="gaff")
```

### Polycaprolactone (PCL)
**SMILES:** `C1CCCCC(=O)O`
**Structure:** Cyclic ester (lactone)
**Notes:** Can be used for ring-opening polymerization

---

## Heteroatom Polymers (NEW)

AutoPoly now supports polymers with heteroatom backbones (containing O, N, etc.) through SMARTS-based mechanism detection.

### Polylactic Acid (PLA) - Polyester
**SMILES:** `CC(C(=O)O)O`
**Structure:** CH₃CH(COOH)OH (Lactic acid)
**Mechanism:** Esterification (condensation)
**Force Field:** GAFF (recommended)
**Notes:** Biodegradable polyester with C-O backbone

```python
from AutoPy.monomer_generator import MonomerGenerator

generator = MonomerGenerator(
    base_name="PLA",
    output_dir="./monomers",
    mechanism='esterification',
    is_gaff=True
)
variants = generator.generate_variants(smiles="[*]CC(C(=O)O)O[*]")
files = generator.generate_lt_files(variants)
```

### Nylon-6,6 - Polyamide
**SMILES:** `NCC(=O)O` (amino acid monomer)
**Structure:** H₂N-CH₂-COOH (glycine as proxy)
**Mechanism:** Amidation (condensation)
**Force Field:** GAFF (recommended)
**Notes:** Polyamide with C-N backbone

```python
generator = MonomerGenerator(
    base_name="Nylon66",
    mechanism='amidation',
    is_gaff=True
)
variants = generator.generate_variants(smiles="[*]NCC(=O)O[*]")
files = generator.generate_lt_files(variants)
```

### Polyethylene Glycol (PEG) - Polyether
**SMILES:** `OCCO`
**Structure:** HO-CH₂-CH₂-OH (ethylene glycol)
**Mechanism:** Etherification (condensation)
**Force Field:** OPLS-AA
**Notes:** Polyether with C-O-C backbone

```python
generator = MonomerGenerator(
    base_name="PEG",
    output_dir="./monomers",
    mechanism='etherification',
    is_gaff=False  # OPLS-AA
)
variants = generator.generate_variants(smiles="[*]OCCO[*]")
files = generator.generate_lt_files(variants)
```

**Key Features:**
- **Element-agnostic**: Works with C, O, N backbones
- **Auto-detection**: Automatically identifies polymerization mechanism
- **SMARTS-based**: Uses chemical pattern matching instead of hard-coded rules
- **DOP-aware**: Respects Degree of Polymerization (see below)

---

## Single Molecule Generation (NEW)

AutoPoly preserves support for generating single molecules and systems of independent molecules (DOP=1).

### Single Ethanol Molecule
**SMILES:** `CCO`
**Use Case:** Solvent molecule, non-polymerizable
**Mechanism:** 'none' (auto-detected for DOP=1)

```python
generator = MonomerGenerator(
    base_name="ethanol",
    output_dir="./monomers",
    mechanism='none',  # Non-polymerizable
    is_gaff=False
)
variants = generator.generate_variants(smiles="[*]CCO[*]")
files = generator.generate_lt_files(variants)
# Generates: ethanol.lt (single molecule, no le/re/i variants)
```

### Water Molecules (100 independent molecules)
```python
from AutoPy import System, Polymer, Polymerization

system = System(out="water_box")

polymer = Polymer(
    ChainNum=100,    # 100 independent molecules
    Sequence=["O"],  # Water SMILES
    DOP=1            # Single molecules (NOT polymerized)
)

poly = Polymerization(
    name="water",
    system=system,
    model=[polymer],
    force_field="oplsaa"
)
# Result: molecule_1, molecule_2, ... molecule_100 (no bonds between them)
```

### Lactic Acid Molecule (with reactive groups)
```python
# Lactic acid has COOH and OH groups, but DOP=1 means keep them intact
generator = MonomerGenerator(
    base_name="lactic_acid",
    output_dir="./monomers",
    mechanism='none',  # Don't polymerize
    is_gaff=True
)
variants = generator.generate_variants(smiles="[*]CC(C(=O)O)O[*]")
files = generator.generate_lt_files(variants)
# Generates: lactic_acid.lt (plain, with all functional groups)
```

**DOP=1 Behavior:**
- Mechanism detection returns 'none' unless obvious vinyl pattern
- No le/re/i variant generation
- No connection point modifications
- Plain .lt file with all functional groups intact
- Systems of multiple independent molecules work correctly

---

## pSMILES (Required for All Monomers)

**pSMILES** (polymer SMILES) with wildcard atoms `[*]` is **required** for all monomer generation in AutoPoly. The wildcards explicitly mark the connection points where polymerization occurs.

**Syntax:** Add `[*]` wildcards at the two connection points in your monomer SMILES

### Why pSMILES is Required
- **Explicit connection points**: Forces users to specify exactly where polymerization occurs
- **Avoids ambiguity**: Prevents incorrect monomer generation from heuristic guessing
- **Works for all polymers**: Linear, ring, vinyl, condensation, single molecules

### How to Convert SMILES to pSMILES
1. Identify the two atoms where polymer bonds will form
2. Replace hydrogen atoms on those positions with `[*]` wildcards
3. The wildcards mark connection points (left and right)

### Examples

#### Polyethylene (Vinyl Polymer)
**SMILES:** `C=C`
**pSMILES:** `[*]C=C[*]`  ← Wildcards mark both ends of double bond

```python
from AutoPoly.monomer_generator import MonomerGenerator

generator = MonomerGenerator(base_name="PE", output_dir="./monomers", mechanism='vinyl_addition')
variants = generator.generate_variants(smiles="[*]C=C[*]")
files = generator.generate_lt_files(variants)
```

#### Polylactic Acid (Condensation Polymer)
**SMILES:** `CC(C(=O)O)O`
**pSMILES:** `[*]CC(C(=O)O)O[*]`  ← Wildcards mark alcohol and carboxyl ends

```python
generator = MonomerGenerator(
    base_name="PLA",
    output_dir="./monomers",
    mechanism='esterification',
    is_gaff=True
)
variants = generator.generate_variants(smiles="[*]CC(C(=O)O)O[*]")
files = generator.generate_lt_files(variants)
```

#### Nylon-6,6 (Polyamide)
**SMILES:** `NCC(=O)O`
**pSMILES:** `[*]NCC(=O)O[*]`  ← Wildcards mark amine and carboxyl ends

```python
generator = MonomerGenerator(
    base_name="Nylon66",
    output_dir="./monomers",
    mechanism='amidation',
    is_gaff=True
)
variants = generator.generate_variants(smiles="[*]NCC(=O)O[*]")
files = generator.generate_lt_files(variants)
```

#### Single Molecules (DOP=1)
Even for non-polymerizable molecules, pSMILES is required:

```python
generator = MonomerGenerator(
    base_name="ethanol",
    output_dir="./monomers",
    mechanism='none',
    is_gaff=False
)
variants = generator.generate_variants(smiles="[*]CCO[*]")
files = generator.generate_lt_files(variants)
```

**Important:** You must provide exactly 2 wildcard atoms (`[*]`) in your pSMILES. If you don't, AutoPoly will raise an error:
```
MonomerGeneratorError: pSMILES must contain exactly 2 wildcard atoms ([*] or *).
Got: C=C. Example correct pSMILES: '[*]C=C[*]'
```

---

## Finding SMILES for Custom Monomers

### Option 1: Online Databases
1. **PubChem** (https://pubchem.ncbi.nlm.nih.gov/)
   - Search for your compound
   - Find the SMILES in the "Chemical and Physical Properties" section

2. **ChemSpider** (http://www.chemspider.com/)
   - Search for your compound
   - Copy the SMILES string

### Option 2: Chemical Structure Tools
- **ChemDraw**: Draw structure → File → Export → SMILES
- **MarvinSketch**: Draw structure → Edit → Copy as → SMILES
- **BKChem**: Open source chemical structure editor

### Option 3: Command Line Tools
```bash
# Using Open Babel
obabel -:"YOUR_STRUCTURE" -osmi

# Example: Get SMILES for isopropanol
obabel -:"CC(O)C" -osmi
```

---

## Examples

### Basic Linear Polymer (Polyethylene)
```python
from AutoPoly import System, Polymer, Polymerization

system = System(out="pe_output")
polymer = Polymer(
    ChainNum=10,
    Sequence=["C=C"] * 100,  # 100 ethylene monomers
    topology="linear"
)

poly = Polymerization(
    name="polyethylene",
    system=system,
    model=[polymer],
    force_field="oplsaa"
)
```

### Functionalized Polymer (PMMA with Tacticity)
```python
system = System(out="pmma_output")
polymer = Polymer(
    ChainNum=5,
    Sequence=["C=C(C)C(=O)OC"] * 50,  # MMA SMILES
    topology="linear",
    tacticity="atactic",  # Random chirality
)

poly = Polymerization(
    name="pmma",
    system=system,
    model=[polymer],
    force_field="gaff"  # Use GAFF for functionalized monomers
)
```

### Ring Polymer
```python
system = System(out="ring_output")
polymer = Polymer(
    ChainNum=5,
    Sequence=["[*]C=C[*]"] * 30,  # pSMILES for ring closure
    topology="ring"
)

poly = Polymerization(
    name="ring_pe",
    system=system,
    model=[polymer],
    force_field="oplsaa"
)
```

### Copolymer (Two Different Monomers)
```python
system = System(out="copolymer_output")
polymer = Polymer(
    ChainNum=10,
    Sequence=["C=C", "C=C(C1=CC=CC=C1)"] * 25,  # Alternating PE and PS
    topology="linear"
)

poly = Polymerization(
    name="pe_ps_copolymer",
    system=system,
    model=[polymer],
    force_field="oplsaa"
)
```

---

## Tips and Best Practices

1. **Validate Your SMILES**
   - Use RDKit to check if your SMILES is valid:
   ```python
   from rdkit import Chem
   mol = Chem.MolFromSmiles("C=C")
   if mol:
       print("Valid SMILES")
   else:
       print("Invalid SMILES")
   ```

2. **Force Field Selection**
   - Simple hydrocarbons (PE, PP, PS): Use `oplsaa`
   - Functionalized monomers (PMMA, PVAc): Use `gaff`
   - Check compatibility with your force field

3. **Tacticity**
   - Monomers with chiral centers: Set `tacticity` parameter
   - `atactic`: Random (most common)
   - `isotactic`: All same configuration
   - `syndiotactic`: Alternating configuration

4. **SMILES Optimization**
   - Use canonical SMILES when possible
   - Avoid unnecessary stereochemistry (`@`, `@@`) unless needed
   - Test your SMILES with small systems first

---

## Troubleshooting

### Problem: "Invalid SMILES string"
**Solution:**
- Check for missing parentheses
- Verify all atoms are properly bonded
- Use an online SMILES validator

### Problem: "Monomer generation failed"
**Solution:**
- Ensure SMILES represents a valid vinyl monomer (C=C double bond)
- Check for unusual atoms or bonds not supported by force field
- Try simplifying the SMILES

### Problem: "Charges are 0.00"
**Solution:**
- For GAFF force field, calculate AM1-BCC charges
- See `example/04_advanced_features/Force_Field_Comparison/` for details
- For OPLS-AA, charges are included automatically

---

## Additional Resources

- **Open Babel**: https://openbabel.org/
- **RDKit**: https://www.rdkit.org/
- **SMILES Tutorial**: https://www.daylight.com/dayhtml_tutorials/languages/smiles/
- **PubChem**: https://pubchem.ncbi.nlm.nih.gov/

---

## Quick Reference Table

| Polymer | Common Name | SMILES | Force Field | Chiral | Mechanism |
|---------|-------------|---------|-------------|--------|-----------|
| PE | Polyethylene | `C=C` | OPLS-AA | No | vinyl_addition |
| PP | Polypropylene | `C=C(C)` | OPLS-AA | Yes | vinyl_addition |
| PS | Polystyrene | `C=C(C1=CC=CC=C1)` | OPLS-AA | No | vinyl_addition |
| PMMA | Poly(methyl methacrylate) | `C=C(C)C(=O)OC` | GAFF | Yes | vinyl_addition |
| PVAc | Poly(vinyl acetate) | `C=C(OC(=O)C)` | GAFF | No | vinyl_addition |
| PAN | Polyacrylonitrile | `C=C(C#N)` | GAFF | No | vinyl_addition |
| PVA | Poly(vinyl alcohol) | `C=C(CO)` | OPLS-AA | No | vinyl_addition |
| **PLA** | **Polylactic acid** | `CC(C(=O)O)O` | **GAFF** | **Yes** | **esterification** |
| **Nylon-6,6** | **Polyamide** | `NCC(=O)O` | **GAFF** | **No** | **amidation** |
| **PEG** | **Polyethylene glycol** | `OCCO` | **OPLS-AA** | **No** | **etherification** |

For more monomers, search online databases or use chemical drawing tools to obtain SMILES strings.

---
