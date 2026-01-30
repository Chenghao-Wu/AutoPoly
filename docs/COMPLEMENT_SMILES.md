# Complement SMILES Format Guide

Complement SMILES is AutoPoly's unique format for defining polymer sequences with precise positional control. This guide explains the concept, syntax, and best practices.

## Table of Contents

- [What is Complement SMILES?](#what-is-complement-smiles)
- [The Three Position Types](#the-three-position-types)
- [Visual Examples](#visual-examples)
- [Common Patterns](#common-patterns)
- [Building Sequences](#building-sequences)
- [Common Monomers Reference](#common-monomers-reference)
- [Best Practices](#best-practices)
- [Troubleshooting](#troubleshooting)

---

## What is Complement SMILES?

**Complement SMILES** extends standard SMILES notation to specify the exact position of each monomer in a polymer chain using **wildcard atoms** (`[*]`).

### Why Complement SMILES?

In polymer chemistry, monomers behave differently depending on their position:

1. **First monomer** - Has a starting end group and one connection point
2. **Middle monomers** - Have two connection points (connecting to neighbors)
3. **Last monomer** - Has one connection point and an ending end group

Standard SMILES cannot distinguish these positions. Complement SMILES solves this by using wildcards to mark connection points.

### The Core Principle

**Number of wildcards = Number of connections**

- First position: **1 wildcard** (right connection only)
- Middle positions: **2 wildcards** (left and right connections)
- Last position: **1 wildcard** (left connection only)

---

## The Three Position Types

### First Position (Chain Start)

**Wildcard pattern:** `CC[*]` (right connection only)

The first monomer starts the chain. It has:
- A free end group on the left (no wildcard)
- One connection point on the right (`[*]`)

**Example - Ethylene (first):**
```
CC[*]
│ │
│ └─ Wildcard: connects to next monomer
└─── Free end group (methyl in this case)
```

### Middle Positions

**Wildcard pattern:** `[*]CC[*]` (both connections)

Middle monomers connect two neighbors. They have:
- A connection point on the left (`[*]`)
- A connection point on the right (`[*]`)

**Example - Ethylene (middle):**
```
[*]CC[*]
 │  │  │
 │  │  └─ Wildcard: connects to next monomer
 │  └──── Carbon backbone
 └─────── Wildcard: connects to previous monomer
```

### Last Position (Chain End)

**Wildcard pattern:** `[*]CC` (left connection only)

The last monomer ends the chain. It has:
- One connection point on the left (`[*]`)
- A free end group on the right (no wildcard)

**Example - Ethylene (last):**
```
[*]CC
 │  │
 │  └─ Free end group (methyl in this case)
 └──── Wildcard: connects to previous monomer
```

---

## Visual Examples

### Example 1: Polyethylene (PE) Chain

**Monomer:** Ethylene (—CH₂—CH₂—)

**Sequence for DOP=5:**
```python
sequence = [
    "CC[*]",        # Position 0: First  (1 wildcard right)
    "[*]CC[*]",     # Position 1: Middle (2 wildcards)
    "[*]CC[*]",     # Position 2: Middle (2 wildcards)
    "[*]CC[*]",     # Position 3: Middle (2 wildcards)
    "[*]CC"         # Position 4: Last   (1 wildcard left)
]

poly = Polymer(chain_num=1, sequence=sequence)
```

**Visual representation:**
```
CH₃—CH₂—CH₂—CH₂—CH₂—CH₂—CH₂—CH₂—CH₂—CH₃
 │      │       │       │       │
First  Mid     Mid     Mid    Last
```

### Example 2: Block Copolymer (PE-PS-PE)

**ABA triblock:** 2 ethylene, 3 styrene, 2 ethylene

```python
sequence = [
    # Block A: Polyethylene (2 units)
    "CC[*]",                    # Position 0: PE (first)
    "[*]CC[*]",                 # Position 1: PE (middle)

    # Block B: Polystyrene (3 units)
    "[*]CC([*])c1ccccc1",       # Position 2: PS (middle)
    "[*]CC([*])c1ccccc1",       # Position 3: PS (middle)
    "[*]CC([*])c1ccccc1",       # Position 4: PS (middle)

    # Block A: Polyethylene (2 units)
    "[*]CC[*]",                 # Position 5: PE (middle)
    "[*]CC"                     # Position 6: PE (last)
]

poly = Polymer(chain_num=5, sequence=sequence)  # DOP = 7
```

**Visual representation:**
```
—CH₂—CH₂—│—CH₂—CH(Ph)—CH₂—CH(Ph)—CH₂—CH(Ph)—│—CH₂—CH₂—CH₂—
  Block A  │           Block B                │  Block A
```

---

## Common Patterns

### Pattern 1: Uniform Homopolymer

**All same monomer, vary only by position:**

```python
# Polyethylene with DOP=50
sequence = ["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"]

poly = Polymer(chain_num=10, sequence=sequence)
```

### Pattern 2: Block Copolymer

**Contiguous blocks:** AAA-BBB-AAA

```python
# PE(10)-PS(20)-PE(10) triblock
sequence = (
    ["CC[*]"] + ["[*]CC[*]"] * 9 +                  # PE block (10 units)
    ["[*]CC([*])c1ccccc1"] * 20 +                   # PS block (20 units)
    ["[*]CC[*]"] * 9 + ["[*]CC"]                    # PE block (10 units)
)

poly = Polymer(chain_num=5, sequence=sequence)  # DOP = 40
```

### Pattern 3: Alternating Copolymer

**Two monomers alternate:** A-B-A-B-A-B...

```python
# Ethylene-Styrene alternating (DOP=6)
sequence = [
    "CC[*]",                      # A (first)
    "[*]CC([*])c1ccccc1",         # B
    "[*]CC[*]",                   # A
    "[*]CC([*])c1ccccc1",         # B
    "[*]CC[*]",                   # A
    "[*]CC(c1ccccc1)"             # B (last)
]
```

---

## Building Sequences

### Step-by-Step Guide

**Step 1: Define your monomer**

Start with middle variant (most common):
```python
middle = "[*]CC[*]"  # Ethylene middle
```

**Step 2: Derive first and last variants**

Remove wildcards appropriately:
```python
first = "CC[*]"      # Remove left wildcard
last = "[*]CC"       # Remove right wildcard
```

**Step 3: Build sequence**

```python
dop = 100
sequence = [first] + [middle] * (dop - 2) + [last]
```

**Step 4: Create polymer**

```python
polymer = Polymer(
    chain_num=10,
    sequence=sequence,
    topology="linear",
    tacticity="atactic"
)
```

---

## Common Monomers Reference

### Commodity Polymers

**Polyethylene (PE)**
```python
first  = "CC[*]"
middle = "[*]CC[*]"
last   = "[*]CC"
```

**Polypropylene (PP)**
```python
first  = "CC(C)[*]"
middle = "[*]CC([*])C"
last   = "[*]CC(C)"
```

**Polystyrene (PS)**
```python
first  = "CC(c1ccccc1)[*]"
middle = "[*]CC([*])c1ccccc1"
last   = "[*]CC(c1ccccc1)"
```

**Poly(vinyl chloride) (PVC)**
```python
first  = "CC(Cl)[*]"
middle = "[*]CC([*])Cl"
last   = "[*]CC(Cl)"
```

### Engineering Polymers

**Poly(methyl methacrylate) (PMMA)**
```python
first  = "CC(C)(C(=O)OC)[*]"
middle = "[*]CC([*])(C)C(=O)OC"
last   = "[*]CC(C)(C(=O)OC)"
```

**Polyacrylonitrile (PAN)**
```python
first  = "CC(C#N)[*]"
middle = "[*]CC([*])C#N"
last   = "[*]CC(C#N)"
```

**Poly(ethylene oxide) (PEO)**
```python
first  = "COC[*]"
middle = "[*]COC[*]"
last   = "[*]COC"
```

### Specialty Polymers

**Poly(tetrafluoroethylene) (PTFE, Teflon)**
```python
first  = "C(F)(F)C(F)(F)[*]"
middle = "[*]C(F)(F)C([*])(F)F"
last   = "[*]C(F)(F)C(F)(F)"
```

**Polyisoprene (Natural Rubber)**
```python
first  = "CC(=C)C[*]"
middle = "[*]CC([*])=CC"
last   = "[*]CC(=C)C"
```

---

## Best Practices

### 1. Always Check Wildcard Count

✓ **Correct:**
```python
first = "CC[*]"          # 1 wildcard ✓
middle = "[*]CC[*]"      # 2 wildcards ✓
last = "[*]CC"           # 1 wildcard ✓
```

❌ **Incorrect:**
```python
first = "[*]CC[*]"       # 2 wildcards ✗ (should be 1)
middle = "CC[*]"         # 1 wildcard ✗ (should be 2)
last = "CC"              # 0 wildcards ✗ (should be 1)
```

### 2. Document Complex Sequences

Add comments for clarity:
```python
sequence = [
    "CC[*]",                      # Position 0: PE (first)
    "[*]CC[*]",                   # Position 1: PE
    "[*]CC[*]",                   # Position 2: PE (end of block A)

    "[*]CC([*])c1ccccc1",         # Position 3: PS (start of block B)
    "[*]CC([*])c1ccccc1",         # Position 4: PS
    "[*]CC([*])c1ccccc1",         # Position 5: PS (end of block B)

    "[*]CC[*]",                   # Position 6: PE (start of block A)
    "[*]CC"                       # Position 7: PE (last)
]
```

### 3. Verify Sequence Length

```python
# Check sequence length matches expected DOP
sequence = [...]  # Your sequence
expected_dop = 50

assert len(sequence) == expected_dop, f"Expected DOP {expected_dop}, got {len(sequence)}"

polymer = Polymer(chain_num=10, sequence=sequence)
print(f"Actual DOP: {polymer.dop}")  # Should match len(sequence)
```

### 4. Use Helper Functions

Create reusable functions:
```python
def uniform_polymer(first: str, middle: str, last: str, dop: int) -> list:
    """Create uniform polymer sequence."""
    if dop == 1:
        return [middle]  # Single unit
    elif dop == 2:
        return [first, last]
    else:
        return [first] + [middle] * (dop - 2) + [last]

# Usage
sequence = uniform_polymer("CC[*]", "[*]CC[*]", "[*]CC", 100)
```

---

## Troubleshooting

### Issue 1: Wrong Chain Length

**Problem:** Generated chain has unexpected DOP

**Solution:**
```python
# DOP = len(sequence) automatically
sequence = [...]
print(f"DOP will be: {len(sequence)}")

polymer = Polymer(chain_num=10, sequence=sequence)
print(f"Actual DOP: {polymer.dop}")
```

### Issue 2: ValidationError - Invalid SMILES

**Error:**
```
ValidationError: Invalid SMILES in sequence
```

**Causes & Solutions:**
1. **Wrong wildcard count** - Check first=1, middle=2, last=1
2. **SMILES syntax error** - Validate using RDKit or online tools
3. **Missing brackets** - Wildcards must be `[*]` not `*`

### Issue 3: Ring Polymer Issues

**Problem:** Ring polymers fail with first/last variants

**Solution:** For ring polymers, ALL positions use middle variant (2 wildcards):

```python
# Linear polymer - needs first/middle/last
linear_sequence = ["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"]

# Ring polymer - ALL middle variants
ring_sequence = ["[*]CC[*]"] * 50  # All positions have 2 wildcards

poly_ring = Polymer(
    chain_num=5,
    sequence=ring_sequence,
    topology="ring"  # Important!
)
```

### Issue 4: Block Boundaries Wrong

**Problem:** Block copolymer has incorrect transitions

**Solution:** All internal positions need 2 wildcards:

```python
# Correct block copolymer
sequence = (
    ["CC[*]"] +                              # First (1 wildcard)
    ["[*]CC[*]"] * 9 +                       # Block A middle
    ["[*]CC([*])c1ccccc1"] * 20 +            # Block B middle
    ["[*]CC[*]"] * 9 +                       # Block A middle
    ["[*]CC"]                                # Last (1 wildcard)
)
```

---

## Ring vs Linear Polymers

### Linear Polymers

**Need three variants** (first, middle, last):
```python
sequence_linear = [
    "CC[*]",        # First (1 wildcard)
    "[*]CC[*]",     # Middle (2 wildcards)
    "[*]CC[*]",     # Middle (2 wildcards)
    "[*]CC"         # Last (1 wildcard)
]

polymer_linear = Polymer(
    chain_num=10,
    sequence=sequence_linear,
    topology="linear"
)
```

### Ring Polymers

**Only use middle variant** (all 2 wildcards):
```python
sequence_ring = [
    "[*]CC[*]",     # Position 0 (2 wildcards)
    "[*]CC[*]",     # Position 1 (2 wildcards)
    "[*]CC[*]",     # Position 2 (2 wildcards)
    "[*]CC[*]"      # Position 3 (2 wildcards)
]  # No first or last!

polymer_ring = Polymer(
    chain_num=10,
    sequence=sequence_ring,
    topology="ring"
)
```

**Why?** Ring polymers have no end groups—every position connects to two neighbors.

---

## Advanced Topics

### Custom End Groups

Modify chain ends:
```python
# Methyl end groups (default)
first = "CC[*]"
last = "[*]CC"

# Hydroxyl end groups
first_OH = "OCC[*]"
last_OH = "[*]CCO"

# Bromine end groups
first_Br = "BrCC[*]"
last_Br = "[*]CCBr"
```

### Explicit Stereochemistry

Use SMILES stereochemistry notation:
```python
# R-configuration (isotactic)
isotactic_middle = "[*]C[C@@H]([*])C"

# S-configuration
syndiotactic_middle = "[*]C[C@H]([*])C"
```

---

## Summary

**Complement SMILES** = Positional SMILES with wildcards

| Position | Wildcards | Example (Ethylene) | Role |
|----------|-----------|-------------------|------|
| **First** | 1 (right) | `CC[*]` | Chain start |
| **Middle** | 2 (both) | `[*]CC[*]` | Internal units |
| **Last** | 1 (left) | `[*]CC` | Chain end |

**Key Rules:**
1. Linear polymers: Use first, middle, and last variants
2. Ring polymers: Use only middle variant (all positions)
3. Wildcard count must match connection points
4. DOP = `len(sequence)` automatically

---

## Next Steps

- Review [API Documentation](API.md) for Polymer class details
- Check [Examples](../examples/) for working code
- See [Force Field Guide](FORCE_FIELDS.md) for simulation setup
- Read [Migration Guide](../MIGRATION.md) if upgrading from v0.x

---

**Questions?** See the [Troubleshooting Guide](TROUBLESHOOTING.md) or open an issue on GitHub.
