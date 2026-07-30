# Complement SMILES

Complement SMILES is AutoPoly's format for defining polymer sequences with precise positional control. It extends standard SMILES notation with **wildcard atoms** (`[*]`) that mark where each monomer connects to its neighbors.

## Why Complement SMILES?

In a polymer chain, monomers behave differently depending on their position:

1. **First monomer** — has a starting end group and one connection point
2. **Middle monomers** — have two connection points (to both neighbors)
3. **Last monomer** — has one connection point and an ending end group

Standard SMILES cannot distinguish these positions. Complement SMILES solves this with one rule:

**Number of wildcards = number of connections**

| Position | Wildcards | Example (ethylene) | Role |
|----------|-----------|--------------------|------|
| **First** | 1 (right) | `"CC[*]"` | Chain start |
| **Middle** | 2 (both sides) | `"[*]CC[*]"` | Interior units |
| **Last** | 1 (left) | `"[*]CC"` | Chain end |

## The Three Position Types

### First position (chain start)

```
CC[*]
│ │
│ └─ Wildcard: connects to the next monomer
└─── Free end group (methyl here)
```

### Middle positions

```
[*]CC[*]
 │  │  │
 │  │  └─ Wildcard: connects to the next monomer
 │  └──── Carbon backbone
 └─────── Wildcard: connects to the previous monomer
```

### Last position (chain end)

```
[*]CC
 │  │
 │  └─ Free end group (methyl here)
 └──── Wildcard: connects to the previous monomer
```

## Building Sequences

A chain is an explicit Python list — one complement SMILES per monomer, in order. The degree of polymerization is always `len(sequence)`.

### Step 1: Define the middle variant

```python
middle = "[*]CC[*]"   # ethylene, interior
```

### Step 2: Derive first and last

Remove the appropriate wildcard:

```python
first = "CC[*]"       # left wildcard removed
last  = "[*]CC"       # right wildcard removed
```

### Step 3: Build the sequence

```python
dop = 100
sequence = [first] + [middle] * (dop - 2) + [last]
```

### Step 4: Create the polymer

```python
from AutoPoly import Polymer

polymer = Polymer(
    chain_num=10,
    sequence=sequence,
    topology="linear",
    tacticity="atactic",
)
print(polymer.dop)   # 100
```

A reusable helper:

```python
def uniform_polymer(first: str, middle: str, last: str, dop: int) -> list:
    """Create a uniform homopolymer sequence."""
    if dop == 1:
        return [middle]
    if dop == 2:
        return [first, last]
    return [first] + [middle] * (dop - 2) + [last]

sequence = uniform_polymer("CC[*]", "[*]CC[*]", "[*]CC", 100)
```

## Common Patterns

### Uniform homopolymer

```python
# Polyethylene, DOP=50
sequence = ["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"]
```

### Block copolymer

Contiguous blocks — just mix monomers in the list:

```python
# PE(10)-PS(20)-PE(10) ABA triblock, DOP=40
sequence = (
    ["CC[*]"] + ["[*]CC[*]"] * 9 +                  # PE block (10)
    ["[*]CC([*])c1ccccc1"] * 20 +                   # PS block (20)
    ["[*]CC[*]"] * 9 + ["[*]CC"]                    # PE block (10)
)
```

### Alternating copolymer

```python
# Ethylene-styrene alternating, DOP=6
sequence = [
    "CC[*]",                      # A (first)
    "[*]CC([*])c1ccccc1",         # B
    "[*]CC[*]",                   # A
    "[*]CC([*])c1ccccc1",         # B
    "[*]CC[*]",                   # A
    "[*]CC(c1ccccc1)",            # B (last)
]
```

!!! tip "Comment complex sequences"
    For long block sequences, annotate block boundaries with inline comments — the sequence is data, and future-you will read it.

## Ring vs Linear Polymers

### Linear polymers

Need all three variants (first, middle, last):

```python
sequence_linear = ["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"]

polymer = Polymer(chain_num=10, sequence=sequence_linear, topology="linear")
```

### Ring polymers

Use **only the middle variant** — every position has two connections:

```python
sequence_ring = ["[*]CC[*]"] * 50    # all middle; no first or last

polymer = Polymer(chain_num=10, sequence=sequence_ring, topology="ring")
```

**Why?** Ring polymers have no end groups — every monomer connects to two neighbors. Using a first or last variant in a ring will fail.

## Common Monomers Reference

### Commodity polymers

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

### Engineering polymers

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

### Specialty polymers

**Poly(tetrafluoroethylene) (PTFE, Teflon)**

```python
first  = "C(F)(F)C(F)(F)[*]"
middle = "[*]C(F)(F)C([*])(F)F"
last   = "[*]C(F)(F)C(F)(F)"
```

**Polyisoprene (natural rubber)**

```python
first  = "CC(=C)C[*]"
middle = "[*]CC([*])=CC"
last   = "[*]CC(=C)C"
```

**Poly(lactic acid) (PLA)** — condensation polymer with an ester backbone

```python
first  = "OC(C)C(=O)[*]"
middle = "[*]OC(C)C(=O)[*]"
last   = "[*]OC(C)C(=O)O"
```

## Best Practices

### 1. Always check the wildcard count

✓ Correct:

```python
first  = "CC[*]"        # 1 wildcard
middle = "[*]CC[*]"     # 2 wildcards
last   = "[*]CC"        # 1 wildcard
```

❌ Incorrect:

```python
first  = "[*]CC[*]"     # 2 wildcards — that's a middle monomer
middle = "CC[*]"        # 1 wildcard — that's a first monomer
last   = "CC"           # 0 wildcards — no connection point at all
```

### 2. Verify the sequence length

```python
sequence = [...]
expected_dop = 50
assert len(sequence) == expected_dop, f"Expected DOP {expected_dop}, got {len(sequence)}"
```

### 3. Validate SMILES before a big run

```python
from rdkit import Chem

def valid(smiles: str) -> bool:
    return Chem.MolFromSmiles(smiles.replace("[*]", "C")) is not None

assert all(valid(s) for s in sequence)
```

## Advanced Topics

### Custom end groups

Change the non-wildcard end of the first/last monomers:

```python
first_OH = "OCC[*]"     # hydroxyl start
last_OH  = "[*]CCO"     # hydroxyl end

first_Br = "BrCC[*]"    # bromine start
last_Br  = "[*]CCBr"    # bromine end
```

### Explicit stereochemistry

Use SMILES chirality markers for tacticity control at the monomer level:

```python
isotactic_middle     = "[*]C[C@@H]([*])C"   # R-configuration
syndiotactic_middle  = "[*]C[C@H]([*])C"    # S-configuration
```

For whole-chain tacticity, prefer the `tacticity` parameter of `Polymer` (`"atactic"`, `"isotactic"`, `"syndiotactic"`).

## Appendix: SMILES Basics

Complement SMILES builds on standard SMILES. The essentials:

- **Atoms** — `C` aliphatic carbon, `c` aromatic carbon, `O`, `N`, halogens `F Cl Br I`
- **Bonds** — single (omitted or `-`), double `=`, triple `#` (e.g. `C=C`, `C#N`)
- **Branches** — parentheses: `CC(C)C` is isobutane
- **Rings** — matching digits close a ring: `C1CCCCC1` cyclohexane; aromatic rings use lowercase: `c1ccccc1` benzene
- **Wildcards** — `[*]` (brackets required; a bare `*` is invalid)

To find a SMILES for a new monomer, look it up in [PubChem](https://pubchem.ncbi.nlm.nih.gov/) or draw it in a chemical editor and export. Then replace the two connection-point hydrogens with `[*]` wildcards, one per side for a middle monomer.

## Next Steps

- [Polymers vs Molecules](polymers-vs-molecules.md) — when *not* to use wildcards
- [Force Fields](force-fields.md) — pick parameters for your chemistry
- [Tutorials](../tutorials/index.md) — sequences in action
- [Polymer API](../reference/polymer.md) — full class reference
