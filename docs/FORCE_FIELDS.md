# Force Field Selection Guide

AutoPoly supports six force fields for polymer and molecule simulations. This guide helps you choose the right force field for your research.

## Quick Selection Table

| Force Field | Best For | Strengths | Limitations |
|------------|----------|-----------|-------------|
| **OPLS-AA** | General organic polymers, biomolecules | Well-tested, broad coverage, good for liquids | May lack parameters for exotic chemistries |
| **LOPLS** | Liquid-phase simulations, improved densities | Optimized liquid properties, better densities | Less extensive than OPLS-AA |
| **GAFF** | Small molecules, drug-like compounds | Comprehensive organic coverage, well-validated | Requires charge calculation (AM1-BCC/RESP) |
| **GAFF2** | Updated GAFF applications | Improved parameters over GAFF, better accuracy | Newer, less extensively tested |
| **DREIDING** | Generic systems, metals, inorganics | Universal coverage, handles diverse atom types | Generic parameters, less accurate for organics |
| **COMPASS** | Commercial polymers, condensed phases | Excellent for polymers, accurate densities | Commercial FF, may need licensing |

## Detailed Force Field Descriptions

### OPLS-AA (Optimized Potentials for Liquid Simulations - All Atom)

**String value:** `"oplsaa"`

**Recommended for:**
- General organic polymers (polyethylene, polypropylene, polystyrene)
- Biomolecules (proteins, nucleic acids)
- Organic liquids
- Systems where broad parameter coverage is needed

**Strengths:**
- Extensively parameterized and validated
- Good reproduction of liquid-phase properties
- Wide coverage of organic functional groups
- Well-tested for polymers and biomolecules
- Large user community and extensive literature

**Limitations:**
- May not have parameters for all exotic chemical groups
- Not optimized for gas-phase properties
- Some parameters may be outdated compared to newer force fields

**Typical Applications:**
```python
# Commodity polymers
polymerization = Polymerization(
    name="polyethylene",
    system=system,
    model=[pe_polymer],
    force_field="oplsaa"
)
```

**Citations:**
- Jorgensen, W. L., et al. *J. Am. Chem. Soc.* **1996**, 118, 11225-11236.

---

### LOPLS (Liquid Optimized OPLS)

**String value:** `"lopls"`

**Recommended for:**
- Liquid-phase simulations
- Systems where accurate densities are critical
- Polymer melts and solutions
- Systems requiring improved thermodynamic properties

**Strengths:**
- Optimized for liquid-phase properties
- Better reproduction of densities and enthalpies of vaporization
- Improved over OPLS-AA for condensed phases
- Good for polymer melt simulations

**Limitations:**
- Less extensively parameterized than OPLS-AA
- Fewer validation studies
- May not have parameters for all OPLS-AA functional groups

**Typical Applications:**
```python
# Polymer melts with accurate density requirements
polymerization = Polymerization(
    name="polymer_melt",
    system=system,
    model=[polymer],
    force_field="lopls"
)
```

**Citations:**
- Siu, S. W. I., et al. *J. Chem. Theory Comput.* **2012**, 8, 1459-1470.

---

### GAFF (General AMBER Force Field)

**String value:** `"gaff"`

**Recommended for:**
- Small organic molecules
- Drug-like compounds
- Solvents
- Systems with diverse functional groups
- Molecule-only systems

**Strengths:**
- Comprehensive coverage of organic molecules
- Well-validated for drug-like molecules
- Designed to handle diverse functional groups
- Good transferability between similar molecules
- Compatible with AMBER protein force fields

**Limitations:**
- **Requires explicit charge calculation** (AM1-BCC or RESP)
- More complex setup than OPLS-based force fields
- May need additional parameterization for polymers
- Not extensively tested for all polymer systems

**Charge Calculation Requirement:**
GAFF requires pre-calculated partial charges using AM1-BCC or RESP methods. AutoPoly generates topology but you must:
1. Calculate charges using Antechamber (AM1-BCC) or Gaussian (RESP)
2. Update the generated charge file

**Typical Applications:**
```python
# Small molecule solvents
water = Molecule(Count=100, Smiles="O", Name="water")
ethanol = Molecule(Count=20, Smiles="CCO", Name="ethanol")

polymerization = Polymerization(
    name="solvent_mixture",
    system=system,
    model=[water, ethanol],
    force_field="gaff"
)
```

**Citations:**
- Wang, J., et al. *J. Comput. Chem.* **2004**, 25, 1157-1174.

---

### GAFF2 (General AMBER Force Field 2)

**String value:** `"gaff2"`

**Recommended for:**
- Updated GAFF applications
- Small molecules requiring improved accuracy
- Systems where GAFF is appropriate but better parameters are desired

**Strengths:**
- Improved parameters over original GAFF
- Better accuracy for small molecules
- Updated torsional parameters
- Maintains GAFF philosophy with better performance

**Limitations:**
- **Requires explicit charge calculation** (same as GAFF)
- Newer force field with less extensive validation
- May have fewer literature examples than GAFF
- Not all tools support GAFF2 yet

**When to choose GAFF2 over GAFF:**
- New projects (no legacy compatibility needed)
- Systems where improved accuracy justifies using newer FF
- When using current AMBER tools

**Typical Applications:**
```python
# Modern small molecule simulations
polymerization = Polymerization(
    name="molecules",
    system=system,
    model=[molecule],
    force_field="gaff2"
)
```

**Citations:**
- Wang, J., et al. *J. Mol. Graphics Modell.* **2006**, 25, 247-260.

---

### DREIDING (Generic Force Field)

**String value:** `"dreiding"`

**Recommended for:**
- Exploratory simulations
- Systems with unusual atom types
- Metal-organic frameworks (MOFs)
- Inorganic materials
- Systems lacking specific force field parameters

**Strengths:**
- Universal coverage - handles most elements
- Simple, transferable parameters
- Good for initial structure generation
- Handles metals and inorganic atoms
- Easy to apply to novel systems

**Limitations:**
- Generic parameters less accurate than specific FFs
- Not optimized for organic polymers
- May give poor quantitative results
- Better for structure generation than production simulations
- Energetics may be less reliable

**When to use DREIDING:**
- Initial exploration of new systems
- When specific force fields unavailable
- Structure generation and equilibration
- Systems with diverse/unusual atom types

**Typical Applications:**
```python
# Initial structure generation for novel polymers
polymerization = Polymerization(
    name="novel_polymer",
    system=system,
    model=[polymer],
    force_field="dreiding"
)
```

**Citations:**
- Mayo, S. L., et al. *J. Phys. Chem.* **1990**, 94, 8897-8909.

---

### COMPASS (Condensed-phase Optimized Molecular Potentials for Atomistic Simulation Studies)

**String value:** `"compass"`

**Recommended for:**
- Commercial polymers
- Industrial polymer applications
- Condensed-phase simulations
- Systems requiring high accuracy for polymer properties

**Strengths:**
- Excellent for polymer simulations
- Very accurate densities and mechanical properties
- Extensive validation for commercial polymers
- Good for polyolefins, polyesters, polyamides
- Includes cross-terms for improved accuracy

**Limitations:**
- Originally commercial (licensing may apply)
- More complex functional form
- Computationally more expensive than simpler FFs
- May require specific LAMMPS compilation

**Typical Applications:**
```python
# Industrial polymer research
polymerization = Polymerization(
    name="commercial_polymer",
    system=system,
    model=[polymer],
    force_field="compass"
)
```

**Citations:**
- Sun, H. *J. Phys. Chem. B* **1998**, 102, 7338-7364.

---

## Selection Workflow

### Step 1: Identify Your System Type

**Pure Polymer Systems:**
- Commodity polymers (PE, PP, PS) → OPLS-AA or LOPLS
- Commercial polymers → COMPASS
- Novel polymers → OPLS-AA or DREIDING (exploratory)

**Small Molecule Systems:**
- Organic solvents → GAFF or GAFF2
- Simple molecules → OPLS-AA or GAFF
- Drug-like molecules → GAFF or GAFF2

**Mixed Systems (Polymer + Solvent):**
- Polymer in organic solvent → GAFF or GAFF2 (handles both)
- Polymer in water → GAFF (with TIP3P water) or OPLS-AA

**Special Cases:**
- MOFs, metal-containing → DREIDING
- Exploratory work → DREIDING
- High-accuracy polymer properties → COMPASS

### Step 2: Consider Your Priorities

**If accuracy is critical:**
- Polymers → COMPASS or LOPLS
- Small molecules → GAFF2 or GAFF

**If speed/simplicity is critical:**
- Use OPLS-AA (simple setup, fast)
- Avoid GAFF (requires charge calculation)

**If universality is critical:**
- Use DREIDING (handles everything)
- May sacrifice accuracy

### Step 3: Check Parameter Availability

Before committing to a force field:
1. Check if your specific monomers/molecules are parameterized
2. Review literature for similar systems
3. Test with small systems first

---

## Force Field Comparison Examples

### Example 1: Polyethylene Density

Different force fields give different densities for PE at 300K, 1 atm:

| Force Field | Density (g/cm³) | Experimental |
|------------|----------------|--------------|
| OPLS-AA | 0.85 | 0.92-0.96 |
| LOPLS | 0.92 | 0.92-0.96 |
| COMPASS | 0.94 | 0.92-0.96 |
| DREIDING | 0.80 | 0.92-0.96 |

**Conclusion:** For accurate PE density, use COMPASS or LOPLS.

### Example 2: Small Molecule Solvation

For drug molecules in water, GAFF/GAFF2 are industry standard:

```python
# Drug-like molecule in water
drug = Molecule(Count=1, Smiles="CC(=O)Oc1ccccc1C(=O)O", Name="aspirin")
water = Molecule(Count=1000, Smiles="O", Name="water")

polymerization = Polymerization(
    name="drug_solvation",
    system=system,
    model=[drug, water],
    force_field="gaff"  # Industry standard for this application
)
```

---

## Common Pitfalls

### Mixing Force Fields Incorrectly

❌ **Don't mix force fields without validation:**
```python
# BAD: Using OPLS-AA polymer with GAFF solvent without combination rules
```

✓ **Do use compatible force fields:**
```python
# GOOD: Use GAFF for both polymer and solvent
polymerization = Polymerization(
    name="compatible_system",
    system=system,
    model=[polymer, solvent],
    force_field="gaff"  # Same FF for all components
)
```

### Forgetting GAFF Charge Calculation

❌ **Don't skip charge calculation for GAFF:**
```python
# BAD: Using GAFF without calculating charges
force_field="gaff"  # Charges will be wrong!
```

✓ **Do calculate charges properly:**
```bash
# After AutoPoly generates topology, calculate charges:
antechamber -i molecule.mol2 -fi mol2 -o charged.mol2 -fo mol2 -c bcc
# Then update system.in.charges file
```

### Using DREIDING for Production

❌ **Don't use DREIDING for final quantitative results:**
```python
# BAD: Using DREIDING for accurate property prediction
force_field="dreiding"  # Generic parameters, poor accuracy
```

✓ **Do use DREIDING for exploration, then switch:**
```python
# GOOD: Explore with DREIDING, produce with specific FF
# Step 1: Explore structure
polymerization = Polymerization(..., force_field="dreiding")

# Step 2: Production with accurate FF
polymerization = Polymerization(..., force_field="compass")
```

---

## Validation and Testing

Before running production simulations, always validate your force field choice:

### 1. Literature Comparison
- Search for similar systems in literature
- Compare your results with published data
- Check if your force field is appropriate

### 2. Property Benchmarking
- Calculate basic properties (density, energy)
- Compare with experimental values
- If off by >10%, consider different FF

### 3. Small Test Systems
- Start with small systems (10-20 chains, DOP 50)
- Verify behavior before scaling up
- Check for anomalies or instabilities

---

## Additional Resources

- **OPLS-AA Parameters:** [GROMACS Force Field](http://www.gromacs.org)
- **GAFF Parameters:** [AmberTools](http://ambermd.org/AmberTools.php)
- **DREIDING Parameters:** [Materials Studio Documentation](https://www.3ds.com/products/biovia/materials-studio)
- **COMPASS:** [LAMMPS COMPASS](https://docs.lammps.org/Howto_bioffe.html)

---

## Summary Decision Tree

```
System Type?
│
├─ Polymer Only
│  ├─ Commodity (PE, PP, PS) → OPLS-AA or LOPLS
│  ├─ Commercial/Industrial → COMPASS
│  └─ Novel/Exploratory → DREIDING or OPLS-AA
│
├─ Small Molecules Only
│  ├─ Organic solvents → GAFF or GAFF2
│  └─ Simple molecules → OPLS-AA
│
├─ Mixed (Polymer + Solvent)
│  ├─ Organic solvent → GAFF/GAFF2
│  └─ Aqueous → OPLS-AA or GAFF
│
└─ Special Systems
   ├─ MOFs, metals → DREIDING
   └─ Inorganics → DREIDING
```

---

## Questions?

If you're unsure which force field to use:

1. Check the [Troubleshooting Guide](TROUBLESHOOTING.md)
2. Review the [API Documentation](API.md)
3. Examine [example files](../examples/) for similar systems
4. Search literature for your specific polymer/molecule
5. Start with OPLS-AA (safe default for organic systems)

---

**Note:** Force field choice significantly impacts simulation results. Always validate your choice with experimental data or high-level calculations before drawing conclusions.
