# AutoPoly Examples

Welcome to the AutoPoly examples directory! This collection contains tutorial examples for generating polymer structures with AutoPoly.

## Overview

AutoPoly is a Python package for generating atomistic and coarse-grained polymer structures for LAMMPS molecular dynamics simulations. These examples demonstrate various features and workflows, including the pSMILES-based approach for defining monomers with connection points.

## Available Examples

### 1. Linear Polymethyl Methacrylate (PMMA) - `example_pmma_linear.py`

**Difficulty:** Beginner
**Estimated Runtime:** 2-5 minutes
**Output Size:** ~5-10 MB

Demonstrates the complete workflow for generating linear PMMA polymers using pSMILES notation.

**Key Features:**
- pSMILES-based monomer definition (`"[*]CC([*])(C)C(=O)OC"` for methyl methacrylate)
- Vinyl addition mechanism (C-C backbone)
- Automatic monomer variant generation via MonomerGenerator
- OPLS-AA force field for vinyl polymers
- Comprehensive explanations of each step

**What You'll Learn:**
- How to define polymers using pSMILES strings with [*] connection points
- Understanding vinyl addition vs condensation mechanisms
- Why OPLS-AA is recommended for vinyl polymers
- How AutoPoly automatically generates first/middle/last monomer variants
- What steps are needed before running LAMMPS

**Run the Example:**
```bash
cd /home/zhenghaowu/mcp_lammps_poly/examples
python example_pmma_linear.py
```

**Output Structure:**
```
pmma_tutorial/
├── input/                    # Moltemplate files (reusable)
│   ├── monomer_0_0le.lt     # Left-end monomer (chain start)
│   ├── monomer_0_1i.lt      # Internal monomer (middle of chain)
│   ├── monomer_0_2re.lt     # Right-end monomer (chain end)
│   ├── monomer_0_0le_T1.lt  # Left-end variant (mirror tacticity)
│   ├── monomer_0_1i_T1.lt   # Internal variant (mirror tacticity)
│   ├── monomer_0_2re_T1.lt  # Right-end variant (mirror tacticity)
│   ├── poly_1.lt            # Polymer chain definition
│   ├── oplsaa.lt            # OPLS-AA force field import
│   └── oplsaa.lt.prm        # OPLS-AA parameters
│
└── output/                   # LAMMPS input files
    ├── system.data          # Atom positions, bonds, angles, etc.
    ├── system.in            # LAMMPS input script
    ├── system.in.settings   # Force field parameters
    └── system.in.charges    # Atomic charges (manual calc. required)
```

## Prerequisites

Before running any examples, ensure you have:

1. **Python 3.7 or higher**
   ```bash
   python --version
   ```

2. **AutoPoly Installed**
   ```bash
   cd /home/zhenghaowu/mcp_lammps_poly/AutoPoly
   pip install -e .
   ```

3. **Required Python Packages**
   - rdkit
   - numpy
   (Installed automatically with AutoPoly)

4. **LAMMPS** (for running simulations)
   - Download from: https://lammps.sandia.gov/
   - Or use conda: `conda install -c conda-forge lammps`

5. **Moltemplate** (included with AutoPoly)
   - AutoPoly handles Moltemplate integration automatically

## Understanding pSMILES Notation

AutoPoly uses **pSMILES** (polymer SMILES) notation, which extends standard SMILES with `[*]` wildcards to mark connection points for polymerization.

### What are [*] Wildcards?

The `[*]` symbols in pSMILES mark the atoms where monomers will connect to form polymer chains:

- **First monomer:** Left `[*]` is replaced with hydrogen (chain start)
- **Middle monomers:** Both `[*]` connect to adjacent monomers
- **Last monomer:** Right `[*]` is replaced with hydrogen (chain end)

AutoPoly **automatically generates all three variants** from a single pSMILES string!

### Example: PMMA pSMILES

```
pSMILES: "[*]CC([*])(C)C(=O)OC"
          ^   ^
          |   |
        Connection points (backbone carbons)
```

This single pSMILES string generates:
- `monomer_0_0le.lt` - Left-end variant (first `[*]` → H)
- `monomer_0_1i.lt` - Internal variant (both `[*]` active)
- `monomer_0_2re.lt` - Right-end variant (second `[*]` → H)
- Plus T1 variants for stereochemistry control

## Polymerization Mechanisms

### Vinyl Addition (Chain-Growth)

**Mechanism:** C=C double bond opens to form C-C single bonds

**Characteristics:**
- All-carbon backbone
- No byproducts eliminated
- Side groups remain as pendant groups
- Common examples: PE, PP, PS, PMMA

**PMMA Example:**
```
Monomer: CH₂=C(CH₃)COOCH₃
Reaction: n CH₂=C(CH₃)COOCH₃ → [-CH₂-C(CH₃)(COOCH₃)-]n
Backbone: ...-CH₂-C(CH₃)(COOCH₃)-CH₂-C(CH₃)(COOCH₃)-...
```

### Condensation (Step-Growth)

**Mechanism:** Functional groups react with elimination of small molecules

**Characteristics:**
- Heteroatom backbone (C-O, C-N, etc.)
- Byproducts eliminated (H₂O, etc.)
- Examples: Polyesters, polyamides

**PLA Example (Esterification):**
```
Monomer: Lactic acid
Reaction: n HO-CH(CH₃)-COOH → [-O-CH(CH₃)-CO-]n + n H₂O
Backbone: ...-O-CH(CH₃)-CO-O-CH(CH₃)-CO-...
```

## Common pSMILES Examples

### Vinyl Polymers (Vinyl Addition)
- **Ethylene (PE):** `"[*]CC[*]"`
- **Propylene (PP):** `"[*]CC(C)[*]"`
- **Styrene (PS):** `"[*]CC(c1ccccc1)[*]"`
- **Methyl methacrylate (PMMA):** `"[*]CC([*])(C)C(=O)OC"`

### Condensation Polymers
- **Lactic acid (PLA):** `"[*]OC(C)C(=O)[*]"` - Esterification
- **Glycine (Nylon):** `"[*]NCC(=O)[*]"` - Amidation

### Ring Polymers
- Use `topology="ring"` parameter
- Example: `"[*]CC[*]"` with `topology="ring"` creates cyclic polyethylene

## Common Parameters Reference

### Polymer Topology
- `"linear"`: Linear chains with distinct start and end
- `"ring"`: Cyclic polymers (all monomers equivalent)

### Tacticity (Stereochemistry)
- `"atactic"`: Random stereochemistry (default, most common)
- `"isotactic"`: Uniform configuration at all chiral centers
- `"syndiotactic"`: Alternating configuration

### Force Fields

#### OPLS-AA (Optimized Potentials for Liquid Simulations - All Atom)
- **Best for:** Vinyl polymers, general organic molecules
- **Examples:** PE, PP, PS, PMMA
- **Advantages:**
  - Well-parameterized for hydrocarbons
  - Good thermodynamic properties
  - Excellent for vinyl polymers with ester side groups
- **Note:** Charges must be calculated (AM1-BCC or RESP)

#### GAFF (General Amber Force Field)
- **Best for:** Polyesters, biomolecules, drug-like molecules
- **Examples:** PLA, PGA, Nylon
- **Advantages:**
  - Well-validated ester parameters
  - Extensive coverage of functional groups
  - Good for condensation polymers
- **Note:** Charges must be calculated (AM1-BCC or RESP)

#### LOPLS (Long-chain optimized OPLS)
- **Best for:** Alkanes, long hydrocarbon chains
- **Examples:** Polyethylene, polypropylene
- **Advantages:**
  - Optimized for long hydrocarbon chains
  - Improved transferability
- **Note:** Less coverage for functional groups

## Troubleshooting

### Issue: "No module named 'AutoPoly'"
**Solution:** Install AutoPoly:
```bash
cd /home/zhenghaowu/mcp_lammps_poly/AutoPoly
pip install -e .
```

### Issue: "Monomer not found"
**Solution:**
- Check pSMILES string is valid (include [*] connection points)
- Verify polymerization mechanism is supported
- Try running with `verbose=True` in Polymerization

### Issue: "OPLS-AA/GAFF charges are 0.00"
**Solution:** This is expected! You need to calculate charges manually:
1. Use Antechamber: `antechamber -c bcc -m molecule.mol2 -cf charges.mol2`
2. Or use RESP with quantum chemistry calculations
3. Update the `system.in.charges` file

### Issue: "Moltemplate execution failed"
**Solution:**
- Check .lt file syntax in the `input/` directory
- Verify Moltemplate is installed: `which moltemplate`
- Check error messages in terminal output

### Issue: Simulation explodes in LAMMPS
**Solution:**
- Ensure proper equilibration (minimize → NVT → NPT)
- Check timestep (typically 0.5-1.0 fs for atomistic)
- Verify charges are calculated and assigned
- Check for bad contacts or overlaps
- Review force field parameters for your polymer

## Best Practices

### 1. Start Small
For testing and development, use small systems:
```python
polymer = Polymer(
    ChainNum=2,              # Start with 2 chains
    Sequence=["[*]CC[*]"],   # Simple pSMILES
    DOP=10,                  # Short chains for testing
    topology="linear"
)
```

### 2. Use pSMILES with [*] Wildcards
Recommended approach for all polymers:
- **Easier:** Define once, AutoPoly generates all variants
- **More flexible:** Works for linear, ring, and branched topologies
- **Less error-prone:** No manual variant creation needed

### 3. Validate Before Production
Before running long production simulations:
- Check system density matches experimental values
  - PMMA: ~1.18 g/cm³ (amorphous)
  - PLA: ~1.24-1.26 g/cm³
- Verify energy conservation in NVE ensemble
- Test with shorter runs first
- Validate structural properties (Rg, end-to-end distance)

### 4. Force Field Selection
Choose force fields based on your polymer:
- **Vinyl polymers (PE, PP, PS, PMMA):** Use OPLS-AA
- **Polyesters (PLA, PGA):** Use GAFF
- **Polyamides (Nylon):** Use GAFF or OPLS-AA
- **Hydrocarbons (long chains):** Use LOPLS

### 5. Charge Calculation
Both OPLS-AA and GAFF require manual charge calculation:
- **AM1-BCC:** Fast, reasonable accuracy (~0.1-0.2 kcal/mol error)
- **RESP:** More accurate (~0.01-0.05 kcal/mol error), but slower
- **Tools:** Antechamber, RESP, Gaussian, ORCA

### 6. Equilibration Protocol
Recommended workflow:
1. **Energy minimization:** 1000-10000 steps
2. **NVT heating:** Heat from 100 K to target temperature
3. **NPT compression:** Apply pressure to reach target density
4. **NPT production:** Run at target T and P

## Polymer-Specific Considerations

### PMMA (Polymethyl Methacrylate)
- **Glass transition:** ~378 K (105 °C)
- **Density:** ~1.18-1.20 g/cm³ (amorphous)
- **Tacticity:** Atactic (most common commercial grade)
- **Properties:**
  - High optical clarity
  - Good weather resistance
  - Moderate thermal stability
- **Applications:** Optical lenses, displays, coatings, bone cement
- **Validation checklist:**
  - [ ] Calculate partial charges (AM1-BCC or RESP)
  - [ ] Update system.in.charges file
  - [ ] Check density (~1.18 g/cm³ for amorphous)
  - [ ] Verify Tg (~378 K)
  - [ ] Test optical properties (refractive index ~1.49)

### PLA (Polylactic Acid)
- **Glass transition:** ~330 K (57 °C)
- **Density:** ~1.24-1.26 g/cm³
- **Biodegradable:** Yes
- **Applications:** Biomedical implants, packaging, 3D printing

### PS (Polystyrene)
- **Glass transition:** ~373 K (100 °C)
- **Density:** ~1.04-1.06 g/cm³
- **Properties:** Rigid, transparent, good electrical insulation

## Additional Resources

### Documentation
- **AutoPoly README:** `/home/zhenghaowu/mcp_lammps_poly/AutoPoly/README.md`
- **API Documentation:** `/home/zhenghaowu/mcp_lammps_poly/AutoPoly/docs/API.md`
- **SMILES Guide:** `/home/zhenghaowu/mcp_lammps_poly/AutoPoly/docs/SMILES_GUIDE.md`

### External Resources
- **LAMMPS Documentation:** https://docs.lammps.org/
- **Moltemplate:** https://moltemplate.org/
- **RDKit (SMILES):** https://www.rdkit.org/
- **OPLS-AA Force Field:** https://doi.org/10.1021/ja9621760
- **GAFF Force Field:** https://ambermd.org/antechamber/gaff.php
- **Antechamber:** https://ambermd.org/antechamber/antechamber.html

## Citation

If you use AutoPoly in your research, please cite:

```bibtex
@software{autopoly2024,
  title={AutoPoly: Automated Polymer Generation for LAMMPS},
  author={Wu, Zhenghao},
  year={2024},
  url={https://github.com/your-repo/autopoly}
}
```

## Contributing

We welcome contributions! If you'd like to add an example:

1. Create a new script in this directory
2. Follow the naming convention: `example_<description>.py`
3. Include comprehensive comments and explanations
4. Use pSMILES notation with [*] wildcards
5. Update this README with your example
6. Test thoroughly before submitting

## Support

For issues, questions, or suggestions:
- Open an issue on GitHub
- Check existing documentation
- Review troubleshooting section above

---

**Happy Simulating!** 🔬
