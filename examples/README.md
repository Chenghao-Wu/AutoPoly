# AutoPoly Examples

Welcome to the AutoPoly examples directory! This collection contains tutorial examples for generating polymer structures with AutoPoly.

## Overview

AutoPoly is a Python package for generating atomistic and coarse-grained polymer structures for LAMMPS molecular dynamics simulations. These examples demonstrate various features and workflows.

## Available Examples

### 1. Linear Polylactic Acid (PLA) - `example_pla_linear.py`

**Difficulty:** Beginner
**Estimated Runtime:** 2-5 minutes
**Output Size:** ~5-10 MB

Demonstrates the complete workflow for generating linear PLA polymers using SMILES notation.

**Key Features:**
- SMILES-based monomer definition (`"CC(C(=O)O)O"` for lactic acid)
- Esterification mechanism (C-O backbone)
- Automatic monomer variant generation
- GAFF force field for polyesters
- Comprehensive explanations of each step

**What You'll Learn:**
- How to define polymers using SMILES strings
- Understanding esterification vs vinyl addition
- Why GAFF is recommended for polyesters
- How to navigate generated output files
- What steps are needed before running LAMMPS

**Run the Example:**
```bash
cd /home/zhenghaowu/mcp_lammps_poly/examples
python example_pla_linear.py
```

**Output Structure:**
```
pla_tutorial/
├── input/                    # Moltemplate files (reusable)
│   ├── monomer_0i.lt        # Internal monomer (middle of chain)
│   ├── monomer_0le.lt       # Left-end monomer (chain start)
│   ├── monomer_0re.lt       # Right-end monomer (chain end)
│   ├── monomer_0i_T1.lt     # Chirality variant
│   ├── monomer_0le_T1.lt    # Chirality variant
│   ├── monomer_0re_T1.lt    # Chirality variant
│   ├── poly_1.lt            # Polymer chain definition
│   └── gaff.lt              # GAFF force field
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

## Common Parameters Reference

### Polymer Topology
- `"linear"`: Linear chains with distinct start and end
- `"ring"`: Cyclic polymers (all monomers equivalent)

### Tacticity (Stereochemistry)
- `"atactic"`: Random stereochemistry (default)
- `"isotactic"`: Uniform configuration at all chiral centers
- `"syndiotactic"`: Alternating configuration

### Force Fields
- `"gaff"`: General Amber Force Field
  - **Best for:** Polyesters, biomolecules, drug-like molecules
  - **Advantages:** Well-validated ester parameters, extensive coverage
  - **Note:** Charges must be calculated (AM1-BCC or RESP)

- `"oplsaa"`: OPLS All-Atom Force Field
  - **Best for:** General organic molecules, liquid crystals
  - **Advantages:** Good thermodynamics properties
  - **Note:** May require validation for specific polymers

- `"lopls"`: Long-chain optimized OPLS
  - **Best for:** Alkanes, polyethylene
  - **Advantages:** Optimized for long hydrocarbon chains
  - **Note:** Less coverage for functional groups

## Common SMILES Examples

### Vinyl Polymers (Vinyl Addition Mechanism)
- **Ethylene (PE):** `"C=C"` or `"[*]C=C[*]"`
- **Propylene (PP):** `"C=C(C)"` or `"[*]C=C(C)[*]"`
- **Styrene (PS):** `"C=C(C1=CC=CC=C1)"` or `"[*]C=C(C1=CC=CC=C1)[*]"`
- **Methyl methacrylate (PMMA):** `"C=C(C)C(=O)OC"` or `"[*]C=C(C)C(=O)OC[*]"`

### Condensation Polymers
- **Lactic acid (PLA):** `"CC(C(=O)O)O"` - Esterification
- **Glycine (Nylon):** `"NCC(=O)O"` - Amidation

### Ring Polymers
- Add `[*]` wildcards for ring closure
- Example: `"[*]C=C[*]"` with `topology="ring"`

## Troubleshooting

### Issue: "No module named 'AutoPoly'"
**Solution:** Install AutoPoly:
```bash
cd /home/zhenghaowu/mcp_lammps_poly/AutoPoly
pip install -e .
```

### Issue: "Monomer not found"
**Solution:**
- Check SMILES string is valid
- Verify polymerization mechanism is supported
- Try running with `verbose=True` in Polymerization

### Issue: "GAFF charges are 0.00"
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

## Best Practices

### 1. Start Small
For testing and development, use small systems:
```python
polymer = Polymer(
    ChainNum=2,     # Start with 2 chains
    Sequence=["C=C"],
    DOP=10,         # Short chains for testing
    topology="linear"
)
```

### 2. Validate Before Production
Before running long production simulations:
- Check system density matches experimental values
- Verify energy conservation in NVE ensemble
- Test with shorter runs first
- Validate structural properties (Rg, end-to-end distance)

### 3. Force Field Selection
Choose force fields based on your polymer:
- **Polyesters (PLA, PGA):** Use GAFF
- **Hydrocarbons (PE, PP):** Use OPLS-AA or LOPLS
- **Polyamides (Nylon):** Use GAFF or OPLS-AA

### 4. Charge Calculation
GAFF requires manual charge calculation:
- **AM1-BCC:** Fast, reasonable accuracy (~0.1-0.2 kcal/mol error)
- **RESP:** More accurate (~0.01-0.05 kcal/mol error), but slower
- **Tools:** Antechamber, RESP, Gaussian, ORCA

### 5. Equilibration Protocol
Recommended workflow:
1. **Energy minimization:** 1000-10000 steps
2. **NVT heating:** Heat from 100 K to target temperature
3. **NPT compression:** Apply pressure to reach target density
4. **NPT production:** Run at target T and P

## Additional Resources

### Documentation
- **AutoPoly README:** `/home/zhenghaowu/mcp_lammps_poly/AutoPoly/README.md`
- **API Documentation:** `/home/zhenghaowu/mcp_lammps_poly/AutoPoly/docs/API.md`
- **SMILES Guide:** `/home/zhenghaowu/mcp_lammps_poly/AutoPoly/docs/SMILES_GUIDE.md`

### External Resources
- **LAMMPS Documentation:** https://docs.lammps.org/
- **Moltemplate:** https://moltemplate.org/
- **RDKit (SMILES):** https://www.rdkit.org/
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
4. Update this README with your example
5. Test thoroughly before submitting

## Support

For issues, questions, or suggestions:
- Open an issue on GitHub
- Check existing documentation
- Review troubleshooting section above

---

**Happy Simulating!** 🧪🔬
