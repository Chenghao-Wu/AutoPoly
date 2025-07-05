# CML to LT File Converter with Tacticity Support

## Overview

`rdlt_avogadro.py` is a Python script that converts CML (Chemical Markup Language) files into LAMMPS Template (`.lt`) files for use in polymer simulations. It supports tacticity (atactic, isotactic, syndiotactic) and generates all required monomer variants for the AutoPoly polymerization system.

---

## Features
- **CML to LT conversion** for polymer monomers
- **Tacticity support**: atactic, isotactic, syndiotactic
- **Stereochemical variants**: normal and T1 (180° rotation)
- **OPLSAA and LOPLS** atom typing
- **AutoPoly compatibility**
- **Command-line interface**

---

## Quickstart

### 1. Install Requirements
- Python 3.7+
- RDKit
- NumPy

Install with pip (RDKit may require conda):
```bash
conda install -c conda-forge rdkit
pip install numpy
```

### 2. Prepare Your CML Files
Place your monomer CML files in the same directory:
- `PEi.cml` (internal)
- `PEl.cml` (left-end)
- `PEr.cml` (right-end)

### 3. Run the Converter
```bash
python rdlt_avogadro.py --internal PEi.cml --left PEl.cml --right PEr.cml \
    --base-name PE --output-dir output --tacticity atactic
```

---

## Command-Line Usage

```bash
python rdlt_avogadro.py \
    --internal PEi.cml \
    --left PEl.cml \
    --right PEr.cml \
    --base-name PE \
    --output-dir output \
    --tacticity {atactic,isotactic,syndiotactic} \
    [--lopls] [--opls-fdef opls.fdef] [--lopls-fdef lopls.fdef]
```

**Arguments:**
- `--internal`   Internal monomer CML file
- `--left`       Left-end monomer CML file
- `--right`      Right-end monomer CML file
- `--base-name`  Base name for output files (e.g. PE)
- `--output-dir` Output directory (default: current)
- `--tacticity`  Tacticity type: atactic, isotactic, syndiotactic
- `--lopls`      Use LOPLS atom typing (optional)
- `--opls-fdef`  Path to OPLS feature definition file (optional)
- `--lopls-fdef` Path to LOPLS feature definition file (optional)

---

## Output
Depending on tacticity, the script generates:
- Atactic: `PEi.lt`, `PEi_T1.lt`, `PEle.lt`, `PEle_T1.lt`, `PEre.lt`, `PEre_T1.lt`
- Isotactic: `PEi.lt`, `PEle.lt`, `PEre.lt`
- Syndiotactic: `PEi.lt`, `PEi_T1.lt`, `PEle.lt`, `PEle_T1.lt`, `PEre.lt`, `PEre_T1.lt`

All files are ready for use with AutoPoly and LAMMPS.

---

## Tutorial: Generating Tacticity-Controlled Monomers

### Step 1: Prepare Monomer CML Files
- Use Avogadro or another tool to create and export your monomer structures as `.cml` files.
- Name them according to their role: `*i.cml` (internal), `*l.cml` (left), `*r.cml` (right).

### Step 2: Run the Script for Atactic Polymer
```bash
python rdlt_avogadro.py --internal PEi.cml --left PEl.cml --right PEr.cml \
    --base-name PE --output-dir output --tacticity atactic
```
- This will generate both normal and T1 variants for each monomer type.

### Step 3: Run for Isotactic or Syndiotactic
```bash
# Isotactic (all same stereochemistry)
python rdlt_avogadro.py --internal PEi.cml --left PEl.cml --right PEr.cml \
    --base-name PE --output-dir output --tacticity isotactic

# Syndiotactic (alternating stereochemistry)
python rdlt_avogadro.py --internal PEi.cml --left PEl.cml --right PEr.cml \
    --base-name PE --output-dir output --tacticity syndiotactic
```

### Step 4: Use in AutoPoly
- The generated `.lt` files can be directly used in the AutoPoly system for polymerization.
- Select the appropriate tacticity for your simulation.

---

## Advanced: LOPLS Support
To use LOPLS atom typing, add the `--lopls` flag:
```bash
python rdlt_avogadro.py --internal PEi.cml --left PEl.cml --right PEr.cml \
    --base-name PE --output-dir output --tacticity atactic --lopls
```

---

## Troubleshooting
- **Missing dependencies?** Install RDKit and NumPy as shown above.
- **CML parsing errors?** Ensure your CML files are valid and contain 3D coordinates.
- **File not generated?** Check the script output for warnings about missing or invalid files.

---

## License
MIT License

---

## Citation
If you use this script in your research, please cite the AutoPoly project and this repository. 