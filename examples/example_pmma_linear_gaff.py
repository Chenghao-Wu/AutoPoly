#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Linear Polymethyl Methacrylate (PMMA) Polymer Tutorial - GAFF Force Field

This script demonstrates how to generate atomistic linear Polymethyl Methacrylate (PMMA)
polymer structures for LAMMPS simulations using AutoPoly's pSMILES-based approach
with the GAFF (General AMBER Force Field).

Key Concepts Demonstrated:
- pSMILES notation for monomer definition (using [*] wildcards for connection points)
- Vinyl addition mechanism for PMMA (C-C backbone)
- Automatic monomer variant generation via MonomerGenerator
- GAFF force field for organic polymers
- LAMMPS data file generation via Moltemplate

Author: AutoPoly Development Team
"""

import sys
import os
from pathlib import Path

# ============================================================================
# SECTION 1: Configuration
# ============================================================================

"""
POLYMETHYL METHACRYLATE (PMMA) CONFIGURATION - GAFF
====================================================

pSMILES: "[*]CC([*])(C)C(=O)OC"
  Chemical Name: Methyl methacrylate (MMA) repeat unit
  Functional Groups: Backbone carbons (connection points), Methyl (C), Ester (COOCH3)
  Connection Points: [*] marks where monomers connect in the polymer chain

Polymerization Mechanism: Vinyl Addition
  Reaction: C=C double bond opens to form C-C single bonds
  Backbone: All carbon backbone (C-C single bonds)
  Connection Atoms: Both backbone carbons marked with [*]

pSMILES NOTATION:
  - [*] wildcards indicate connection points for polymerization
  - For PMMA: [*]CC([*])(C)C(=O)OC
    - First [*]: connects to previous monomer (or H at chain start)
    - Second [*]: connects to next monomer (or H at chain end)
  - AutoPoly generates first/middle/last variants automatically

Force Field: GAFF (General AMBER Force Field)
  - Developed for organic molecules and drug-like compounds
  - Uses letter-based atom types (c3, hc, os, c, o, etc.)
  - Well-suited for polymers with diverse functional groups
  - Requires AM1-BCC or RESP charges for accurate electrostatics
  - Reference: Wang et al., J. Comput. Chem. 2004, 25, 1157-1174
"""

PMMA_SMILES = "[*]CC([*])(C)C(=O)OC"  # pSMILES with [*] connection points for backbone carbons
CHAIN_NUM = 10               # Number of polymer chains
DOP = 50                     # Degree of polymerization (monomers per chain)
TOPOLOGY = "linear"          # Linear chain topology
TACTICITY = "atactic"        # Random stereochemistry
FORCE_FIELD = "gaff"         # GAFF force field (General AMBER Force Field)

# ============================================================================
# SECTION 2: Imports and Setup
# ============================================================================

def setup_imports():
    """
    Add AutoPoly to Python path and import required modules.

    AutoPoly should be installed via: pip install -e /path/to/AutoPoly
    """
    # Add parent directory to path to import AutoPoly
    project_root = Path(__file__).parent.parent
    sys.path.insert(0, str(project_root / "AutoPoly"))

    try:
        from AutoPoly import System, Polymer, Polymerization
        return System, Polymer, Polymerization
    except ImportError as e:
        print(f"\nError: Could not import AutoPoly modules: {e}")
        print("\nPlease ensure AutoPoly is installed:")
        print("  cd /home/zhenghaowu/mcp_lammps_poly/AutoPoly")
        print("  pip install -e .")
        sys.exit(1)

System, Polymer, Polymerization = setup_imports()

# ============================================================================
# SECTION 3: Create System for Output Management
# ============================================================================

def create_system():
    """
    Step 1: Create System for Output Management

    The System class manages file paths and creates output directories.

    Output Structure:
        pmma_tutorial_gaff/
        ├── moltemplate/    # Intermediate Moltemplate files
        ├── input/          # Generated .lt monomer files
        └── output/         # Final LAMMPS input files
    """
    print("\n" + "=" * 70)
    print("STEP 1: Creating System for Output Management")
    print("=" * 70)

    try:
        # Create System object with GAFF-specific output directory
        system = System(out="pmma_tutorial_gaff")

        output_path = system.get_folder_path()
        print(f"\nOutput directory: {output_path}")
        print("System created successfully!")

        return system

    except Exception as e:
        print(f"\nError creating system: {e}")
        sys.exit(1)

# ============================================================================
# SECTION 4: Define Linear PMMA Polymer
# ============================================================================

def define_pmmapolymer():
    """
    Step 2: Define Linear PMMA Polymer

    The Polymer class defines polymer structure using pSMILES notation.

    Key Parameters:
    - ChainNum: Number of polymer chains
    - Sequence: List of monomer pSMILES strings (with [*] connection points)
    - DOP: Degree of polymerization (monomers per chain)
    - topology: Chain topology ("linear" or "ring")
    - tacticity: Stereochemistry arrangement

    For PMMA:
    - pSMILES "[*]CC([*])(C)C(=O)OC" represents MMA repeat unit
    - [*] marks backbone connection points
    - Atactic tacticity = random configuration (most common commercially)
    """
    print("\n" + "=" * 70)
    print("STEP 2: Defining Linear PMMA Polymer")
    print("=" * 70)

    print("\nPolymer Configuration:")
    print(f"  Monomer pSMILES: {PMMA_SMILES}")
    print(f"  Chemical Name: Methyl methacrylate (MMA)")
    print(f"  Number of Chains: {CHAIN_NUM}")
    print(f"  Degree of Polymerization: {DOP}")
    print(f"  Total Monomers: {CHAIN_NUM * DOP}")
    print(f"  Topology: {TOPOLOGY}")
    print(f"  Tacticity: {TACTICITY}")

    try:
        # Create Polymer object
        polymer = Polymer(
            ChainNum=CHAIN_NUM,
            Sequence=[PMMA_SMILES],  # MMA pSMILES with [*] connection points
            DOP=DOP,
            topology=TOPOLOGY,
            tacticity=TACTICITY
        )

        print("\nPolymer defined successfully!")
        print(f"  Sequence set generated: {len(polymer.sequenceSet)} chains")
        if len(polymer.sequenceSet) > 0:
            print(f"  Monomer variants per chain: {len(polymer.sequenceSet[0])}")

        return polymer

    except Exception as e:
        print(f"\nError defining polymer: {e}")
        sys.exit(1)

# ============================================================================
# SECTION 5: Explain GAFF Force Field
# ============================================================================

def explain_gaff_force_field():
    """
    Step 3: GAFF Force Field Overview

    This section explains the GAFF force field and its advantages
    for organic polymer simulations.
    """
    print("\n" + "=" * 70)
    print("STEP 3: Understanding GAFF Force Field")
    print("=" * 70)

    print("""
WHAT IS GAFF?
-------------
GAFF (General AMBER Force Field) is a force field designed for:
  - Organic molecules and drug-like compounds
  - Diverse functional groups (esters, ethers, amides, etc.)
  - Compatible with AMBER protein force fields
  - Widely used in pharmaceutical and materials science

GAFF ATOM TYPES:
----------------
GAFF uses letter-based atom type naming:
  - c3: sp3 carbon (tetrahedral)
  - c : sp2 carbonyl carbon (C=O)
  - hc: hydrogen on sp3 carbon
  - os: ether/ester oxygen (-O-)
  - o : carbonyl oxygen (C=O)

For PMMA, typical GAFF atom types:
  - Backbone CH2: c3 + hc
  - Quaternary C: c3
  - Methyl CH3: c3 + hc
  - Ester C=O: c + o
  - Ester O-C: os
  - Methoxy CH3: c3 + hc

GAFF VS OPLS-AA:
----------------
| Feature          | GAFF              | OPLS-AA           |
|------------------|-------------------|-------------------|
| Atom types       | Letter-based      | Numeric           |
| Origin           | AMBER (AmberTools)| Jorgensen group   |
| Focus            | Drug molecules    | Proteins/liquids  |
| Compatibility    | AMBER FF          | CHARMM-like       |
| Charges          | AM1-BCC/RESP      | OPLS charges      |

WHY USE GAFF FOR PMMA?
----------------------
1. Excellent coverage of ester functional groups
2. Well-parameterized for organic molecules
3. Compatible with common charge methods (AM1-BCC)
4. Widely validated for polymer simulations
5. Good transferability to similar compounds
    """)

# ============================================================================
# SECTION 6: Explain Monomer Generation with GAFF
# ============================================================================

def explain_monomer_generation():
    """
    Step 4: Automatic Monomer Generation with GAFF

    This section explains how AutoPoly automatically generates monomer
    variants with GAFF atom types from pSMILES strings.
    """
    print("\n" + "=" * 70)
    print("STEP 4: Automatic Monomer Generation (GAFF)")
    print("=" * 70)

    print("""
For PMMA (pSMILES: "[*]CC([*])(C)C(=O)OC"), AutoPoly will automatically:

1. PARSE pSMILES AND BUILD CHAIN:
   - [*] wildcards mark the backbone connection points
   - MonomerGenerator builds a short chain (3+ monomers)
   - GAFF atom types assigned based on chemical environment

2. GENERATE 6 MONOMER VARIANT FILES:
   monomer_0_0le.lt     - Left-end monomer (first, chain start)
   monomer_0_1i.lt      - Internal monomer (middle of chain)
   monomer_0_49re.lt    - Right-end monomer (last, chain end)
   monomer_0_0le_T1.lt  - Left-end variant (mirror tacticity)
   monomer_0_1i_T1.lt   - Internal variant (mirror tacticity)
   monomer_0_49re_T1.lt - Right-end variant (mirror tacticity)

3. ASSIGN GAFF ATOM TYPES:
   - Uses SMARTS patterns from gaff_lt.fdefn
   - Letter-based types: c3, hc, c, o, os
   - Partial charges set to 0.00 (AM1-BCC calculation required)

4. LT FILE FORMAT (GAFF):
   - Imports "gaff.lt" (GAFF force field definition)
   - Class inherits from GAFF
   - Atom types use @atom:c3, @atom:hc, etc.

EXAMPLE GAFF LT FILE:
---------------------
import "gaff.lt"
monomer_0_1i inherits GAFF {
  write("Data Atoms") {
    $atom:C1 $mol:... @atom:c3 0.00   0.000  0.000  0.000
    $atom:C2 $mol:... @atom:c3 0.00   1.540  0.000  0.000
    $atom:H3 $mol:... @atom:hc 0.00  -0.350  1.000  0.000
    ...
  }
}

NOTE: This process happens AUTOMATICALLY in Polymerization!
You don't need to manually generate monomers.
    """)

# ============================================================================
# SECTION 7: Run Polymerization
# ============================================================================

def run_polymerization(system, polymer):
    """
    Step 5: Run Polymerization to Generate LAMMPS Files

    The Polymerization class orchestrates the complete workflow:
    1. Generates monomer variants from SMILES with GAFF typing
    2. Creates polymer chains
    3. Generates GAFF force field parameters
    4. Runs Moltemplate to create LAMMPS data files

    Key Features:
    - Automatic SMILES to monomer conversion
    - GAFF force field integration
    - Moltemplate automation
    - Complete LAMMPS input generation
    """
    print("\n" + "=" * 70)
    print("STEP 5: Running Polymerization (GAFF)")
    print("=" * 70)

    print("\nStarting polymerization process...")
    print(f"  Force field: {FORCE_FIELD.upper()}")
    print(f"  Output directory: {system.get_folder_path()}")
    print("\nThis may take 2-5 minutes...")

    try:
        # Create Polymerization object
        # This automatically triggers the full workflow
        polymerization = Polymerization(
            name="pmma_linear_gaff",
            system=system,
            model=[polymer],
            force_field=FORCE_FIELD,
            run=True  # Automatically execute polymerization
        )

        print("\nPolymerization completed successfully!")
        return polymerization

    except FileNotFoundError as e:
        print(f"\nError: Required file not found: {e}")
        print("Solution: Ensure AutoPoly is properly installed")
        sys.exit(1)

    except ValueError as e:
        print(f"\nError: Invalid parameter: {e}")
        print("Solution: Check pSMILES string and force field choice")
        sys.exit(1)

    except RuntimeError as e:
        print(f"\nError: Moltemplate execution failed: {e}")
        print("Solution: Check .lt file syntax and Moltemplate installation")
        sys.exit(1)

    except Exception as e:
        print(f"\nUnexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

# ============================================================================
# SECTION 8: Explain Generated Output Files
# ============================================================================

def explain_outputs(system):
    """
    Step 6: Explain Generated Output Files (GAFF)

    AutoPoly generates a complete set of LAMMPS input files using GAFF.
    This section explains what files are created and what they contain.
    """
    print("\n" + "=" * 70)
    print("STEP 6: Generated Output Files (GAFF)")
    print("=" * 70)

    output_path = Path(system.get_folder_path()) / "pmma_linear_gaff"

    print("\nOutput directory structure:")
    print(f"""
{output_path}/
├── input/                   # Moltemplate input files (reusable)
│   ├── monomer_0_0le.lt     # Left-end monomer (inherits GAFF)
│   ├── monomer_0_1i.lt      # Internal monomer (inherits GAFF)
│   ├── monomer_0_49re.lt    # Right-end monomer (inherits GAFF)
│   ├── monomer_0_0le_T1.lt  # Left-end (mirror tacticity)
│   ├── monomer_0_1i_T1.lt   # Internal (mirror tacticity)
│   ├── monomer_0_49re_T1.lt # Right-end (mirror tacticity)
│   ├── poly_1.lt            # Polymer chain definitions
│   └── gaff.lt              # GAFF force field import
│
└── output/                  # LAMMPS input files (ready for simulation)
    ├── system.data          # Atom positions, topology
    ├── system.in            # LAMMPS input script
    ├── system.in.settings   # GAFF force field parameters
    └── system.in.charges    # Atomic charges (0.00 - needs AM1-BCC)
    """)

    # Check if files exist and display info
    if output_path.exists():
        input_dir = output_path / "input"
        output_dir = output_path / "output"

        print("\nGenerated files:")

        if input_dir.exists():
            lt_files = list(input_dir.glob("monomer_*.lt"))
            print(f"  Monomer .lt files: {len(lt_files)}")

            # Display file sizes
            for lt_file in sorted(lt_files)[:3]:  # Show first 3
                size = lt_file.stat().st_size / 1024  # KB
                print(f"    - {lt_file.name}: {size:.1f} KB")

            if len(lt_files) > 3:
                print(f"    ... and {len(lt_files) - 3} more files")

        if output_dir.exists():
            data_files = ["system.data", "system.in", "system.in.settings"]
            for f in data_files:
                fpath = output_dir / f
                if fpath.exists():
                    size = fpath.stat().st_size / 1024  # KB
                    print(f"  {f}: {size:.1f} KB")
    else:
        print("\nNote: Output directory not yet created.")
        print("Files will be generated when the script completes.")

# ============================================================================
# SECTION 9: Next Steps and Best Practices
# ============================================================================

def explain_next_steps():
    """
    Step 7: Next Steps for Running LAMMPS Simulations with GAFF

    This section provides guidance on what to do after generating
    the PMMA polymer structure with GAFF.
    """
    print("\n" + "=" * 70)
    print("STEP 7: Next Steps for LAMMPS Simulations (GAFF)")
    print("=" * 70)

    print("""
IMPORTANT NOTES FOR GAFF:
-------------------------
1. GAFF CHARGES:
   - All atomic charges are currently set to 0.00
   - You MUST calculate partial charges using:
     * AM1-BCC method (recommended for GAFF)
     * RESP method (more accurate, requires Gaussian)
   - Update pmma_tutorial_gaff/pmma_linear_gaff/output/system.in.charges

2. AM1-BCC CHARGE CALCULATION:
   Using Antechamber (from AmberTools):
   $ antechamber -i molecule.mol2 -fi mol2 -o molecule_bcc.mol2 \\
                 -fo mol2 -c bcc -s 2

   Or using Open Babel + Antechamber:
   $ obabel -ismi "CC(C)(C(=O)OC)C" -omol2 -O monomer.mol2 --gen3d
   $ antechamber -i monomer.mol2 -fi mol2 -o monomer_bcc.mol2 \\
                 -fo mol2 -c bcc -s 2

3. LAMMPS SIMULATION:
   - Review pmma_tutorial_gaff/pmma_linear_gaff/output/system.in
   - Modify simulation parameters as needed
   - Run: lmp -in pmma_tutorial_gaff/pmma_linear_gaff/output/system.in

4. RECOMMENDED EQUILIBRATION PROTOCOL:
   a) Energy minimization: minimize 1.0e-4 1000 10000
   b) NVT heating: 100 K to target temperature
   c) NPT compression: Apply pressure to reach target density
   d) Production run: 10-100 ns depending on properties

GAFF-SPECIFIC CONSIDERATIONS:
-----------------------------
- GAFF uses harmonic bond/angle potentials
- Dihedral parameters from AMBER parameterization
- 1-4 scaling: 1/2 for electrostatics, 1/2 for vdW
- Combine with AM1-BCC charges for best results
- Compatible with TIP3P water model

PMMA-SPECIFIC CONSIDERATIONS:
------------------------------
- Glass transition: ~378 K (105 °C)
- Density: ~1.18-1.20 g/cm³ (amorphous)
- Tacticity: Atactic (most common commercial grade)
- Optical clarity: Highly transparent
- Weather resistance: Good outdoor durability
- Applications: Optical lenses, displays, coatings

VALIDATION CHECKLIST:
--------------------
- [ ] Calculate AM1-BCC partial charges
- [ ] Update system.in.charges file
- [ ] Check density (~1.18 g/cm³ for amorphous PMMA)
- [ ] Verify bond lengths and angles
- [ ] Test energy conservation in NVE ensemble
- [ ] Check glass transition temperature (~378 K)
- [ ] Compare with OPLS-AA results if available
    """)

# ============================================================================
# SECTION 10: Main Execution Function
# ============================================================================

def main():
    """
    Execute Complete PMMA Polymer Generation Workflow (GAFF)

    This function demonstrates the complete workflow:
    1. Create System for output management
    2. Define PMMA polymer using pSMILES
    3. Explain GAFF force field
    4. Explain automatic monomer generation
    5. Run polymerization (automatic monomer generation)
    6. Explain output files
    7. Provide next steps guidance

    Total execution time: ~2-5 minutes for 10 chains of 50 monomers
    """
    print("\n" + "=" * 70)
    print("AUTOPOLY LINEAR PMMA TUTORIAL (GAFF FORCE FIELD)")
    print("=" * 70)
    print("\nGenerating linear Polymethyl Methacrylate (PMMA) polymer structure")
    print("for LAMMPS molecular dynamics simulations using GAFF.")
    print("\nConfiguration:")
    print(f"  Monomer: Methyl methacrylate (pSMILES: {PMMA_SMILES})")
    print(f"  Chains: {CHAIN_NUM}")
    print(f"  DOP: {DOP}")
    print(f"  Force Field: {FORCE_FIELD.upper()}")

    try:
        # Step 1: Create System
        system = create_system()

        # Step 2: Define Polymer
        polymer = define_pmmapolymer()

        # Step 3: Explain GAFF force field
        explain_gaff_force_field()

        # Step 4: Explain monomer generation
        explain_monomer_generation()

        # Step 5: Run polymerization
        polymerization = run_polymerization(system, polymer)

        # Step 6: Explain outputs
        explain_outputs(system)

        # Step 7: Next steps
        explain_next_steps()

        print("\n" + "=" * 70)
        print("TUTORIAL COMPLETED SUCCESSFULLY!")
        print("=" * 70)
        print("\nGenerated files are in: pmma_tutorial_gaff/")
        print("Thank you for using AutoPoly with GAFF!")

    except KeyboardInterrupt:
        print("\n\nTutorial interrupted by user.")
        sys.exit(0)
    except Exception as e:
        print(f"\n\nUnexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    main()

