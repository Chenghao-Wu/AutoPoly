#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Linear Polymethyl Methacrylate (PMMA) Polymer Tutorial

This script demonstrates how to generate atomistic linear Polymethyl Methacrylate (PMMA)
polymer structures for LAMMPS simulations using AutoPoly's pSMILES-based approach.

Key Concepts Demonstrated:
- pSMILES notation for monomer definition (using [*] wildcards for connection points)
- Vinyl addition mechanism for PMMA (C-C backbone)
- Automatic monomer variant generation via MonomerGenerator
- OPLS-AA force field for vinyl polymers
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
POLYMETHYL METHACRYLATE (PMMA) CONFIGURATION
=============================================

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

Force Field: OPLS-AA (Optimized Potentials for Liquid Simulations - All Atom)
  Well-parameterized for vinyl polymers with ester side groups
  Note: Charges should be calculated using AM1-BCC or RESP methods
         for accurate electrostatic interactions
"""

PMMA_SMILES = "[*]CC([*])(C)C(=O)OC"  # pSMILES with [*] connection points for backbone carbons
CHAIN_NUM = 10               # Number of polymer chains
DOP = 50                     # Degree of polymerization (monomers per chain)
# With new API, DOP is derived from sequence length
TOPOLOGY = "linear"          # Linear chain topology
TACTICITY = "atactic"        # Random stereochemistry
FORCE_FIELD = "oplsaa"       # OPLS-AA force field (recommended for vinyl polymers)
# or gaff force field
#FORCE_FIELD = "gaff"



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
        pmma_tutorial/
        ├── moltemplate/    # Intermediate Moltemplate files
        ├── input/          # Generated .lt monomer files
        └── output/         # Final LAMMPS input files
    """
    print("\n" + "=" * 70)
    print("STEP 1: Creating System for Output Management")
    print("=" * 70)

    try:
        # Create System object
        system = System(out="pmma_tutorial")

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
    print(f"\nNote: Using new explicit sequence API")
    print(f"  DOP is automatically derived from len(sequence)")

    try:
        # Create explicit sequence (uniform for this example)
        # With new API, we specify each monomer position explicitly
        sequence = [PMMA_SMILES] * DOP

        # Create Polymer object with Pythonic naming
        polymer = Polymer(
            chain_num=CHAIN_NUM,
            sequence=sequence,  # Explicit sequence
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
# SECTION 5: Explain Vinyl Addition Mechanism
# ============================================================================

def explain_vinyl_addition():
    """
    Step 3: Vinyl Addition Mechanism for PMMA

    This section explains how PMMA polymerization works and why it differs
    from condensation polymers like PLA.
    """
    print("\n" + "=" * 70)
    print("STEP 3: Understanding PMMA Vinyl Addition Polymerization")
    print("=" * 70)

    print("""
WHAT IS VINYL ADDITION POLYMERIZATION?
---------------------------------------
Vinyl addition is a chain-growth polymerization where:
  - A carbon-carbon double bond (C=C) opens
  - Forms single bonds to adjacent monomers
  - Creates polymer backbone with C-C single bonds
  - No byproducts are eliminated

PMMA POLYMERIZATION:
--------------------
Monomer: Methyl methacrylate (CH₂=C(CH₃)COOCH₃)

Reaction:
  n CH₂=C(CH₃)COOCH₃ -> [-CH₂-C(CH₃)(COOCH₃)-]n

Backbone Structure:
  ...-CH₂-C(CH₃)(COOCH₃)-CH₂-C(CH₃)(COOCH₃)-...
        ^    ^
        C    C
        (all carbon backbone)

Connection Atoms:
  - Both carbons from C=C double bond
  - Double bond opens to form two single bonds
  - No atoms are lost during polymerization

WHY THIS MATTERS:
-----------------
1. All-carbon backbone: Different from condensation polymers
2. No elimination: No water or other byproducts formed
3. Side groups: Ester groups remain as pendant groups
4. Properties: Affects thermal, mechanical, and optical properties

pSMILES FOR VINYL POLYMERS:
----------------------------
AutoPoly uses pSMILES notation with [*] wildcards:
  - pSMILES: "[*]CC([*])(C)C(=O)OC"
  - [*] marks the backbone connection points
  - First [*]: connects to previous monomer
  - Second [*]: connects to next monomer

AUTOPOLY VARIANT GENERATION:
-----------------------------
From a single pSMILES, AutoPoly generates:
  - First variant: left connection capped with H
  - Middle variant: both connections active
  - Last variant: right connection capped with H
  - T1 variants: mirror stereochemistry for tacticity
    """)

# ============================================================================
# SECTION 6: Explain Monomer Generation
# ============================================================================

def explain_monomer_generation():
    """
    Step 4: Automatic Monomer Generation

    This section explains how AutoPoly automatically generates monomer
    variants from pSMILES strings. This happens internally during Polymerization.
    """
    print("\n" + "=" * 70)
    print("STEP 4: Automatic Monomer Generation")
    print("=" * 70)

    print("""
For PMMA (pSMILES: "[*]CC([*])(C)C(=O)OC"), AutoPoly will automatically:

1. PARSE pSMILES AND BUILD CHAIN:
   - [*] wildcards mark the backbone connection points
   - MonomerGenerator builds a short chain (3+ monomers)
   - Atom types assigned on the chain for correct chemical environment

2. GENERATE 6 MONOMER VARIANT FILES:
   monomer_0_0le.lt     - Left-end monomer (first, chain start)
   monomer_0_1i.lt      - Internal monomer (middle of chain)
   monomer_0_2re.lt     - Right-end monomer (last, chain end)
   monomer_0_0le_T1.lt  - Left-end variant (mirror tacticity)
   monomer_0_1i_T1.lt   - Internal variant (mirror tacticity)
   monomer_0_2re_T1.lt  - Right-end variant (mirror tacticity)

3. ASSIGN FORCE FIELD ATOM TYPES:
   - Uses SMARTS patterns from .fdefn files
   - C, H, O atoms typed according to OPLS-AA/GAFF definitions
   - Partial charges set to 0.00 (manual calculation required)

4. GENERATE 3D CONFORMERS:
   - Generate 3D coordinates using ETKDG method
   - Align backbone along X-axis for chain assembly
   - Connection atoms placed FIRST in LT file for AutoPoly compatibility

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
    1. Generates monomer variants from SMILES
    2. Creates polymer chains
    3. Generates force field parameters
    4. Runs Moltemplate to create LAMMPS data files

    Key Features:
    - Automatic SMILES to monomer conversion
    - Force field integration (OPLS-AA/GAFF)
    - Moltemplate automation
    - Complete LAMMPS input generation
    """
    print("\n" + "=" * 70)
    print("STEP 5: Running Polymerization")
    print("=" * 70)

    print("\nStarting polymerization process...")
    print(f"  Force field: {FORCE_FIELD.upper()}")
    print(f"  Output directory: {system.get_folder_path()}")
    print("\nThis may take 2-5 minutes...")

    try:
        # Create Polymerization object
        # This automatically triggers the full workflow
        polymerization = Polymerization(
            name="pmma_linear",
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
    Step 6: Explain Generated Output Files

    AutoPoly generates a complete set of LAMMPS input files.
    This section explains what files are created and what they contain.
    """
    print("\n" + "=" * 70)
    print("STEP 6: Generated Output Files")
    print("=" * 70)

    output_path = Path(system.get_folder_path()) / "pmma_linear"

    print("\nOutput directory structure:")
    print(f"""
{output_path}/
├── input/                   # Moltemplate input files (reusable)
│   ├── monomer_0_0le.lt     # Left-end monomer (first, chain start)
│   ├── monomer_0_1i.lt      # Internal monomer (middle of chain)
│   ├── monomer_0_2re.lt     # Right-end monomer (last, chain end)
│   ├── monomer_0_0le_T1.lt  # Left-end (mirror tacticity)
│   ├── monomer_0_1i_T1.lt   # Internal (mirror tacticity)
│   ├── monomer_0_2re_T1.lt  # Right-end (mirror tacticity)
│   ├── poly_1.lt            # Polymer chain definitions
│   ├── oplsaa.lt            # OPLS-AA force field import
│   └── oplsaa.lt.prm        # OPLS-AA parameters
│
└── output/                  # LAMMPS input files (ready for simulation)
    ├── system.data          # Atom positions, topology
    ├── system.in            # LAMMPS input script
    ├── system.in.settings   # Force field parameters
    └── system.in.charges    # Atomic charges (0.00 for OPLS-AA)
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
    Step 7: Next Steps for Running LAMMPS Simulations

    This section provides guidance on what to do after generating
    the PMMA polymer structure.
    """
    print("\n" + "=" * 70)
    print("STEP 7: Next Steps for LAMMPS Simulations")
    print("=" * 70)

    print("""
IMPORTANT NOTES:
----------------
1. OPLS-AA CHARGES:
   - All atomic charges are currently set to 0.00
   - You MUST calculate partial charges using:
     * AM1-BCC method (faster, reasonable accuracy)
     * RESP method (more accurate, requires quantum calculations)
   - Update pmma_tutorial/pmma_linear/output/system.in.charges

2. CHARGE CALCULATION TOOLS:
   - Antechamber (from AmberTools): antechamber -c bcc -m molecule.mol2
   - Open Babel: obabel -ismi -h -o mol2
   - RESP: Gaussian quantum chemistry calculations

3. LAMMPS SIMULATION:
   - Review pmma_tutorial/pmma_linear/output/system.in
   - Modify simulation parameters as needed
   - Run: lmp -in pmma_tutorial/pmma_linear/output/system.in

4. RECOMMENDED EQUILIBRATION PROTOCOL:
   a) Energy minimization: minimize 1.0e-4 1000 10000
   b) NVT heating: 100 K to target temperature
   c) NPT compression: Apply pressure to reach target density
   d) Production run: 10-100 ns depending on properties

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
- [ ] Calculate partial charges (AM1-BCC or RESP)
- [ ] Update system.in.charges file
- [ ] Check density (~1.18 g/cm³ for amorphous PMMA)
- [ ] Verify bond lengths and angles
- [ ] Test energy conservation in NVE ensemble
- [ ] Check glass transition temperature (~378 K)
- [ ] Verify optical properties (refractive index ~1.49)
    """)

# ============================================================================
# SECTION 10: Main Execution Function
# ============================================================================

def main():
    """
    Execute Complete PMMA Polymer Generation Workflow

    This function demonstrates the complete workflow:
    1. Create System for output management
    2. Define PMMA polymer using pSMILES
    3. Explain vinyl addition mechanism
    4. Explain automatic monomer generation
    5. Run polymerization (automatic monomer generation)
    6. Explain output files
    7. Provide next steps guidance

    Total execution time: ~2-5 minutes for 10 chains of 50 monomers
    """
    print("\n" + "=" * 70)
    print("AUTOPOLY LINEAR PMMA TUTORIAL")
    print("=" * 70)
    print("\nGenerating linear Polymethyl Methacrylate (PMMA) polymer structure")
    print("for LAMMPS molecular dynamics simulations.")
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

        # Step 3: Explain vinyl addition mechanism
        explain_vinyl_addition()

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
        print("\nGenerated files are in: pmma_tutorial/")
        print("Thank you for using AutoPoly!")

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
