#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Linear Polylactic Acid (PLA) Polymer Tutorial

This script demonstrates how to generate atomistic linear Polylactic Acid (PLA)
polymer structures for LAMMPS simulations using AutoPoly's SMILES-based approach.

Key Concepts Demonstrated:
- SMILES notation for monomer definition
- Esterification mechanism for PLA (C-O backbone)
- Automatic monomer variant generation
- GAFF force field for polyesters
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
POLYLACTIC ACID (PLA) CONFIGURATION
====================================

SMILES: "CC(C(=O)O)O"
  Chemical Name: Lactic acid (2-hydroxypropanoic acid)
  Functional Groups: Hydroxyl (-OH), Carboxyl (-COOH)
  Chirality: One chiral center at C2

Polymerization Mechanism: Esterification
  Reaction: -COOH + -OH -> -COO- + H2O
  Backbone: Alternating C-O bonds (heteroatom backbone)
  Connection Atoms: Carbon (from carboxyl) and Oxygen (from hydroxyl)

Force Field: GAFF (General Amber Force Field)
  Recommended for polyesters due to better ester group parameters
  Note: GAFF charges are set to 0.00 and require manual calculation
         using AM1-BCC or RESP methods
"""

PLA_SMILES = "CC(C(=O)O)O"   # Lactic acid SMILES notation
CHAIN_NUM = 10               # Number of polymer chains
DOP = 50                     # Degree of polymerization (monomers per chain)
TOPOLOGY = "linear"          # Linear chain topology
TACTICITY = "atactic"        # Random stereochemistry
FORCE_FIELD = "gaff"         # GAFF force field (recommended for polyesters)

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
        pla_tutorial/
        ├── moltemplate/    # Intermediate Moltemplate files
        ├── input/          # Generated .lt monomer files
        └── output/         # Final LAMMPS input files
    """
    print("\n" + "=" * 70)
    print("STEP 1: Creating System for Output Management")
    print("=" * 70)

    try:
        # Create System object
        system = System(out="pla_tutorial")

        output_path = system.get_folder_path()
        print(f"\nOutput directory: {output_path}")
        print("System created successfully!")

        return system

    except Exception as e:
        print(f"\nError creating system: {e}")
        sys.exit(1)

# ============================================================================
# SECTION 4: Define Linear PLA Polymer
# ============================================================================

def define_pla_polymer():
    """
    Step 2: Define Linear PLA Polymer

    The Polymer class defines polymer structure using SMILES notation.

    Key Parameters:
    - ChainNum: Number of polymer chains
    - Sequence: List of monomer SMILES strings
    - DOP: Degree of polymerization (monomers per chain)
    - topology: Chain topology ("linear" or "ring")
    - tacticity: Stereochemistry arrangement

    For PLA:
    - SMILES "CC(C(=O)O)O" represents lactic acid
    - Esterification creates C-O backbone bonds
    - Atactic tacticity = random chiral configuration
    """
    print("\n" + "=" * 70)
    print("STEP 2: Defining Linear PLA Polymer")
    print("=" * 70)

    print("\nPolymer Configuration:")
    print(f"  Monomer SMILES: {PLA_SMILES}")
    print(f"  Chemical Name: Lactic acid (2-hydroxypropanoic acid)")
    print(f"  Number of Chains: {CHAIN_NUM}")
    print(f"  Degree of Polymerization: {DOP}")
    print(f"  Total Monomers: {CHAIN_NUM * DOP}")
    print(f"  Topology: {TOPOLOGY}")
    print(f"  Tacticity: {TACTICITY}")

    try:
        # Create Polymer object
        polymer = Polymer(
            ChainNum=CHAIN_NUM,
            Sequence=[PLA_SMILES],  # Lactic acid SMILES
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
# SECTION 5: Explain Esterification Mechanism
# ============================================================================

def explain_esterification():
    """
    Step 3: Esterification Mechanism for PLA

    This section explains how PLA polymerization works and why it differs
    from vinyl addition polymers like polyethylene.
    """
    print("\n" + "=" * 70)
    print("STEP 3: Understanding PLA Esterification")
    print("=" * 70)

    print("""
WHAT IS ESTERIFICATION?
-----------------------
Esterification is a condensation reaction where:
  - Carboxylic acid (-COOH) reacts with alcohol (-OH)
  - Produces ester linkage (-COO-) + water (H2O)
  - Forms polymer backbone with alternating C-O bonds

PLA POLYMERIZATION:
-------------------
Monomer: Lactic acid (CH3-CH(OH)-COOH)

Reaction:
  n HO-CH(CH3)-COOH -> [-O-CH(CH3)-CO-]n + n H2O

Backbone Structure:
  ...-O-CH(CH3)-C(=O)-O-CH(CH3)-C(=O)-O-...
       ^       ^      ^       ^
       O       C      O       C
       (alternating C-O bonds)

Connection Atoms:
  - Carbon (C): From carboxyl group
  - Oxygen (O): From hydroxyl group
  - Bond formed: C-O single bond (ester linkage)

WHY THIS MATTERS:
-----------------
1. Heteroatom backbone: Different from vinyl polymers (C-C backbone)
2. Polarity: Ester groups are polar, affecting material properties
3. Biodegradability: Ester bonds are hydrolytically cleavable
4. Crystallinity: Affects thermal and mechanical properties

AUTOPOLY AUTO-DETECTION:
------------------------
AutoPoly automatically detects esterification when:
  - SMILES contains both -COOH and -OH groups
  - Distance between groups allows cyclization
  - DOP > 1 (polymer mode, not single molecule)
    """)

# ============================================================================
# SECTION 6: Explain Monomer Generation
# ============================================================================

def explain_monomer_generation():
    """
    Step 4: Automatic Monomer Generation

    This section explains how AutoPoly automatically generates monomer
    variants from SMILES strings. This happens internally during Polymerization.
    """
    print("\n" + "=" * 70)
    print("STEP 4: Automatic Monomer Generation")
    print("=" * 70)

    print("""
For PLA (SMILES: "CC(C(=O)O)O"), AutoPoly will automatically:

1. DETECT ESTERIFICATION MECHANISM:
   - Carboxyl group (-COOH) provides Carbon connection
   - Hydroxyl group (-OH) provides Oxygen connection
   - Forms ester linkage (-COO-) during polymerization

2. GENERATE 6 MONOMER VARIANT FILES:
   monomer_0i.lt       - Internal monomer (middle of chain)
   monomer_0le.lt      - Left-end monomer (chain start)
   monomer_0re.lt      - Right-end monomer (chain end)
   monomer_0i_T1.lt    - Internal variant (mirror chirality)
   monomer_0le_T1.lt   - Left-end variant (mirror chirality)
   monomer_0re_T1.lt   - Right-end variant (mirror chirality)

3. ASSIGN GAFF ATOM TYPES:
   - C, H, O atoms typed according to GAFF definitions
   - Connection points modified for polymerization
   - Partial charges set to 0.00 (manual calculation required)

4. OPTIMIZE GEOMETRY:
   - Generate 3D coordinates using ETKDG method
   - Align backbone along X-axis
   - Prepare for chain assembly

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
    - Force field integration (GAFF/OPLS-AA)
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
            name="pla_linear",
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
        print("Solution: Check SMILES string and force field choice")
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

    output_path = Path(system.get_folder_path()) / "pla_linear"

    print("\nOutput directory structure:")
    print(f"""
{output_path}/
├── input/               # Moltemplate input files (reusable)
│   ├── monomer_0i.lt   # Internal monomer (middle of chain)
│   ├── monomer_0le.lt  # Left-end monomer (chain start)
│   ├── monomer_0re.lt  # Right-end monomer (chain end)
│   ├── monomer_0i_T1.lt  # Internal (mirror chirality)
│   ├── monomer_0le_T1.lt # Left-end (mirror)
│   ├── monomer_0re_T1.lt # Right-end (mirror)
│   ├── poly_1.lt       # Polymer chain definitions
│   ├── gaff.lt         # GAFF force field import
│   └── gaff.lt.prm     # GAFF parameters
│
└── output/             # LAMMPS input files (ready for simulation)
    ├── system.data     # Atom positions, topology
    ├── system.in       # LAMMPS input script
    ├── system.in.settings  # Force field parameters
    └── system.in.charges    # Atomic charges (0.00 for GAFF)
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
    the PLA polymer structure.
    """
    print("\n" + "=" * 70)
    print("STEP 7: Next Steps for LAMMPS Simulations")
    print("=" * 70)

    print("""
IMPORTANT NOTES:
----------------
1. GAFF CHARGES:
   - All atomic charges are currently set to 0.00
   - You MUST calculate partial charges using:
     * AM1-BCC method (faster, reasonable accuracy)
     * RESP method (more accurate, requires quantum calculations)
   - Update pla_tutorial/pla_linear/output/system.in.charges

2. CHARGE CALCULATION TOOLS:
   - Antechamber (from AmberTools): antechamber -c bcc -m molecule.mol2
   - Open Babel: obabel -ismi -h -o mol2
   - RESP: Gaussian quantum chemistry calculations

3. LAMMPS SIMULATION:
   - Review pla_tutorial/pla_linear/output/system.in
   - Modify simulation parameters as needed
   - Run: lmp -in pla_tutorial/pla_linear/output/system.in

4. RECOMMENDED EQUILIBRATION PROTOCOL:
   a) Energy minimization: minimize 1.0e-4 1000 10000
   b) NVT heating: 100 K to target temperature
   c) NPT compression: Apply pressure to reach target density
   d) Production run: 10-100 ns depending on properties

PLA-SPECIFIC CONSIDERATIONS:
----------------------------
- Glass transition: ~330 K (depends on tacticity)
- Density: ~1.24-1.25 g/cm³ (amorphous)
- Crystallinity: Isotactic > syndiotactic > atactic
- Degradation: Ester bonds can hydrolyze at high T

VALIDATION CHECKLIST:
--------------------
- [ ] Calculate partial charges (AM1-BCC or RESP)
- [ ] Update system.in.charges file
- [ ] Check density (~1.24 g/cm³ for amorphous PLA)
- [ ] Verify bond lengths and angles
- [ ] Test energy conservation in NVE ensemble
- [ ] Check glass transition temperature
    """)

# ============================================================================
# SECTION 10: Main Execution Function
# ============================================================================

def main():
    """
    Execute Complete PLA Polymer Generation Workflow

    This function demonstrates the complete workflow:
    1. Create System for output management
    2. Define PLA polymer using SMILES
    3. Explain esterification mechanism
    4. Explain automatic monomer generation
    5. Run polymerization (automatic monomer generation)
    6. Explain output files
    7. Provide next steps guidance

    Total execution time: ~2-5 minutes for 10 chains of 50 monomers
    """
    print("\n" + "=" * 70)
    print("AUTOPOLY LINEAR PLA TUTORIAL")
    print("=" * 70)
    print("\nGenerating linear Polylactic Acid (PLA) polymer structure")
    print("for LAMMPS molecular dynamics simulations.")
    print("\nConfiguration:")
    print(f"  Monomer: Lactic acid (SMILES: {PLA_SMILES})")
    print(f"  Chains: {CHAIN_NUM}")
    print(f"  DOP: {DOP}")
    print(f"  Force Field: {FORCE_FIELD.upper()}")

    try:
        # Step 1: Create System
        system = create_system()

        # Step 2: Define Polymer
        polymer = define_pla_polymer()

        # Step 3: Explain esterification mechanism
        explain_esterification()

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
        print("\nGenerated files are in: pla_tutorial/")
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
