#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Small Molecule Generation Tutorial

This script demonstrates how to generate LAMMPS data files for small molecules
(e.g., water, benzene, ethanol) using AutoPoly's Molecule class.

Key Concepts Demonstrated:
- Regular SMILES notation (without [*] wildcards)
- Single molecule generation
- Molecule mixtures
- Polymer + molecule mixtures
- Force field selection (GAFF recommended for small molecules)
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
MOLECULE GENERATION CONFIGURATION
==================================

Unlike polymers that use pSMILES with [*] wildcards for connection points,
molecules use regular SMILES notation without any wildcards.

Common SMILES Examples:
  Water: "O"
  Ethanol: "CCO"
  Benzene: "c1ccccc1"
  Methane: "C"
  Acetone: "CC(=O)C"
  Toluene: "Cc1ccccc1"

Force Field Recommendations:
  - GAFF (General AMBER Force Field): Recommended for small organic molecules
  - OPLS-AA: Can be used for some molecules, but GAFF is more comprehensive

Key Differences from Polymer:
  - No [*] wildcards in SMILES
  - Count parameter instead of ChainNum
  - DOP is always 1 (single molecule)
"""

# Example 1: Water molecules
WATER_SMILES = "O"
WATER_COUNT = 100
FORCE_FIELD = "gaff"  # GAFF recommended for small molecules

# Example 2: Molecule mixture (water + ethanol)
ETHANOL_SMILES = "CCO"
ETHANOL_COUNT = 20

# Example 3: Polymer + molecule mixture
# Complement SMILES for Polyethylene
PE_FIRST = "CC[*]"
PE_MIDDLE = "[*]CC[*]"
PE_LAST = "[*]CC"
POLYMER_CHAIN_NUM = 5
POLYMER_DOP = 50


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
        from AutoPoly import System, Molecule, Polymer, Polymerization
        return System, Molecule, Polymer, Polymerization
    except ImportError as e:
        print(f"\nError: Could not import AutoPoly modules: {e}")
        print("\nPlease ensure AutoPoly is installed:")
        print("  cd /home/zhenghaowu/mcp_lammps_poly/AutoPoly")
        print("  pip install -e .")
        sys.exit(1)

System, Molecule, Polymer, Polymerization = setup_imports()


# ============================================================================
# SECTION 3: Example 1 - Single Molecule Type
# ============================================================================

def example_single_molecule():
    """
    Example 1: Generate a system with a single molecule type (water).

    This demonstrates the simplest use case - generating multiple copies
    of the same molecule type.
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 1: Single Molecule Type (Water)")
    print("=" * 70)

    print("\nConfiguration:")
    print(f"  Molecule: Water (SMILES: {WATER_SMILES})")
    print(f"  Count: {WATER_COUNT}")
    print(f"  Force Field: {FORCE_FIELD.upper()}")

    try:
        # Step 1: Create System
        print("\n--- Step 1: Creating System ---")
        system = System(out="water_box")
        print(f"Output directory: {system.get_folder_path()}")

        # Step 2: Define Molecule
        print("\n--- Step 2: Defining Water Molecule ---")
        water = Molecule(
            Count=WATER_COUNT,
            Smiles=WATER_SMILES,
            Name="water"
        )
        print(f"Created {water.Count} water molecules")
        print(f"Molecule name: {water.molecule_name}")
        print(f"SMILES: {water.Smiles}")

        # Step 3: Run Polymerization (generates LAMMPS files)
        print("\n--- Step 3: Generating LAMMPS Files ---")
        polymerization = Polymerization(
            name="water_box",
            system=system,
            model=[water],
            force_field=FORCE_FIELD
        )

        print("\nExample 1 completed successfully!")
        return system, water

    except Exception as e:
        print(f"\nError in Example 1: {e}")
        import traceback
        traceback.print_exc()
        return None, None


# ============================================================================
# SECTION 4: Example 2 - Molecule Mixture
# ============================================================================

def example_molecule_mixture():
    """
    Example 2: Generate a system with multiple molecule types.

    This demonstrates how to create mixtures of different molecules,
    useful for solvent systems or binary mixtures.
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 2: Molecule Mixture (Water + Ethanol)")
    print("=" * 70)

    print("\nConfiguration:")
    print(f"  Water: SMILES={WATER_SMILES}, Count={WATER_COUNT}")
    print(f"  Ethanol: SMILES={ETHANOL_SMILES}, Count={ETHANOL_COUNT}")
    print(f"  Force Field: {FORCE_FIELD.upper()}")

    try:
        # Step 1: Create System
        print("\n--- Step 1: Creating System ---")
        system = System(out="water_ethanol_mixture")
        print(f"Output directory: {system.get_folder_path()}")

        # Step 2: Define Molecules
        print("\n--- Step 2: Defining Molecules ---")

        # Water molecules
        water = Molecule(
            Count=WATER_COUNT,
            Smiles=WATER_SMILES,
            Name="water"
        )
        print(f"Created {water.Count} water molecules")

        # Ethanol molecules
        ethanol = Molecule(
            Count=ETHANOL_COUNT,
            Smiles=ETHANOL_SMILES,
            Name="ethanol"
        )
        print(f"Created {ethanol.Count} ethanol molecules")

        print(f"\nMixture composition:")
        print(f"  Water: {WATER_COUNT} molecules ({100*WATER_COUNT/(WATER_COUNT+ETHANOL_COUNT):.1f}%)")
        print(f"  Ethanol: {ETHANOL_COUNT} molecules ({100*ETHANOL_COUNT/(WATER_COUNT+ETHANOL_COUNT):.1f}%)")

        # Step 3: Run Polymerization
        print("\n--- Step 3: Generating LAMMPS Files ---")
        polymerization = Polymerization(
            name="water_ethanol",
            system=system,
            model=[water, ethanol],  # List of molecules
            force_field=FORCE_FIELD
        )

        print("\nExample 2 completed successfully!")
        return system, water, ethanol

    except Exception as e:
        print(f"\nError in Example 2: {e}")
        import traceback
        traceback.print_exc()
        return None, None, None


# ============================================================================
# SECTION 5: Example 3 - Polymer + Molecule Mixture
# ============================================================================

def example_polymer_molecule_mixture():
    """
    Example 3: Generate a system with both polymers and molecules.

    This demonstrates how to create a polymer solution, e.g.,
    polymer chains solvated in explicit solvent molecules.
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 3: Polymer + Molecule Mixture (PE in Water)")
    print("=" * 70)

    print("\nConfiguration:")
    print(f"  Polymer: Polyethylene (complement SMILES)")
    print(f"  Chains: {POLYMER_CHAIN_NUM}")
    print(f"  DOP: {POLYMER_DOP}")
    print(f"  Solvent: Water (SMILES: {WATER_SMILES}, Count={WATER_COUNT})")
    print(f"  Force Field: {FORCE_FIELD.upper()}")

    try:
        # Step 1: Create System
        print("\n--- Step 1: Creating System ---")
        system = System(out="pe_in_water")
        print(f"Output directory: {system.get_folder_path()}")

        # Step 2: Define Polymer
        print("\n--- Step 2: Defining Polyethylene Polymer ---")
        # Build complement SMILES sequence: first + (DOP-2)*middle + last
        sequence = [PE_FIRST] + [PE_MIDDLE] * (POLYMER_DOP - 2) + [PE_LAST]
        pe = Polymer(
            chain_num=POLYMER_CHAIN_NUM,
            sequence=sequence,
            topology="linear"
        )
        print(f"Created {pe.chain_num} PE chains")
        print(f"  DOP: {pe.dop}")
        print(f"  Total monomers: {pe.chain_num * pe.dop}")

        # Step 3: Define Solvent Molecules
        print("\n--- Step 3: Defining Water Solvent ---")
        water = Molecule(
            Count=WATER_COUNT,
            Smiles=WATER_SMILES,
            Name="water"
        )
        print(f"Created {water.Count} water molecules")

        # Step 4: Run Polymerization
        print("\n--- Step 4: Generating LAMMPS Files ---")
        polymerization = Polymerization(
            name="pe_water",
            system=system,
            model=[pe, water],  # List containing both Polymer and Molecule
            force_field=FORCE_FIELD
        )

        print("\nExample 3 completed successfully!")
        return system, pe, water

    except Exception as e:
        print(f"\nError in Example 3: {e}")
        import traceback
        traceback.print_exc()
        return None, None, None


# ============================================================================
# SECTION 6: Explain Key Differences
# ============================================================================

def explain_key_differences():
    """
    Explain key differences between Polymer and Molecule classes.
    """
    print("\n" + "=" * 70)
    print("KEY DIFFERENCES: Polymer vs Molecule")
    print("=" * 70)

    print("""
POLYMER CLASS:
  - Purpose: Define polymer chains with repeat units
  - SMILES: Complement SMILES format with explicit first/middle/last variants
  - Parameters: chain_num, sequence, topology, tacticity
  - DOP: Automatically derived from sequence length
  - Variants: User provides explicit first/middle/last monomer SMILES

MOLECULE CLASS:
  - Purpose: Define small molecules (solvents, additives)
  - SMILES: Regular SMILES without wildcards (e.g., "O", "CCO")
  - Parameters: Count, Smiles, Name
  - DOP: Always 1 (single molecule)
  - Variants: Single .lt file (no variants needed)

SIMILARITIES:
  - Both use the same Polymerization workflow
  - Both can be mixed in the same system
  - Both support the same force fields (OPLS-AA, GAFF, LOPLS)
  - Both generate LAMMPS data files via Moltemplate

USAGE PATTERN:
  # Polymer with complement SMILES
  PE_FIRST = "CC[*]"
  PE_MIDDLE = "[*]CC[*]"
  PE_LAST = "[*]CC"
  sequence = [PE_FIRST] + [PE_MIDDLE] * 98 + [PE_LAST]
  polymer = Polymer(chain_num=5, sequence=sequence)

  # Molecule
  water = Molecule(Count=100, Smiles="O", Name="water")

  # Mixed system
  polymerization = Polymerization(
      name="polymer_solution",
      system=system,
      model=[polymer, water],  # Both types work together
      force_field="gaff"
  )
    """)


# ============================================================================
# SECTION 7: Common SMILES Examples
# ============================================================================

def show_common_smiles():
    """
    Show common SMILES strings for small molecules.
    """
    print("\n" + "=" * 70)
    print("COMMON SMILES FOR SMALL MOLECULES")
    print("=" * 70)

    print("""
WATER AND SOLVENTS:
  Water: "O"
  Methanol: "CO"
  Ethanol: "CCO"
  Isopropanol: "CC(C)O"
  Acetone: "CC(=O)C"
  DMSO: "CS(C)=O"
  Chloroform: "ClC(Cl)Cl"
  Benzene: "c1ccccc1"
  Toluene: "Cc1ccccc1"

HYDROCARBONS:
  Methane: "C"
  Ethane: "CC"
  Propane: "CCC"
  Butane: "CCCC"
  Cyclohexane: "C1CCCCC1"

ACIDS AND BASES:
  Acetic acid: "CC(=O)O"
  Formic acid: "C(=O)O"
  Ammonia: "N"

ESTERS:
  Methyl acetate: "CC(=O)OC"
  Ethyl acetate: "CC(=O)OCC"

AMINES:
  Methylamine: "CN"
  Ethylamine: "CCN"
  Aniline: "c1ccccc1N"

TIPS FOR FINDING SMILES:
  1. Use online databases: PubChem, ChemSpider
  2. Search: "[molecule name] SMILES"
  3. Use chemical drawing software: ChemDraw, MarvinSketch
  4. Verify SMILES using tools: Open Babel, RDKit
    """)


# ============================================================================
# SECTION 8: Explain Force Field Choice
# ============================================================================

def explain_force_fields():
    """
    Explain force field selection for small molecules.
    """
    print("\n" + "=" * 70)
    print("FORCE FIELD SELECTION FOR MOLECULES")
    print("=" * 70)

    print("""
GAFF (General AMBER Force Field):
  Recommended for: Most small organic molecules
  Strengths:
    - Comprehensive atom types for organic molecules
    - Well-tested for drug-like molecules
    - Good for diverse functional groups
  Limitations:
    - Requires charge calculation (AM1-BCC or RESP)
    - May not cover all inorganic/metal-organic compounds

OPLS-AA (Optimized Potentials for Liquid Simulations - All Atom):
  Recommended for: Some organic molecules, polymers
  Strengths:
    - Good for liquids and biomolecules
    - Well-parameterized for many organic compounds
  Limitations:
    - Less comprehensive than GAFF for diverse molecules
    - May not have parameters for all functional groups

LOPLS (Liquid Optimized OPLS):
  Recommended for: Liquid-phase simulations
  Strengths:
    - Optimized for liquid properties
    - Better density and enthalpy of vaporization
  Limitations:
    - Less extensively parameterized than OPLS-AA

RECOMMENDATION:
  For small molecules, use GAFF as the default choice.
  Only use OPLS-AA/LOPLS if you have specific parameterization needs.

CHARGE CALCULATION (CRITICAL FOR GAFF):
  1. AM1-BCC method (faster, reasonable accuracy):
     - Use Antechamber: antechamber -c bcc -m molecule.mol2
  2. RESP method (more accurate, requires QM):
     - Use Gaussian + RESP fitting
  3. Update system.in.charges file with calculated charges
    """)


# ============================================================================
# SECTION 9: Main Execution Function
# ============================================================================

def main():
    """
    Execute molecule generation examples.

    This function runs three examples:
    1. Single molecule type (water)
    2. Molecule mixture (water + ethanol)
    3. Polymer + molecule mixture (PE in water)
    """
    print("\n" + "=" * 70)
    print("AUTOPOLY SMALL MOLECULE GENERATION TUTORIAL")
    print("=" * 70)
    print("\nThis tutorial demonstrates how to generate LAMMPS data files")
    print("for small molecules using AutoPoly's Molecule class.")

    try:
        # Show key differences
        explain_key_differences()

        # Show common SMILES
        show_common_smiles()

        # Explain force fields
        explain_force_fields()

        # Ask user which example to run
        print("\n" + "=" * 70)
        print("SELECT EXAMPLE TO RUN")
        print("=" * 70)
        print("\n1. Single molecule type (100 water molecules)")
        print("2. Molecule mixture (80 water + 20 ethanol)")
        print("3. Polymer + molecule mixture (5 PE chains + 100 water)")
        print("4. Run all examples")
        print("5. Exit")

        choice = input("\nEnter choice (1-5): ").strip()

        if choice == "1":
            example_single_molecule()
        elif choice == "2":
            example_molecule_mixture()
        elif choice == "3":
            example_polymer_molecule_mixture()
        elif choice == "4":
            print("\n" + "=" * 70)
            print("RUNNING ALL EXAMPLES")
            print("=" * 70)
            example_single_molecule()
            example_molecule_mixture()
            example_polymer_molecule_mixture()
        elif choice == "5":
            print("\nExiting tutorial.")
            return
        else:
            print("\nInvalid choice. Exiting.")
            return

        print("\n" + "=" * 70)
        print("TUTORIAL COMPLETED SUCCESSFULLY!")
        print("=" * 70)
        print("\nThank you for using AutoPoly!")

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
