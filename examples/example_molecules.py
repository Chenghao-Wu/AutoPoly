#!/usr/bin/env python3
"""
Small Molecule Generation Examples

Demonstrates AutoPoly's Molecule class for small molecules (solvents,
additives) using regular SMILES — no [*] wildcards, which are only
used for polymer connection points.

Examples:
  1. Single molecule type (100 water molecules)
  2. Molecule mixture (water + ethanol)
  3. Polymer + molecule mixture (PE chains in water)

Usage:
  python example_molecules.py            # run all examples
  python example_molecules.py --example 2

Key differences from the Polymer class:
  - SMILES: regular SMILES (e.g. "O", "CCO"), no wildcards
  - Parameters: Count, Smiles, Name
  - DOP is always 1 (single molecule)

Force field: GAFF is recommended for small organic molecules.
Gasteiger charges are assigned automatically; for production runs,
replace them with AM1-BCC or RESP charges in system.in.charges.

Requires: pip install -e .  (from the AutoPoly repo root)
"""

import argparse

from AutoPoly import System, Molecule, Polymer, Polymerization

FORCE_FIELD = "gaff"  # GAFF recommended for small molecules

# Complement SMILES for polyethylene (used in example 3)
PE_FIRST = "CC[*]"
PE_MIDDLE = "[*]CC[*]"
PE_LAST = "[*]CC"


def example_single_molecule():
    """Example 1: a box of 100 water molecules."""
    print("\n" + "=" * 70)
    print("EXAMPLE 1: Single Molecule Type (Water)")
    print("=" * 70)

    system = System(out="water_box")

    water = Molecule(Count=100, Smiles="O", Name="water")

    Polymerization(
        name="water_box",
        system=system,
        model=[water],
        force_field=FORCE_FIELD,
    )

    print("\nExample 1 completed: water_box/water_box/")


def example_molecule_mixture():
    """Example 2: a water + ethanol binary mixture."""
    print("\n" + "=" * 70)
    print("EXAMPLE 2: Molecule Mixture (Water + Ethanol)")
    print("=" * 70)

    system = System(out="water_ethanol_mixture")

    water = Molecule(Count=100, Smiles="O", Name="water")
    ethanol = Molecule(Count=20, Smiles="CCO", Name="ethanol")

    Polymerization(
        name="water_ethanol",
        system=system,
        model=[water, ethanol],  # multiple components in one system
        force_field=FORCE_FIELD,
    )

    print("\nExample 2 completed: water_ethanol_mixture/water_ethanol/")
    print("Composition: 100 water (83.3%) + 20 ethanol (16.7%)")


def example_polymer_molecule_mixture():
    """Example 3: PE chains solvated in explicit water."""
    print("\n" + "=" * 70)
    print("EXAMPLE 3: Polymer + Molecule Mixture (PE in Water)")
    print("=" * 70)

    system = System(out="pe_in_water")

    # Complement SMILES sequence: first + (DOP-2)*middle + last
    dop = 50
    sequence = [PE_FIRST] + [PE_MIDDLE] * (dop - 2) + [PE_LAST]
    pe = Polymer(chain_num=5, sequence=sequence, topology="linear")

    water = Molecule(Count=100, Smiles="O", Name="water")

    Polymerization(
        name="pe_water",
        system=system,
        model=[pe, water],  # Polymer and Molecule mix freely
        force_field=FORCE_FIELD,
    )

    print("\nExample 3 completed: pe_in_water/pe_water/")
    print(f"System: 5 PE chains (DOP={dop}) + 100 water molecules")


EXAMPLES = {
    1: example_single_molecule,
    2: example_molecule_mixture,
    3: example_polymer_molecule_mixture,
}


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--example", "-e",
        type=int,
        choices=list(EXAMPLES),
        default=None,
        help="Run a single example (default: run all)",
    )
    args = parser.parse_args()

    selected = [EXAMPLES[args.example]] if args.example else EXAMPLES.values()
    for run in selected:
        run()

    print("\nAll requested examples completed.")


if __name__ == "__main__":
    main()
