#!/usr/bin/env python3
"""
PEO Force Field Test Example

Tests all supported force fields with polyethylene oxide (PEO):
- OPLS-AA: Standard force field for vinyl polymers
- GAFF: General Amber Force Field
- GAFF2: Extended GAFF
- L-OPLS: Long-chain optimized OPLS
- DREIDING: Generic force field (requires external charges)
- COMPASS: Class2 force field (requires LAMMPS CLASS2 package)

Configuration:
- Monomer: Ethylene oxide (-CH2-CH2-O-)
- Chain length: 10 monomers
- Number of chains: 10
"""

import sys
from pathlib import Path

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from AutoPoly import System, Polymer, Polymerization

# Force fields to test
FORCE_FIELDS = ["oplsaa", "gaff", "gaff2", "lopls", "dreiding", "compass"]

# PEO configuration
PEO_CONFIG = {
    "first_smiles": "CCO[*]",
    "middle_smiles": "[*]CCO[*]",
    "last_smiles": "[*]CCO",
    "chain_num": 10,
    "dop": 10,
    "topology": "linear",
    "tacticity": "atactic"
}


def build_sequence(config):
    """Build complement SMILES sequence."""
    dop = config["dop"]
    if dop == 1:
        return ["CCO"]  # Single monomer, no wildcards
    elif dop == 2:
        return [config["first_smiles"], config["last_smiles"]]
    else:
        return (
            [config["first_smiles"]] +
            [config["middle_smiles"]] * (dop - 2) +
            [config["last_smiles"]]
        )


def test_force_field(force_field, output_base="peo_test"):
    """Test PEO with a specific force field."""
    output_name = f"{output_base}_{force_field}"
    print(f"\n{'='*60}")
    print(f"Testing {force_field.upper()} force field")
    print(f"{'='*60}")

    try:
        # Create system
        system = System(out=output_name)

        # Build sequence
        sequence = build_sequence(PEO_CONFIG)
        print(f"Sequence length: {len(sequence)} monomers")

        # Create polymer
        polymer = Polymer(
            chain_num=PEO_CONFIG["chain_num"],
            sequence=sequence,
            topology=PEO_CONFIG["topology"],
            tacticity=PEO_CONFIG["tacticity"]
        )

        # Run polymerization
        poly = Polymerization(
            name=f"peo_{force_field}",
            system=system,
            model=[polymer],
            force_field=force_field,
            run=True
        )

        print(f"SUCCESS: {force_field} completed")
        return True

    except Exception as e:
        print(f"FAILED: {force_field} - {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all force field tests."""
    print("PEO Force Field Comparison Test")
    print("================================")
    print(f"Configuration:")
    print(f"  - Chains: {PEO_CONFIG['chain_num']}")
    print(f"  - DOP: {PEO_CONFIG['dop']}")
    print(f"  - Topology: {PEO_CONFIG['topology']}")
    print(f"  - Tacticity: {PEO_CONFIG['tacticity']}")

    results = {}
    for ff in FORCE_FIELDS:
        results[ff] = test_force_field(ff)

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for ff, success in results.items():
        status = "PASS" if success else "FAIL"
        print(f"  {ff:12s}: {status}")

    passed = sum(results.values())
    total = len(results)
    print(f"\nTotal: {passed}/{total} passed")


if __name__ == "__main__":
    main()
