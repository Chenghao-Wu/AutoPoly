#!/usr/bin/env python3
"""
PEO Monte Carlo Placement Example

Demonstrates the Monte Carlo placement features in AutoPoly:
1. MC Random Placement: Random positions/orientations with collision detection
2. MC Chain Growth: Self-avoiding random walk for polymer chains

Compares:
- Grid placement (default): Deterministic grid arrangement
- MC random placement: Random placement with collision avoidance
- MC chain growth: Self-avoiding random walk for chain conformations

Configuration:
- Monomer: Ethylene oxide (-CH2-CH2-O-)
- Chain length: 20 monomers
- Number of chains: 5

Requires: pip install -e .  (from the AutoPoly repo root)
"""

from AutoPoly import System, Polymer, GeometryConfig, generate

# PEO configuration
PEO_CONFIG = {
    "first_smiles": "CCO[*]",
    "middle_smiles": "[*]CCO[*]",
    "last_smiles": "[*]CCO",
    "chain_num": 5,
    "dop": 20,
    "topology": "linear",
    "tacticity": "atactic"
}

# MC configuration
MC_CONFIG = {
    "mc_max_attempts": 10000,
    "monomer_density": 0.03  # Low density for easy placement
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


def test_grid_placement(system, polymer):
    """
    Standard grid placement (default).

    Grid placement arranges polymers in a deterministic grid pattern.
    This is fast and predictable, but may not represent realistic
    initial configurations for equilibration.
    """
    print("\n[1/3] Testing grid placement...")

    try:
        generate(
            system,
            "peo_grid",
            [polymer],
            force_field="oplsaa",
            strategy="grid",
        )
        print("  SUCCESS: peo_grid completed")
        return True
    except Exception as e:
        print(f"  FAILED: peo_grid - {e}")
        import traceback
        traceback.print_exc()
        return False


def test_mc_random_placement(system, polymer):
    """
    MC random placement with collision detection.

    Monte Carlo random placement positions polymers at random
    locations and orientations within the simulation box, while
    avoiding overlaps between molecules.
    """
    print("\n[2/3] Testing MC random placement...")

    try:
        generate(
            system,
            "peo_mc_random",
            [polymer],
            force_field="oplsaa",
            strategy="mc_random",
            mc_max_attempts=MC_CONFIG["mc_max_attempts"],
            monomer_density=MC_CONFIG["monomer_density"],
        )
        print("  SUCCESS: peo_mc_random completed")
        return True
    except Exception as e:
        print(f"  FAILED: peo_mc_random - {e}")
        import traceback
        traceback.print_exc()
        return False


def test_mc_chain_growth(system, polymer):
    """
    MC chain growth with random placement.

    Monte Carlo chain growth generates polymer conformations using
    a self-avoiding random walk algorithm. Combined with MC random
    placement, this produces more realistic initial polymer structures
    with proper chain entanglement.
    """
    print("\n[3/3] Testing MC chain growth...")

    try:
        generate(
            system,
            "peo_mc_chain",
            [polymer],
            force_field="oplsaa",
            strategy="mc_random",
            geometry_config=GeometryConfig(use_mc_chain_growth=True),
            mc_max_attempts=MC_CONFIG["mc_max_attempts"],
            monomer_density=MC_CONFIG["monomer_density"],
        )
        print("  SUCCESS: peo_mc_chain completed")
        return True
    except Exception as e:
        print(f"  FAILED: peo_mc_chain - {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all placement method tests."""
    print("PEO Monte Carlo Placement Example")
    print("=================================")
    print("Configuration:")
    print(f"  - Chains: {PEO_CONFIG['chain_num']}")
    print(f"  - DOP: {PEO_CONFIG['dop']}")
    print(f"  - Total monomers: {PEO_CONFIG['chain_num'] * PEO_CONFIG['dop']}")
    print(f"  - Topology: {PEO_CONFIG['topology']}")
    print(f"  - Tacticity: {PEO_CONFIG['tacticity']}")
    print(f"\nMC Settings:")
    print(f"  - Max attempts: {MC_CONFIG['mc_max_attempts']}")
    print(f"  - Monomer density: {MC_CONFIG['monomer_density']} monomers/A^3")

    # Build sequence
    sequence = build_sequence(PEO_CONFIG)

    # Track results
    results = {}

    # Test 1: Grid placement
    system1 = System(out="peo_mc_example")
    polymer1 = Polymer(
        chain_num=PEO_CONFIG["chain_num"],
        sequence=sequence,
        topology=PEO_CONFIG["topology"],
        tacticity=PEO_CONFIG["tacticity"]
    )
    results["grid"] = test_grid_placement(system1, polymer1)

    # Test 2: MC random placement
    system2 = System(out="peo_mc_example")
    polymer2 = Polymer(
        chain_num=PEO_CONFIG["chain_num"],
        sequence=sequence,
        topology=PEO_CONFIG["topology"],
        tacticity=PEO_CONFIG["tacticity"]
    )
    results["mc_random"] = test_mc_random_placement(system2, polymer2)

    # Test 3: MC chain growth
    system3 = System(out="peo_mc_example")
    polymer3 = Polymer(
        chain_num=PEO_CONFIG["chain_num"],
        sequence=sequence,
        topology=PEO_CONFIG["topology"],
        tacticity=PEO_CONFIG["tacticity"]
    )
    results["mc_chain"] = test_mc_chain_growth(system3, polymer3)

    # Summary
    print(f"\n{'='*40}")
    print("SUMMARY")
    print(f"{'='*40}")
    for method, success in results.items():
        status = "PASS" if success else "FAIL"
        print(f"  {method:12s}: {status}")

    passed = sum(results.values())
    total = len(results)
    print(f"\nTotal: {passed}/{total} passed")
    print(f"\nOutput: peo_mc_example/")
    print("  - peo_grid/        : Grid placement output")
    print("  - peo_mc_random/   : MC random placement output")
    print("  - peo_mc_chain/    : MC chain growth output")


if __name__ == "__main__":
    main()
