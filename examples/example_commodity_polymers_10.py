#!/usr/bin/env python3
"""Generate 10 commodity polymers for LAMMPS using AutoPoly.

Polymers: PE, PP, PS, PVC, PVAc, PMMA, PAN, PB, PI, PEO
Force field: OPLS-AA (except PEO uses GAFF)
"""

import sys
import argparse
from pathlib import Path
from typing import Dict, List

POLYMER_CONFIGS: Dict[str, Dict] = {
    "polyethylene": {"smiles": "[*]CC[*]", "chain_num": 10, "dop": 50,
                     "topology": "linear", "tacticity": "atactic"},
    "polypropylene": {"smiles": "[*]CC([*])(C)", "chain_num": 10, "dop": 50,
                      "topology": "linear", "tacticity": "atactic"},
    "polystyrene": {"smiles": "[*]CC([*])c1ccccc1", "chain_num": 10, "dop": 10,
                    "topology": "linear", "tacticity": "atactic"},
    "polyvinyl_chloride": {"smiles": "[*]CC([*])Cl", "chain_num": 10, "dop": 50,
                           "topology": "linear", "tacticity": "atactic"},
    "polyvinyl_acetate": {"smiles": "[*]CC([*])OC(=O)C", "chain_num": 10, "dop": 50,
                          "topology": "linear", "tacticity": "atactic"},
    "pmma": {"smiles": "[*]CC([*])(C)C(=O)OC", "chain_num": 10, "dop": 50,
             "topology": "linear", "tacticity": "atactic"},
    "polyacrylonitrile": {"smiles": "[*]CC([*])C#N", "chain_num": 10, "dop": 50,
                          "topology": "linear", "tacticity": "atactic"},
    "polybutadiene": {"smiles": "[*]C=CC[*]", "chain_num": 10, "dop": 50,
                      "topology": "linear", "tacticity": "atactic"},
    "polyisoprene": {"smiles": "[*]C=CC(C)[*]", "chain_num": 10, "dop": 50,
                     "topology": "linear", "tacticity": "atactic"},
    "polyethylene_oxide": {"smiles": "[*]CCO[*]", "chain_num": 10, "dop": 50,
                           "topology": "linear", "tacticity": "atactic"}
}

FORCE_FIELD = "oplsaa"


def setup_imports():
    """Add AutoPoly to Python path and import modules."""
    project_root = Path(__file__).parent.parent
    sys.path.insert(0, str(project_root))

    try:
        from AutoPoly import System, Polymer, Polymerization
        return System, Polymer, Polymerization
    except ImportError as e:
        print(f"Error importing AutoPoly: {e}")
        print("Install with: pip install -e /path/to/AutoPoly")
        sys.exit(1)


System, Polymer, Polymerization = setup_imports()


def generate_single_polymer(polymer_id: str, config: Dict, base_output_dir: str = "commodity_polymers"):
    """Generate a single polymer structure."""
    print(f"\n{'='*60}\nGenerating: {polymer_id}")
    print(f"pSMILES: {config['smiles']}, Chains: {config['chain_num']}, DOP: {config['dop']}")

    try:
        output_path = f"{base_output_dir}/{polymer_id}"
        system = System(out=output_path)

        sequence = [config["smiles"]] * config["dop"]
        polymer = Polymer(
            chain_num=config["chain_num"],
            sequence=sequence,
            topology=config["topology"],
            tacticity=config["tacticity"]
        )

        force_field = config.get("force_field", FORCE_FIELD)
        print(f"Running polymerization ({force_field.upper()})...")

        Polymerization(name=polymer_id, system=system, model=[polymer], force_field=force_field, run=True)
        print(f"✓ {polymer_id} completed")
        return True

    except Exception as e:
        print(f"✗ {polymer_id} failed: {e}")
        return False


def generate_commodity_polymers(polymer_ids: List[str] = None, base_output_dir: str = "commodity_polymers"):
    """Generate multiple polymers in batch."""
    if polymer_ids is None:
        polymer_ids = list(POLYMER_CONFIGS.keys())

    print(f"\nGenerating {len(polymer_ids)} polymers ({FORCE_FIELD.upper()})")
    results = {"successful": [], "failed": [], "total": len(polymer_ids)}

    for i, polymer_id in enumerate(polymer_ids, 1):
        print(f"\n[{i}/{len(polymer_ids)}] {polymer_id}")
        config = POLYMER_CONFIGS[polymer_id]
        if generate_single_polymer(polymer_id, config, base_output_dir):
            results["successful"].append(polymer_id)
        else:
            results["failed"].append(polymer_id)

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate commodity polymers with AutoPoly")
    parser.add_argument('--polymers', '-p', nargs='+',
                        choices=list(POLYMER_CONFIGS.keys()) + ['all'], default=['all'],
                        help='Polymers to generate (default: all)')
    args = parser.parse_args()

    polymer_ids = None if 'all' in args.polymers else args.polymers
    if polymer_ids:
        print(f"\nGenerating: {', '.join(polymer_ids)}")

    try:
        results = generate_commodity_polymers(polymer_ids)

        print(f"\n{'='*60}\nSUMMARY")
        print(f"Total: {results['total']}, Success: {len(results['successful'])}, Failed: {len(results['failed'])}")
        if results['successful']:
            print("\nSuccess:", ', '.join(results['successful']))
        if results['failed']:
            print("\nFailed:", ', '.join(results['failed']))

        print(f"\nOutput: commodity_polymers/")

    except KeyboardInterrupt:
        print("\nInterrupted by user")
        sys.exit(0)
