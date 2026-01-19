"""
Example: Create a system of 100 benzene molecules using AutoPoly

This script demonstrates how to build a molecular system containing
100 benzene molecules for LAMMPS simulations.

Benzene (C6H6) is created using the Molecule class with SMILES notation.
"""
import sys
import os
from pathlib import Path
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

# Create output directory structure
system = System(out="benzene_system_100")

# Define benzene molecule
# SMILES "c1ccccc1" represents the benzene ring structure
benzene = Molecule(
    Count=100,              # Number of benzene molecules
    Smiles="c1ccccc1",     # Benzene SMILES notation (aromatic ring)
    Name="benzene"         # Molecule identifier
)

# Run polymerization workflow to generate LAMMPS files
# GAFF force field is recommended for organic molecules with aromatic rings
polymerization = Polymerization(
    name="benzene_100",
    system=system,
    model=[benzene],
    force_field="gaff"     # GAFF (General AMBER Force Field)
)

print("Benzene system created successfully!")
print(f"Output directory: benzene_system_100/benzene_100/")
print(f"Number of benzene molecules: 100")
print(f"Total atoms: 100 x 12 = 1200 atoms (C6H6)")
