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
    """Add AutoPoly to Python path and import modules."""
    project_root = Path(__file__).parent.parent
    sys.path.insert(0, str(project_root))

    try:
        from AutoPoly import System, Molecule, Polymerization
        return System, Molecule, Polymerization
    except ImportError as e:
        print(f"Error importing AutoPoly: {e}")
        print("Install with: pip install -e /path/to/AutoPoly")
        sys.exit(1)


System, Molecule, Polymerization = setup_imports()
# Create output directory structure
system = System(out="d4ppd")

# Define benzene molecule
# SMILES "c1ccccc1" represents the benzene ring structure
benzene = Molecule(
    Count=1,              # Number of benzene molecules
    Smiles="CC(C)Nc1c(C)cc(Nc2ccc(C)cc2)cc1",     # Benzene SMILES notation (aromatic ring)
    Name="d4ppd"         # Molecule identifier
)

# Run polymerization workflow to generate LAMMPS files
# GAFF force field is recommended for organic molecules with aromatic rings
polymerization = Polymerization(
    name="d4ppd",
    system=system,
    model=[benzene],
    force_field="gaff2"    # GAFF2 (includes nq and other extended atom types)
)

print("Benzene system created successfully!")
print(f"Output directory: benzene_system_100/benzene_100/")
print(f"Number of benzene molecules: 100")
print(f"Total atoms: 100 x 12 = 1200 atoms (C6H6)")
