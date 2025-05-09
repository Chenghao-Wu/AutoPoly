from pathlib import Path
from AutoPoly.extern.rdlt import RDlt

MONOMER_BANK = Path(__file__).parent / 'extern' / 'Monomer_bank'


def create_monomer_from_smiles(smiles: str, name: str, monomer_type: str = 'middle',
                               opls_fdef: str = None, lopls_fdef: str = None, charge: bool = True, refresh: bool = False, loplsflag: bool = False):
    """
    Generate a monomer .lt file from a SMILES string using RDlt.
    Args:
        smiles (str): The SMILES string for the monomer.
        name (str): The base name for the monomer (e.g., 'PEAA').
        monomer_type (str): 'middle', 'left', or 'right' (default: 'middle').
        opls_fdef (str): Path to OPLS feature definition file (optional).
        lopls_fdef (str): Path to LOPLS feature definition file (optional).
        charge (bool): Whether to assign charges (default: True).
        refresh (bool): Whether to refresh feature definitions (default: False).
        loplsflag (bool): Whether to use LOPLS atom typing (default: False).
    """
    # Determine suffix for monomer type
    if monomer_type == 'middle':
        suffix = 'i'
    elif monomer_type == 'left':
        suffix = 'le'
    elif monomer_type == 'right':
        suffix = 're'
    else:
        raise ValueError(f"Unknown monomer_type: {monomer_type}")

    lt_name = f"{name}{suffix}"
    lt_file = MONOMER_BANK / f"{lt_name}.lt"

    rdlt = RDlt(smiles=smiles)
    rdlt.run(
        to_file=str(lt_file),
        name=lt_name,
        fdef=opls_fdef if opls_fdef else str(Path(__file__).parent / 'extern' / 'rdlt_data' / 'opls_lt.fdefn'),
        lfdef=lopls_fdef if lopls_fdef else str(Path(__file__).parent / 'extern' / 'rdlt_data' / 'lopls_lt.fdefn'),
        charge=charge,
        refresh=refresh,
        loplsflag=loplsflag
    )
    # Optionally store in bank (overwrites if exists)
    rdlt.store_bank(flag=True)

def convert_monomer_to_reacted_monomer(middle):


# Example usage:
# create_monomer_from_smiles('CC', 'PEAA_test', monomer_type='middle')
