from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import AllChem
import os
import xml.etree.ElementTree as ET

# Path to your OPLS feature definition file
fdef_path = "../AutoPoly/extern/rdlt_data/opls_lt.fdefn"
cml_path = "PDMS_3monomer.cml"

# Check if files exist
if not os.path.exists(fdef_path):
    raise FileNotFoundError(f"OPLS feature definition file not found: {fdef_path}")
if not os.path.exists(cml_path):
    raise FileNotFoundError(f"CML file not found: {cml_path}")

def parse_cml(filename):
    tree = ET.parse(filename)
    root = tree.getroot()
    atoms = []
    atom_id_map = {}
    for atom in root.find('.//atomArray'):
        aid = atom.attrib['id']
        element = atom.attrib['elementType']
        x = float(atom.attrib['x3'])
        y = float(atom.attrib['y3'])
        z = float(atom.attrib['z3'])
        charge = float(atom.attrib.get('formalCharge', 0))
        atoms.append({'id': aid, 'element': element, 'x': x, 'y': y, 'z': z, 'charge': charge})
        atom_id_map[aid] = len(atoms) - 1
    bonds = []
    for bond in root.find('.//bondArray'):
        refs = bond.attrib['atomRefs2'].split()
        order = bond.attrib.get('order', '1')
        bonds.append({'a1': refs[0], 'a2': refs[1], 'order': order})
    return atoms, bonds, atom_id_map

def build_rdkit_mol(atoms, bonds, atom_id_map):
    mol = Chem.RWMol()
    atom_idx_map = {}
    for i, atom in enumerate(atoms):
        a = Chem.Atom(atom['element'])
        idx = mol.AddAtom(a)
        atom_idx_map[atom['id']] = idx
    for bond in bonds:
        a1 = atom_idx_map[bond['a1']]
        a2 = atom_idx_map[bond['a2']]
        order = bond['order']
        if order == '1':
            btype = Chem.BondType.SINGLE
        elif order == '2':
            btype = Chem.BondType.DOUBLE
        elif order == '3':
            btype = Chem.BondType.TRIPLE
        else:
            btype = Chem.BondType.SINGLE
        mol.AddBond(a1, a2, btype)
    m = mol.GetMol()
    conf = Chem.Conformer(m.GetNumAtoms())
    for i, atom in enumerate(atoms):
        conf.SetAtomPosition(i, (atom['x'], atom['y'], atom['z']))
    m.AddConformer(conf, assignId=True)
    Chem.SanitizeMol(m)
    return m

# 1. Read the CML file as an RDKit molecule
atoms, bonds, atom_id_map = parse_cml(cml_path)
mol = build_rdkit_mol(atoms, bonds, atom_id_map)

# 2. Add hydrogens and generate 3D coordinates if needed
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())

# 3. Build the feature factory and assign atom types
factory = ChemicalFeatures.BuildFeatureFactory(fdef_path)
features = factory.GetFeaturesForMol(mol)

# 4. Assign atom types as properties
for f in features:
    atom_idx = f.GetAtomIds()[0]
    atom_type = f.GetType()
    mol.GetAtomWithIdx(atom_idx).SetProp('AtomType', atom_type)

# 5. Print atom types for each atom
for atom in mol.GetAtoms():
    idx = atom.GetIdx()
    symbol = atom.GetSymbol()
    try:
        atype = atom.GetProp('AtomType')
    except KeyError:
        atype = "UNASSIGNED"
    print(f"Atom {idx+1} ({symbol}): {atype}")
