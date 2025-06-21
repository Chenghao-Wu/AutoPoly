import os
import pathlib
import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures

# Get path to opls_lt.fdefn using AutoPoly package
import AutoPoly
FDEF = str(pathlib.Path(AutoPoly.__file__).parent / 'extern' / 'rdlt_data' / 'opls_lt.fdefn')

if not os.path.exists(FDEF):
    raise FileNotFoundError(f"Could not find opls_lt.fdefn at {FDEF}. Please ensure the file exists in the correct location.")

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

def assign_atom_types(m, fdef=FDEF):
    factory = ChemicalFeatures.BuildFeatureFactory(fdef)
    features = factory.GetFeaturesForMol(m)
    for f in features:
        m.GetAtomWithIdx(f.GetAtomIds()[0]).SetProp('AtomType', f.GetType())
    # Check all atoms have types
    for at in m.GetAtoms():
        try:
            at.GetProp('AtomType')
        except KeyError:
            raise RuntimeError(f"Atom {at.GetIdx()} does not have an assigned atom type!")
    return m

def write_lt(m, atoms, bonds, atom_id_map, outname):
    molname = pathlib.Path(outname).stem
    with open(outname, 'w') as f:
        f.write('import "oplsaa.lt"\n\n')
        f.write(f'{molname} inherits OPLSAA {{\n\n')
        f.write('    # atom-id mol-id atom-type charge X Y Z # comments\n\n')
        f.write('   write("Data Atoms") {\n')
        conf = m.GetConformer(0)
        for i, at in enumerate(m.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            symbol = at.GetSymbol()
            atom_id = f'{symbol}{i+1}'
            atom_type = at.GetProp('AtomType')
            f.write(f'        $atom:{atom_id}  $mol:...  {atom_type}  0.00  {pos.x:8.6f}  {pos.y:8.6f}  {pos.z:8.6f}  # {symbol}\n')
        f.write('    }\n\n')
        f.write('    write("Data Bond List") {\n')
        for j, bond in enumerate(bonds):
            a1_idx = atom_id_map[bond['a1']]
            a2_idx = atom_id_map[bond['a2']]
            a1 = m.GetAtomWithIdx(a1_idx)
            a2 = m.GetAtomWithIdx(a2_idx)
            a1id = f'{a1.GetSymbol()}{a1_idx+1}'
            a2id = f'{a2.GetSymbol()}{a2_idx+1}'
            f.write(f'        $bond:{j+1}  $atom:{a1id}  $atom:{a2id}  # order:{bond["order"]}\n')
        f.write('    }\n')
        f.write('}\n')

if __name__ == '__main__':
    cml_file = 'PLA.cml'  # Change this path if needed
    outname = os.path.splitext(cml_file)[0] + '.lt'
    atoms, bonds, atom_id_map = parse_cml(cml_file)
    m = build_rdkit_mol(atoms, bonds, atom_id_map)
    m = assign_atom_types(m)
    write_lt(m, atoms, bonds, atom_id_map, outname)
    print(f'Wrote {outname}')
