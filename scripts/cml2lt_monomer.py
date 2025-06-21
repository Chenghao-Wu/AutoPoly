import os
import glob
import pathlib
import xml.etree.ElementTree as ET
import argparse
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


def reorder_atoms(atoms, bonds, m=None, left_type=None, right_type=None):
    # Typical valence for common elements
    valence = {'C': 4, 'N': 3, 'O': 2, 'S': 2, 'P': 3, 'F': 1, 'Cl': 1, 'Br': 1, 'I': 1}
    # Build bond info for each atom
    bond_info = {a['id']: {'count': 0, 'order_sum': 0, 'has_double_or_triple': False} for a in atoms}
    for b in bonds:
        order = int(float(b['order']))
        bond_info[b['a1']]['count'] += 1
        bond_info[b['a2']]['count'] += 1
        bond_info[b['a1']]['order_sum'] += order
        bond_info[b['a2']]['order_sum'] += order
        if order > 1:
            bond_info[b['a1']]['has_double_or_triple'] = True
            bond_info[b['a2']]['has_double_or_triple'] = True
    # Identify heavy atoms (non-H)
    heavy_atoms = [(i, a) for i, a in enumerate(atoms) if a['element'] != 'H']
    if len(heavy_atoms) < 2:
        return atoms, bonds, {a['id']: i for i, a in enumerate(atoms)}
    # If left_type/right_type are given, use them to select ends
    def get_atom_type(idx):
        if m is not None:
            at = m.GetAtomWithIdx(idx)
            return (at.GetSymbol(), at.GetProp('AtomType') if at.HasProp('AtomType') else None)
        else:
            return (atoms[idx]['element'], None)
    if left_type and right_type:
        # Find all heavy atoms with matching type
        left_candidates = [i for i, a in heavy_atoms if get_atom_type(i) == left_type]
        right_candidates = [i for i, a in heavy_atoms if get_atom_type(i) == right_type]
        # Use x as tiebreaker
        left_idx = min(left_candidates, key=lambda i: atoms[i]['x']) if left_candidates else heavy_atoms[0][0]
        right_idx = max(right_candidates, key=lambda i: atoms[i]['x']) if right_candidates else heavy_atoms[-1][0]
    else:
        # Mark unsaturated heavy atoms
        unsat_heavy = []
        for i, a in heavy_atoms:
            v = valence.get(a['element'], 4)
            binfo = bond_info[a['id']]
            if binfo['order_sum'] < v:
                unsat_heavy.append((i, a))
        if unsat_heavy:
            left_idx, _ = min(unsat_heavy, key=lambda x: x[1]['x'])
            right_idx, _ = max(unsat_heavy, key=lambda x: x[1]['x'])
        else:
            left_idx, _ = min(heavy_atoms, key=lambda x: x[1]['x'])
            right_idx, _ = max(heavy_atoms, key=lambda x: x[1]['x'])
        print(left_idx, right_idx)
    if left_idx == right_idx:
        return atoms, bonds, {a['id']: i for i, a in enumerate(atoms)}
    
    new_order = [left_idx, right_idx] + [i for i in range(len(atoms)) if i not in (left_idx, right_idx)]
    new_atoms = [atoms[i] for i in new_order]
    new_atom_id_map = {a['id']: i for i, a in enumerate(new_atoms)}
    new_bonds = []
    for b in bonds:
        new_bonds.append({'a1': b['a1'], 'a2': b['a2'], 'order': b['order']})
    return new_atoms, new_bonds, new_atom_id_map


def get_end_atom_types(m, atoms, left_idx,right_idx):
    # Returns (left_type, right_type) as (element, atom_type) tuples
    heavy_idxs = [i for i, a in enumerate(atoms) if a['element'] != 'H']
    if not heavy_idxs:
        return None, None
    left_atom = m.GetAtomWithIdx(left_idx)
    right_atom = m.GetAtomWithIdx(right_idx)
    left_type = (left_atom.GetSymbol(), left_atom.GetProp('AtomType') if left_atom.HasProp('AtomType') else None)
    right_type = (right_atom.GetSymbol(), right_atom.GetProp('AtomType') if right_atom.HasProp('AtomType') else None)
    return left_type, right_type


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Convert CML files to LT files with consistent atom ordering')
    parser.add_argument('--internal', '-i', required=True, help='Internal monomer CML file (e.g., *i.cml)')
    parser.add_argument('--monomers', '-m', nargs='+', required=True, help='List of monomer CML files to process')
    args = parser.parse_args()

    # Process internal monomer first
    if not os.path.exists(args.internal):
        print(f'Error: Internal monomer file {args.internal} not found!')
        return
        
    if not args.internal.endswith('i.cml'):
        print(f'Error: Internal monomer file {args.internal} must end with "i.cml"!')
        return
        
    print(f'Processing {args.internal} (internal) ...')
    atoms, bonds, atom_id_map = parse_cml(args.internal)
    atoms, bonds, atom_id_map = reorder_atoms(atoms, bonds)
    m = build_rdkit_mol(atoms, bonds, atom_id_map)
    m = assign_atom_types(m)
    left_type, right_type = get_end_atom_types(m, atoms, 0, 1)
    outname = os.path.splitext(args.internal)[0] + '.lt'
    write_lt(m, atoms, bonds, atom_id_map, outname)
    print(f'Wrote {outname}')
    
    # Process other files using end types from *i.cml
    for cml in args.monomers:
        if not os.path.exists(cml):
            print(f'Warning: Monomer file {cml} not found, skipping...')
            continue
            
        print(f'Processing {cml} ...')
        atoms, bonds, atom_id_map = parse_cml(cml)
        m_tmp = build_rdkit_mol(atoms, bonds, atom_id_map)
        m_tmp = assign_atom_types(m_tmp)
        atoms, bonds, atom_id_map = reorder_atoms(atoms, bonds, m=m_tmp, left_type=left_type, right_type=right_type)
        m = build_rdkit_mol(atoms, bonds, atom_id_map)
        m = assign_atom_types(m)
        outname = os.path.splitext(cml)[0] + '.lt'
        write_lt(m, atoms, bonds, atom_id_map, outname)
        print(f'Wrote {outname}')

if __name__ == '__main__':
    main()
