#!/usr/bin/env python3
"""
CML to LT file converter for polymer monomers.

This module converts CML (Chemical Markup Language) files to LT (LAMMPS Template) files
for use with the AutoPoly polymerization system. It generates three monomer variants:
- Internal monomer (suffix 'i'): Two connection points
- Left-end monomer (suffix 'le'): One terminal end, one connection point
- Right-end monomer (suffix 're'): One connection point, one terminal end
"""

import sys
import os
import xml.etree.ElementTree as ET
from pathlib import Path
import argparse
from typing import Dict, List, Tuple, Optional

# RDKit imports
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures


class CMLToLTConverter:
    """Convert CML files to LT files for polymer monomers."""
    
    def __init__(self, opls_fdef_path: Optional[str] = None, lopls_fdef_path: Optional[str] = None):
        """
        Initialize the converter.
        
        Args:
            opls_fdef_path: Path to OPLS feature definition file
            lopls_fdef_path: Path to LOPLS feature definition file
        """
        self.opls_fdef_path = opls_fdef_path
        self.lopls_fdef_path = lopls_fdef_path
        
        # Atom type mapping for polymer monomers
        self.atom_type_mapping = {
            'internal': {'C1': '@atom:81', 'C2': '@atom:82'},  # CH2, CH2
            'left': {'C1': '@atom:80', 'C2': '@atom:82'},      # CH3, CH2
            'right': {'C1': '@atom:81', 'C2': '@atom:81'}      # CH2, CH2
        }
    
    def parse_cml_file(self, cml_file_path: str) -> Tuple[List[Dict], List[Dict], Dict]:
        """
        Parse CML file and extract atom and bond information.
        
        Args:
            cml_file_path: Path to the CML file
            
        Returns:
            Tuple of (atoms, bonds, atom_id_map)
            - atoms: List of atom dictionaries with id, element, x, y, z
            - bonds: List of bond dictionaries with a1, a2, order
            - atom_id_map: Mapping from atom id to index
        """
        try:
            tree = ET.parse(cml_file_path)
            root = tree.getroot()
            
            # Parse atoms
            atoms = []
            atom_id_map = {}
            
            # Find atomArray - try without namespace first
            atom_array = root.find('.//atomArray')
            if atom_array is None:
                raise ValueError("No atomArray found in CML file")
            
            for i, atom_elem in enumerate(atom_array.findall('.//atom')):
                atom_id = atom_elem.get('id')
                element = atom_elem.get('elementType')
                x = float(atom_elem.get('x3', 0))
                y = float(atom_elem.get('y3', 0))
                z = float(atom_elem.get('z3', 0))
                
                atom_data = {
                    'id': atom_id,
                    'element': element,
                    'x': x,
                    'y': y,
                    'z': z,
                    'index': i
                }
                atoms.append(atom_data)
                atom_id_map[atom_id] = i
            
            # Parse bonds
            bonds = []
            bond_array = root.find('.//bondArray')
            if bond_array is not None:
                for bond_elem in bond_array.findall('.//bond'):
                    atom_refs = bond_elem.get('atomRefs2', '').split()
                    if len(atom_refs) == 2:
                        bond_data = {
                            'a1': atom_refs[0],
                            'a2': atom_refs[1],
                            'order': bond_elem.get('order', '1')
                        }
                        bonds.append(bond_data)
            
            return atoms, bonds, atom_id_map
            
        except Exception as e:
            raise ValueError(f"Error parsing CML file {cml_file_path}: {str(e)}")
    
    def build_rdkit_molecule(self, atoms: List[Dict], bonds: List[Dict], atom_id_map: Dict) -> Chem.Mol:
        """
        Build RDKit molecule from parsed CML data.
        
        Args:
            atoms: List of atom dictionaries
            bonds: List of bond dictionaries
            atom_id_map: Mapping from atom id to index
            
        Returns:
            RDKit molecule with 3D coordinates
        """
        mol = Chem.RWMol()
        
        # Add atoms
        for atom in atoms:
            rdkit_atom = Chem.Atom(atom['element'])
            mol.AddAtom(rdkit_atom)
        
        # Add bonds
        for bond in bonds:
            a1_idx = atom_id_map[bond['a1']]
            a2_idx = atom_id_map[bond['a2']]
            order = int(bond['order'])
            
            if order == 1:
                bond_type = Chem.BondType.SINGLE
            elif order == 2:
                bond_type = Chem.BondType.DOUBLE
            elif order == 3:
                bond_type = Chem.BondType.TRIPLE
            else:
                bond_type = Chem.BondType.SINGLE
            
            mol.AddBond(a1_idx, a2_idx, bond_type)
        
        # Convert to molecule and add conformer
        m = mol.GetMol()
        conf = Chem.Conformer(m.GetNumAtoms())
        
        for atom in atoms:
            conf.SetAtomPosition(atom['index'], (atom['x'], atom['y'], atom['z']))
        
        m.AddConformer(conf, assignId=True)
        
        # Sanitize molecule
        Chem.SanitizeMol(m)
        
        return m
    
    def identify_connection_carbons(self, mol: Chem.Mol) -> Tuple[int, int]:
        """
        Identify the two carbons that will serve as connection points.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Tuple of (c1_idx, c2_idx) - indices of connection carbons
        """
        carbon_atoms = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C':
                carbon_atoms.append(atom.GetIdx())
        
        if len(carbon_atoms) < 2:
            raise ValueError("Molecule must have at least 2 carbon atoms for polymer building")
        
        # Sort by atom index for consistent ordering
        carbon_atoms.sort()
        return carbon_atoms[0], carbon_atoms[1]
    
    def assign_atom_types(self, mol: Chem.Mol, monomer_type: str, c1_idx: int, c2_idx: int) -> Chem.Mol:
        """
        Assign atom types based on monomer variant.
        
        Args:
            mol: RDKit molecule
            monomer_type: 'internal', 'left', or 'right'
            c1_idx: Index of first connection carbon
            c2_idx: Index of second connection carbon
            
        Returns:
            Molecule with updated atom types
        """
        # First assign default atom types using feature factory
        if self.opls_fdef_path and os.path.exists(self.opls_fdef_path):
            factory = ChemicalFeatures.BuildFeatureFactory(self.opls_fdef_path)
            features = factory.GetFeaturesForMol(mol)
            
            # Assign atom types from features
            for feature in features:
                for atom_id in feature.GetAtomIds():
                    mol.GetAtomWithIdx(atom_id).SetProp('AtomType', feature.GetType())
        
        # Override connection carbons based on monomer type
        if monomer_type not in self.atom_type_mapping:
            raise ValueError(f"Unknown monomer_type: {monomer_type}")
        
        type_map = self.atom_type_mapping[monomer_type]
        mol.GetAtomWithIdx(c1_idx).SetProp('AtomType', type_map['C1'])
        mol.GetAtomWithIdx(c2_idx).SetProp('AtomType', type_map['C2'])
        
        # Assign default types to other atoms if not already assigned
        for atom in mol.GetAtoms():
            if not atom.HasProp('AtomType'):
                # Simple default assignment based on element and connectivity
                if atom.GetSymbol() == 'H':
                    atom.SetProp('AtomType', '@atom:85')  # H
                elif atom.GetSymbol() == 'C':
                    # Count hydrogens for carbon typing
                    h_count = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == 'H')
                    if h_count == 3:
                        atom.SetProp('AtomType', '@atom:80')  # CH3
                    elif h_count == 2:
                        atom.SetProp('AtomType', '@atom:81')  # CH2
                    elif h_count == 1:
                        atom.SetProp('AtomType', '@atom:82')  # CH
                    else:
                        atom.SetProp('AtomType', '@atom:83')  # C
        
        return mol
    
    def write_lt_header(self, monomer_name: str, loplsflag: bool = False) -> str:
        """Generate LT file header."""
        header = f'import "oplsaa.lt"    # <-- defines the standard "OPLSAA" force field\n'
        if loplsflag:
            header += '''import "loplsaa.lt"   # <-- custom parameters for long alkane chains taken from
                      #     Sui et al. J.Chem.Theory.Comp (2012), 8, 1459
                      #     To use the ordinary OPLSAA force field parameters,
                      #     (instead of the Sui et al. parameters), change the
                      #     atom types below from "@atom:81L","@atom:85LCH2" to
                      #     "@atom:81" and "@atom:85"  (defined in "oplsaa.lt")\n'''
        header += f'{monomer_name} inherits OPLSAA {{\n'
        return header
    
    def write_lt_atoms(self, mol: Chem.Mol) -> str:
        """Generate LT file atoms section."""
        atoms_section = '\n# atom-id  mol-id  atom-type charge      X         Y        Z\n'
        atoms_section += '  write("Data Atoms") {\n'
        
        conf = mol.GetConformer(0)
        for atom in mol.GetAtoms():
            point = conf.GetAtomPosition(atom.GetIdx())
            atom_type = atom.GetProp('AtomType')
            atoms_section += f'\t$atom:{atom.GetSymbol()}{atom.GetIdx()+1} $mol:... {atom_type} 0.00 {point.x:8.3f} {point.y:8.3f} {point.z:8.3f}\n'
        
        atoms_section += '  }\n'
        return atoms_section
    
    def write_lt_bonds(self, mol: Chem.Mol) -> str:
        """Generate LT file bonds section."""
        bonds_section = '\n  write(\'Data Bond List\') {\n'
        
        for bond in mol.GetBonds():
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            begin_name = f"{begin_atom.GetSymbol()}{begin_atom.GetIdx()+1}"
            end_name = f"{end_atom.GetSymbol()}{end_atom.GetIdx()+1}"
            bonds_section += f'\t$bond:{begin_name}{end_name}\t$atom:{begin_name}\t$atom:{end_name}\n'
        
        bonds_section += '  }\n'
        return bonds_section
    
    def write_lt_footer(self, monomer_name: str) -> str:
        """Generate LT file footer."""
        footer = f'}} # {monomer_name}\n\n'
        footer += '# Note: You don\'t need to supply the partial partial charges of the atoms.\n'
        footer += '#       If you like, just fill the fourth column with zeros ("0.000").\n'
        footer += '#       Moltemplate and LAMMPS will automatically assign the charge later\n'
        return footer
    
    def generate_lt_file(self, mol: Chem.Mol, monomer_name: str, output_file: str, loplsflag: bool = False) -> str:
        """
        Generate complete LT file.
        
        Args:
            mol: RDKit molecule with atom types assigned
            monomer_name: Name for the monomer
            output_file: Output file path
            loplsflag: Whether to use LOPLS atom typing
            
        Returns:
            Path to generated file
        """
        # Validate that all atoms have types
        for atom in mol.GetAtoms():
            if not atom.HasProp('AtomType'):
                raise ValueError(f"Atom {atom.GetIdx()} does not have an assigned atom type!")
        
        # Generate LT file content
        content = self.write_lt_header(monomer_name, loplsflag)
        content += self.write_lt_atoms(mol)
        content += self.write_lt_bonds(mol)
        content += self.write_lt_footer(monomer_name)
        
        # Write to file
        with open(output_file, 'w') as f:
            f.write(content)
        
        return output_file
    
    def process_cml_to_lt(self, cml_file_path: str, monomer_name: str, monomer_type: str, 
                         output_file: str, loplsflag: bool = False) -> str:
        """
        Process a single CML file to LT format.
        
        Args:
            cml_file_path: Path to input CML file
            monomer_name: Name for the monomer
            monomer_type: 'internal', 'left', or 'right'
            output_file: Path to output LT file
            loplsflag: Whether to use LOPLS atom typing
            
        Returns:
            Path to generated LT file
        """
        # Parse CML file
        atoms, bonds, atom_id_map = self.parse_cml_file(cml_file_path)
        
        # Build RDKit molecule
        mol = self.build_rdkit_molecule(atoms, bonds, atom_id_map)
        
        # Identify connection carbons
        c1_idx, c2_idx = self.identify_connection_carbons(mol)
        
        # Assign atom types
        mol = self.assign_atom_types(mol, monomer_type, c1_idx, c2_idx)
        
        # Generate LT file
        return self.generate_lt_file(mol, monomer_name, output_file, loplsflag)
    
    def process_polymer_monomers(self, cml_files: Dict[str, str], base_name: str, 
                                output_dir: str = ".", loplsflag: bool = False) -> Dict[str, str]:
        """
        Process all three CML files for polymer monomers.
        
        Args:
            cml_files: Dictionary mapping monomer types to CML file paths
            base_name: Base name for the monomer (e.g., 'PE', 'PS')
            output_dir: Output directory
            loplsflag: Whether to use LOPLS atom typing
            
        Returns:
            Dictionary mapping monomer types to generated LT file paths
        """
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Define monomer variants
        variants = {
            'internal': f"{base_name}i",
            'left': f"{base_name}le",
            'right': f"{base_name}re"
        }
        
        generated_files = {}
        
        for monomer_type, cml_file in cml_files.items():
            if monomer_type not in variants:
                print(f"Warning: Unknown monomer type '{monomer_type}', skipping...")
                continue
            
            if not os.path.exists(cml_file):
                print(f"Warning: CML file '{cml_file}' not found, skipping...")
                generated_files[monomer_type] = None
                continue
            
            monomer_name = variants[monomer_type]
            output_file = os.path.join(output_dir, f"{monomer_name}.lt")
            
            try:
                generated_file = self.process_cml_to_lt(
                    cml_file, monomer_name, monomer_type, output_file, loplsflag
                )
                generated_files[monomer_type] = generated_file
                print(f"Generated {monomer_type} monomer: {generated_file}")
                
            except Exception as e:
                print(f"Error generating {monomer_type} monomer: {e}")
                generated_files[monomer_type] = None
        
        return generated_files


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(description='Convert CML files to LT files for polymer monomers')
    parser.add_argument('--internal', '-i', required=True, help='Internal monomer CML file (e.g., *i.cml)')
    parser.add_argument('--left', '-l', required=True, help='Left-end monomer CML file (e.g., *l.cml)')
    parser.add_argument('--right', '-r', required=True, help='Right-end monomer CML file (e.g., *r.cml)')
    parser.add_argument('--base-name', '-b', required=True, help='Base name for the monomer (e.g., PE, PS)')
    parser.add_argument('--output-dir', '-o', default='.', help='Output directory for LT files')
    parser.add_argument('--opls-fdef', help='Path to OPLS feature definition file')
    parser.add_argument('--lopls-fdef', help='Path to LOPLS feature definition file')
    parser.add_argument('--lopls', action='store_true', help='Use LOPLS atom typing')
    
    args = parser.parse_args()
    
    # Set up converter
    converter = CMLToLTConverter(args.opls_fdef, args.lopls_fdef)
    
    # Define CML files
    cml_files = {
        'internal': args.internal,
        'left': args.left,
        'right': args.right
    }
    
    # Process files
    try:
        generated_files = converter.process_polymer_monomers(
            cml_files, args.base_name, args.output_dir, args.lopls
        )
        
        print("\nProcessing completed!")
        print("Generated files:")
        for monomer_type, file_path in generated_files.items():
            if file_path:
                print(f"  {monomer_type}: {file_path}")
            else:
                print(f"  {monomer_type}: FAILED")
                
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
