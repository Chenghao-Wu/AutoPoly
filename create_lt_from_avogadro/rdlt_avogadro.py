#!/usr/bin/env python3
"""
CML to LT file converter for polymer monomers with tacticity support.

This module converts CML (Chemical Markup Language) files to LT (LAMMPS Template) files
for use with the AutoPoly polymerization system. It generates monomer variants with
tacticity support for different stereochemical configurations.

Tacticity Types:
- Atactic: Random stereochemistry
- Isotactic: All monomers have the same stereochemistry  
- Syndiotactic: Alternating stereochemistry

Monomer Variants:
- Internal monomer (suffix 'i'): Two connection points
- Left-end monomer (suffix 'le'): One terminal end, one connection point
- Right-end monomer (suffix 're'): One connection point, one terminal end

Stereochemical Variants:
- Normal: Standard configuration
- T1: Stereochemical variant (180° rotation around backbone)
"""

import sys
import os
import random
import xml.etree.ElementTree as ET
from pathlib import Path
import argparse
from typing import Dict, List, Tuple, Optional

# RDKit imports
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
import numpy as np


class CMLToLTConverterTacticity:
    """Convert CML files to LT files for polymer monomers with tacticity support."""
    
    def __init__(self, opls_fdef_path: Optional[str] = None, lopls_fdef_path: Optional[str] = None,
                 atom_type_mapping: Optional[Dict] = None):
        """
        Initialize the converter.
        
        Args:
            opls_fdef_path: Path to OPLS feature definition file
            lopls_fdef_path: Path to LOPLS feature definition file
            atom_type_mapping: Custom atom type mapping dictionary
        """
        self.opls_fdef_path = opls_fdef_path
        self.lopls_fdef_path = lopls_fdef_path
        
        # Default atom type mapping for polymer monomers
        default_atom_type_mapping = {
            'internal': {'atom1': '@atom:81', 'atom2': '@atom:82'},  # CH2, CH2
            'left': {'atom1': '@atom:80', 'atom2': '@atom:82'},      # CH3, CH2
            'right': {'atom1': '@atom:81', 'atom2': '@atom:81'}      # CH2, CH2
        }
        
        # Use custom mapping if provided, otherwise use default
        if atom_type_mapping is not None:
            # Validate custom mapping
            self.validate_atom_type_mapping(atom_type_mapping)
            self.atom_type_mapping = atom_type_mapping
        else:
            self.atom_type_mapping = default_atom_type_mapping
        
        # Tacticity options
        self.tacticity_options = ['atactic', 'isotactic', 'syndiotactic']
        
        # Stereochemical transformation matrix for T1 variant (180° rotation around backbone)
        self.t1_transformation = np.array([
            [-1, 0, 0],  # Flip x coordinates
            [0, -1, 0],  # Flip y coordinates  
            [0, 0, 1]    # Keep z coordinates
        ])
    
    def validate_atom_type_mapping(self, atom_type_mapping: Dict) -> None:
        """
        Validate the provided atom type mapping.
        
        Args:
            atom_type_mapping: Dictionary containing atom type mappings
            
        Raises:
            ValueError: If the mapping is invalid
        """
        required_keys = ['internal', 'left', 'right']
        
        # Check if all required monomer types are present
        for key in required_keys:
            if key not in atom_type_mapping:
                raise ValueError(f"Missing required monomer type '{key}' in atom_type_mapping")
            
            # Check that each monomer type has at least 2 atoms defined
            atoms = atom_type_mapping[key]
            if not isinstance(atoms, dict) or len(atoms) < 2:
                raise ValueError(f"Monomer type '{key}' must have at least 2 atoms defined")
            
            # Validate atom type format for all atoms
            for atom_name, atom_type in atoms.items():
                if not isinstance(atom_type, str) or not atom_type.startswith('@atom:'):
                    raise ValueError(f"Invalid atom type format '{atom_type}' for {key}.{atom_name}. Must start with '@atom:'")
    
    def set_atom_type_mapping(self, atom_type_mapping: Dict) -> None:
        """
        Set custom atom type mapping after initialization.
        
        Args:
            atom_type_mapping: Dictionary containing atom type mappings
        """
        self.validate_atom_type_mapping(atom_type_mapping)
        self.atom_type_mapping = atom_type_mapping
    
    def get_atom_type_mapping(self) -> Dict:
        """
        Get the current atom type mapping.
        
        Returns:
            Dictionary containing current atom type mappings
        """
        return self.atom_type_mapping.copy()
    
    def print_atom_type_mapping(self) -> None:
        """Print the current atom type mapping in a readable format."""
        print("Current Atom Type Mapping:")
        print("=" * 50)
        for monomer_type, atoms in self.atom_type_mapping.items():
            print(f"{monomer_type.upper()} monomer:")
            # Sort atoms by name for consistent output
            for atom_name in sorted(atoms.keys()):
                print(f"  {atom_name}: {atoms[atom_name]}")
            print()
    
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
    
    def apply_stereochemical_transformation(self, atoms: List[Dict], variant: str = 'normal') -> List[Dict]:
        """
        Apply stereochemical transformation to atom coordinates.
        
        Args:
            atoms: List of atom dictionaries
            variant: 'normal' or 'T1' for stereochemical variant
            
        Returns:
            List of transformed atom dictionaries
        """
        if variant == 'normal':
            return atoms
        
        elif variant == 'T1':
            # Apply T1 transformation (180° rotation around backbone)
            transformed_atoms = []
            
            for atom in atoms:
                # Get original coordinates
                coords = np.array([atom['x'], atom['y'], atom['z']])
                
                # Apply transformation
                transformed_coords = np.dot(self.t1_transformation, coords)
                
                # Create new atom data with transformed coordinates
                transformed_atom = atom.copy()
                transformed_atom['x'] = transformed_coords[0]
                transformed_atom['y'] = transformed_coords[1]
                transformed_atom['z'] = transformed_coords[2]
                
                transformed_atoms.append(transformed_atom)
            
            return transformed_atoms
        
        else:
            raise ValueError(f"Unknown stereochemical variant: {variant}")
    
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
    
    def identify_connection_carbons(self, mol: Chem.Mol) -> List[int]:
        """
        Identify all carbon atoms that will serve as connection points.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            List of indices of connection carbon atoms
        """
        carbon_atoms = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C':
                carbon_atoms.append(atom.GetIdx())
        
        if len(carbon_atoms) < 2:
            raise ValueError("Molecule must have at least 2 carbon atoms for polymer building")
        
        # Sort by atom index for consistent ordering
        carbon_atoms.sort()
        return carbon_atoms
    
    def assign_atom_types(self, mol: Chem.Mol, monomer_type: str, connection_indices: List[int]) -> Chem.Mol:
        """
        Assign atom types based on monomer variant.
        
        Args:
            mol: RDKit molecule
            monomer_type: 'internal', 'left', or 'right'
            connection_indices: List of indices of connection atoms (backbone atoms)
            
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
        
        # Override connection atoms based on monomer type
        if monomer_type not in self.atom_type_mapping:
            raise ValueError(f"Unknown monomer_type: {monomer_type}")
        
        type_map = self.atom_type_mapping[monomer_type]
        
        # Get the atom names from the type map (sorted for consistent ordering)
        atom_names = sorted(type_map.keys())
        
        # Assign types to connection atoms
        for i, atom_idx in enumerate(connection_indices):
            if i < len(atom_names):
                atom_name = atom_names[i]
                mol.GetAtomWithIdx(atom_idx).SetProp('AtomType', type_map[atom_name])
            else:
                # If we have more connection atoms than defined types, use the last defined type
                last_atom_name = atom_names[-1]
                mol.GetAtomWithIdx(atom_idx).SetProp('AtomType', type_map[last_atom_name])
        
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
    
    def process_cml_to_lt_with_tacticity(self, cml_file_path: str, monomer_name: str, monomer_type: str, 
                                       output_file: str, tacticity: str = 'atactic', 
                                       stereochemical_variant: str = 'normal', loplsflag: bool = False) -> str:
        """
        Process a single CML file to LT format with tacticity support.
        
        Args:
            cml_file_path: Path to input CML file
            monomer_name: Name for the monomer
            monomer_type: 'internal', 'left', or 'right'
            output_file: Path to output LT file
            tacticity: 'atactic', 'isotactic', or 'syndiotactic'
            stereochemical_variant: 'normal' or 'T1'
            loplsflag: Whether to use LOPLS atom typing
            
        Returns:
            Path to generated LT file
        """
        # Parse CML file
        atoms, bonds, atom_id_map = self.parse_cml_file(cml_file_path)
        
        # Apply stereochemical transformation
        atoms = self.apply_stereochemical_transformation(atoms, stereochemical_variant)
        
        # Build RDKit molecule
        mol = self.build_rdkit_molecule(atoms, bonds, atom_id_map)
        
        # Identify connection carbons
        connection_carbons = self.identify_connection_carbons(mol)
        
        # Assign atom types
        mol = self.assign_atom_types(mol, monomer_type, connection_carbons)
        
        # Generate LT file
        return self.generate_lt_file(mol, monomer_name, output_file, loplsflag)
    
    def generate_tacticity_variants(self, cml_files: Dict[str, str], base_name: str, 
                                  output_dir: str = ".", tacticity: str = 'atactic', 
                                  loplsflag: bool = False) -> Dict[str, List[str]]:
        """
        Generate monomer variants with tacticity support.
        
        Args:
            cml_files: Dictionary mapping monomer types to CML file paths
            base_name: Base name for the monomer (e.g., 'PE', 'PS')
            output_dir: Output directory
            tacticity: 'atactic', 'isotactic', or 'syndiotactic'
            loplsflag: Whether to use LOPLS atom typing
            
        Returns:
            Dictionary mapping monomer types to lists of generated file paths
        """
        if tacticity not in self.tacticity_options:
            raise ValueError(f"Unknown tacticity: {tacticity}. Must be one of {self.tacticity_options}")
        
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
                generated_files[monomer_type] = []
                continue
            
            monomer_name = variants[monomer_type]
            generated_files[monomer_type] = []
            
            # Generate variants based on tacticity
            if tacticity == 'atactic':
                # For atactic, generate both normal and T1 variants
                for variant in ['normal', 'T1']:
                    suffix = '_T1' if variant == 'T1' else ''
                    output_file = os.path.join(output_dir, f"{monomer_name}{suffix}.lt")
                    
                    try:
                        generated_file = self.process_cml_to_lt_with_tacticity(
                            cml_file, f"{monomer_name}{suffix}", monomer_type, 
                            output_file, tacticity, variant, loplsflag
                        )
                        generated_files[monomer_type].append(generated_file)
                        print(f"Generated {monomer_type} monomer ({variant}): {generated_file}")
                    except Exception as e:
                        print(f"Error generating {monomer_type} monomer ({variant}): {e}")
            
            elif tacticity == 'isotactic':
                # For isotactic, generate only normal variant
                output_file = os.path.join(output_dir, f"{monomer_name}.lt")
                
                try:
                    generated_file = self.process_cml_to_lt_with_tacticity(
                        cml_file, monomer_name, monomer_type, 
                        output_file, tacticity, 'normal', loplsflag
                    )
                    generated_files[monomer_type].append(generated_file)
                    print(f"Generated {monomer_type} monomer (isotactic): {generated_file}")
                except Exception as e:
                    print(f"Error generating {monomer_type} monomer (isotactic): {e}")
            
            elif tacticity == 'syndiotactic':
                # For syndiotactic, generate both variants (will be used alternately)
                for variant in ['normal', 'T1']:
                    suffix = '_T1' if variant == 'T1' else ''
                    output_file = os.path.join(output_dir, f"{monomer_name}{suffix}.lt")
                    
                    try:
                        generated_file = self.process_cml_to_lt_with_tacticity(
                            cml_file, f"{monomer_name}{suffix}", monomer_type, 
                            output_file, tacticity, variant, loplsflag
                        )
                        generated_files[monomer_type].append(generated_file)
                        print(f"Generated {monomer_type} monomer ({variant}): {generated_file}")
                    except Exception as e:
                        print(f"Error generating {monomer_type} monomer ({variant}): {e}")
        
        return generated_files


def parse_atom_type_mapping(mapping_string: str) -> Dict:
    """
    Parse atom type mapping from command line string.
    
    Expected format: "internal:atom1@atom:81,atom2@atom:82;left:atom1@atom:80,atom2@atom:82;right:atom1@atom:81,atom2@atom:81"
    Atom names can be any identifier (e.g., C1, C2, Si1, Si2, etc.)
    
    Args:
        mapping_string: String containing atom type mapping
        
    Returns:
        Dictionary containing parsed atom type mapping
        
    Raises:
        ValueError: If the mapping string format is invalid
    """
    atom_type_mapping = {}
    
    try:
        # Split by monomer types
        monomer_types = mapping_string.split(';')
        
        for monomer_type_str in monomer_types:
            if not monomer_type_str.strip():
                continue
            # DEBUG: print the string being parsed
            # print(f"Parsing monomer_type_str: '{monomer_type_str}'")
            # Split monomer type and atoms (only at the first colon)
            parts = monomer_type_str.split(':', 1)
            if len(parts) != 2:
                raise ValueError(f"Invalid format in '{monomer_type_str}'. Expected 'type:atom1@atom:XX,atom2@atom:XX'")
            monomer_type = parts[0].strip()
            atoms_str = parts[1].strip()
            # Parse atoms
            atoms = {}
            atom_pairs = atoms_str.split(',')
            for atom_pair in atom_pairs:
                if not atom_pair.strip():
                    continue
                if '@' not in atom_pair:
                    raise ValueError(f"Invalid atom format in '{atom_pair}'. Expected 'atomname@atom:XX'")
                atom_name, atom_type = atom_pair.split('@', 1)
                atom_name = atom_name.strip()
                atom_type = f"@{atom_type.strip()}"
                atoms[atom_name] = atom_type
            atom_type_mapping[monomer_type] = atoms
    
    except Exception as e:
        raise ValueError(f"Error parsing atom type mapping: {str(e)}")
    
    return atom_type_mapping


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(description='Convert CML files to LT files for polymer monomers with tacticity support')
    parser.add_argument('--internal', '-i', help='Internal monomer CML file (e.g., *i.cml)')
    parser.add_argument('--left', '-l', help='Left-end monomer CML file (e.g., *l.cml)')
    parser.add_argument('--right', '-r', help='Right-end monomer CML file (e.g., *r.cml)')
    parser.add_argument('--base-name', '-b', help='Base name for the monomer (e.g., PE, PS)')
    parser.add_argument('--output-dir', '-o', default='.', help='Output directory for LT files')
    parser.add_argument('--tacticity', '-t', choices=['atactic', 'isotactic', 'syndiotactic'], 
                       default='atactic', help='Polymer tacticity')
    parser.add_argument('--opls-fdef', help='Path to OPLS feature definition file')
    parser.add_argument('--lopls-fdef', help='Path to LOPLS feature definition file')
    parser.add_argument('--lopls', action='store_true', help='Use LOPLS atom typing')
    parser.add_argument('--atom-type-mapping', '-m', help='Custom atom type mapping. Format: "internal:atom1@atom:81,atom2@atom:82;left:atom1@atom:80,atom2@atom:82;right:atom1@atom:81,atom2@atom:81" (atom names can be any identifier)')
    parser.add_argument('--show-mapping', action='store_true', help='Show current atom type mapping and exit')
    
    args = parser.parse_args()
    
    # Show mapping if requested (no other arguments needed)
    if args.show_mapping:
        # Parse custom atom type mapping if provided
        atom_type_mapping = None
        if args.atom_type_mapping:
            try:
                atom_type_mapping = parse_atom_type_mapping(args.atom_type_mapping)
                print("Using custom atom type mapping:")
                for monomer_type, atoms in atom_type_mapping.items():
                    atom_list = [f"{name}={atom_type}" for name, atom_type in sorted(atoms.items())]
                    print(f"  {monomer_type}: {', '.join(atom_list)}")
                print()
            except ValueError as e:
                print(f"Error: {e}")
                sys.exit(1)
        
        converter = CMLToLTConverterTacticity(atom_type_mapping=atom_type_mapping)
        converter.print_atom_type_mapping()
        sys.exit(0)
    
    # Check required arguments for normal operation
    if not all([args.internal, args.left, args.right, args.base_name]):
        parser.error("--internal, --left, --right, and --base-name are required for conversion (unless using --show-mapping)")
    
    # Parse custom atom type mapping if provided
    atom_type_mapping = None
    if args.atom_type_mapping:
        try:
            atom_type_mapping = parse_atom_type_mapping(args.atom_type_mapping)
            print("Using custom atom type mapping:")
            for monomer_type, atoms in atom_type_mapping.items():
                atom_list = [f"{name}={atom_type}" for name, atom_type in sorted(atoms.items())]
                print(f"  {monomer_type}: {', '.join(atom_list)}")
            print()
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)
    
    # Set up converter
    converter = CMLToLTConverterTacticity(args.opls_fdef, args.lopls_fdef, atom_type_mapping)
    
    # Define CML files
    cml_files = {
        'internal': args.internal,
        'left': args.left,
        'right': args.right
    }
    
    # Process files with tacticity
    try:
        generated_files = converter.generate_tacticity_variants(
            cml_files, args.base_name, args.output_dir, args.tacticity, args.lopls
        )
        
        print(f"\nProcessing completed for {args.tacticity} tacticity!")
        print("Generated files:")
        for monomer_type, file_list in generated_files.items():
            if file_list:
                print(f"  {monomer_type}:")
                for file_path in file_list:
                    print(f"    {file_path}")
            else:
                print(f"  {monomer_type}: FAILED")
                
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
