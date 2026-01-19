#!/usr/bin/env python3
"""
Standalone Monomer Generator for AutoPoly

Generates polymer monomer templates from SMILES/chain molecules with
automatic force field typing and LT file generation.

This is a COMPLETELY STANDALONE implementation that:
- Uses only RDKit for chemistry operations
- Parses RDlt .fdefn files directly for SMARTS patterns
- Does NOT depend on any AutoPoly internal classes

Key design decisions:
- Atom typing happens on the CHAIN (before splitting) for correct chemical environment
- Uses atom map numbers to track connection points through all operations
- Connection atoms are placed FIRST in LT file for AutoPoly compatibility
"""

import os
import re
import pickle
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, ChemicalFeatures

# Setup logging
logger = logging.getLogger(__name__)

# =============================================================================
# Constants for atom map tracking
# =============================================================================
LEFT_CONN_MAP = 1      # Map number for left connection atom
RIGHT_CONN_MAP = 2     # Map number for right connection atom
INTER_BOND_MAP_START = 100  # Starting map number for inter-monomer bond markers
CHARGE_MAP_START = 10000   # Starting map number for Gasteiger charge tracking (avoid conflicts)
CAP_H_MAP_START = 20000    # Starting map number for cap H atoms (chain ends)
TERMINAL_DUMMY_ISOTOPE = 99  # Isotope to mark terminal dummy atoms (chain ends)

# =============================================================================
# Exceptions
# =============================================================================

class MonomerGeneratorError(Exception):
    """Base exception for monomer_generator module."""
    pass


class AtomTypingError(MonomerGeneratorError):
    """Raised when atom typing fails."""
    pass


class ChainBuildingError(MonomerGeneratorError):
    """Raised when chain building fails."""
    pass


class ChainSplittingError(MonomerGeneratorError):
    """Raised when chain splitting fails."""
    pass


class LTWritingError(MonomerGeneratorError):
    """Raised when LT file generation fails."""
    pass


class ValidationError(MonomerGeneratorError):
    """Raised when validation fails."""
    pass


# =============================================================================
# MonomerVariant Dataclass
# =============================================================================

@dataclass
class MonomerVariant:
    """
    Immutable container for monomer variant data.
    
    Attributes:
        base_name: Base name for monomer (e.g., "PE", "PMMA")
        variant_type: 'first', 'middle', 'last', 'single', or 'ring'
        mol: RDKit Mol object with conformer and atom types
        smiles: SMILES representation
        connection_atoms: Tuple of (left_idx, right_idx) for polymerization
        force_field: 'oplsaa' or 'gaff'
        position: Position in original chain (0-based)
        atom_ids: Mapping from atom index to LT atom ID (e.g., {0: 'C1', 1: 'H2'})
    """
    base_name: str
    variant_type: str
    mol: Chem.Mol
    smiles: str
    connection_atoms: Tuple[int, int]
    force_field: str
    position: int = 0
    atom_ids: Dict[int, str] = field(default_factory=dict)
    is_t1: bool = False
    
    def get_atom_type_at_connection(self, side: str) -> Optional[str]:
        """Get atom type at left or right connection point."""
        side_idx = 0 if side == 'left' else 1
        idx = self.connection_atoms[side_idx]
        atom = self.mol.GetAtomWithIdx(idx)
        try:
            return atom.GetProp('AtomType')
        except KeyError:
            return None


# =============================================================================
# AtomMapTracker - Track atoms through operations using map numbers
# =============================================================================

class AtomMapTracker:
    """
    Track atoms through operations using RDKit atom map numbers.
    
    Atom map numbers survive atom removal, reordering, and fragmentation,
    making them ideal for tracking connection points.
    """
    
    @staticmethod
    def mark_connection_atoms(mol: Chem.Mol, left_idx: int, right_idx: int) -> None:
        """
        Mark connection atoms with map numbers.
        
        Args:
            mol: RDKit Mol object (modified in-place)
            left_idx: Index of left connection atom
            right_idx: Index of right connection atom
        """
        mol.GetAtomWithIdx(left_idx).SetAtomMapNum(LEFT_CONN_MAP)
        mol.GetAtomWithIdx(right_idx).SetAtomMapNum(RIGHT_CONN_MAP)
    
    @staticmethod
    def find_by_map_num(mol: Chem.Mol, map_num: int) -> Optional[int]:
        """
        Find atom index by map number.
        
        Args:
            mol: RDKit Mol object
            map_num: Map number to search for
            
        Returns:
            Atom index if found, None otherwise
        """
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == map_num:
                return atom.GetIdx()
        return None
    
    @staticmethod
    def get_connection_atoms(mol: Chem.Mol) -> Tuple[int, int]:
        """
        Get current indices of connection atoms.
        
        Args:
            mol: RDKit Mol object with marked connection atoms
            
        Returns:
            Tuple of (left_idx, right_idx)
            
        Raises:
            MonomerGeneratorError: If connection atoms not found
        """
        left = AtomMapTracker.find_by_map_num(mol, LEFT_CONN_MAP)
        right = AtomMapTracker.find_by_map_num(mol, RIGHT_CONN_MAP)
        if left is None or right is None:
            raise MonomerGeneratorError(
                f"Connection atoms not found. Left={left}, Right={right}"
            )
        return (left, right)
    
    @staticmethod
    def clear_map_nums(mol: Chem.Mol) -> None:
        """Clear all atom map numbers from molecule."""
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
    
    @staticmethod
    def get_all_mapped_atoms(mol: Chem.Mol) -> Dict[int, int]:
        """
        Get all atoms with non-zero map numbers.
        
        Returns:
            Dict mapping map_num -> atom_idx
        """
        result = {}
        for atom in mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num > 0:
                result[map_num] = atom.GetIdx()
        return result


# =============================================================================
# SMARTSTyper - Atom typing from .fdefn files
# =============================================================================

class SMARTSTyper:
    """
    Assigns force field atom types by parsing RDlt .fdefn files directly.
    
    This standalone implementation:
    - Parses .fdefn files to extract SMARTS patterns
    - Uses RDKit's ChemicalFeatures.BuildFeatureFactory()
    - Matches patterns to assign atom types
    - Supports OPLS-AA and GAFF force fields
    """
    
    def __init__(self, force_field: str = 'oplsaa', verbose: bool = True):
        """
        Initialize SMARTS typer.
        
        Args:
            force_field: 'oplsaa' or 'gaff'
            verbose: Enable logging
        """
        self.force_field = force_field
        self.verbose = verbose
        
        # Locate .fdefn files
        module_dir = Path(__file__).parent
        
        if force_field == 'gaff':
            self.fdef_path = str(module_dir / 'extern' / 'rdlt_data' / 'gaff_lt.fdefn')
        elif force_field == 'lopls':
            self.fdef_path = str(module_dir / 'extern' / 'rdlt_data' / 'lopls_lt.fdefn')
        else:  # oplsaa
            self.fdef_path = str(module_dir / 'extern' / 'rdlt_data' / 'opls_lt.fdefn')
        
        # Load charge dictionary
        self.charge_dict = self._load_charges()
        
        # Build feature factory from .fdefn file
        try:
            self.factory = ChemicalFeatures.BuildFeatureFactory(self.fdef_path)
            if verbose:
                logger.info(f"Loaded feature factory from {self.fdef_path}")
        except Exception as e:
            raise AtomTypingError(f"Failed to load feature factory from {self.fdef_path}: {e}")
    
    def _load_charges(self) -> dict:
        """Load charge dictionary from pickle file."""
        dict_path = self.fdef_path.replace('.fdefn', '_dict.pkl')
        try:
            with open(dict_path, 'rb') as f:
                return pickle.load(f)
        except FileNotFoundError:
            logger.warning(f"Charge dictionary not found: {dict_path}")
            return {}
    
    def assign_atom_types(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Assign force field atom types using SMARTS patterns.
        
        IMPORTANT: This should be called on the CHAIN molecule before splitting,
        so atoms have the correct chemical environment for typing.
        
        Args:
            mol: RDKit Mol object (will be modified in-place)
            
        Returns:
            Mol with AtomType properties set on each atom
            
        Raises:
            AtomTypingError: If typing fails
        """
        try:
            # Get features from feature factory
            features = self.factory.GetFeaturesForMol(mol)
            
            # Assign atom types based on feature matches
            typed_atoms = set()
            for feature in features:
                atom_ids = feature.GetAtomIds()
                if len(atom_ids) == 0:
                    continue
                    
                atom_idx = atom_ids[0]
                atom = mol.GetAtomWithIdx(atom_idx)
                
                # Skip dummy atoms
                if atom.GetAtomicNum() == 0:
                    continue
                
                # Get atom type from feature
                atom_type = feature.GetType()
                atom.SetProp('AtomType', atom_type)
                typed_atoms.add(atom_idx)
                
                if self.verbose:
                    logger.debug(
                        f"Atom {atom_idx} ({atom.GetSymbol()}): {atom_type}"
                    )
            
            # Check for untyped atoms (excluding dummy atoms and hydrogens we might add later)
            untyped = []
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 0:  # Skip dummy atoms
                    continue
                if atom.GetIdx() not in typed_atoms:
                    untyped.append((atom.GetIdx(), atom.GetSymbol()))
            
            if untyped:
                logger.warning(
                    f"Untyped atoms found: {untyped}. "
                    f"This chemistry may not be fully supported."
                )
            
            if self.verbose:
                logger.info(f"Successfully typed {len(typed_atoms)} atoms")
            
            return mol
            
        except Exception as e:
            raise AtomTypingError(f"Atom typing failed: {e}")


# =============================================================================
# ConformerGenerator - ETKDG embedding
# =============================================================================

class ConformerGenerator:
    """Generates simple conformers using ETKDG (no optimization)."""
    
    def __init__(self, random_seed: int = 42, verbose: bool = True):
        """
        Initialize conformer generator.
        
        Args:
            random_seed: Seed for reproducible conformer generation
            verbose: Enable logging
        """
        self.random_seed = random_seed
        self.verbose = verbose
    
    def generate_conformer(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Generate simple conformer using ETKDG (no optimization).
        
        Args:
            mol: RDKit Mol object (with explicit H)
            
        Returns:
            Mol with single conformer
            
        Note:
            Uses ETKDG embedding without optimization per requirements.
            Conformers generated AFTER chain splitting.
            No hydrogen capping is performed - monomers keep their connection points.
        """
        # Generate conformer using ETKDG
        params = AllChem.ETKDGv3()
        params.randomSeed = self.random_seed
        
        result = AllChem.EmbedMolecule(mol, params)
        
        if result == -1:
            # Try with less strict parameters
            params.useRandomCoords = True
            result = AllChem.EmbedMolecule(mol, params)
            
            if result == -1:
                raise MonomerGeneratorError(
                    f"Failed to generate conformer for molecule with {mol.GetNumAtoms()} atoms"
                )
        
        # NO optimization (per requirements)
        
        if self.verbose:
            logger.info(f"Generated conformer for {mol.GetNumAtoms()} atoms")
        
        return mol


# =============================================================================
# ChainBuilder - Connect monomers with map number tracking
# =============================================================================

class ChainBuilder:
    """
    Builds polymer chains by connecting monomer units.
    
    Uses atom map numbers to track inter-monomer bonds for later splitting.
    """
    
    def __init__(self, verbose: bool = True):
        """
        Initialize chain builder.
        
        Args:
            verbose: Enable logging
        """
        self.verbose = verbose
    
    def build_chain(
        self,
        smiles: str,
        n_monomers: int
    ) -> Tuple[Chem.Mol, List[Tuple[int, int]]]:
        """
        Build a polymer chain by connecting monomer units.
        
        Args:
            smiles: SMILES string with wildcards (e.g., "[*]CC[*]")
            n_monomers: Number of monomer units to connect
            
        Returns:
            Tuple of (chain_mol, inter_bond_markers) where:
                - chain_mol: RDKit Mol object representing the polymer chain
                - inter_bond_markers: List of (map_num1, map_num2) tuples identifying
                  atoms at inter-monomer bonds (using map numbers, not indices!)
                  
        Raises:
            ChainBuildingError: If SMILES doesn't have exactly 2 wildcards
        """
        if self.verbose:
            logger.info(f"Building chain with {n_monomers} monomers from SMILES: {smiles}")
        
        # Step 1: Validate SMILES
        mol_template = Chem.MolFromSmiles(smiles)
        if mol_template is None:
            raise ChainBuildingError(f"Invalid SMILES: {smiles}")
        
        # Add explicit hydrogens
        mol_template = Chem.AddHs(mol_template)
        
        # Check for exactly 2 wildcards (dummy atoms)
        dummy_indices = [atom.GetIdx() for atom in mol_template.GetAtoms() if atom.GetAtomicNum() == 0]
        if len(dummy_indices) != 2:
            raise ChainBuildingError(
                f"SMILES must contain exactly 2 wildcard atoms ([*] or *), got {len(dummy_indices)}"
            )
        
        # Step 2: Mark terminal dummy atoms with special isotope
        # These are the chain end dummies that should be capped with H after splitting
        for idx in dummy_indices:
            mol_template.GetAtomWithIdx(idx).SetIsotope(TERMINAL_DUMMY_ISOTOPE)
        
        # Step 3: Build chain
        chain_mol, inter_bond_markers = self._connect_monomer_units(mol_template, n_monomers)
        
        if self.verbose:
            logger.info(f"Built chain with {chain_mol.GetNumAtoms()} atoms")
            logger.info(f"Inter-monomer bond markers: {inter_bond_markers}")
        
        return chain_mol, inter_bond_markers
    
    def _find_dummy_neighbors(self, mol: Chem.Mol) -> Tuple[int, int, int, int]:
        """
        Find dummy atoms and their neighbors (connection points).
        
        Returns:
            Tuple of (left_dummy, left_neighbor, right_dummy, right_neighbor)
        """
        dummy_info = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                dummy_idx = atom.GetIdx()
                # Find non-dummy neighbor
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() > 0:
                        dummy_info.append((dummy_idx, neighbor.GetIdx()))
                        break
        
        if len(dummy_info) != 2:
            raise ChainBuildingError(f"Expected 2 dummy atoms, found {len(dummy_info)}")
        
        # Sort by dummy index to get consistent left/right ordering
        dummy_info.sort(key=lambda x: x[0])
        return (dummy_info[0][0], dummy_info[0][1], dummy_info[1][0], dummy_info[1][1])
    
    def _connect_monomer_units(
        self,
        mol_template: Chem.Mol,
        n_monomers: int
    ) -> Tuple[Chem.Mol, List[Tuple[int, int]]]:
        """
        Connect multiple monomer units into a chain.
        
        Args:
            mol_template: Monomer Mol object with dummy atoms (wildcards)
            n_monomers: Number of units to connect
            
        Returns:
            Tuple of (chain_mol, inter_bond_markers)
        """
        inter_bond_markers = []
        map_counter = INTER_BOND_MAP_START
        
        # Start with first monomer (make a copy)
        chain_mol = Chem.RWMol(mol_template)
        
        # Find dummy atoms in chain so far
        left_dummy, left_ne, right_dummy, right_ne = self._find_dummy_neighbors(chain_mol)
        
        # Connect remaining monomers
        for i in range(1, n_monomers):
            if self.verbose:
                logger.debug(f"Adding monomer {i+1}/{n_monomers}")
            
            # Create a copy of template
            new_monomer = Chem.RWMol(mol_template)
            new_left_dummy, new_left_ne, new_right_dummy, new_right_ne = self._find_dummy_neighbors(new_monomer)
            
            # Offset for new monomer atoms
            offset = chain_mol.GetNumAtoms()
            
            # Add atoms from new monomer (preserve isotope for terminal tracking)
            for atom in new_monomer.GetAtoms():
                new_atom = Chem.Atom(atom.GetAtomicNum())
                new_atom.SetFormalCharge(atom.GetFormalCharge())
                new_atom.SetNumExplicitHs(atom.GetNumExplicitHs())
                new_atom.SetIsotope(atom.GetIsotope())  # Preserve isotope marker
                chain_mol.AddAtom(new_atom)
            
            # Add bonds from new monomer
            for bond in new_monomer.GetBonds():
                begin_idx = bond.GetBeginAtomIdx() + offset
                end_idx = bond.GetEndAtomIdx() + offset
                bond_type = bond.GetBondType()
                chain_mol.AddBond(begin_idx, end_idx, bond_type)
            
            # Mark the atoms that will form the inter-monomer bond
            # right_ne of current chain connects to new_left_ne of new monomer
            chain_mol.GetAtomWithIdx(right_ne).SetAtomMapNum(map_counter)
            chain_mol.GetAtomWithIdx(new_left_ne + offset).SetAtomMapNum(map_counter + 1)
            inter_bond_markers.append((map_counter, map_counter + 1))
            map_counter += 2
            
            # Form bond between chain's right neighbor and new monomer's left neighbor
            chain_mol.AddBond(right_ne, new_left_ne + offset, Chem.BondType.SINGLE)
            
            # Remove dummy atoms (in reverse order by index)
            # right_dummy from chain, new_left_dummy from new monomer
            dummies_to_remove = sorted([right_dummy, new_left_dummy + offset], reverse=True)
            for dummy_idx in dummies_to_remove:
                chain_mol.RemoveAtom(dummy_idx)
            
            # Update right_dummy and right_ne for next iteration
            # After removal, indices shift - need to recalculate
            chain_mol_temp = chain_mol.GetMol()
            
            # Find the rightmost dummy atom in the updated chain
            for atom in chain_mol_temp.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    right_dummy = atom.GetIdx()
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() > 0:
                            right_ne = neighbor.GetIdx()
                            break
            
            # Convert back to RWMol for next iteration
            chain_mol = Chem.RWMol(chain_mol_temp)
        
        # Final molecule
        final_mol = chain_mol.GetMol()
        Chem.SanitizeMol(final_mol)
        
        return final_mol, inter_bond_markers


# =============================================================================
# ChainSplitter - Fragment and extract monomers
# =============================================================================

class ChainSplitter:
    """
    Splits polymer chains into monomer fragments.
    
    Uses correct FragmentOnBonds with BOND indices (not atom indices).
    """
    
    def __init__(self, conformer_gen: ConformerGenerator, verbose: bool = True):
        """
        Initialize chain splitter.
        
        Args:
            conformer_gen: ConformerGenerator instance
            verbose: Enable logging
        """
        self.conformer_gen = conformer_gen
        self.verbose = verbose
    
    def split_chain(
        self,
        chain_mol: Chem.Mol,
        inter_bond_markers: List[Tuple[int, int]],
        base_name: str,
        force_field: str
    ) -> List[MonomerVariant]:
        """
        Split chain into monomer fragments.
        
        Args:
            chain_mol: Polymer chain RDKit Mol object (with atom types assigned!)
            inter_bond_markers: List of (map_num1, map_num2) tuples for inter-monomer bonds
            base_name: Base name for monomers
            force_field: Force field name
            
        Returns:
            List of MonomerVariant objects
        """
        if self.verbose:
            logger.info(f"Splitting chain into {len(inter_bond_markers)+1} monomers")

        # Extract cap H charge mapping from chain_mol
        # This was stored by from_smiles() or from_chain() methods
        cap_h_charge_map = {}
        if chain_mol.HasProp('_CapHChargeMap'):
            import json
            try:
                cap_h_json = chain_mol.GetProp('_CapHChargeMap')
                parsed = json.loads(cap_h_json)
                cap_h_charge_map = {int(k): tuple(v) for k, v in parsed.items()}
                if self.verbose:
                    logger.debug(f"Extracted cap_h_charge_map with {len(cap_h_charge_map)} entries")
            except Exception as e:
                logger.warning(f"Failed to parse cap H charge map: {e}")

        n_monomers = len(inter_bond_markers) + 1
        
        # Convert map numbers to bond indices
        bond_indices_to_break = []
        for map1, map2 in inter_bond_markers:
            idx1 = AtomMapTracker.find_by_map_num(chain_mol, map1)
            idx2 = AtomMapTracker.find_by_map_num(chain_mol, map2)
            
            if idx1 is None or idx2 is None:
                raise ChainSplittingError(
                    f"Could not find atoms for map numbers {map1}, {map2}"
                )
            
            bond = chain_mol.GetBondBetweenAtoms(idx1, idx2)
            if bond is None:
                raise ChainSplittingError(
                    f"No bond between atoms {idx1} and {idx2}"
                )
            
            bond_indices_to_break.append(bond.GetIdx())
        
        if self.verbose:
            logger.debug(f"Bond indices to break: {bond_indices_to_break}")

        # Preserve Gasteiger charges through fragmentation
        # FragmentOnBonds() loses atom properties, so we need to manually track them
        # Charges are already stored with map numbers from MonomerGenerator
        gasteiger_charges_by_map_num = {}
        if force_field.lower() == 'gaff':
            # Store Gasteiger charges keyed by map number (using existing map numbers)
            for atom in chain_mol.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if map_num > 0:  # Only store charges for atoms with map numbers
                    try:
                        charge = atom.GetProp('_GasteigerCharge')
                        gasteiger_charges_by_map_num[map_num] = charge
                    except KeyError:
                        pass  # No charge stored

            if self.verbose:
                logger.debug(f"Stored {len(gasteiger_charges_by_map_num)} Gasteiger charges")
                # Debug: print first few charges
                for i, (map_num, charge) in enumerate(list(gasteiger_charges_by_map_num.items())[:5]):
                    logger.debug(f"  Map {map_num} -> charge {charge}")

        # Fragment on BOND indices (not atom indices!)
        fragmented = Chem.FragmentOnBonds(
            chain_mol,
            bond_indices_to_break,
            addDummies=True
        )

        # Get fragments as separate molecules
        fragments = Chem.GetMolFrags(fragmented, asMols=True, sanitizeFrags=False)

        # Restore Gasteiger charges to fragments
        if gasteiger_charges_by_map_num:
            restored_count = 0
            for frag in fragments:
                for atom in frag.GetAtoms():
                    map_num = atom.GetAtomMapNum()
                    if map_num in gasteiger_charges_by_map_num:
                        atom.SetProp('_GasteigerCharge', str(gasteiger_charges_by_map_num[map_num]))
                        restored_count += 1

            if self.verbose:
                logger.debug(f"Restored {restored_count} Gasteiger charges to {len(fragments)} fragments")
                # Debug: verify restoration by checking first fragment
                if fragments:
                    first_frag = fragments[0]
                    for atom in list(first_frag.GetAtoms())[:5]:
                        map_num = atom.GetAtomMapNum()
                        try:
                            charge = atom.GetProp('_GasteigerCharge')
                            logger.debug(f"  Verified: atom {atom.GetIdx()} (map={map_num}, elem={atom.GetSymbol()}) has charge {charge}")
                        except KeyError:
                            logger.debug(f"  Verified: atom {atom.GetIdx()} (map={map_num}, elem={atom.GetSymbol()}) has NO CHARGE")

            # Clear temporary CHARGE map numbers from original chain
            # Don't clear INTER_BOND_MAP_START range - those are needed for connection tracking
            for atom in chain_mol.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if map_num >= CHARGE_MAP_START:
                    atom.SetAtomMapNum(0)
        
        if len(fragments) != n_monomers:
            raise ChainSplittingError(
                f"Expected {n_monomers} fragments, got {len(fragments)}"
            )
        
        # Process each fragment
        variants = []
        for i, frag in enumerate(fragments):
            # Determine variant type
            if n_monomers == 1:
                variant_type = 'single'
            elif i == 0:
                variant_type = 'first'
            elif i == n_monomers - 1:
                variant_type = 'last'
            else:
                variant_type = 'middle'
            
            # Process fragment
            processed_mol, conn_atoms = self._process_fragment(frag, i, variant_type, force_field, cap_h_charge_map)
            
            # Create variant
            variant = MonomerVariant(
                base_name=base_name,
                variant_type=variant_type,
                mol=processed_mol,
                smiles=Chem.MolToSmiles(processed_mol),
                connection_atoms=conn_atoms,
                force_field=force_field,
                position=i
            )
            variants.append(variant)
            
            if self.verbose:
                logger.info(f"Created {variant_type} variant at position {i}")
        
        return variants
    
    def _process_fragment(
        self,
        frag: Chem.Mol,
        position: int,
        variant_type: str,
        force_field: str,
        cap_h_charge_map: Dict[int, Tuple[int, str]]
    ) -> Tuple[Chem.Mol, Tuple[int, int]]:
        """
        Process fragment: find connections, handle terminal/connection dummies, generate conformer.

        Terminal dummies (isotope 99): Chain ends - cap with H
        Connection dummies (no isotope): Polymerization points - identify as connection atoms

        Args:
            frag: Fragment molecule with dummy atoms
            position: Position in chain
            variant_type: 'first', 'middle', 'last', or 'single'
            force_field: Force field type for atom typing ('oplsaa', 'lopls', 'gaff')
            cap_h_charge_map: Mapping of neighbor_map_num -> (cap_h_map_num, charge) from full chain

        Returns:
            Tuple of (processed_mol, (left_conn_idx, right_conn_idx))
        """
        rw_mol = Chem.RWMol(frag)
        
        # 1. Classify dummy atoms as terminal (isotope 99) or connection (no isotope)
        terminal_dummies = []  # [(dummy_idx, neighbor_idx), ...] - will be capped with H
        connection_dummies = []  # [(dummy_idx, neighbor_idx), ...] - connection points
        
        for atom in rw_mol.GetAtoms():
            if atom.GetAtomicNum() == 0:  # Dummy atom
                dummy_idx = atom.GetIdx()
                neighbor_idx = None
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() > 0:
                        neighbor_idx = neighbor.GetIdx()
                        break
                
                if neighbor_idx is not None:
                    if atom.GetIsotope() == TERMINAL_DUMMY_ISOTOPE:
                        terminal_dummies.append((dummy_idx, neighbor_idx))
                    else:
                        connection_dummies.append((dummy_idx, neighbor_idx))
        
        if self.verbose:
            logger.debug(f"Terminal dummies (cap with H): {terminal_dummies}")
            logger.debug(f"Connection dummies (connection points): {connection_dummies}")
        
        # Sort by dummy index for consistent left/right ordering
        terminal_dummies.sort(key=lambda x: x[0])
        connection_dummies.sort(key=lambda x: x[0])
        
        # 2. Mark connection atoms with map numbers (only for connection dummies)
        for i, (dummy_idx, neighbor_idx) in enumerate(connection_dummies):
            if i == 0:
                rw_mol.GetAtomWithIdx(neighbor_idx).SetAtomMapNum(LEFT_CONN_MAP)
            elif i == 1:
                rw_mol.GetAtomWithIdx(neighbor_idx).SetAtomMapNum(RIGHT_CONN_MAP)
        
        # 3. Replace terminal dummies with H atoms (cap chain ends)
        # This maintains proper valence for terminal atoms
        # Track added H atoms, their neighbors, and neighbor's map number for charge restoration
        added_cap_h_atoms = []  # [(h_idx, neighbor_idx, neighbor_map_num), ...]

        for dummy_idx, neighbor_idx in terminal_dummies:
            # Get the bond type between dummy and neighbor
            bond = rw_mol.GetBondBetweenAtoms(dummy_idx, neighbor_idx)
            bond_type = bond.GetBondType() if bond else Chem.BondType.SINGLE

            # Add new H atom
            h_idx = rw_mol.AddAtom(Chem.Atom(1))  # Add hydrogen

            # Get the neighbor's map number (for charge lookup)
            neighbor_map_num = rw_mol.GetAtomWithIdx(neighbor_idx).GetAtomMapNum()
            added_cap_h_atoms.append((h_idx, neighbor_idx, neighbor_map_num))

            # Add bond from neighbor to new H
            rw_mol.AddBond(neighbor_idx, h_idx, bond_type)
        
        # 3b. Assign atom types to cap H atoms based on force field
        # This ensures the cap atoms have valid force field types
        for h_idx, neighbor_idx, _ in added_cap_h_atoms:
            h_atom = rw_mol.GetAtomWithIdx(h_idx)
            neighbor_atom = rw_mol.GetAtomWithIdx(neighbor_idx)
            
            if force_field in ('oplsaa', 'lopls'):
                # OPLS-AA: H on sp3 carbon = type 85 (HC)
                h_atom.SetProp('AtomType', '@atom:85')
                
                # Check if neighbor carbon needs its type updated
                # This happens when the original typing couldn't match due to dummy atoms
                if neighbor_atom.GetAtomicNum() == 6:
                    # Check if the existing type is invalid (e.g., '@atom:c' fallback)
                    existing_type = neighbor_atom.GetProp('AtomType') if neighbor_atom.HasProp('AtomType') else ''
                    # Invalid if it's just element symbol or doesn't have numeric part
                    is_invalid = not existing_type or existing_type in ('@atom:c', '@atom:C', 'c', 'C')
                    
                    if is_invalid:
                        # Count H neighbors to determine sp3 carbon type
                        h_count = sum(1 for n in neighbor_atom.GetNeighbors() if n.GetAtomicNum() == 1)
                        
                        # OPLS-AA sp3 carbon types:
                        # 80: CH3 (methyl) - 3 H neighbors
                        # 81: CH2 (methylene) - 2 H neighbors  
                        # 82: CH (methine) - 1 H neighbor
                        # 84: C (quaternary) - 0 H neighbors
                        if h_count >= 3:
                            neighbor_atom.SetProp('AtomType', '@atom:80')
                        elif h_count == 2:
                            neighbor_atom.SetProp('AtomType', '@atom:81')
                        elif h_count == 1:
                            neighbor_atom.SetProp('AtomType', '@atom:82')
                        else:
                            neighbor_atom.SetProp('AtomType', '@atom:84')
            else:
                # GAFF: H on sp3 carbon = hc
                h_atom.SetProp('AtomType', '@atom:hc')
                
                # For terminal carbons in GAFF, use c3 if not already typed
                if neighbor_atom.GetAtomicNum() == 6:
                    existing_type = neighbor_atom.GetProp('AtomType') if neighbor_atom.HasProp('AtomType') else ''
                    is_invalid = not existing_type or existing_type in ('@atom:c', '@atom:C', 'c', 'C')
                    if is_invalid:
                        neighbor_atom.SetProp('AtomType', '@atom:c3')
        
        # 4. Remove dummy atoms (terminal dummies replaced with H, connection dummies just removed)
        # We remove ALL dummies, but terminal positions already have H added
        all_dummy_indices = [d[0] for d in terminal_dummies] + [d[0] for d in connection_dummies]
        for idx in sorted(all_dummy_indices, reverse=True):
            rw_mol.RemoveAtom(idx)
        
        # 5. Get molecule
        mol = rw_mol.GetMol()

        # 5b. For GAFF force field, restore cap H charges from full chain
        # Cap H charges were calculated on the full chain and preserved via cap_h_charge_map
        if force_field.lower() == 'gaff':
            restored_count = 0
            for h_idx, neighbor_idx, neighbor_map_num in added_cap_h_atoms:
                if neighbor_map_num in cap_h_charge_map:
                    _, charge = cap_h_charge_map[neighbor_map_num]
                    # Find the cap H atom in the final molecule by looking for an H atom
                    # bonded to the atom with the matching neighbor_map_num
                    neighbor_atom_idx = AtomMapTracker.find_by_map_num(mol, neighbor_map_num)
                    if neighbor_atom_idx is not None:
                        # Find the H atom bonded to this neighbor
                        cap_h_found = False
                        for neighbor in mol.GetAtomWithIdx(neighbor_atom_idx).GetNeighbors():
                            if neighbor.GetAtomicNum() == 1:  # Hydrogen
                                try:
                                    neighbor.GetProp('_GasteigerCharge')
                                    # This H already has a charge, skip it
                                    continue
                                except KeyError:
                                    # This is the cap H atom - set its charge
                                    neighbor.SetProp('_GasteigerCharge', charge)
                                    cap_h_found = True
                                    restored_count += 1
                                    break
                        if not cap_h_found:
                            logger.warning(f"Could not find cap H atom for neighbor map {neighbor_map_num}")
                    else:
                        logger.warning(f"Could not find neighbor atom with map {neighbor_map_num}")
                else:
                    logger.error(
                        f"Missing cap H charge for neighbor map {neighbor_map_num}. "
                        f"This should not happen - charges should be preserved from full chain."
                    )

            if self.verbose:
                logger.debug(f"Restored {restored_count}/{len(added_cap_h_atoms)} cap H charges from full chain")

        # 6. Sanitize molecule
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            logger.warning(f"Sanitization warning: {e}")
        
        # 7. Generate conformer
        mol = self.conformer_gen.generate_conformer(mol)
        
        # 8. Find connection atoms via map numbers
        left_conn = AtomMapTracker.find_by_map_num(mol, LEFT_CONN_MAP)
        right_conn = AtomMapTracker.find_by_map_num(mol, RIGHT_CONN_MAP)
        
        # Handle edge cases based on variant type
        if variant_type == 'first':
            # First monomer: has terminal on left, connection on right
            # The connection point should be on the right side
            if left_conn is not None and right_conn is None:
                # Only one connection marked, it's the right connection
                right_conn = left_conn
                left_conn = 0  # Will be reassigned
        elif variant_type == 'last':
            # Last monomer: has connection on left, terminal on right
            if left_conn is not None and right_conn is None:
                # left_conn is correct, right_conn is terminal (capped)
                right_conn = mol.GetNumAtoms() - 1  # Placeholder
        elif variant_type == 'single':
            # Single monomer: both ends are terminal (capped)
            left_conn = 0
            right_conn = mol.GetNumAtoms() - 1
        
        # Default if None
        if left_conn is None:
            left_conn = 0
        if right_conn is None:
            right_conn = mol.GetNumAtoms() - 1
        
        if self.verbose:
            logger.debug(f"Connection atoms: left={left_conn}, right={right_conn}")
            logger.debug(f"Total atoms after processing: {mol.GetNumAtoms()}")
        
        return mol, (left_conn, right_conn)


# =============================================================================
# BackboneAligner - Align to X-axis, reorder atoms
# =============================================================================

class BackboneAligner:
    """
    Aligns backbone and reorders atoms for AutoPoly compatibility.
    
    CRITICAL: Ensures connection atoms are FIRST in Data Atoms block
    for AutoPoly's read_lt_end_atoms() compatibility.
    """
    
    def __init__(self, verbose: bool = True):
        """
        Initialize backbone aligner.
        
        Args:
            verbose: Enable logging
        """
        self.verbose = verbose
    
    def align_for_lt_file(self, variant: MonomerVariant) -> MonomerVariant:
        """
        Align monomer variant for LT file generation.
        
        Performs:
        1. Align backbone to X-axis
        2. Reorder atoms (connection atoms FIRST!)
        
        Args:
            variant: MonomerVariant to align
            
        Returns:
            New MonomerVariant with aligned mol
        """
        mol = variant.mol
        conn_left, conn_right = variant.connection_atoms
        
        # Step 1: Align backbone to X-axis
        mol_aligned = self._align_backbone_to_x_axis(mol, conn_left, conn_right)
        
        # Step 2: Generate atom IDs with connection atoms first
        atom_order, atom_ids = self._generate_atom_order_and_ids(
            mol_aligned, conn_left, conn_right, variant.variant_type
        )
        
        # Create new variant with aligned mol and atom IDs
        aligned_variant = MonomerVariant(
            base_name=variant.base_name,
            variant_type=variant.variant_type,
            mol=mol_aligned,
            smiles=Chem.MolToSmiles(mol_aligned),
            connection_atoms=(conn_left, conn_right),
            force_field=variant.force_field,
            position=variant.position,
            atom_ids=atom_ids
        )
        
        return aligned_variant
    
    def _align_backbone_to_x_axis(
        self,
        mol: Chem.Mol,
        left_idx: int,
        right_idx: int
    ) -> Chem.Mol:
        """
        Rotate and translate molecule to align backbone with X-axis.
        
        Args:
            mol: RDKit Mol object
            left_idx: Left connection atom index
            right_idx: Right connection atom index
            
        Returns:
            Aligned molecule
        """
        if mol.GetNumConformers() == 0:
            return mol
        
        conf = mol.GetConformer(0)
        
        # Get positions
        left_pos = np.array([
            conf.GetAtomPosition(left_idx).x,
            conf.GetAtomPosition(left_idx).y,
            conf.GetAtomPosition(left_idx).z
        ])
        right_pos = np.array([
            conf.GetAtomPosition(right_idx).x,
            conf.GetAtomPosition(right_idx).y,
            conf.GetAtomPosition(right_idx).z
        ])
        
        # Calculate midpoint and backbone vector
        midpoint = (right_pos + left_pos) / 2
        backbone_vec = right_pos - left_pos
        backbone_length = np.linalg.norm(backbone_vec)
        
        if backbone_length < 0.001:
            if self.verbose:
                logger.warning("Connection atoms too close, skipping alignment")
            return mol
        
        # Normalize backbone vector
        backbone_unit = backbone_vec / backbone_length
        
        # Target vector (X-axis)
        target_vec = np.array([1.0, 0.0, 0.0])
        
        # Calculate rotation
        rotation_axis = np.cross(backbone_unit, target_vec)
        axis_norm = np.linalg.norm(rotation_axis)
        
        new_mol = Chem.Mol(mol)
        new_conf = new_mol.GetConformer(0)
        
        if axis_norm < 0.001:
            # Vectors are parallel or anti-parallel
            if np.dot(backbone_unit, target_vec) < 0:
                # Anti-parallel - flip
                for i in range(mol.GetNumAtoms()):
                    pos = np.array([
                        conf.GetAtomPosition(i).x,
                        conf.GetAtomPosition(i).y,
                        conf.GetAtomPosition(i).z
                    ])
                    pos_centered = pos - midpoint
                    pos_flipped = np.array([-pos_centered[0], pos_centered[1], pos_centered[2]])
                    new_conf.SetAtomPosition(i, pos_flipped)
            else:
                # Already aligned, just center
                for i in range(mol.GetNumAtoms()):
                    pos = np.array([
                        conf.GetAtomPosition(i).x,
                        conf.GetAtomPosition(i).y,
                        conf.GetAtomPosition(i).z
                    ])
                    pos_centered = pos - midpoint
                    new_conf.SetAtomPosition(i, pos_centered)
        else:
            # Perform rotation
            rotation_axis = rotation_axis / axis_norm
            cos_theta = np.dot(backbone_unit, target_vec)
            sin_theta = axis_norm
            
            for i in range(mol.GetNumAtoms()):
                pos = np.array([
                    conf.GetAtomPosition(i).x,
                    conf.GetAtomPosition(i).y,
                    conf.GetAtomPosition(i).z
                ])
                pos_centered = pos - midpoint
                pos_rotated = self._rodrigues_rotation(
                    pos_centered, rotation_axis, cos_theta, sin_theta
                )
                new_conf.SetAtomPosition(i, pos_rotated)
        
        return new_mol
    
    def _rodrigues_rotation(
        self,
        v: np.ndarray,
        k: np.ndarray,
        cos_t: float,
        sin_t: float
    ) -> np.ndarray:
        """Rotate vector v around unit axis k using Rodrigues' formula."""
        return (v * cos_t +
                np.cross(k, v) * sin_t +
                k * np.dot(k, v) * (1 - cos_t))
    
    def _generate_atom_order_and_ids(
        self,
        mol: Chem.Mol,
        conn_left: int,
        conn_right: int,
        variant_type: str
    ) -> Tuple[List[int], Dict[int, str]]:
        """
        Generate atom ordering and IDs with connection atoms first.
        
        Args:
            mol: RDKit Mol object
            conn_left: Left connection atom index
            conn_right: Right connection atom index
            variant_type: Type of variant
            
        Returns:
            Tuple of (atom_order, atom_ids) where:
                - atom_order: List of atom indices in desired order
                - atom_ids: Dict mapping atom_idx -> LT atom ID (e.g., 'C1', 'H2')
        """
        # Build atom order with connection atoms first
        # For first/last monomers, only ONE side is a real connection point
        atom_order = []
        
        # Determine which atoms are real connection points
        if variant_type == 'first':
            # First monomer: only RIGHT side is a connection point (left side is terminal/capped)
            # Put the right connection atom first
            atom_order = [conn_right]
        elif variant_type == 'last':
            # Last monomer: only LEFT side is a connection point (right side is terminal/capped)
            # Put the left connection atom first
            atom_order = [conn_left]
        elif variant_type == 'single':
            # Single monomer: no real connections (both ends are capped)
            atom_order = []
        else:
            # Middle: both sides are connection points
            atom_order = [conn_left, conn_right] if conn_left != conn_right else [conn_left]
        
        # Add remaining atoms (heavy atoms first, then hydrogens for better readability)
        conn_set = set(atom_order)
        heavy_atoms = []
        h_atoms = []
        for i in range(mol.GetNumAtoms()):
            if i not in conn_set:
                atom = mol.GetAtomWithIdx(i)
                if atom.GetAtomicNum() == 1:  # Hydrogen
                    h_atoms.append(i)
                else:
                    heavy_atoms.append(i)
        atom_order.extend(heavy_atoms)
        atom_order.extend(h_atoms)
        
        # Generate sequential IDs: C1, H2, O3, etc.
        atom_ids = {}
        for seq_num, atom_idx in enumerate(atom_order, 1):
            atom = mol.GetAtomWithIdx(atom_idx)
            element = atom.GetSymbol()
            atom_ids[atom_idx] = f"{element}{seq_num}"
        
        return atom_order, atom_ids


# =============================================================================
# LTWriter - AutoPoly-compatible format
# =============================================================================

class LTWriter:
    """
    Writes Moltemplate .lt files for monomer variants.
    
    Matches exact AutoPoly format from example files.
    """
    
    def __init__(
        self,
        force_field: str,
        charge_dict: Dict[str, float],
        verbose: bool = True
    ):
        """
        Initialize LT writer.
        
        Args:
            force_field: 'oplsaa' or 'gaff'
            charge_dict: Atom type to charge mapping
            verbose: Enable logging
        """
        self.force_field = force_field
        self.charge_dict = charge_dict
        self.verbose = verbose

    def write_variant(
        self,
        variant: MonomerVariant,
        output_dir: str,
        generate_t1: bool = True
    ) -> List[str]:
        """
        Write .lt file(s) for a monomer variant.
        
        Args:
            variant: MonomerVariant object
            output_dir: Output directory path
            generate_t1: Whether to generate T1 chirality variant
            
        Returns:
            List of generated file paths
        """
        files = []
        
        # Generate filename based on variant type
        if variant.variant_type == 'first':
            suffix = 'le'
        elif variant.variant_type == 'last':
            suffix = 're'
        elif variant.variant_type == 'single':
            suffix = 'single'
        else:
            suffix = 'i'
        
        # Use position in name for unique identification
        filename = f"{variant.base_name}_{variant.position}{suffix}.lt"
        filepath = os.path.join(output_dir, filename)
        
        # Write main variant
        self._write_lt_file(variant, filepath)
        files.append(filepath)
        
        if self.verbose:
            logger.info(f"Generated: {filepath}")
        
        # Generate T1 variant if requested
        if generate_t1:
            variant_t1 = self._create_t1_variant(variant)
            filename_t1 = f"{variant.base_name}_{variant.position}{suffix}_T1.lt"
            filepath_t1 = os.path.join(output_dir, filename_t1)
            
            self._write_lt_file(variant_t1, filepath_t1)
            files.append(filepath_t1)
            
            if self.verbose:
                logger.info(f"Generated: {filepath_t1}")
        
        return files
    
    def _create_t1_variant(self, variant: MonomerVariant) -> MonomerVariant:
        """
        Create T1 variant with opposite chirality.
        
        Inverts Z-coordinates to create mirror image.
        """
        mol_t1 = Chem.Mol(variant.mol)
        
        if mol_t1.GetNumConformers() > 0:
            conf = variant.mol.GetConformer(0)
            conf_t1 = mol_t1.GetConformer(0)
            
            # Invert Z coordinates
            for i in range(variant.mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                conf_t1.SetAtomPosition(i, (pos.x, pos.y, -pos.z))
        
        # Create new variant with is_t1 flag
        variant_t1 = MonomerVariant(
            base_name=variant.base_name,
            variant_type=variant.variant_type,
            mol=mol_t1,
            smiles=Chem.MolToSmiles(mol_t1),
            connection_atoms=variant.connection_atoms,
            force_field=variant.force_field,
            position=variant.position,
            atom_ids=variant.atom_ids,
            is_t1=True
        )
        
        return variant_t1
    
    def _write_lt_file(self, variant: MonomerVariant, filepath: str) -> None:
        """Write a single .lt file."""
        with open(filepath, 'w') as f:
            # Write header
            self._write_header(f, variant)
            
            # Write atoms block (connection atoms FIRST!)
            self._write_atoms_block(f, variant)
            
            # Write bonds block
            self._write_bonds_block(f, variant)
            
            # Write footer
            self._write_footer(f, variant)
    
    def _write_header(self, f, variant: MonomerVariant) -> None:
        """Write LT file header."""
        # Use consistent suffix with filename: le=first, re=last, i=middle, single=single
        suffix_map = {'first': 'le', 'last': 're', 'middle': 'i', 'single': 'single'}
        suffix = suffix_map.get(variant.variant_type, variant.variant_type[0])
        class_name = f"{variant.base_name}_{variant.position}{suffix}"
        if variant.is_t1:
            class_name += "_T1"
        
        if variant.force_field == 'gaff':
            f.write('import "gaff.lt"    # <-- defines the GAFF (General Amber Force Field)\n')
            f.write('# NOTE: GAFF requires user-supplied charges (AM1-BCC or RESP recommended)\n')
            f.write('# See: http://ambermd.org/antechamber/gaff.pdf\n')
            f.write(f'{class_name} inherits GAFF {{\n\n')
        else:
            f.write('import "oplsaa.lt"    # <-- defines the OPLS-AA force field\n')
            f.write(f'{class_name} inherits OPLSAA {{\n\n')
        
        f.write('# atom-id  mol-id  atom-type charge      X         Y        Z\n\n')
    
    def _write_atoms_block(self, f, variant: MonomerVariant) -> None:
        """
        Write Data Atoms block with connection atoms FIRST.
        
        CRITICAL: AutoPoly's read_lt_end_atoms() reads the first two atoms
        from this block to identify connection points.
        
        For first/last monomers, only ONE side is a real connection point.
        """
        f.write('  write("Data Atoms") {\n')
        
        mol = variant.mol
        atom_ids = variant.atom_ids
        conn_left, conn_right = variant.connection_atoms
        
        # Build atom order with connection atoms first (based on variant type)
        if variant.variant_type == 'first':
            # First monomer: only RIGHT side is a connection point
            atom_order = [conn_right]
        elif variant.variant_type == 'last':
            # Last monomer: only LEFT side is a connection point
            atom_order = [conn_left]
        elif variant.variant_type == 'single':
            # Single monomer: no real connections
            atom_order = []
        else:
            # Middle: both sides are connection points
            atom_order = [conn_left]
            if conn_right != conn_left:
                atom_order.append(conn_right)
        
        # Add remaining atoms (heavy atoms first, then hydrogens)
        conn_set = set(atom_order)
        heavy_atoms = []
        h_atoms = []
        for i in range(mol.GetNumAtoms()):
            if i not in conn_set:
                atom = mol.GetAtomWithIdx(i)
                if atom.GetAtomicNum() == 1:  # Hydrogen
                    h_atoms.append(i)
                else:
                    heavy_atoms.append(i)
        atom_order.extend(heavy_atoms)
        atom_order.extend(h_atoms)
        
        # Regenerate atom IDs based on output order
        output_atom_ids = {}
        for seq_num, atom_idx in enumerate(atom_order, 1):
            atom = mol.GetAtomWithIdx(atom_idx)
            element = atom.GetSymbol()
            output_atom_ids[atom_idx] = f"{element}{seq_num}"
        
        # Get conformer for coordinates
        if mol.GetNumConformers() > 0:
            conf = mol.GetConformer(0)
        else:
            conf = None

        for atom_idx in atom_order:
            atom = mol.GetAtomWithIdx(atom_idx)
            atom_id = output_atom_ids[atom_idx]

            # Get atom type
            try:
                atom_type = atom.GetProp('AtomType')
            except KeyError:
                # Fallback atom type based on element
                atom_type = f"@atom:{atom.GetSymbol().lower()}"

            # Get charge - try Gasteiger property first (calculated on full chain),
            # then fallback to charge_dict for OPLS or legacy behavior
            charge = 0.0
            if self.force_field.lower() == 'gaff':
                try:
                    raw_charge = atom.GetProp('_GasteigerCharge')
                    charge = float(raw_charge)
                    # Handle NaN/Inf values (Gasteiger can produce these)
                    if (charge == float('inf') or
                        charge == float('-inf') or
                        charge != charge):  # NaN check
                        charge = 0.0
                except (KeyError, ValueError):
                    charge = 0.0
            else:
                charge = self.charge_dict.get(atom_type, 0.0)
            
            # Get coordinates
            if conf is not None:
                pos = conf.GetAtomPosition(atom_idx)
                x, y, z = pos.x, pos.y, pos.z
            else:
                x, y, z = 0.0, 0.0, 0.0
            
            f.write(f'\t$atom:{atom_id} $mol:... {atom_type} {charge:.4f}')
            f.write(f'    {x:.3f}   {y:.3f}   {z:.3f}\n')
        
        f.write('  }\n\n')
        
        # Store the output atom IDs for bonds block
        variant.atom_ids = output_atom_ids
    
    def _write_bonds_block(self, f, variant: MonomerVariant) -> None:
        """Write Data Bond List using atom IDs."""
        f.write("  write('Data Bond List') {\n")
        
        mol = variant.mol
        atom_ids = variant.atom_ids
        
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            
            id1 = atom_ids.get(a1, f"X{a1}")
            id2 = atom_ids.get(a2, f"X{a2}")
            
            f.write(f'\t$bond:{id1}{id2}\t$atom:{id1}\t$atom:{id2}\n')
        
        f.write('  }\n')
    
    def _write_footer(self, f, variant: MonomerVariant) -> None:
        """Write LT file footer."""
        # Use consistent suffix with filename: le=first, re=last, i=middle, single=single
        suffix_map = {'first': 'le', 'last': 're', 'middle': 'i', 'single': 'single'}
        suffix = suffix_map.get(variant.variant_type, variant.variant_type[0])
        class_name = f"{variant.base_name}_{variant.position}{suffix}"
        if variant.is_t1:
            class_name += "_T1"
        f.write(f'}}  # {class_name}\n\n')
        f.write("# Note: You don't need to supply the partial partial charges of the atoms.\n")
        f.write("#       If you like, just fill the fourth column with zeros (\"0.000\").\n")
        f.write("#       Moltemplate and LAMMPS will automatically assign the charge later\n\n")


# =============================================================================
# MonomerGenerator - Main API class
# =============================================================================

class MonomerGenerator:
    """
    Main class for generating monomer variants from SMILES or chain molecules.
    
    Example:
        >>> generator = MonomerGenerator(
        ...     base_name="PE",
        ...     force_field="gaff",
        ...     output_dir="./monomers"
        ... )
        >>> variants = generator.from_smiles("[*]CC[*]", n_monomers=5)
        >>> files = generator.write_lt_files(variants)
        >>> print(files)
        ['./monomers/PE_0le.lt', './monomers/PE_1i.lt', ...]
    """
    
    def __init__(
        self,
        base_name: str,
        force_field: str = 'gaff',
        output_dir: str = './monomers',
        verbose: bool = True
    ):
        """
        Initialize monomer generator.
        
        Args:
            base_name: Base name for monomers (e.g., "PE", "PMMA")
            force_field: 'oplsaa' or 'gaff'
            output_dir: Directory for .lt files
            verbose: Enable verbose logging
        """
        self.base_name = base_name
        self.force_field = force_field
        self.output_dir = output_dir
        self.verbose = verbose
        
        # Validate force field
        if force_field not in ['oplsaa', 'gaff', 'lopls']:
            raise MonomerGeneratorError(
                f"Unknown force field: '{force_field}'. Use 'oplsaa', 'gaff', or 'lopls'"
            )
        
        # Initialize components
        self.atom_typer = SMARTSTyper(force_field, verbose)
        self.conformer_gen = ConformerGenerator(verbose=verbose)
        self.chain_builder = ChainBuilder(verbose=verbose)
        self.chain_splitter = ChainSplitter(self.conformer_gen, verbose)
        self.backbone_aligner = BackboneAligner(verbose=verbose)
        self.lt_writer = LTWriter(force_field, self.atom_typer.charge_dict, verbose)
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        if verbose:
            logger.info(f"MonomerGenerator initialized for '{base_name}'")
            logger.info(f"Force field: {force_field}")
            logger.info(f"Output directory: {output_dir}")
    
    def from_smiles(
        self,
        smiles: str,
        n_monomers: int = 3
    ) -> List[MonomerVariant]:
        """
        Generate monomer variants from SMILES string.
        
        Args:
            smiles: SMILES string with wildcards (e.g., "[*]CC[*]")
            n_monomers: Number of monomers to build chain with (default 3)
            
        Returns:
            List of MonomerVariant objects
            
        Note:
            For proper atom typing, we need to build a chain first so atoms
            have the correct chemical environment. The chain is then split
            back into individual monomers.
        """
        if self.verbose:
            logger.info(f"Generating monomers from SMILES: {smiles}")
        
        if n_monomers < 3:
            logger.warning("Using at least 3 monomers for proper chemical environment")
            n_monomers = 3
        
        # 1. Build chain
        chain_mol, inter_bond_markers = self.chain_builder.build_chain(smiles, n_monomers)

        # 1b. Cap terminal dummies BEFORE typing
        # This ensures terminal atoms have correct chemistry (e.g., alcohol OH, not ether O-dummy)
        # for SMARTS pattern matching. Without this, terminal atoms bonded to dummies fail to match
        # patterns that require specific neighbor atoms (e.g., ether O needs C-O-C, not dummy-O-C)
        rw_mol = Chem.RWMol(chain_mol)

        # First collect terminal dummies and their neighbors (don't modify while iterating)
        terminal_dummies_to_cap = []  # [(dummy_idx, neighbor_idx, bond_type), ...]
        for atom in rw_mol.GetAtoms():
            if atom.GetAtomicNum() == 0 and atom.GetIsotope() == TERMINAL_DUMMY_ISOTOPE:
                dummy_idx = atom.GetIdx()
                # Find the real neighbor atom
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() > 0:
                        neighbor_idx = neighbor.GetIdx()
                        bond = rw_mol.GetBondBetweenAtoms(dummy_idx, neighbor_idx)
                        bond_type = bond.GetBondType() if bond else Chem.BondType.SINGLE
                        terminal_dummies_to_cap.append((dummy_idx, neighbor_idx, bond_type))
                        break

        # Now add H atoms to cap terminal positions
        for dummy_idx, neighbor_idx, bond_type in terminal_dummies_to_cap:
            h_idx = rw_mol.AddAtom(Chem.Atom(1))
            rw_mol.AddBond(neighbor_idx, h_idx, bond_type)

        # Remove terminal dummy atoms (reverse order to preserve indices)
        dummy_indices = [d[0] for d in terminal_dummies_to_cap]
        for idx in sorted(dummy_indices, reverse=True):
            rw_mol.RemoveAtom(idx)

        chain_mol = rw_mol.GetMol()

        # 2. Assign atom types on CHAIN (now terminal atoms have correct chemistry!)
        chain_mol = self.atom_typer.assign_atom_types(chain_mol)

        # 3. Calculate Gasteiger charges on FULL chain (before splitting!)
        # This ensures atoms have full chemical environment for charge calculation
        # NOTE: Gasteiger doesn't work with dummy atoms, so we use map numbers to track atoms
        if self.force_field.lower() == 'gaff':
            try:
                # Assign map numbers to all atoms for tracking
                for i, atom in enumerate(chain_mol.GetAtoms()):
                    if atom.GetAtomMapNum() == 0:
                        atom.SetAtomMapNum(CHARGE_MAP_START + i)

                # Record which atoms are adjacent to terminal dummies (before removing them)
                # These atoms will later have cap H atoms added by AddHs
                terminal_neighbor_map_nums = set()
                for atom in chain_mol.GetAtoms():
                    if atom.GetIsotope() == TERMINAL_DUMMY_ISOTOPE:
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetAtomicNum() > 0:
                                map_num = neighbor.GetAtomMapNum()
                                if map_num > 0:
                                    terminal_neighbor_map_nums.add(map_num)
                                break

                # Create a version for charge calculation that preserves chemical environment
                # Strategy: Remove terminal dummies, then add explicit H atoms
                # This preserves the chemical environment better than just removing dummies
                rw_mol = Chem.RWMol(chain_mol)

                # Process dummy atoms: remove all dummies (both terminal and connection)
                # We'll add explicit H atoms afterwards to satisfy valence
                dummy_indices = []
                for atom in rw_mol.GetAtoms():
                    if atom.GetAtomicNum() == 0:  # Dummy atom
                        dummy_indices.append(atom.GetIdx())

                # Remove dummies in reverse order to preserve indices
                for idx in sorted(dummy_indices, reverse=True):
                    rw_mol.RemoveAtom(idx)

                # Get the molecule without dummies
                chain_for_charges = rw_mol.GetMol()

                # Sanitize molecule BEFORE AddHs to properly recognize unsatisfied valence
                try:
                    Chem.SanitizeMol(chain_for_charges)
                except Exception as sanitize_error:
                    logger.warning(f"Sanitization warning before AddHs: {sanitize_error}")

                # Add explicit hydrogens to satisfy valence
                # This ensures proper chemical environment for charge calculation
                chain_for_charges = Chem.AddHs(chain_for_charges, addCoords=True)

                # Assign map numbers to cap H atoms added by AddHs
                # Cap H atoms are those added to atoms that were adjacent to terminal dummies
                cap_h_charge_map = {}  # {neighbor_map_num: (cap_h_map_num, charge)}

                for atom in chain_for_charges.GetAtoms():
                    if atom.GetAtomicNum() == 1:  # Hydrogen
                        # Check if this H was added by AddHs (has no map number yet)
                        if atom.GetAtomMapNum() == 0:
                            # Check if this H is bonded to an atom that was adjacent to a terminal dummy
                            for neighbor in atom.GetNeighbors():
                                neighbor_map_num = neighbor.GetAtomMapNum()
                                if neighbor_map_num in terminal_neighbor_map_nums:
                                    # This is a cap H atom - assign it a map number
                                    cap_h_map_num = CAP_H_MAP_START + len(cap_h_charge_map)
                                    atom.SetAtomMapNum(cap_h_map_num)
                                    cap_h_charge_map[neighbor_map_num] = (cap_h_map_num, None)  # Charge added later
                                    break

                # Sanitize molecule before calculating Gasteiger charges
                try:
                    Chem.SanitizeMol(chain_for_charges)
                except Exception as sanitize_error:
                    logger.warning(f"Sanitization warning before charge calc: {sanitize_error}")

                # Calculate Gasteiger charges on the chain with explicit H's
                AllChem.ComputeGasteigerCharges(chain_for_charges)

                # Store cap H charges with neighbor map number as key
                for neighbor_map_num, (cap_h_map_num, _) in list(cap_h_charge_map.items()):
                    cap_h_atom_idx = AtomMapTracker.find_by_map_num(chain_for_charges, cap_h_map_num)
                    if cap_h_atom_idx is not None:
                        cap_h_atom = chain_for_charges.GetAtomWithIdx(cap_h_atom_idx)
                        try:
                            charge = cap_h_atom.GetProp('_GasteigerCharge')
                            cap_h_charge_map[neighbor_map_num] = (cap_h_map_num, charge)
                        except KeyError:
                            logger.warning(f"Cap H atom map={cap_h_map_num} has no charge")

                # Store charges by map number
                charges_by_map = {}
                for atom in chain_for_charges.GetAtoms():
                    map_num = atom.GetAtomMapNum()
                    if map_num > 0:
                        try:
                            charge = atom.GetProp('_GasteigerCharge')
                            charges_by_map[map_num] = charge
                        except KeyError:
                            pass

                # Transfer charges back to original chain using map numbers
                for atom in chain_mol.GetAtoms():
                    map_num = atom.GetAtomMapNum()
                    if map_num in charges_by_map:
                        atom.SetProp('_GasteigerCharge', charges_by_map[map_num])

                # Store cap_h_charge_map as molecule property for split_chain to use
                import json
                cap_h_json = json.dumps({str(k): v for k, v in cap_h_charge_map.items()})
                chain_mol.SetProp('_CapHChargeMap', cap_h_json)

                if self.verbose:
                    logger.info("Calculated Gasteiger charges on full polymer chain")
                    logger.info(f"Preserved {len(cap_h_charge_map)} cap H charges from full chain")
            except Exception as e:
                logger.error(f"Failed to compute Gasteiger charges on chain: {e}")

        # 4. Split into monomers (charges preserved in fragments)
        variants = self.chain_splitter.split_chain(
            chain_mol,
            inter_bond_markers,
            self.base_name,
            self.force_field
        )

        # Clean up cap_h_charge_map property from chain_mol
        if chain_mol.HasProp('_CapHChargeMap'):
            chain_mol.ClearProp('_CapHChargeMap')

        # 5. Align for LT file generation
        aligned_variants = [
            self.backbone_aligner.align_for_lt_file(v)
            for v in variants
        ]
        
        return aligned_variants
    
    def from_chain(
        self,
        chain_mol: Chem.Mol,
        inter_bond_markers: List[Tuple[int, int]]
    ) -> List[MonomerVariant]:
        """
        Split pre-assembled chain into monomer variants.
        
        Use this for chains from external sources (e.g., RadonPy's connect_mols).
        
        Args:
            chain_mol: Polymer chain RDKit Mol object
            inter_bond_markers: Map number pairs for inter-monomer bonds
            
        Returns:
            List of MonomerVariant objects with conformers
        """
        if self.verbose:
            logger.info(f"Processing chain with {len(inter_bond_markers)+1} monomers")
        
        # 1. Assign atom types on chain
        chain_mol = self.atom_typer.assign_atom_types(chain_mol)

        # 2. Calculate Gasteiger charges on FULL chain (before splitting!)
        # This ensures atoms have full chemical environment for charge calculation
        # NOTE: Gasteiger doesn't work with dummy atoms, so we use map numbers to track atoms
        if self.force_field.lower() == 'gaff':
            try:
                # Assign map numbers to all atoms for tracking
                for i, atom in enumerate(chain_mol.GetAtoms()):
                    if atom.GetAtomMapNum() == 0:
                        atom.SetAtomMapNum(CHARGE_MAP_START + i)

                # Record which atoms are adjacent to terminal dummies (before removing them)
                # These atoms will later have cap H atoms added by AddHs
                terminal_neighbor_map_nums = set()
                for atom in chain_mol.GetAtoms():
                    if atom.GetIsotope() == TERMINAL_DUMMY_ISOTOPE:
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetAtomicNum() > 0:
                                map_num = neighbor.GetAtomMapNum()
                                if map_num > 0:
                                    terminal_neighbor_map_nums.add(map_num)
                                break

                # Create a version for charge calculation that preserves chemical environment
                # Strategy: Remove terminal dummies, then add explicit H atoms
                # This preserves the chemical environment better than just removing dummies
                rw_mol = Chem.RWMol(chain_mol)

                # Process dummy atoms: remove all dummies (both terminal and connection)
                # We'll add explicit H atoms afterwards to satisfy valence
                dummy_indices = []
                for atom in rw_mol.GetAtoms():
                    if atom.GetAtomicNum() == 0:  # Dummy atom
                        dummy_indices.append(atom.GetIdx())

                # Remove dummies in reverse order to preserve indices
                for idx in sorted(dummy_indices, reverse=True):
                    rw_mol.RemoveAtom(idx)

                # Get the molecule without dummies
                chain_for_charges = rw_mol.GetMol()

                # Sanitize molecule BEFORE AddHs to properly recognize unsatisfied valence
                try:
                    Chem.SanitizeMol(chain_for_charges)
                except Exception as sanitize_error:
                    logger.warning(f"Sanitization warning before AddHs: {sanitize_error}")

                # Add explicit hydrogens to satisfy valence
                # This ensures proper chemical environment for charge calculation
                chain_for_charges = Chem.AddHs(chain_for_charges, addCoords=True)

                # Assign map numbers to cap H atoms added by AddHs
                # Cap H atoms are those added to atoms that were adjacent to terminal dummies
                cap_h_charge_map = {}  # {neighbor_map_num: (cap_h_map_num, charge)}

                for atom in chain_for_charges.GetAtoms():
                    if atom.GetAtomicNum() == 1:  # Hydrogen
                        # Check if this H was added by AddHs (has no map number yet)
                        if atom.GetAtomMapNum() == 0:
                            # Check if this H is bonded to an atom that was adjacent to a terminal dummy
                            for neighbor in atom.GetNeighbors():
                                neighbor_map_num = neighbor.GetAtomMapNum()
                                if neighbor_map_num in terminal_neighbor_map_nums:
                                    # This is a cap H atom - assign it a map number
                                    cap_h_map_num = CAP_H_MAP_START + len(cap_h_charge_map)
                                    atom.SetAtomMapNum(cap_h_map_num)
                                    cap_h_charge_map[neighbor_map_num] = (cap_h_map_num, None)  # Charge added later
                                    break

                # Sanitize molecule before calculating Gasteiger charges
                try:
                    Chem.SanitizeMol(chain_for_charges)
                except Exception as sanitize_error:
                    logger.warning(f"Sanitization warning before charge calc: {sanitize_error}")

                # Calculate Gasteiger charges on the chain with explicit H's
                AllChem.ComputeGasteigerCharges(chain_for_charges)

                # Store cap H charges with neighbor map number as key
                for neighbor_map_num, (cap_h_map_num, _) in list(cap_h_charge_map.items()):
                    cap_h_atom_idx = AtomMapTracker.find_by_map_num(chain_for_charges, cap_h_map_num)
                    if cap_h_atom_idx is not None:
                        cap_h_atom = chain_for_charges.GetAtomWithIdx(cap_h_atom_idx)
                        try:
                            charge = cap_h_atom.GetProp('_GasteigerCharge')
                            cap_h_charge_map[neighbor_map_num] = (cap_h_map_num, charge)
                        except KeyError:
                            logger.warning(f"Cap H atom map={cap_h_map_num} has no charge")

                # Store charges by map number
                charges_by_map = {}
                for atom in chain_for_charges.GetAtoms():
                    map_num = atom.GetAtomMapNum()
                    if map_num > 0:
                        try:
                            charge = atom.GetProp('_GasteigerCharge')
                            charges_by_map[map_num] = charge
                        except KeyError:
                            pass

                # Transfer charges back to original chain using map numbers
                for atom in chain_mol.GetAtoms():
                    map_num = atom.GetAtomMapNum()
                    if map_num in charges_by_map:
                        atom.SetProp('_GasteigerCharge', charges_by_map[map_num])

                # Store cap_h_charge_map as molecule property for split_chain to use
                import json
                cap_h_json = json.dumps({str(k): v for k, v in cap_h_charge_map.items()})
                chain_mol.SetProp('_CapHChargeMap', cap_h_json)

                if self.verbose:
                    logger.info("Calculated Gasteiger charges on full polymer chain")
                    logger.info(f"Preserved {len(cap_h_charge_map)} cap H charges from full chain")
            except Exception as e:
                logger.error(f"Failed to compute Gasteiger charges on chain: {e}")

        # 3. Split chain
        variants = self.chain_splitter.split_chain(
            chain_mol,
            inter_bond_markers,
            self.base_name,
            self.force_field
        )

        # Clean up cap_h_charge_map property from chain_mol
        if chain_mol.HasProp('_CapHChargeMap'):
            chain_mol.ClearProp('_CapHChargeMap')

        # 4. Align for LT file generation
        aligned_variants = [
            self.backbone_aligner.align_for_lt_file(v)
            for v in variants
        ]
        
        return aligned_variants
    
    def write_lt_files(
        self,
        variants: List[MonomerVariant],
        generate_t1: bool = True
    ) -> List[str]:
        """
        Write .lt files for all variants.
        
        Args:
            variants: List of MonomerVariant objects
            generate_t1: Generate T1 chirality variants
            
        Returns:
            List of generated file paths
        """
        files = []
        for variant in variants:
            variant_files = self.lt_writer.write_variant(
                variant, self.output_dir, generate_t1
            )
            files.extend(variant_files)
        
        if self.verbose:
            logger.info(f"Generated {len(files)} LT files in {self.output_dir}")
        
        return files
    
    def from_single_molecule(
        self,
        smiles: str,
        molecule_name: Optional[str] = None
    ) -> MonomerVariant:
        """
        Generate LT file for a single molecule (no wildcards).
        
        This is for non-polymer molecules like solvents, additives, or
        small molecules that need to be included in LAMMPS simulations.
        
        Args:
            smiles: SMILES string WITHOUT wildcards (e.g., "CCO" for ethanol)
            molecule_name: Optional name for the molecule (defaults to base_name)
            
        Returns:
            MonomerVariant object representing the molecule
            
        Example:
            >>> generator = MonomerGenerator(base_name="ethanol", force_field="gaff")
            >>> variant = generator.from_single_molecule("CCO")
            >>> generator.write_single_molecule(variant)
        """
        if molecule_name is None:
            molecule_name = self.base_name
        
        if self.verbose:
            logger.info(f"Processing single molecule: {smiles}")
        
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise MonomerGeneratorError(f"Invalid SMILES: {smiles}")
        
        # Check for wildcards - this function is for non-polymer molecules
        has_wildcard = any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms())
        if has_wildcard:
            raise MonomerGeneratorError(
                f"SMILES contains wildcards. Use from_smiles() for polymer monomers. "
                f"Got: {smiles}"
            )
        
        # Add explicit hydrogens
        mol = Chem.AddHs(mol)
        
        # Assign atom types
        mol = self.atom_typer.assign_atom_types(mol)
        
        # Generate conformer
        mol = self.conformer_gen.generate_conformer(mol)
        
        # Create MonomerVariant (using 'single' type for standalone molecules)
        variant = MonomerVariant(
            base_name=molecule_name,
            variant_type='molecule',  # Special type for non-polymer molecules
            mol=mol,
            smiles=Chem.MolToSmiles(mol),
            connection_atoms=(0, 0),  # No connection points
            force_field=self.force_field,
            position=0
        )
        
        if self.verbose:
            logger.info(f"Created molecule variant: {molecule_name} with {mol.GetNumAtoms()} atoms")
        
        return variant
    
    def write_single_molecule(
        self,
        variant: MonomerVariant,
        generate_t1: bool = False
    ) -> List[str]:
        """
        Write LT file for a single molecule.
        
        Args:
            variant: MonomerVariant object from from_single_molecule()
            generate_t1: Generate T1 chirality variant (usually False for molecules)
            
        Returns:
            List of generated file paths
        """
        files = []
        
        # Generate filename
        filename = f"{variant.base_name}.lt"
        filepath = os.path.join(self.output_dir, filename)
        
        # Write the LT file using a specialized writer for molecules
        self._write_molecule_lt_file(variant, filepath)
        files.append(filepath)
        
        if self.verbose:
            logger.info(f"Generated: {filepath}")
        
        if generate_t1:
            variant_t1 = self.lt_writer._create_t1_variant(variant)
            filepath_t1 = os.path.join(self.output_dir, f"{variant.base_name}_T1.lt")
            self._write_molecule_lt_file(variant_t1, filepath_t1)
            files.append(filepath_t1)
            
            if self.verbose:
                logger.info(f"Generated: {filepath_t1}")
        
        return files
    
    def _write_molecule_lt_file(self, variant: MonomerVariant, filepath: str) -> None:
        """Write LT file for a single molecule (no connection atoms)."""
        mol = variant.mol
        
        with open(filepath, 'w') as f:
            # Write header
            if variant.force_field == 'gaff':
                f.write('import "gaff.lt"    # <-- defines the GAFF (General Amber Force Field)\n')
                f.write('# NOTE: GAFF requires user-supplied charges (AM1-BCC or RESP recommended)\n')
                f.write('# See: http://ambermd.org/antechamber/gaff.pdf\n')
                f.write(f'{variant.base_name} inherits GAFF {{\n\n')
            else:
                f.write('import "oplsaa.lt"    # <-- defines the OPLS-AA force field\n')
                f.write(f'{variant.base_name} inherits OPLSAA {{\n\n')
            
            f.write('# atom-id  mol-id  atom-type charge      X         Y        Z\n\n')
            
            # Write atoms block (heavy atoms first, then hydrogens)
            f.write('  write("Data Atoms") {\n')
            
            # Order atoms: heavy atoms first, then hydrogens
            atom_order = []
            h_atoms = []
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                if atom.GetAtomicNum() == 1:
                    h_atoms.append(i)
                else:
                    atom_order.append(i)
            atom_order.extend(h_atoms)
            
            # Generate atom IDs
            atom_ids = {}
            for seq_num, atom_idx in enumerate(atom_order, 1):
                atom = mol.GetAtomWithIdx(atom_idx)
                element = atom.GetSymbol()
                atom_ids[atom_idx] = f"{element}{seq_num}"
            
            # Get conformer
            conf = mol.GetConformer(0) if mol.GetNumConformers() > 0 else None
            
            for atom_idx in atom_order:
                atom = mol.GetAtomWithIdx(atom_idx)
                atom_id = atom_ids[atom_idx]
                
                try:
                    atom_type = atom.GetProp('AtomType')
                except KeyError:
                    atom_type = f"@atom:{atom.GetSymbol().lower()}"
                
                charge = self.lt_writer.charge_dict.get(atom_type, 0.0)
                
                if conf is not None:
                    pos = conf.GetAtomPosition(atom_idx)
                    x, y, z = pos.x, pos.y, pos.z
                else:
                    x, y, z = 0.0, 0.0, 0.0
                
                f.write(f'\t$atom:{atom_id} $mol:... {atom_type} {charge:.4f}')
                f.write(f'    {x:.3f}   {y:.3f}   {z:.3f}\n')
            
            f.write('  }\n\n')
            
            # Write bonds block
            f.write("  write('Data Bond List') {\n")
            
            for bond in mol.GetBonds():
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                id1 = atom_ids.get(a1, f"X{a1}")
                id2 = atom_ids.get(a2, f"X{a2}")
                f.write(f'\t$bond:{id1}{id2}\t$atom:{id1}\t$atom:{id2}\n')
            
            f.write('  }\n')
            
            # Write footer
            f.write(f'}}  # {variant.base_name}\n\n')
            f.write("# Note: You don't need to supply the partial partial charges of the atoms.\n")
            f.write("#       If you like, just fill the fourth column with zeros (\"0.000\").\n")
            f.write("#       Moltemplate and LAMMPS will automatically assign the charge later\n\n")


# =============================================================================
# Convenience functions
# =============================================================================

def generate_monomers(
    smiles: str,
    base_name: str,
    force_field: str = 'gaff',
    output_dir: str = './monomers',
    n_monomers: int = 3,
    generate_t1: bool = True,
    verbose: bool = True
) -> List[str]:
    """
    Convenience function to generate monomer LT files from SMILES.
    
    Args:
        smiles: SMILES string with wildcards (e.g., "[*]CC[*]")
        base_name: Base name for monomers (e.g., "PE")
        force_field: 'oplsaa' or 'gaff'
        output_dir: Directory for output files
        n_monomers: Number of monomers for chain building
        generate_t1: Generate T1 chirality variants
        verbose: Enable verbose output
        
    Returns:
        List of generated file paths
        
    Example:
        >>> files = generate_monomers(
        ...     smiles="[*]CC[*]",
        ...     base_name="PE",
        ...     force_field="gaff",
        ...     output_dir="./pe_monomers"
        ... )
    """
    generator = MonomerGenerator(
        base_name=base_name,
        force_field=force_field,
        output_dir=output_dir,
        verbose=verbose
    )
    
    variants = generator.from_smiles(smiles, n_monomers=n_monomers)
    files = generator.write_lt_files(variants, generate_t1=generate_t1)
    
    return files


def generate_molecule_lt(
    smiles: str,
    molecule_name: str,
    force_field: str = 'gaff',
    output_dir: str = './molecules',
    generate_t1: bool = False,
    verbose: bool = True
) -> List[str]:
    """
    Generate LT file for a single molecule (no wildcards, non-polymer).
    
    Use this for solvents, additives, or any small molecule that needs
    to be included in LAMMPS simulations via AutoPoly/Moltemplate.
    
    Args:
        smiles: SMILES string WITHOUT wildcards (e.g., "CCO" for ethanol)
        molecule_name: Name for the molecule (used in LT file)
        force_field: 'oplsaa' or 'gaff'
        output_dir: Directory for output files
        generate_t1: Generate T1 chirality variant
        verbose: Enable verbose output
        
    Returns:
        List of generated file paths
        
    Example:
        >>> # Generate ethanol LT file
        >>> files = generate_molecule_lt(
        ...     smiles="CCO",
        ...     molecule_name="ethanol",
        ...     force_field="gaff",
        ...     output_dir="./molecules"
        ... )
        
        >>> # Generate water (TIP3P compatible structure)
        >>> files = generate_molecule_lt(
        ...     smiles="O",
        ...     molecule_name="water",
        ...     force_field="gaff"
        ... )
        
        >>> # Generate acetone
        >>> files = generate_molecule_lt(
        ...     smiles="CC(=O)C",
        ...     molecule_name="acetone",
        ...     force_field="gaff"
        ... )
    """
    generator = MonomerGenerator(
        base_name=molecule_name,
        force_field=force_field,
        output_dir=output_dir,
        verbose=verbose
    )
    
    variant = generator.from_single_molecule(smiles, molecule_name)
    files = generator.write_single_molecule(variant, generate_t1=generate_t1)
    
    return files


# =============================================================================
# Main entry point
# =============================================================================

if __name__ == '__main__':
    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Example usage
    print("Monomer Generator for AutoPoly")
    print("=" * 50)
    
    # Example: Polyethylene
    try:
        files = generate_monomers(
            smiles="[*]CC[*]",
            base_name="PE",
            force_field="gaff",
            output_dir="./test_monomers",
            n_monomers=3,
            verbose=True
        )
        print(f"\nGenerated files: {files}")
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

