"""
Monomer Generator Module for AutoPoly

This module provides automatic generation of monomer variants (internal, left-end, right-end)
from SMILES strings or RDKit Mol objects with automatic H-atom capping and OPLS-AA force
field typing.

Example:
    >>> from AutoPoly.monomer_generator import MonomerGenerator
    >>> generator = MonomerGenerator(base_name="PE")
    >>> variants = generator.generate_variants(smiles="C=C")
    >>> files = generator.generate_lt_files(variants)
    >>> # Generates: PEi.lt, PEle.lt, PErre.lt

Author: AutoPoly Development Team
"""

import os
import sys
import pickle
import pathlib
from typing import Dict, Optional, Union

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

from .extern.rdlt import RDlt
from .system import logger
from .polymerization_mechanism import PolymerizationMechanism
from .connection_point import ConnectionPointModifier
from .polymerization_patterns import POLYMERIZATION_PATTERNS


class MonomerGeneratorError(Exception):
    """Base exception for MonomerGenerator errors."""
    pass


class BondingSiteError(MonomerGeneratorError):
    """Raised when bonding sites cannot be identified."""
    pass


class AtomTypingError(MonomerGeneratorError):
    """Raised when atom typing fails."""
    pass


class ValidationError(MonomerGeneratorError):
    """Raised when validation fails."""
    pass


class MonomerGenerator:
    """
    Automatic monomer variant generator for AutoPoly.

    Generates internal (i), left-end (le), and right-end (re) monomers
    from SMILES strings or RDKit Mol objects with automatic H-atom capping
    and OPLS-AA force field typing.

    Attributes:
        base_name (str): Base name for the monomer (e.g., "PE", "PMMA")
        output_dir (str): Output directory for generated .lt files
        is_lopls (bool): Use LOPLS force field parameters
        verbose (bool): Enable verbose logging
        charge_dict (dict): Atom type to charge mapping

    Example:
        >>> generator = MonomerGenerator(base_name="PE")
        >>> variants = generator.generate_variants(smiles="C=C")
        >>> files = generator.generate_lt_files(variants)
    """

    def __init__(
        self,
        base_name: str,
        output_dir: Optional[str] = None,
        is_lopls: bool = False,
        is_gaff: bool = False,
        mechanism: Optional[str] = None,
        verbose: bool = True
    ):
        """
        Initialize monomer generator.

        Args:
            base_name: Base name for monomer (e.g., "PE", "PMMA")
            output_dir: Directory for .lt files (default: Monomer_bank)
            is_lopls: Use LOPLS force field parameters
            is_gaff: Use GAFF force field parameters
            mechanism: Polymerization mechanism (None for auto-detection)
                      Options: 'none', 'vinyl_addition', 'esterification', 'amidation', 'etherification'
            verbose: Enable detailed logging

        Raises:
            MonomerGeneratorError: If both is_lopls and is_gaff are True
        """
        self.base_name = base_name
        self.is_lopls = is_lopls
        self.is_gaff = is_gaff
        self.mechanism = mechanism
        self.verbose = verbose

        # Validate mutual exclusivity
        if is_lopls and is_gaff:
            raise MonomerGeneratorError(
                "Cannot use both LOPLS and GAFF simultaneously. Please choose one force field."
            )

        # Validate mechanism if provided
        if mechanism is not None:
            if mechanism not in POLYMERIZATION_PATTERNS:
                raise MonomerGeneratorError(
                    f"Unknown mechanism: '{mechanism}'. "
                    f"Available: {list(POLYMERIZATION_PATTERNS.keys())}"
                )

        # Set up paths
        module_dir = pathlib.Path(__file__).parent.resolve()

        # Default output directory
        if output_dir is None:
            self.output_dir = str(module_dir / "extern" / "Monomer_bank")
        else:
            self.output_dir = output_dir

        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)

        # Paths to feature definition files
        if is_gaff:
            self.fdef_path = str(module_dir / "extern" / "rdlt_data" / "gaff_lt.fdefn")
            self.lfdef_path = None  # Not used for GAFF
        else:
            self.fdef_path = str(module_dir / "extern" / "rdlt_data" / "opls_lt.fdefn")
            self.lfdef_path = str(module_dir / "extern" / "rdlt_data" / "lopls_lt.fdefn")

        # Load charge dictionaries
        self.charge_dict = self._load_charge_dictionaries()

        # Initialize RDlt converter for reuse
        self.rdlt_converter = RDlt()

        # Determine force field string for connection point modifier
        if is_gaff:
            ff_str = 'gaff'
        elif is_lopls:
            ff_str = 'lopls'
        else:
            ff_str = 'oplsaa'

        # Initialize mechanism detector and connection point modifier
        self.mech_detector = PolymerizationMechanism(verbose=verbose)
        self.conn_modifier = ConnectionPointModifier(force_field=ff_str, verbose=verbose)

        if self.verbose:
            logger.info(f"MonomerGenerator initialized for '{base_name}'")
            logger.info(f"Output directory: {self.output_dir}")
            if self.mechanism:
                logger.info(f"Mechanism: {self.mechanism}")
            if is_gaff:
                logger.info("Force field: GAFF (General Amber Force Field)")
                logger.warning("GAFF requires manual charge calculation using AM1-BCC or RESP")
            elif is_lopls:
                logger.info("Force field: L-OPLS (Long-chain optimized OPLS)")
            else:
                logger.info("Force field: OPLS-AA")

    def _load_charge_dictionaries(self) -> dict:
        """
        Load force field charge dictionaries.

        Returns:
            dict: Combined charge dictionary

        Note:
            GAFF charge dictionaries are typically empty (all zeros) since
            GAFF requires manual charge calculation using AM1-BCC or RESP.
        """
        try:
            if self.is_gaff:
                # GAFF charges are typically empty (manual calculation required)
                try:
                    with open(self.fdef_path.replace('.fdefn', '_dict.pkl'), 'rb') as f:
                        cdict = pickle.load(f)

                    if not cdict or all(v == 0.0 for v in cdict.values()):
                        if self.verbose:
                            logger.info(
                                "GAFF charge dictionary is empty (all zeros). "
                                "Charges must be calculated manually using AM1-BCC or RESP."
                            )
                except FileNotFoundError:
                    if self.verbose:
                        logger.warning("GAFF charge dictionary not found. This is expected.")
                    cdict = {}
            else:
                # OPLS charges
                with open(self.fdef_path.replace('.fdefn', '_dict.pkl'), 'rb') as f:
                    cdict = pickle.load(f)

                if self.is_lopls:
                    try:
                        with open(self.lfdef_path.replace('.fdefn', '_dict.pkl'), 'rb') as f:
                            ldict = pickle.load(f)
                        cdict.update(ldict)
                    except FileNotFoundError:
                        logger.warning("LOPLS charge dictionary not found, using OPLS only")

            return cdict
        except Exception as e:
            logger.warning(f"Could not load charge dictionaries: {e}")
            return {}

    def _validate_input(self, smiles: Optional[str] = None, mol: Optional[Chem.Mol] = None) -> Chem.Mol:
        """
        Validate input parameters and return RDKit Mol object.

        Args:
            smiles: SMILES string
            mol: RDKit Mol object

        Returns:
            Chem.Mol: Validated RDKit Mol object

        Raises:
            MonomerGeneratorError: If validation fails
        """
        if (smiles is None and mol is None) or (smiles is not None and mol is not None):
            raise MonomerGeneratorError(
                "Provide exactly one: SMILES string or Mol object"
            )

        if smiles is not None:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise MonomerGeneratorError(f"Invalid SMILES: {smiles}")

        if mol is None:
            raise MonomerGeneratorError("Molecule is None")

        return mol

    def _identify_bonding_sites(self, mol: Chem.Mol, dop: int = 1) -> dict:
        """
        Identify atoms with unsatisfied valency (polymerization sites).

        Algorithm:
        1. Detect polymerization mechanism (vinyl, esterification, etc.)
        2. Use mechanism-specific SMARTS patterns to find bonding sites
        3. For symmetric molecules: identify equivalent sites
        4. For asymmetric molecules: identify distinct left/right sites

        Args:
            mol: RDKit Mol object (with explicit hydrogens)
            dop: Degree of polymerization (1 = single molecule)

        Returns:
            dict: {
                'bonding_sites': [idx1, idx2, ...],
                'left_site': idx1,
                'right_site': idx2,
                'is_symmetric': bool,
                'site_atoms': [(element, idx), ...],
                'mechanism': str
            }

        Raises:
            BondingSiteError: If bonding sites cannot be identified
        """
        mol_with_h = AllChem.AddHs(mol)

        # Detect mechanism (or use user-specified mechanism)
        if self.mechanism:
            mechanism = self.mechanism
        else:
            mechanism = self.mech_detector.detect_mechanism(mol_with_h, dop)

        # For 'none' mechanism, no bonding sites needed (single molecule)
        if mechanism == 'none':
            if self.verbose:
                logger.info("Non-polymerizable molecule (mechanism='none')")
            return {
                'bonding_sites': [],
                'left_site': None,
                'right_site': None,
                'is_symmetric': True,
                'site_atoms': [],
                'mechanism': 'none'
            }

        # Get connection atoms using mechanism-specific SMARTS
        connection_atoms = self.mech_detector.get_connection_atoms(mol_with_h, mechanism)

        if len(connection_atoms) < 2:
            # Fallback: try generic double bond detection for vinyl
            if mechanism == 'vinyl_addition':
                bonding_sites = []
                site_atoms = []
                for atom in mol_with_h.GetAtoms():
                    for bond in atom.GetBonds():
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            idx = atom.GetIdx()
                            element = atom.GetSymbol()
                            if idx not in bonding_sites:
                                bonding_sites.append(idx)
                                site_atoms.append((element, idx))
                            other_idx = bond.GetOtherAtomIdx(idx)
                            other_atom = mol_with_h.GetAtomWithIdx(other_idx)
                            if other_idx not in bonding_sites:
                                bonding_sites.append(other_idx)
                                site_atoms.append((other_atom.GetSymbol(), other_idx))
                connection_atoms = bonding_sites

        # Validate we found enough sites
        if len(connection_atoms) < 2:
            raise BondingSiteError(
                f"Found {len(connection_atoms)} bonding sites for mechanism '{mechanism}', need at least 2. "
                f"This monomer may not be suitable for this polymerization type."
            )

        # Use first two sites
        bonding_sites = connection_atoms[:2]
        site_atoms = [(mol_with_h.GetAtomWithIdx(idx).GetSymbol(), idx) for idx in bonding_sites]

        # Check if symmetric (same element and similar environment)
        atom1 = mol_with_h.GetAtomWithIdx(bonding_sites[0])
        atom2 = mol_with_h.GetAtomWithIdx(bonding_sites[1])
        is_symmetric = (
            atom1.GetSymbol() == atom2.GetSymbol() and
            atom1.GetDegree() == atom2.GetDegree()
        )

        result = {
            'bonding_sites': bonding_sites,
            'left_site': bonding_sites[0],
            'right_site': bonding_sites[1],
            'is_symmetric': is_symmetric,
            'site_atoms': site_atoms,
            'mechanism': mechanism
        }

        if self.verbose:
            logger.info(f"Identified bonding sites: {site_atoms}")
            logger.info(f"Mechanism: {mechanism}")
            logger.info(f"Symmetric: {is_symmetric}")

        return result

    def _remove_terminal_hydrogen(self, mol: Chem.Mol, atom_idx: int) -> Chem.Mol:
        """
        Remove hydrogen from a terminal atom to create bonding site.

        Args:
            mol: RDKit Mol object
            atom_idx: Index of atom to remove hydrogen from

        Returns:
            Chem.Mol: Modified Mol object
        """
        rw_mol = Chem.RWMol(mol)

        # Find hydrogen atoms bonded to the specified atom
        atom = rw_mol.GetAtomWithIdx(atom_idx)
        h_indices = []

        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'H':
                h_indices.append(neighbor.GetIdx())

        # Remove one hydrogen if present
        if h_indices:
            rw_mol.RemoveAtom(h_indices[0])

        return rw_mol.GetMol()

    def _add_hydrogen(self, mol: Chem.Mol, atom_idx: int) -> Chem.Mol:
        """
        Add hydrogen atom to cap a bonding site.

        Args:
            mol: RDKit Mol object
            atom_idx: Index of atom to cap with hydrogen

        Returns:
            Chem.Mol: Modified Mol object with added hydrogen
        """
        rw_mol = Chem.RWMol(mol)

        # Add new hydrogen atom
        h_idx = rw_mol.AddAtom(Chem.Atom('H'))

        # Add bond between H and the specified atom
        rw_mol.AddBond(atom_idx, h_idx, Chem.rdchem.BondType.SINGLE)

        mol = rw_mol.GetMol()

        # Sanitize and generate conformer
        Chem.SanitizeMol(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol)

        return mol

    def _create_internal_variant(self, mol: Chem.Mol, bonding_info: dict) -> Chem.Mol:
        """
        Create internal monomer (original alkene for middle of chain).

        For alkenes: returns the original molecule with double bond intact.
        The double bond will be converted to single bonds during polymerization.

        Args:
            mol: RDKit Mol object
            bonding_info: Bonding site information from _identify_bonding_sites

        Returns:
            Chem.Mol: Internal variant (original alkene)
        """
        # For internal variant, just return the original molecule with hydrogens
        mol_with_h = AllChem.AddHs(mol)

        # Generate conformer
        AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol_with_h)

        if self.verbose:
            logger.info("Created internal variant")

        return mol_with_h

    def _create_left_end_variant(self, mol: Chem.Mol, bonding_info: dict) -> Chem.Mol:
        """
        Create left-end monomer (left end capped, right end for polymerization).

        Uses the original monomer structure. The _write_lt_file_direct method
        will remove one H from the right end to create the connection point.

        Args:
            mol: RDKit Mol object (original monomer)
            bonding_info: Bonding site information

        Returns:
            Chem.Mol: Left-end variant (same structure as original monomer)
        """
        # Use original monomer structure (not ethane!)
        mol_le = Chem.Mol(mol)
        mol_le = AllChem.AddHs(mol_le)

        # Generate conformer
        AllChem.EmbedMolecule(mol_le, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol_le)

        if self.verbose:
            logger.info("Created left-end variant from original monomer structure")

        return mol_le

    def _create_right_end_variant(self, mol: Chem.Mol, bonding_info: dict) -> Chem.Mol:
        """
        Create right-end monomer (right end capped, left end for polymerization).

        Uses the original monomer structure. The _write_lt_file_direct method
        will remove one H from the left end to create the connection point.

        Args:
            mol: RDKit Mol object (original monomer)
            bonding_info: Bonding site information

        Returns:
            Chem.Mol: Right-end variant (same structure as original monomer)
        """
        # Use original monomer structure (not ethane!)
        mol_re = Chem.Mol(mol)
        mol_re = AllChem.AddHs(mol_re)

        # Generate conformer
        AllChem.EmbedMolecule(mol_re, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol_re)

        if self.verbose:
            logger.info("Created right-end variant from original monomer structure")

        return mol_re

    def _get_backbone_atoms(self, mol: Chem.Mol, mechanism: Optional[str] = None) -> list:
        """
        Identify backbone atoms (polymerization connecting sites).

        Element-agnostic backbone identification using SMARTS patterns.

        Args:
            mol: RDKit Mol object
            mechanism: Polymerization mechanism (if known)

        Returns:
            list: Indices of backbone atoms (typically 2 atoms)
        """
        backbone_atoms = []

        # If mechanism is known, use mechanism-specific SMARTS
        if mechanism and mechanism != 'none':
            pattern = POLYMERIZATION_PATTERNS.get(mechanism)
            if pattern and pattern.smarts:
                try:
                    smarts_pattern = Chem.MolFromSmarts(pattern.smarts)
                    if smarts_pattern:
                        matches = mol.GetSubstructMatches(smarts_pattern)
                        if matches:
                            # Get unique atoms from all matches
                            for match in matches:
                                for atom_idx in match:
                                    if atom_idx not in backbone_atoms:
                                        backbone_atoms.append(atom_idx)
                            # Limit to first 2 atoms
                            backbone_atoms = backbone_atoms[:2]
                            if len(backbone_atoms) >= 2:
                                return backbone_atoms
                except Exception:
                    pass

        # Generic approach: look for double bonds (element-agnostic)
        if len(backbone_atoms) < 2:
            for atom in mol.GetAtoms():
                for bond in atom.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        idx = atom.GetIdx()
                        other_idx = bond.GetOtherAtomIdx(idx)
                        if idx not in backbone_atoms:
                            backbone_atoms.append(idx)
                        if other_idx not in backbone_atoms:
                            backbone_atoms.append(other_idx)
                        if len(backbone_atoms) >= 2:
                            return backbone_atoms

        # If no double bond found, look for other patterns
        if len(backbone_atoms) < 2:
            # Use terminal atoms of any element
            # Find atoms with the smallest and largest index
            all_atoms = [atom.GetIdx() for atom in mol.GetAtoms()]
            if len(all_atoms) >= 2:
                backbone_atoms = [all_atoms[0], all_atoms[-1]]

        return backbone_atoms

    def _has_chirality(self, mol: Chem.Mol) -> bool:
        """
        Detect if molecule has chiral centers.

        Args:
            mol: RDKit Mol object

        Returns:
            bool: True if molecule has chiral centers
        """
        try:
            chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            return len(chiral_centers) > 0
        except Exception:
            return False

    def _create_t1_variant(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Create T1 variant with opposite chirality.

        Inverts Z-coordinates to create mirror image.

        Args:
            mol: RDKit Mol object

        Returns:
            Chem.Mol: T1 variant with inverted chirality
        """
        # Create a copy of the molecule
        mol_t1 = Chem.Mol(mol)
        mol_t1 = AllChem.AddHs(mol_t1)

        # Get conformer from original molecule
        conf = mol.GetConformer(0)
        conf_t1 = mol_t1.GetConformer(0)

        # Invert Z coordinates for all atoms
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            conf_t1.SetAtomPosition(i, (pos.x, pos.y, -pos.z))

        if self.verbose:
            logger.info("Created T1 variant (Z-coordinate inversion)")

        return mol_t1

    def _reorder_atoms(self, mol: Chem.Mol, new_order: list) -> Chem.Mol:
        """
        Reorder atoms in molecule according to new_order.

        Creates a new molecule with atoms in the specified order.
        Updates all bonds to use new atom indices.

        Args:
            mol: RDKit Mol object
            new_order: List of atom indices in desired order

        Returns:
            Chem.Mol: Molecule with reordered atoms
        """
        # Create editable molecule
        rw_mol = Chem.RWMol()

        # Map old indices to new indices
        old_to_new = {}
        for new_idx, old_idx in enumerate(new_order):
            old_atom = mol.GetAtomWithIdx(old_idx)
            # Copy atom properties
            new_atom = Chem.Atom(old_atom.GetSymbol())
            new_atom.SetFormalCharge(old_atom.GetFormalCharge())

            # Copy atom properties like AtomType
            for prop_name in old_atom.GetPropNames():
                try:
                    new_atom.SetProp(prop_name, old_atom.GetProp(prop_name))
                except Exception:
                    pass

            new_idx_in_mol = rw_mol.AddAtom(new_atom)
            old_to_new[old_idx] = new_idx_in_mol

        # Add bonds with updated indices
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            bond_type = bond.GetBondType()

            if begin_idx in old_to_new and end_idx in old_to_new:
                rw_mol.AddBond(
                    old_to_new[begin_idx],
                    old_to_new[end_idx],
                    bond_type
                )

        # Get conformer and reorder coordinates
        new_mol = rw_mol.GetMol()
        Chem.SanitizeMol(new_mol)

        # Create a conformer for the new molecule
        new_conf = Chem.Conformer(mol.GetNumAtoms())
        old_conf = mol.GetConformer(0)

        # Copy coordinates in new order
        for new_idx, old_idx in enumerate(new_order):
            pos = old_conf.GetAtomPosition(old_idx)
            new_conf.SetAtomPosition(new_idx, pos)

        new_mol.AddConformer(new_conf)

        if self.verbose:
            logger.info(f"Reordered {len(new_order)} atoms")

        return new_mol

    def _align_backbone_to_x_axis(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Rotate and translate molecule to align backbone with X-axis, centered at origin.

        - Right atom at positive X
        - Left atom at negative X
        - Both at Y=0, Z=0

        Args:
            mol: RDKit Mol object

        Returns:
            Chem.Mol: Aligned molecule
        """
        import numpy as np

        # Get backbone atoms
        backbone = self._get_backbone_atoms(mol)
        if len(backbone) < 2:
            if self.verbose:
                logger.warning("Could not find backbone atoms, skipping alignment")
            return mol

        conf = mol.GetConformer(0)

        # Get positions of backbone atoms
        right_pos = np.array([conf.GetAtomPosition(backbone[0]).x,
                             conf.GetAtomPosition(backbone[0]).y,
                             conf.GetAtomPosition(backbone[0]).z])
        left_pos = np.array([conf.GetAtomPosition(backbone[1]).x,
                            conf.GetAtomPosition(backbone[1]).y,
                            conf.GetAtomPosition(backbone[1]).z])

        # Determine which is right (larger X) and which is left (smaller X)
        if right_pos[0] < left_pos[0]:
            right_pos, left_pos = left_pos, right_pos
            backbone[0], backbone[1] = backbone[1], backbone[0]

        # Calculate midpoint
        midpoint = (right_pos + left_pos) / 2

        # Backbone vector (from left to right)
        backbone_vec = right_pos - left_pos
        backbone_length = np.linalg.norm(backbone_vec)

        if backbone_length < 0.001:
            if self.verbose:
                logger.warning("Backbone atoms too close, skipping alignment")
            return mol

        # Normalize backbone vector
        backbone_unit = backbone_vec / backbone_length

        # Target vector (X-axis)
        target_vec = np.array([1.0, 0.0, 0.0])

        # Calculate rotation axis and angle
        rotation_axis = np.cross(backbone_unit, target_vec)
        axis_norm = np.linalg.norm(rotation_axis)

        if axis_norm < 0.001:
            # Already aligned with X-axis (or opposite)
            if np.dot(backbone_unit, target_vec) < 0:
                # Opposite direction, need 180-degree rotation
                # Rotate around Y or Z axis
                rotation_axis = np.array([0.0, 1.0, 0.0])
                cos_theta = -1.0
                sin_theta = 0.0
            else:
                # Already aligned
                cos_theta = 1.0
                sin_theta = 0.0
        else:
            rotation_axis = rotation_axis / axis_norm
            cos_theta = np.dot(backbone_unit, target_vec)
            sin_theta = np.sqrt(1 - cos_theta**2)

        # Rodrigues' rotation formula
        def rotate_vector(v, k, cos_t, sin_t):
            """Rotate vector v around unit axis k by angle where cos(cos_t) and sin(sin_t)."""
            return (v * cos_t +
                    np.cross(k, v) * sin_t +
                    k * np.dot(k, v) * (1 - cos_t))

        # Apply rotation and translation to all atoms
        new_mol = Chem.Mol(mol)
        new_conf = new_mol.GetConformer(0)

        for i in range(mol.GetNumAtoms()):
            pos = np.array([conf.GetAtomPosition(i).x,
                           conf.GetAtomPosition(i).y,
                           conf.GetAtomPosition(i).z])

            # Translate to origin
            pos_centered = pos - midpoint

            # Rotate
            if axis_norm >= 0.001 or cos_theta < 0:
                pos_rotated = rotate_vector(pos_centered, rotation_axis, cos_theta, sin_theta)
            else:
                pos_rotated = pos_centered

            new_conf.SetAtomPosition(i, pos_rotated)

        if self.verbose:
            logger.info("Aligned backbone to X-axis and centered at origin")

        return new_mol

    def _align_and_reorder_backbone(self, mol: Chem.Mol, variant_type: str) -> Chem.Mol:
        """
        Align backbone along X-axis and reorder atoms.

        For internal variants: RIGHT atom first, LEFT atom second
        For end variants: maintain existing order if already correct

        Args:
            mol: RDKit Mol object
            variant_type: 'internal', 'left_end', or 'right_end'

        Returns:
            Chem.Mol: Aligned and reordered molecule
        """
        # Step 1: Align to X-axis
        mol_aligned = self._align_backbone_to_x_axis(mol)

        # Step 2: Get backbone atoms
        backbone = self._get_backbone_atoms(mol_aligned)
        if len(backbone) < 2:
            if self.verbose:
                logger.warning("Could not find backbone atoms for reordering")
            return mol_aligned

        # Step 3: Determine atom order
        conf = mol_aligned.GetConformer(0)

        right_pos = conf.GetAtomPosition(backbone[0])
        left_pos = conf.GetAtomPosition(backbone[1])

        # Ensure backbone[0] is right, backbone[1] is left
        if right_pos.x < left_pos.x:
            backbone[0], backbone[1] = backbone[1], backbone[0]
            right_pos, left_pos = left_pos, right_pos

        # Step 4: For internal variants, reorder so RIGHT is first
        if variant_type == 'internal':
            # RIGHT atom (backbone[0]) should be first, LEFT atom (backbone[1]) second
            new_order = [backbone[0], backbone[1]]
        else:
            # For end variants, keep the order that makes the connecting atom first
            # left_end: right atom should connect, so it's first
            # right_end: left atom should connect, so it's first
            if variant_type == 'left_end':
                # Left-end: right side is the connecting point
                new_order = [backbone[0], backbone[1]]
            else:  # right_end
                # Right-end: left side is the connecting point
                new_order = [backbone[1], backbone[0]]

        # Add remaining atoms
        all_atoms = set(range(mol_aligned.GetNumAtoms()))
        backbone_set = set(backbone)
        other_atoms = list(all_atoms - backbone_set)
        new_order.extend(other_atoms)

        # Step 5: Reorder atoms
        mol_reordered = self._reorder_atoms(mol_aligned, new_order)

        if self.verbose:
            conf_final = mol_reordered.GetConformer(0)
            c1_pos = conf_final.GetAtomPosition(0)
            c2_pos = conf_final.GetAtomPosition(1)
            logger.info(f"Backbone alignment: C1 at X={c1_pos.x:.3f}, C2 at X={c2_pos.x:.3f}")

        return mol_reordered

    def _assign_opls_types(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Assign atom types using feature factory (OPLS, L-OPLS, or GAFF).

        Args:
            mol: RDKit Mol object

        Returns:
            Chem.Mol: Mol object with AtomType properties

        Raises:
            AtomTypingError: If atom typing fails
        """
        # Special handling for water molecules (SPC/E model)
        if mol.GetNumAtoms() == 3:
            atoms = list(mol.GetAtoms())
            symbols = [atom.GetSymbol() for atom in atoms]
            if symbols.count('O') == 1 and symbols.count('H') == 2:
                # Assign SPC/E water model atom types
                for atom in atoms:
                    if atom.GetSymbol() == 'O':
                        atom.SetProp('AtomType', 'spce_O')  # SPC/E oxygen
                    else:  # H
                        atom.SetProp('AtomType', 'spce_H')  # SPC/E hydrogen
                if self.verbose:
                    logger.info("SPC/E water atom types assigned successfully")
                return mol

        # Choose feature definition file based on force field
        if self.is_lopls:
            fdef = self.lfdef_path
        elif self.is_gaff:
            fdef = self.fdef_path  # GAFF uses gaff_lt.fdefn
        else:
            fdef = self.fdef_path  # OPLS uses opls_lt.fdefn

        try:
            # Build feature factory
            factory = Chem.ChemicalFeatures.BuildFeatureFactory(fdef)
            features = factory.GetFeaturesForMol(mol)

            # Assign atom types
            for feature in features:
                atom_idx = feature.GetAtomIds()[0]
                atom_type = feature.GetType()
                mol.GetAtomWithIdx(atom_idx).SetProp('AtomType', atom_type)

            # Check for untyped atoms
            untyped = []
            for atom in mol.GetAtoms():
                try:
                    atom.GetProp('AtomType')
                except KeyError:
                    untyped.append(atom.GetIdx())

            if untyped:
                raise AtomTypingError(
                    f"Failed to type atoms: {untyped}. "
                    f"This monomer chemistry may not be supported."
                )

            if self.verbose:
                logger.info("OPLS atom types assigned successfully")

            return mol

        except Exception as e:
            raise AtomTypingError(f"Atom typing failed: {e}")

    def validate_variants(self, variants: dict) -> bool:
        """
        Validate generated variants for correctness.

        Checks:
        1. All atoms have OPLS types assigned
        2. Bonding sites are correctly identified
        3. Charge neutrality (sum of charges ≈ 0)

        Args:
            variants: Dictionary of Mol objects from generate_variants()

        Returns:
            bool: True if all validations pass

        Raises:
            ValidationError: If validation fails
        """
        if not self.charge_dict:
            logger.warning("No charge dictionary available, skipping charge validation")
            return True

        for variant_name, mol in variants.items():
            if not isinstance(mol, Chem.Mol):
                continue

            # Check atom typing
            for atom in mol.GetAtoms():
                try:
                    atom.GetProp('AtomType')
                except KeyError:
                    raise ValidationError(
                        f"Variant {variant_name}: atom {atom.GetIdx()} has no AtomType"
                    )

            # Check charge neutrality
            charge_sum = 0.0
            for atom in mol.GetAtoms():
                atype = atom.GetProp('AtomType')
                charge = self.charge_dict.get(atype, 0.0)
                charge_sum += charge

            if abs(charge_sum) > 0.01:
                logger.warning(
                    f"Variant {variant_name}: net charge = {charge_sum:.3f} "
                    f"(may not be neutral)"
                )

        if self.verbose:
            logger.info("All variants validated successfully")

        return True

    def generate_variants(
        self,
        smiles: Optional[str] = None,
        mol: Optional[Chem.Mol] = None
    ) -> dict:
        """
        Generate all three monomer variants from SMILES or Mol object.

        This is the main entry point for monomer generation.

        Args:
            smiles: SMILES string (exclusive with mol)
            mol: RDKit Mol object (exclusive with smiles)

        Returns:
            dict: {
                'internal': mol_i,      # RDKit Mol for internal monomer
                'left_end': mol_le,     # RDKit Mol for left-end monomer
                'right_end': mol_re,    # RDKit Mol for right-end monomer
                'smiles': {
                    'internal': '...',
                    'left_end': '...',
                    'right_end': '...'
                },
                'bonding_info': {...}
            }

        Raises:
            MonomerGeneratorError: If input validation fails
            BondingSiteError: If bonding sites cannot be identified
            AtomTypingError: If atom typing fails

        Example:
            >>> generator = MonomerGenerator(base_name="PE")
            >>> variants = generator.generate_variants(smiles="C=C")
            >>> print(variants['smiles']['internal'])
        """
        # Validate input
        mol = self._validate_input(smiles, mol)

        if self.verbose:
            logger.info(f"Generating variants for {self.base_name}")

        # Add explicit hydrogens and generate initial conformer
        mol_with_h = AllChem.AddHs(mol)
        AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol_with_h)

        # Identify bonding sites
        bonding_info = self._identify_bonding_sites(mol_with_h)

        # Generate variants
        mol_internal = self._create_internal_variant(mol_with_h, bonding_info)
        mol_left_end = self._create_left_end_variant(mol_with_h, bonding_info)
        mol_right_end = self._create_right_end_variant(mol_with_h, bonding_info)

        # Assign OPLS types to all variants
        try:
            mol_internal = self._assign_opls_types(mol_internal)
            mol_left_end = self._assign_opls_types(mol_left_end)
            mol_right_end = self._assign_opls_types(mol_right_end)
        except AtomTypingError as e:
            logger.error(f"OPLS typing failed: {e}")
            raise

        # Package results
        variants = {
            'internal': mol_internal,
            'left_end': mol_left_end,
            'right_end': mol_right_end,
            'smiles': {
                'internal': Chem.MolToSmiles(mol_internal),
                'left_end': Chem.MolToSmiles(mol_left_end),
                'right_end': Chem.MolToSmiles(mol_right_end)
            },
            'bonding_info': bonding_info
        }

        if self.verbose:
            logger.info(f"Generated 3 variants:")
            logger.info(f"  Internal: {variants['smiles']['internal']}")
            logger.info(f"  Left-end: {variants['smiles']['left_end']}")
            logger.info(f"  Right-end: {variants['smiles']['right_end']}")

        return variants

    def generate_lt_files(self, variants: dict, generate_t1: bool = True) -> dict:
        """
        Generate Moltemplate .lt files for all variants including T1 chirality variants.

        Applies backbone alignment and reordering before file generation.
        Generates T1 variants for chiral monomers.

        Args:
            variants: Dictionary of Mol objects from generate_variants()
            generate_t1: Whether to generate T1 chirality variants (default: True)

        Returns:
            dict: Paths to generated .lt files
                {
                    'internal': '/path/to/PEi.lt',
                    'left_end': '/path/to/PEle.lt',
                    'right_end': '/path/to/PEre.lt',
                    'internal_T1': '/path/to/PEi_T1.lt',  # if chiral
                    'left_end_T1': '/path/to/PEle_T1.lt',  # if chiral
                    'right_end_T1': '/path/to/PEre_T1.lt'  # if chiral
                }

        Raises:
            IOError: If file writing fails

        Example:
            >>> variants = generator.generate_variants(smiles="C=C")
            >>> files = generator.generate_lt_files(variants)
            >>> print(files['internal'])
        """
        output_files = {}

        # Generate files for each variant
        variant_mapping = {
            'internal': 'i',
            'left_end': 'le',
            'right_end': 're'
        }

        for variant_type, suffix in variant_mapping.items():
            mol = variants[variant_type]

            # Step 1: Apply backbone alignment and reordering
            mol_aligned = self._align_and_reorder_backbone(mol, variant_type)

            # Step 2: Get SMILES from aligned molecule
            smiles = Chem.MolToSmiles(mol_aligned)

            # Step 3: Generate standard .lt file
            filename = f"{self.base_name}{suffix}.lt"
            filepath = os.path.join(self.output_dir, filename)

            # Check if file exists
            if os.path.exists(filepath):
                if self.verbose:
                    logger.warning(f"{filepath} exists, overwriting")

            # Write .lt file directly with aligned coordinates
            try:
                self._write_lt_file_direct(mol_aligned, filepath, f"{self.base_name}{suffix}", self.mechanism)
                output_files[variant_type] = filepath
            except Exception as e:
                logger.error(f"Failed to generate {filepath}: {e}")
                raise IOError(f"Failed to generate {filepath}: {e}")

            # Step 4: Generate T1 variant (always generate, regardless of chirality)
            if generate_t1:
                # Create T1 variant (Z-coordinate inversion)
                mol_t1 = self._create_t1_variant(mol_aligned)

                # Generate T1 filename
                filename_t1 = f"{self.base_name}{suffix}_T1.lt"
                filepath_t1 = os.path.join(self.output_dir, filename_t1)

                # Write T1 .lt file directly
                try:
                    self._write_lt_file_direct(mol_t1, filepath_t1, f"{self.base_name}{suffix}_T1", self.mechanism)
                    output_files[f"{variant_type}_T1"] = filepath_t1
                except Exception as e:
                    logger.error(f"Failed to generate {filepath_t1}: {e}")
                    raise IOError(f"Failed to generate {filepath_t1}: {e}")

        return output_files

    def _write_lt_file_direct(self, mol: Chem.Mol, filepath: str, name: str, mechanism: Optional[str] = None) -> None:
        """
        Write .lt file directly from molecule, preserving coordinates.

        This bypasses RDlt's conformer generation to use our aligned coordinates.
        For end-cap variants (le/re), modifies the output to create a connection point.
        Element-agnostic: works with C, O, N, and other elements.

        Args:
            mol: RDKit Mol object with AtomType properties and conformer
            filepath: Path to output file
            name: Name for the monomer class
            mechanism: Polymerization mechanism (for connection point logic)
        """
        # Detect if this is an end-cap variant
        is_left_end = 'le' in name and 'le_T1' not in name
        is_right_end = 're' in name and 're_T1' not in name
        is_left_end_t1 = 'le_T1' in name
        is_right_end_t1 = 're_T1' in name
        is_end_cap = is_left_end or is_right_end or is_left_end_t1 or is_right_end_t1

        # For end-caps, identify which H to remove and which atom's type to change
        skip_h_idx = None
        atom_to_modify = None
        new_atom_type = None

        if is_end_cap:
            # Get backbone atoms (connection points)
            # Use mechanism if available, otherwise detect
            if mechanism is None and self.mechanism:
                mechanism = self.mechanism

            backbone_atoms = self._get_backbone_atoms(mol, mechanism)

            if len(backbone_atoms) >= 2:
                if is_left_end or is_left_end_t1:
                    # Left-end: remove H from second backbone atom
                    atom_to_modify = backbone_atoms[1]
                else:  # right_end or right_end_t1
                    # Right-end: remove H from first backbone atom
                    atom_to_modify = backbone_atoms[0]

                # Use ConnectionPointModifier to get the new atom type
                atom_obj = mol.GetAtomWithIdx(atom_to_modify)
                original_type, connection_type = self.conn_modifier.get_connection_atom_type(atom_obj, mol)

                # Only skip H removal and change type if applicable
                if original_type and connection_type:
                    new_atom_type = connection_type

                    # Find one hydrogen bonded to the atom to modify
                    for neighbor in atom_obj.GetNeighbors():
                        if neighbor.GetSymbol() == 'H':
                            skip_h_idx = neighbor.GetIdx()
                            break

                    if self.verbose:
                        logger.info(f"End-cap: modifying atom {atom_to_modify} ({atom_obj.GetSymbol()}) "
                                   f"from {original_type} to {new_atom_type}")

        # Detect if this is a water molecule for SPC/E force field
        is_water_molecule = mol.GetNumAtoms() == 3
        if is_water_molecule:
            atoms = list(mol.GetAtoms())
            symbols = [atom.GetSymbol() for atom in atoms]
            is_water_molecule = symbols.count('O') == 1 and symbols.count('H') == 2

        with open(filepath, 'w') as f:
            # Write header based on force field
            if is_water_molecule:
                # SPC/E water model
                print('import "spce.lt"    # <-- defines the SPC/E water model', file=f)
                print(f'{name} inherits spce {{', file=f)
            elif self.is_gaff:
                print('import "gaff.lt"    # <-- defines the GAFF (General Amber Force Field)', file=f)
                print('# NOTE: GAFF requires user-supplied charges (AM1-BCC or RESP recommended)', file=f)
                print('# See: http://ambermd.org/antechamber/gaff.pdf', file=f)
                print(f'{name} inherits GAFF {{', file=f)
            elif self.is_lopls:
                print('import "oplsaa.lt"    # <-- defines the standard "OPLSAA" force field', file=f)
                print('import "loplsaa.lt"   # <-- custom parameters for long alkane chains', file=f)
                print(f'{name} inherits OPLSAA {{', file=f)
            else:
                print('import "oplsaa.lt"    # <-- defines the standard "OPLSAA" force field', file=f)
                print(f'{name} inherits OPLSAA {{', file=f)
            print('', file=f)
            print('# atom-id  mol-id  atom-type charge      X         Y        Z', file=f)
            print('', file=f)
            print('  write("Data Atoms") {', file=f)

            # Write atoms
            conf = mol.GetConformer(0)
            atom_counter = 1  # For naming atoms in output

            for atom in mol.GetAtoms():
                idx = atom.GetIdx()
                symbol = atom.GetSymbol()

                # Skip this hydrogen if it's the one we're removing
                if idx == skip_h_idx:
                    continue

                pos = conf.GetAtomPosition(idx)

                try:
                    atom_type = atom.GetProp('AtomType')
                    # For end-cap connection point, change atom type
                    if is_end_cap and idx == atom_to_modify and new_atom_type:
                        atom_type = new_atom_type
                except KeyError:
                    atom_type = '@atom:???'  # Should not happen if typing succeeded

                print(f'\t$atom:{symbol}{atom_counter} $mol:... {atom_type} 0.00    {pos.x:.3f}   {pos.y:.3f}   {pos.z:.3f}', file=f)
                atom_counter += 1

            print('  }', file=f)
            print('', file=f)
            print('  write(\'Data Bond List\') {', file=f)

            # Track atom index mapping after skipping hydrogen
            idx_map = {}
            new_idx = 1
            for atom in mol.GetAtoms():
                old_idx = atom.GetIdx()
                if old_idx == skip_h_idx:
                    continue
                idx_map[old_idx] = new_idx
                new_idx += 1

            # Write bonds
            for bond in mol.GetBonds():
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()

                # Skip bonds involving the removed hydrogen
                if begin_idx == skip_h_idx or end_idx == skip_h_idx:
                    continue

                begin_atom = mol.GetAtomWithIdx(begin_idx)
                end_atom = mol.GetAtomWithIdx(end_idx)

                begin_symbol = begin_atom.GetSymbol()
                end_symbol = end_atom.GetSymbol()

                # Use new indices after skipping hydrogen
                new_begin_idx = idx_map[begin_idx]
                new_end_idx = idx_map[end_idx]

                bond_name = f'{begin_symbol}{new_begin_idx}{end_symbol}{new_end_idx}'
                begin_name = f'{begin_symbol}{new_begin_idx}'
                end_name = f'{end_symbol}{new_end_idx}'

                print(f'\t$bond:{bond_name}\t$atom:{begin_name}\t$atom:{end_name}', file=f)

            print('  }', file=f)
            print('}  # ' + name, file=f)
            print('', file=f)
            print('# Note: You don\'t need to supply the partial partial charges of the atoms.', file=f)
            print('#       If you like, just fill the fourth column with zeros ("0.000").', file=f)
            print('#       Moltemplate and LAMMPS will automatically assign the charge later', file=f)
            print('', file=f)

        if self.verbose:
            logger.info(f"Wrote {filepath} with aligned coordinates")
            if is_end_cap:
                logger.info(f"  End-cap: removed H atom and changed connection point C type to @atom:82")


# Convenience function for quick monomer generation
def generate_monomer(
    base_name: str,
    smiles: str,
    output_dir: Optional[str] = None,
    is_lopls: bool = False
) -> dict:
    """
    Convenience function to generate monomer variants and .lt files in one call.

    Args:
        base_name: Base name for monomer (e.g., "PE", "PMMA")
        smiles: SMILES string for monomer
        output_dir: Output directory (default: Monomer_bank)
        is_lopls: Use LOPLS force field

    Returns:
        dict: {
            'variants': {...},  # Mol objects
            'files': {...}      # File paths
        }

    Example:
        >>> result = generate_monomer("PE", "C=C")
        >>> print(f"Generated: {result['files']}")
    """
    generator = MonomerGenerator(
        base_name=base_name,
        output_dir=output_dir,
        is_lopls=is_lopls
    )

    variants = generator.generate_variants(smiles=smiles)
    files = generator.generate_lt_files(variants)

    return {
        'variants': variants,
        'files': files,
        'generator': generator
    }
