"""Contract tests for external dependencies."""

import pytest
from rdkit import Chem


@pytest.mark.integration
class TestExternalContracts:
    """Test external dependency contracts."""

    def test_rdkit_mol_from_smiles_returns_valid_mol(self):
        """Test that RDKit MolFromSmiles returns valid Mol object."""
        mol = Chem.MolFromSmiles("O")
        assert mol is not None
        # RDKit doesn't add implicit H by default, so we only count explicit atoms
        # Water SMILES "O" has 1 atom (oxygen), hydrogens are implicit
        assert mol.GetNumAtoms() == 1

    def test_rdkit_add_hs_preserves_heavy_atoms(self):
        """Test that AddHs preserves heavy atom count."""
        mol = Chem.MolFromSmiles("C")  # Methane
        mol_with_h = Chem.AddHs(mol)
        assert mol_with_h.GetNumAtoms() == 5  # 1 C + 4 H

    def test_rdkit_has_property_returns_false_for_missing(self):
        """Test that HasProp returns False for missing properties."""
        mol = Chem.MolFromSmiles("O")
        assert not mol.HasProp("missing_property")

    def test_rdkit_smiles_with_wildcards(self):
        """Test RDKit handles wildcards in SMILES."""
        mol = Chem.MolFromSmiles("[*]C[*]")
        assert mol is not None
        # Should have 3 atoms: wildcard, C, wildcard
        assert mol.GetNumAtoms() == 3

    def test_rdkit_aromatic_smiles(self):
        """Test RDKit handles aromatic SMILES correctly."""
        mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene
        assert mol is not None
        # RDKit doesn't add implicit H by default
        # Benzene SMILES only has 6 explicit C atoms
        assert mol.GetNumAtoms() == 6
        # With AddHs, we get 12 atoms (6 C + 6 H)
        mol_with_h = Chem.AddHs(mol)
        assert mol_with_h.GetNumAtoms() == 12

    def test_rdkit_bond_iteration(self):
        """Test that RDKit bond iteration works."""
        mol = Chem.MolFromSmiles("CCO")  # Ethanol
        assert mol is not None
        num_bonds = mol.GetNumBonds()
        assert num_bonds == 2  # C-C and C-O bonds

    def test_rdkit_atom_properties(self):
        """Test that RDKit atom properties are accessible."""
        mol = Chem.MolFromSmiles("CCO")
        assert mol is not None

        atoms = mol.GetAtoms()
        assert len(list(atoms)) == 3  # 2 C + 1 O
