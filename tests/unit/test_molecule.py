"""Tests for the Molecule class."""

import pytest
from AutoPoly.molecule import Molecule


class TestMoleculeInitialization:
    """Test Molecule class initialization."""

    def test_molecule_init_with_valid_parameters(self):
        """Test Molecule initialization with water SMILES."""
        mol = Molecule(Count=10, Smiles="O")
        assert mol.Count == 10
        assert mol.Smiles == "O"
        assert mol.molecule_name == "molecule_O"  # Auto-generated from SMILES

    def test_molecule_init_with_custom_name(self):
        """Test Molecule initialization with custom name."""
        mol = Molecule(Count=10, Smiles="O", Name="my_water")
        assert mol.molecule_name == "my_water"

    def test_molecule_init_with_benzene(self):
        """Test Molecule initialization with benzene SMILES."""
        mol = Molecule(Count=5, Smiles="c1ccccc1")
        assert mol.Count == 5
        assert mol.Smiles == "c1ccccc1"
        # Auto-generated name should clean special characters
        assert "molecule" in mol.molecule_name

    def test_molecule_init_with_ethanol(self):
        """Test Molecule initialization with ethanol SMILES."""
        mol = Molecule(Count=20, Smiles="CCO")
        assert mol.Count == 20
        assert mol.Smiles == "CCO"
        assert mol.molecule_name == "molecule_CCO"


class TestMoleculeValidation:
    """Test Molecule parameter validation."""

    def test_molecule_count_none_raises_value_error(self):
        """Test that Count=None raises ValueError."""
        with pytest.raises(ValueError, match="Count must be a positive integer"):
            Molecule(Count=None, Smiles="O")

    def test_molecule_count_zero_raises_value_error(self):
        """Test that Count=0 raises ValueError."""
        with pytest.raises(ValueError, match="Count must be a positive integer"):
            Molecule(Count=0, Smiles="O")

    def test_molecule_count_negative_raises_value_error(self):
        """Test that negative Count raises ValueError."""
        with pytest.raises(ValueError, match="Count must be a positive integer"):
            Molecule(Count=-5, Smiles="O")

    def test_molecule_smiles_none_raises_value_error(self):
        """Test that Smiles=None raises ValueError."""
        with pytest.raises(ValueError, match="Smiles cannot be None or empty"):
            Molecule(Count=10, Smiles=None)

    def test_molecule_smiles_empty_raises_value_error(self):
        """Test that empty Smiles raises ValueError."""
        with pytest.raises(ValueError, match="Smiles cannot be None or empty"):
            Molecule(Count=10, Smiles="")

    def test_molecule_smiles_whitespace_only_raises_value_error(self):
        """Test that whitespace-only Smiles raises ValueError."""
        with pytest.raises(ValueError, match="Smiles cannot be None or empty"):
            Molecule(Count=10, Smiles="   ")

    def test_molecule_smiles_with_wildcard_raises_value_error(self):
        """Test that SMILES with wildcard [*] raises ValueError."""
        # Molecules should NOT have connection points (wildcards)
        with pytest.raises(ValueError, match="should not contain wildcards"):
            Molecule(Count=10, Smiles="[*]C[*]")

    def test_molecule_smiles_with_single_wildcard_raises_value_error(self):
        """Test that SMILES with single wildcard * raises ValueError."""
        with pytest.raises(ValueError, match="should not contain wildcards"):
            Molecule(Count=10, Smiles="C*")


class TestMoleculeStructure:
    """Test molecule structure setup."""

    def test_molecule_sequence_set_structure(self):
        """Test that sequenceSet has correct structure."""
        mol = Molecule(Count=10, Smiles="O")
        assert len(mol.sequenceSet) == mol.Count
        assert all(len(item) == 1 for item in mol.sequenceSet)  # Each is a single-element list
        assert all(item == ["molecule_O"] for item in mol.sequenceSet)

    def test_molecule_sequence_name_structure(self):
        """Test that sequenceName has correct structure."""
        mol = Molecule(Count=5, Smiles="CCO")
        assert len(mol.sequenceName) == mol.Count
        assert all(len(item) == 1 for item in mol.sequenceName)
        assert all(item == ["molecule_CCO"] for item in mol.sequenceName)

    def test_molecule_mer_set_contains_single_molecule(self):
        """Test that mer contains the molecule name."""
        mol = Molecule(Count=10, Smiles="O")
        assert mol.merSet == ["molecule_O"]
        assert len(mol.merSet) == 1

    def test_molecule_dop_is_always_one(self):
        """Test that DOP is always 1 for molecules."""
        mol = Molecule(Count=10, Smiles="O")
        assert mol.DOP == 1

    def test_molecule_is_molecule_flag(self):
        """Test that _is_molecule flag is set to True."""
        mol = Molecule(Count=10, Smiles="O")
        assert mol._is_molecule is True


class TestMoleculeGetters:
    """Test Molecule getter methods."""

    def test_get_count(self):
        """Test get_count returns correct value."""
        mol = Molecule(Count=42, Smiles="O")
        assert mol.get_count() == 42

    def test_get_smiles(self):
        """Test get_smiles returns correct value."""
        mol = Molecule(Count=10, Smiles="CCO")
        assert mol.get_smiles() == "CCO"

    def test_get_name(self):
        """Test get_name returns correct value."""
        mol = Molecule(Count=10, Smiles="O", Name="water")
        assert mol.get_name() == "water"

    def test_get_sequence_set(self):
        """Test get_sequence_set returns correct structure."""
        mol = Molecule(Count=3, Smiles="O")
        seq_set = mol.get_sequence_set()
        assert len(seq_set) == 3
        assert all(isinstance(item, list) for item in seq_set)

    def test_get_sequence_names(self):
        """Test get_sequence_names returns correct structure."""
        mol = Molecule(Count=3, Smiles="O", Name="water")
        seq_names = mol.get_sequence_names()
        assert len(seq_names) == 3
        assert all(item == ["water"] for item in seq_names)

    def test_get_mer_set(self):
        """Test get_mer_set returns correct value."""
        mol = Molecule(Count=10, Smiles="c1ccccc1", Name="benzene")
        mer_set = mol.get_mer_set()
        assert mer_set == ["benzene"]

    def test_get_molecule_info(self):
        """Test get_molecule_info returns complete dictionary."""
        mol = Molecule(Count=5, Smiles="O", Name="water")
        info = mol.get_molecule_info()

        assert 'count' in info
        assert 'smiles' in info
        assert 'name' in info
        assert 'dop' in info
        assert 'mer_set' in info
        assert 'sequence_set' in info
        assert 'sequence_names' in info

        assert info['count'] == 5
        assert info['smiles'] == "O"
        assert info['name'] == "water"
        assert info['dop'] == 1
        assert info['mer_set'] == ["water"]


class TestMoleculeRepr:
    """Test Molecule string representation."""

    def test_repr(self):
        """Test __repr__ returns correct string."""
        mol = Molecule(Count=10, Smiles="O", Name="water")
        repr_str = repr(mol)
        assert "Molecule" in repr_str
        assert "Count=10" in repr_str
        assert "Smiles='O'" in repr_str
        assert "Name='water'" in repr_str
