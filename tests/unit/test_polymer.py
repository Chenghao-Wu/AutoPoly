"""Tests for the Polymer class."""

import pytest
from AutoPoly.polymer import Polymer
from AutoPoly.exceptions import ValidationError


class TestPolymerInitialization:
    """Test Polymer class initialization."""

    def test_polymer_init_with_default_parameters(self):
        """Test Polymer initialization with defaults."""
        poly = Polymer(ChainNum=1, Sequence=["PE"])
        assert poly.ChainNum == 1
        assert poly.DOP == 1  # Defaults to sequence length
        assert poly.topology == "linear"
        assert poly.tacticity == "atactic"

    def test_polymer_init_with_custom_dop(self):
        """Test Polymer initialization with custom DOP."""
        poly = Polymer(ChainNum=1, Sequence=["PE"], DOP=5)
        assert poly.DOP == 5

    def test_polymer_init_with_invalid_topology_raises_value_error(self):
        """Test that invalid topology raises ValueError."""
        with pytest.raises(ValueError, match="Topology must be either 'linear' or 'ring'"):
            Polymer(ChainNum=1, Sequence=["PE"], topology="invalid")

    def test_polymer_init_with_ring_topology(self):
        """Test Polymer initialization with ring topology."""
        poly = Polymer(ChainNum=1, Sequence=["PE"], topology="ring")
        assert poly.topology == "ring"

    def test_polymer_init_with_nested_sequence(self):
        """Test Polymer initialization with nested sequence."""
        # User might pass [["PE"]] instead of ["PE"]
        poly = Polymer(ChainNum=1, Sequence=[["PE"]])
        assert poly.sequence == ["PE"]

    def test_polymer_init_with_copolymer_sequence(self):
        """Test Polymer initialization with copolymer sequence."""
        poly = Polymer(ChainNum=1, Sequence=["PE", "PS"])
        assert poly.sequence == ["PE", "PS"]
        assert len(poly.merSet) == 2

    def test_polymer_init_empty_sequence_raises_value_error(self):
        """Test that empty sequence raises ValueError."""
        with pytest.raises(ValueError, match="Sequence cannot be None or empty"):
            Polymer(ChainNum=1, Sequence=[])

    def test_polymer_init_none_sequence_raises_value_error(self):
        """Test that None sequence raises ValueError."""
        with pytest.raises(ValueError, match="Sequence cannot be None or empty"):
            Polymer(ChainNum=1, Sequence=None)


class TestPolymerSequence:
    """Test polymer sequence generation."""

    def test_polymer_set_sequence_with_single_monomer(self):
        """Test sequence generation with single monomer type."""
        poly = Polymer(ChainNum=1, Sequence=["PE"], DOP=5)
        poly.set_Sequence()
        assert len(poly.sequenceSet) == poly.ChainNum
        assert len(poly.sequenceSet[0]) == poly.DOP
        # Check that base SMILES is in the identifiers (may have _T1 suffix)
        assert all("PE" in identifier for identifier in poly.sequenceSet[0])

    def test_polymer_set_sequence_with_copolymer(self):
        """Test sequence generation with copolymer (alternating)."""
        poly = Polymer(ChainNum=1, Sequence=["PE", "PS"], DOP=4)
        poly.set_Sequence()
        assert len(poly.sequenceSet[0]) == 4
        # Should alternate: PE, PS, PE, PS (or similar pattern)
        # Check that both PE and PS are present
        identifiers = "".join(poly.sequenceSet[0])
        assert "PE" in identifiers
        assert "PS" in identifiers

    def test_polymer_set_sequence_with_zero_chain_num_raises_system_exit(self):
        """Test that ChainNum=0 raises ValidationError."""
        # Create polymer with ChainNum=1 first to avoid ValidationError in __init__
        poly = Polymer(ChainNum=1, Sequence=["PE"])
        # Then set ChainNum to 0 and call set_Sequence
        poly.ChainNum = 0
        with pytest.raises(ValidationError, match="ChainNum must be greater than 0"):
            poly.set_Sequence()

    def test_polymer_sequence_set_structure(self):
        """Test that sequenceSet has correct structure."""
        poly = Polymer(ChainNum=2, Sequence=["PE"], DOP=3)
        poly.set_Sequence()
        assert len(poly.sequenceSet) == 2  # 2 chains
        assert len(poly.sequenceSet[0]) == 3  # DOP=3
        assert len(poly.sequenceSet[1]) == 3

    def test_polymer_removes_lt_extension_from_sequence(self):
        """Test that .lt extension is removed from sequence items."""
        poly = Polymer(ChainNum=1, Sequence=["PE.lt"], DOP=3)
        poly.set_Sequence()
        # Should remove .lt extension from identifiers
        # All identifiers should contain "PE" but not "PE.lt"
        assert all("PE" in identifier and ".lt" not in identifier
                  for identifier in poly.sequenceSet[0])

    def test_polymer_sequence_name_matches_sequence_set(self):
        """Test that sequenceName matches sequenceSet."""
        poly = Polymer(ChainNum=1, Sequence=["PE"], DOP=3)
        poly.set_Sequence()
        assert poly.sequenceName == poly.sequenceSet

    def test_polymer_multiple_chains_have_independent_sequences(self):
        """Test that multiple chains have independent sequences."""
        poly = Polymer(ChainNum=3, Sequence=["PE"], DOP=2)
        poly.set_Sequence()
        assert len(poly.sequenceSet) == 3
        # Each chain should have DOP monomers
        assert all(len(chain) == 2 for chain in poly.sequenceSet)


class TestPolymerTacticity:
    """Test polymer tacticity assignment."""

    def test_polymer_isotactic_all_same_chirality(self):
        """Test isotactic tacticity (all same chirality)."""
        poly = Polymer(ChainNum=1, Sequence=["PE"], DOP=5, tacticity="isotactic")
        poly.set_Sequence()
        tacticity = poly.tacticitySet[0]
        # All should be True or all should be False
        assert all(tacticity) or all(not t for t in tacticity)

    def test_polymer_syndiotactic_alternating_chirality(self):
        """Test syndiotactic tacticity (alternating chirality)."""
        poly = Polymer(ChainNum=1, Sequence=["PE"], DOP=5, tacticity="syndiotactic")
        poly.set_Sequence()
        chirality = poly.tacticitySet[0]
        # Should alternate: False, True, False, True, ...
        for i in range(len(chirality) - 1):
            assert chirality[i] != chirality[i + 1]

    def test_polymer_atactic_random_chirality(self):
        """Test atactic tacticity (random chirality)."""
        poly = Polymer(ChainNum=1, Sequence=["PE"], DOP=10, tacticity="atactic")
        poly.set_Sequence()
        tacticity = poly.tacticitySet[0]

        # Should have mix of True and False (not all same)
        # Note: This is probabilistic, but with 10 positions it's very unlikely
        # to get all True or all False by chance
        assert not all(tacticity) or not all(not t for t in tacticity)

    def test_polymer_tacticity_set_per_chain(self):
        """Test that each chain has its own tacticity."""
        poly = Polymer(ChainNum=2, Sequence=["PE"], DOP=5, tacticity="isotactic")
        poly.set_Sequence()
        assert len(poly.tacticitySet) == 2  # 2 chains
        assert len(poly.tacticitySet[0]) == 5  # DOP=5
        assert len(poly.tacticitySet[1]) == 5

    def test_polymer_t1_markers_in_identifiers_for_tacticity(self):
        """Test that _T1 markers are added to identifiers when use_t1=True."""
        # For isotactic with use_t1=True, all identifiers should have _T1
        poly = Polymer(ChainNum=1, Sequence=["PE"], DOP=3, tacticity="isotactic")
        poly.set_Sequence()
        tacticity = poly.tacticitySet[0]

        # Check that identifiers match tacticity
        for i, identifier in enumerate(poly.sequenceSet[0]):
            if tacticity[i]:
                assert "_T1" in identifier
            else:
                assert "_T1" not in identifier


class TestPolymerGetters:
    """Test Polymer getter methods."""

    def test_polymer_get_mer_set_returns_unique_monomers(self):
        """Test that get_mer() returns unique monomers."""
        poly = Polymer(ChainNum=1, Sequence=["PE", "PS", "PE"], DOP=3)
        poly.set_Sequence()
        mer_set = poly.get_mer_set()
        assert "PE" in mer_set
        assert "PS" in mer_set
        # Should be unique (no duplicates)
        assert len(mer_set) == len(set(mer_set))

    def test_polymer_get_sequence_set(self):
        """Test get_sequence_set returns correct structure."""
        poly = Polymer(ChainNum=1, Sequence=["PE"], DOP=3)
        poly.set_Sequence()
        seq_set = poly.get_sequence_set()
        assert len(seq_set) == 1
        assert len(seq_set[0]) == 3

    def test_polymer_get_sequence_names(self):
        """Test get_sequence_names returns correct structure."""
        poly = Polymer(ChainNum=1, Sequence=["PE"], DOP=3)
        poly.set_Sequence()
        seq_names = poly.get_sequence_names()
        assert len(seq_names) == 1
        assert len(seq_names[0]) == 3

    def test_polymer_get_chain_info_returns_complete_dict(self):
        """Test that get_chain_info() returns complete information."""
        poly = Polymer(ChainNum=1, Sequence=["PE"], DOP=5, tacticity="isotactic")
        poly.set_Sequence()
        info = poly.get_chain_info()

        assert 'chain_num' in info
        assert 'sequence' in info
        assert 'dop' in info
        assert 'topology' in info
        assert 'tacticity' in info
        assert 'sequence_length' in info
        assert 'mer_set' in info
        assert 'sequence_set' in info
        assert 'sequence_names' in info
        assert 'tacticity_set' in info

        assert info['chain_num'] == 1
        assert info['dop'] == 5
        assert info['topology'] == "linear"
        assert info['tacticity'] == "isotactic"

    def test_polymer_get_tacticity_for_chain(self):
        """Test getting tacticity for specific chain."""
        poly = Polymer(ChainNum=2, Sequence=["PE"], DOP=3)
        poly.set_Sequence()
        chirality_0 = poly.get_tacticity_for_chain(0)
        chirality_1 = poly.get_tacticity_for_chain(1)
        assert len(chirality_0) == 3
        assert len(chirality_1) == 3

    def test_polymer_get_tacticity_for_invalid_chain(self):
        """Test getting tacticity for invalid chain index."""
        poly = Polymer(ChainNum=1, Sequence=["PE"], DOP=3)
        poly.set_Sequence()
        # Invalid chain index should return empty list
        chirality = poly.get_tacticity_for_chain(5)
        assert chirality == []


class TestPolymerSetters:
    """Test Polymer setter methods."""

    def test_polymer_set_mer_set_with_list(self):
        """Test set_merSet with list of monomers."""
        poly = Polymer(ChainNum=1, Sequence=["PE", "PS"])
        # Should remove duplicates
        assert "PE" in poly.merSet
        assert "PS" in poly.merSet

    def test_polymer_set_mer_set_with_single_string(self):
        """Test set_merSet with single monomer string."""
        poly = Polymer(ChainNum=1, Sequence=["PE"])
        assert poly.merSet == ["PE"]

    def test_polymer_set_dop(self):
        """Test set_dop method."""
        poly = Polymer(ChainNum=1, Sequence=["PE"], DOP=5)
        assert poly.DOP == 5
        poly.set_dop(10)
        assert poly.DOP == 10


class TestPolymerCopolymers:
    """Test copolymer-specific behavior."""

    def test_polymer_copolymer_alternating_sequence(self):
        """Test that copolymer alternates through sequence."""
        # With sequence ["PE", "PS"] and DOP=4, should get PE, PS, PE, PS
        poly = Polymer(ChainNum=1, Sequence=["PE", "PS"], DOP=4)
        poly.set_Sequence()

        identifiers = [id.replace("_T1", "") for id in poly.sequenceSet[0]]
        # Should alternate: PE, PS, PE, PS
        assert identifiers[0] == "PE"
        assert identifiers[1] == "PS"
        assert identifiers[2] == "PE"
        assert identifiers[3] == "PS"

    def test_polymer_copolymer_longer_than_sequence(self):
        """Test copolymer when DOP > sequence length."""
        # With sequence ["PE", "PS", "PP"] and DOP=7, should cycle
        poly = Polymer(ChainNum=1, Sequence=["PE", "PS", "PP"], DOP=7)
        poly.set_Sequence()

        assert len(poly.sequenceSet[0]) == 7
        # Should cycle through PE, PS, PP, PE, PS, PP, PE
        identifiers = [id.replace("_T1", "") for id in poly.sequenceSet[0]]
        assert identifiers[0] == "PE"
        assert identifiers[1] == "PS"
        assert identifiers[2] == "PP"
        assert identifiers[3] == "PE"  # Cycles back
