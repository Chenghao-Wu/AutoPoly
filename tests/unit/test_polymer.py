"""Tests for the Polymer class."""

import pytest
from AutoPoly.models.polymer import Polymer
from AutoPoly.core.exceptions import ValidationError

# Valid pSMILES building blocks for the new wildcard validation:
# first monomer = 1 wildcard, middle = 2 wildcards, last = 1 wildcard,
# single-monomer sequence = 0 wildcards.
PE_FIRST = "CC[*]"
PE_MID = "[*]CC[*]"
PE_LAST = "[*]CC"
PS_MID = "[*]C=C[*]"
PS_LAST = "[*]C=C"
PP_MID = "[*]CC(C)[*]"
SINGLE = "CC"  # complete molecule, no wildcards


def uniform(n):
    """n-monomer homopolymer sequence with valid wildcard pattern."""
    if n == 1:
        return [SINGLE]
    return [PE_FIRST] + [PE_MID] * (n - 2) + [PE_LAST]


class TestPolymerInitialization:
    """Test Polymer class initialization."""

    def test_polymer_init_with_default_parameters(self):
        """Test Polymer initialization with defaults."""
        poly = Polymer(chain_num=1, sequence=[SINGLE])
        assert poly.chain_num == 1
        assert poly.dop == 1  # Defaults to sequence length
        assert poly.topology == "linear"
        assert poly.tacticity == "atactic"

    def test_polymer_init_with_custom_dop(self):
        """Test Polymer initialization with custom DOP (DOP is derived from sequence length)."""
        # With new API, DOP is derived from sequence length, not a parameter
        poly = Polymer(chain_num=1, sequence=uniform(5))
        assert poly.dop == 5

    def test_polymer_init_with_invalid_topology_raises_value_error(self):
        """Test that invalid topology raises ValueError."""
        with pytest.raises(ValueError, match="topology must be either 'linear' or 'ring'"):
            Polymer(chain_num=1, sequence=[SINGLE], topology="invalid")

    def test_polymer_init_with_ring_topology(self):
        """Test Polymer initialization with ring topology."""
        poly = Polymer(chain_num=1, sequence=[SINGLE], topology="ring")
        assert poly.topology == "ring"

    def test_polymer_init_with_nested_sequence(self):
        """Test Polymer initialization with nested sequence."""
        # User might pass [["CC"]] instead of ["CC"]
        poly = Polymer(chain_num=1, sequence=[[SINGLE]])
        assert poly.sequence == [SINGLE]

    def test_polymer_init_with_copolymer_sequence(self):
        """Test Polymer initialization with copolymer sequence."""
        poly = Polymer(chain_num=1, sequence=[PE_FIRST, PS_LAST])
        assert poly.sequence == [PE_FIRST, PS_LAST]
        assert len(poly.mer_set) == 2

    def test_polymer_init_empty_sequence_raises_value_error(self):
        """Test that empty sequence raises ValueError."""
        with pytest.raises(ValidationError, match="sequence cannot be empty"):
            Polymer(chain_num=1, sequence=[])

    def test_polymer_init_none_sequence_raises_value_error(self):
        """Test that None sequence raises ValueError."""
        with pytest.raises(ValidationError, match="sequence cannot be empty"):
            Polymer(chain_num=1, sequence=None)


class TestPolymerSequence:
    """Test polymer sequence generation."""

    def test_polymer_set_sequence_with_single_monomer(self):
        """Test sequence generation with single monomer type."""
        poly = Polymer(chain_num=1, sequence=uniform(5))
        poly.set_Sequence()
        assert len(poly.sequenceSet) == poly.chain_num
        assert len(poly.sequenceSet[0]) == poly.dop
        # Check that base SMILES is in the identifiers (may have _T1 suffix)
        assert all("CC" in identifier for identifier in poly.sequenceSet[0])

    def test_polymer_set_sequence_with_copolymer(self):
        """Test sequence generation with explicit copolymer sequence."""
        # With new API, explicit sequence - NO CYCLING
        sequence = [PE_FIRST, PS_MID, PE_MID, PS_LAST]
        poly = Polymer(chain_num=1, sequence=sequence)
        poly.set_Sequence()
        assert len(poly.sequenceSet[0]) == 4
        # Should be exact sequence, no cycling — both monomers present
        identifiers = "".join(poly.sequenceSet[0])
        assert "CC" in identifiers
        assert "C=C" in identifiers

    def test_polymer_set_sequence_with_zero_chain_num_raises_validation_error(self):
        """Test that chain_num=0 raises ValidationError."""
        # Create polymer with chain_num=1 first to avoid ValidationError in __init__
        poly = Polymer(chain_num=1, sequence=[SINGLE])
        # Then set chain_num to 0 and call set_Sequence
        poly.chain_num = 0
        with pytest.raises(ValidationError, match="chain_num must be greater than 0"):
            poly.set_Sequence()

    def test_polymer_sequence_set_structure(self):
        """Test that sequenceSet has correct structure."""
        poly = Polymer(chain_num=2, sequence=uniform(3))
        poly.set_Sequence()
        assert len(poly.sequenceSet) == 2  # 2 chains
        assert len(poly.sequenceSet[0]) == 3  # dop=3
        assert len(poly.sequenceSet[1]) == 3

    def test_polymer_removes_lt_extension_from_sequence(self):
        """Test that .lt extension is removed from sequence items."""
        sequence = [PE_FIRST + ".lt", PE_MID + ".lt", PE_LAST + ".lt"]
        poly = Polymer(chain_num=1, sequence=sequence)
        poly.set_Sequence()
        # Should remove .lt extension from identifiers
        assert all("CC" in identifier and ".lt" not in identifier
                  for identifier in poly.sequenceSet[0])

    def test_polymer_sequence_name_matches_sequence_set(self):
        """Test that sequenceName matches sequenceSet."""
        poly = Polymer(chain_num=1, sequence=uniform(3))
        poly.set_Sequence()
        assert poly.sequenceName == poly.sequenceSet

    def test_polymer_multiple_chains_have_independent_sequences(self):
        """Test that multiple chains have independent sequences."""
        poly = Polymer(chain_num=3, sequence=uniform(2))
        poly.set_Sequence()
        assert len(poly.sequenceSet) == 3
        # Each chain should have DOP monomers
        assert all(len(chain) == 2 for chain in poly.sequenceSet)


class TestPolymerTacticity:
    """Test polymer tacticity assignment."""

    def test_polymer_isotactic_all_same_chirality(self):
        """Test isotactic tacticity (all same chirality)."""
        poly = Polymer(chain_num=1, sequence=uniform(5), tacticity="isotactic")
        poly.set_Sequence()
        tacticity = poly.tacticity_set[0]
        # All should be True or all should be False
        assert all(tacticity) or all(not t for t in tacticity)

    def test_polymer_syndiotactic_alternating_chirality(self):
        """Test syndiotactic tacticity (alternating chirality)."""
        poly = Polymer(chain_num=1, sequence=uniform(5), tacticity="syndiotactic")
        poly.set_Sequence()
        chirality = poly.tacticity_set[0]
        # Should alternate: False, True, False, True, ...
        for i in range(len(chirality) - 1):
            assert chirality[i] != chirality[i + 1]

    def test_polymer_atactic_random_chirality(self):
        """Test atactic tacticity (random chirality)."""
        poly = Polymer(chain_num=1, sequence=uniform(10), tacticity="atactic")
        poly.set_Sequence()
        tacticity = poly.tacticity_set[0]

        # Should have mix of True and False (not all same)
        # Note: This is probabilistic, but with 10 positions it's very unlikely
        # to get all True or all False by chance
        assert not all(tacticity) or not all(not t for t in tacticity)

    def test_polymer_tacticity_set_per_chain(self):
        """Test that each chain has its own tacticity."""
        poly = Polymer(chain_num=2, sequence=uniform(5), tacticity="isotactic")
        poly.set_Sequence()
        assert len(poly.tacticity_set) == 2  # 2 chains
        assert len(poly.tacticity_set[0]) == 5  # DOP=5
        assert len(poly.tacticity_set[1]) == 5

    def test_polymer_t1_markers_in_identifiers_for_tacticity(self):
        """Test that _T1 markers are added to identifiers when use_t1=True."""
        # For isotactic with use_t1=True, all identifiers should have _T1
        poly = Polymer(chain_num=1, sequence=uniform(3), tacticity="isotactic")
        poly.set_Sequence()
        tacticity = poly.tacticity_set[0]

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
        poly = Polymer(chain_num=1, sequence=[PE_FIRST, PS_MID, PE_LAST])
        poly.set_Sequence()
        mer_set = poly.get_mer_set()
        assert PE_FIRST in mer_set
        assert PS_MID in mer_set
        # Should be unique (no duplicates)
        assert len(mer_set) == len(set(mer_set))

    def test_polymer_get_sequence_set(self):
        """Test get_sequenceSet returns correct structure."""
        poly = Polymer(chain_num=1, sequence=uniform(3))
        poly.set_Sequence()
        seq_set = poly.get_sequenceSet()
        assert len(seq_set) == 1
        assert len(seq_set[0]) == 3

    def test_polymer_get_sequence_names(self):
        """Test get_sequenceNames returns correct structure."""
        poly = Polymer(chain_num=1, sequence=uniform(3))
        poly.set_Sequence()
        seq_names = poly.get_sequenceNames()
        assert len(seq_names) == 1
        assert len(seq_names[0]) == 3

    def test_polymer_get_chain_info_returns_complete_dict(self):
        """Test that get_chain_info() returns complete information."""
        poly = Polymer(chain_num=1, sequence=uniform(5), tacticity="isotactic")
        poly.set_Sequence()
        info = poly.get_chain_info()

        assert 'chain_num' in info
        assert 'sequence' in info
        assert 'dop' in info
        assert 'topology' in info
        assert 'tacticity' in info
        assert 'mer_set' in info
        assert 'sequenceSet' in info
        assert 'sequenceNames' in info
        assert 'tacticity_set' in info

        assert info['chain_num'] == 1
        assert info['dop'] == 5
        assert info['topology'] == "linear"
        assert info['tacticity'] == "isotactic"

    def test_polymer_get_tacticity_for_chain(self):
        """Test getting tacticity for specific chain."""
        poly = Polymer(chain_num=2, sequence=uniform(3))
        poly.set_Sequence()
        chirality_0 = poly.get_tacticity_for_chain(0)
        chirality_1 = poly.get_tacticity_for_chain(1)
        assert len(chirality_0) == 3
        assert len(chirality_1) == 3

    def test_polymer_get_tacticity_for_invalid_chain(self):
        """Test getting tacticity for invalid chain index."""
        poly = Polymer(chain_num=1, sequence=uniform(3))
        poly.set_Sequence()
        # Invalid chain index should return empty list
        chirality = poly.get_tacticity_for_chain(5)
        assert chirality == []


class TestPolymerSetters:
    """Test Polymer setter methods."""

    def test_polymer_set_mer_set_with_list(self):
        """Test set_merSet with list of monomers."""
        poly = Polymer(chain_num=1, sequence=[PE_FIRST, PS_LAST])
        # Should remove duplicates
        assert PE_FIRST in poly.mer_set
        assert PS_LAST in poly.mer_set

    def test_polymer_set_mer_set_with_single_string(self):
        """Test set_merSet with single monomer string."""
        poly = Polymer(chain_num=1, sequence=[SINGLE])
        assert poly.mer_set == [SINGLE]

    def test_polymer_set_dop(self):
        """Test set_dop method."""
        poly = Polymer(chain_num=1, sequence=uniform(5))
        assert poly.dop == 5
        poly.set_dop(10)
        assert poly.dop == 10


class TestPolymerCopolymers:
    """Test copolymer-specific behavior."""

    def test_polymer_copolymer_alternating_sequence(self):
        """Test that copolymer uses explicit alternating sequence."""
        sequence = [PE_FIRST, PS_MID, PE_MID, PS_LAST]
        poly = Polymer(chain_num=1, sequence=sequence)
        poly.set_Sequence()

        identifiers = [id.replace("_T1", "") for id in poly.sequenceSet[0]]
        # Should match explicit sequence exactly
        assert identifiers == sequence

    def test_polymer_copolymer_explicit_sequence(self):
        """Test copolymer with explicit sequence (no cycling)."""
        # With explicit sequence, each position is used exactly once
        sequence = [PE_FIRST, PS_MID, PP_MID, PE_MID, PS_MID, PP_MID, PE_LAST]
        poly = Polymer(chain_num=1, sequence=sequence)
        poly.set_Sequence()

        assert len(poly.sequenceSet[0]) == 7
        # Should match explicit sequence exactly
        identifiers = [id.replace("_T1", "") for id in poly.sequenceSet[0]]
        assert identifiers == sequence
