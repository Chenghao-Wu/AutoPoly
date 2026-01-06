"""
Unit tests for the Polymer class.

This module tests the Polymer class functionality including:
- Initialization with various parameters
- Topology validation (linear, ring)
- Tacticity modes (atactic, isotactic, syndiotactic)
- Sequence generation
- Chain information retrieval
"""

import os
import sys

# Add the AutoPoly package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from AutoPoly.polymer import Polymer


class TestPolymerInitialization:
    """Test Polymer class initialization."""

    def test_linear_polymer_initialization(self):
        """Test initialization of a linear polymer."""
        polymer = Polymer(
            ChainNum=1,
            Sequence=["PE", "PE"],
            DOP=2,
            topology="linear",
            tacticity="atactic"
        )
        assert polymer.ChainNum == 1
        assert polymer.sequence == ["PE", "PE"]
        assert polymer.DOP == 2
        assert polymer.topology == "linear"
        assert polymer.tacticity == "atactic"

    def test_ring_polymer_initialization(self, sample_ring_polymer):
        """Test initialization of a ring polymer."""
        assert sample_ring_polymer.topology == "ring"
        assert sample_ring_polymer.ChainNum == 1
        assert sample_ring_polymer.DOP == 3

    def test_polymer_with_dop_zero(self):
        """Test polymer initialization with DOP=0 (uses sequence length)."""
        polymer = Polymer(
            ChainNum=1,
            Sequence=["PE", "PE", "PE"],
            DOP=0,
            topology="linear"
        )
        assert polymer.DOP == 3

    def test_polymer_with_explicit_dop(self):
        """Test polymer initialization with explicit DOP."""
        polymer = Polymer(
            ChainNum=1,
            Sequence=["PE", "PE", "PE"],
            DOP=5,
            topology="linear"
        )
        assert polymer.DOP == 5

    def test_invalid_topology_raises_error(self):
        """Test that invalid topology raises ValueError."""
        with pytest.raises(ValueError, match="Topology must be either 'linear' or 'ring'"):
            Polymer(
                ChainNum=1,
                Sequence=["PE"],
                topology="invalid"
            )

    def test_none_sequence_raises_error(self):
        """Test that None sequence raises ValueError."""
        with pytest.raises(ValueError, match="Sequence cannot be None or empty"):
            Polymer(
                ChainNum=1,
                Sequence=None
            )

    def test_empty_sequence_raises_error(self):
        """Test that empty sequence raises ValueError."""
        with pytest.raises(ValueError, match="Sequence cannot be None or empty"):
            Polymer(
                ChainNum=1,
                Sequence=[]
            )


class TestPolymerTopology:
    """Test Polymer topology handling."""

    def test_linear_topology(self, sample_polymer):
        """Test linear polymer topology."""
        assert sample_polymer.topology == "linear"
        sequence_set = sample_polymer.get_sequence_set()
        assert len(sequence_set) == 1
        # Linear polymers should have different monomers for ends
        assert len(sequence_set[0]) == 2

    def test_ring_topology(self, sample_ring_polymer):
        """Test ring polymer topology."""
        assert sample_ring_polymer.topology == "ring"
        sequence_set = sample_ring_polymer.get_sequence_set()
        assert len(sequence_set) == 1
        # Ring polymers should have all internal monomers
        assert all("i.lt" in monomer for monomer in sequence_set[0])


class TestPolymerTacticity:
    """Test Polymer tacticity modes."""

    def test_atactic_tacticity(self, sample_polymer):
        """Test atactic polymer generation."""
        assert sample_polymer.tacticity == "atactic"
        sequence_set = sample_polymer.get_sequence_set()
        assert len(sequence_set) == 1

    def test_isotactic_tacticity(self, sample_isotactic_polymer):
        """Test isotactic polymer generation."""
        assert sample_isotactic_polymer.tacticity == "isotactic"
        sequence_set = sample_isotactic_polymer.get_sequence_set()
        assert len(sequence_set) == 1

    def test_syndiotactic_tacticity(self, sample_syndiotactic_polymer):
        """Test syndiotactic polymer generation."""
        assert sample_syndiotactic_polymer.tacticity == "syndiotactic"
        sequence_set = sample_syndiotactic_polymer.get_sequence_set()
        assert len(sequence_set) == 1


class TestPolymerSequence:
    """Test Polymer sequence generation and management."""

    def test_sequence_generation(self, sample_polymer):
        """Test that sequences are generated correctly."""
        sequence_set = sample_polymer.get_sequence_set()
        sequence_names = sample_polymer.get_sequence_names()

        assert len(sequence_set) == sample_polymer.ChainNum
        assert len(sequence_names) == sample_polymer.ChainNum

    def test_get_sequence_set(self, sample_polymer):
        """Test get_sequence_set returns list of lists."""
        sequence_set = sample_polymer.get_sequence_set()
        assert isinstance(sequence_set, list)
        assert len(sequence_set) > 0
        assert isinstance(sequence_set[0], list)

    def test_get_sequence_names(self, sample_polymer):
        """Test get_sequence_names returns list of lists."""
        sequence_names = sample_polymer.get_sequence_names()
        assert isinstance(sequence_names, list)
        assert len(sequence_names) > 0
        assert isinstance(sequence_names[0], list)

    def test_get_mer_set(self, sample_polymer):
        """Test get_mer_set returns unique monomers."""
        mer_set = sample_polymer.get_mer_set()
        assert isinstance(mer_set, list)
        assert len(mer_set) > 0

    def test_set_mer_set(self):
        """Test set_mer_set method."""
        polymer = Polymer(
            ChainNum=1,
            Sequence=["PE", "PE", "PS"]
        )
        polymer.set_merSet(["PE", "PS"])
        assert polymer.get_mer_set() == ["PE", "PS"]

    def test_set_dop(self):
        """Test set_dop method."""
        polymer = Polymer(
            ChainNum=1,
            Sequence=["PE", "PE"]
        )
        polymer.set_dop(10)
        assert polymer.DOP == 10


class TestPolymerChainInfo:
    """Test Polymer chain information retrieval."""

    def test_get_chain_info(self, sample_polymer):
        """Test get_chain_info returns comprehensive information."""
        info = sample_polymer.get_chain_info()

        assert isinstance(info, dict)
        assert 'chain_num' in info
        assert 'sequence' in info
        assert 'dop' in info
        assert 'topology' in info
        assert 'tacticity' in info
        assert 'sequence_length' in info
        assert 'mer_set' in info
        assert 'sequence_set' in info
        assert 'sequence_names' in info

        assert info['chain_num'] == 1
        assert info['dop'] == 2
        assert info['topology'] == "linear"
        assert info['tacticity'] == "atactic"

    def test_sequence_len_attribute_fix(self, sample_polymer):
        """Test that SequenceLen typo has been fixed."""
        # The old typo 'SequnceLen' should not exist
        assert not hasattr(sample_polymer, 'SequnceLen')
        # The correct name 'SequenceLen' should exist
        assert hasattr(sample_polymer, 'SequenceLen')
        assert sample_polymer.SequenceLen == 2


class TestPolymerMultipleChains:
    """Test Polymer with multiple chains."""

    def test_multiple_chains(self):
        """Test polymer with multiple chains."""
        polymer = Polymer(
            ChainNum=3,
            Sequence=["PE", "PE"],
            topology="linear"
        )
        sequence_set = polymer.get_sequence_set()
        sequence_names = polymer.get_sequence_names()

        assert len(sequence_set) == 3
        assert len(sequence_names) == 3
        assert polymer.ChainNum == 3


class TestPolymerEdgeCases:
    """Test edge cases in Polymer class."""

    def test_polymer_dop1_ring_topology(self):
        """Test edge case: DOP=1 with ring topology."""
        # DOP=1 with ring should be handled gracefully
        polymer = Polymer(
            ChainNum=1,
            Sequence=["PE"],
            DOP=1,
            topology="ring"
        )
        
        # Should create polymer without errors
        assert polymer.DOP == 1
        assert polymer.topology == "ring"
        
    def test_polymer_dop1_isotactic(self):
        """Test edge case: DOP=1 with tacticity."""
        # DOP=1 with isotactic tacticity
        polymer = Polymer(
            ChainNum=1,
            Sequence=["PP"],
            DOP=1,
            tacticity="isotactic"
        )
        
        # Should create polymer with tacticity set
        assert polymer.DOP == 1
        assert polymer.tacticity == "isotactic"

    def test_polymer_dop1_syndiotactic(self):
        """Test edge case: DOP=1 with syndiotactic tacticity."""
        polymer = Polymer(
            ChainNum=1,
            Sequence=["PP"],
            DOP=1,
            tacticity="syndiotactic"
        )
        
        assert polymer.DOP == 1
        assert polymer.tacticity == "syndiotactic"

    def test_polymer_sequence_cycling(self):
        """Test sequence repetition when DOP > sequence length."""
        # Create a polymer with DOP > sequence length
        # to test sequence cycling behavior
        polymer = Polymer(
            ChainNum=1,
            Sequence=["PE", "PP"],  # 2 monomers
            DOP=5,  # DOP > sequence length
            topology="linear"
        )
        
        # Should generate sequences by cycling through the monomers
        sequences = polymer.sequenceSet
        assert len(sequences) == 1  # One chain
        # Each chain should have 5 monomers (DOP=5)
        assert len(sequences[0]) == 5

    def test_polymer_single_monomer_copolymers(self):
        """Test copolymer with single monomer type."""
        polymer = Polymer(
            ChainNum=2,
            Sequence=["PE"],  # Only one monomer type
            DOP=3,
            topology="linear"
        )
        
        # Should create 2 chains, each with 3 PE monomers
        sequences = polymer.sequenceSet
        assert len(sequences) == 2
        assert len(sequences[0]) == 3
        assert len(sequences[1]) == 3

    def test_polymer_large_dop(self):
        """Test polymer with large DOP value."""
        polymer = Polymer(
            ChainNum=1,
            Sequence=["PE"],
            DOP=100,
            topology="linear"
        )
        
        # Should handle large DOP
        assert polymer.DOP == 100
        sequences = polymer.sequenceSet
        assert len(sequences[0]) == 100

    def test_polymer_multiple_chains_different_topologies(self):
        """Test multiple chains with different topologies are handled."""
        # Create first polymer with linear topology
        polymer1 = Polymer(
            ChainNum=1,
            Sequence=["PE"],
            DOP=5,
            topology="linear"
        )
        assert polymer1.topology == "linear"
        
        # Create second polymer with ring topology
        polymer2 = Polymer(
            ChainNum=1,
            Sequence=["PE"],
            DOP=5,
            topology="ring"
        )
        assert polymer2.topology == "ring"
        
        # Both should be valid
        assert polymer1.DOP == polymer2.DOP
