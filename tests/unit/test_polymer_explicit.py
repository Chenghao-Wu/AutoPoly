"""
Tests for explicit monomer sequence functionality.

This module tests the new explicit sequence API where users specify
the exact monomer at each position in the polymer chain.
"""

import pytest
from AutoPoly import Polymer
from AutoPoly.exceptions import ValidationError


class TestExplicitSequences:
    """Test explicit monomer sequence functionality."""

    def test_explicit_sequence_basic(self):
        """Test basic explicit sequence creation."""
        sequence = ["[*]CC[*]", "[*]C=C[*]", "[*]CC[*]"]
        poly = Polymer(
            chain_num=1,
            sequence=sequence,
            topology="linear",
            tacticity="atactic"
        )
        assert poly.dop == 3
        assert len(poly.sequence_set) == 1
        assert len(poly.sequence_set[0]) == 3

    def test_block_copolymer(self):
        """Test ABA triblock copolymer creation."""
        # ABA: 2 PE-like, 3 PS-like, 2 PE-like
        sequence = [
            "[*]CC[*]",  # Position 0: Ethylene
            "[*]CC[*]",  # Position 1: Ethylene
            "[*]C=C[*]",  # Position 2: Styrene
            "[*]C=C[*]",  # Position 3: Styrene
            "[*]C=C[*]",  # Position 4: Styrene
            "[*]CC[*]",  # Position 5: Ethylene
            "[*]CC[*]"   # Position 6: Ethylene
        ]
        poly = Polymer(chain_num=5, sequence=sequence, tacticity="isotactic")
        assert poly.dop == 7
        assert poly.chain_num == 5
        assert len(poly.sequence_set) == 5

    def test_dop_derived_from_sequence(self):
        """Test that DOP is automatically derived from sequence length."""
        poly = Polymer(
            chain_num=1,
            sequence=["[*]CC[*]"] * 50
        )
        assert poly.dop == 50

    def test_pythonic_naming(self):
        """Test that Pythonic parameter names are used."""
        poly = Polymer(
            chain_num=3,  # Not ChainNum
            sequence=["[*]CC[*]"],  # Not Sequence
            topology="linear",
            tacticity="atactic"
        )
        assert hasattr(poly, 'chain_num')
        assert hasattr(poly, 'sequence')
        assert hasattr(poly, 'topology')
        assert hasattr(poly, 'tacticity')
        assert hasattr(poly, 'dop')
        assert hasattr(poly, 'sequence_set')
        assert hasattr(poly, 'mer_set')

    def test_global_tacticity_isotactic(self):
        """Test that isotactic tacticity applies uniformly."""
        sequence = ["[*]CC[*]", "[*]C=C[*]", "[*]CC[*]"]
        poly = Polymer(
            chain_num=2,
            sequence=sequence,
            tacticity="isotactic"
        )
        # All positions in a chain should have same chirality
        assert len(set(poly.tacticity_set[0])) == 1

    def test_global_tacticity_syndiotactic(self):
        """Test that syndiotactic tacticity alternates."""
        sequence = ["[*]CC[*]"] * 4
        poly = Polymer(
            chain_num=1,
            sequence=sequence,
            tacticity="syndiotactic"
        )
        # Should alternate: False, True, False, True (position % 2 == 1)
        # Position 0: 0 % 2 = 0 = False
        # Position 1: 1 % 2 = 1 = True
        # Position 2: 2 % 2 = 0 = False
        # Position 3: 3 % 2 = 1 = True
        expected = [False, True, False, True]
        assert poly.tacticity_set[0] == expected

    def test_no_cycling_mode(self):
        """Test that explicit sequences don't cycle (old behavior removed)."""
        # Create a sequence with 2 elements, DOP should be 2 (not cycled)
        sequence = ["[*]CC[*]", "[*]C=C[*]"]
        poly = Polymer(chain_num=1, sequence=sequence)

        # DOP is 2, not some larger number
        assert poly.dop == 2

        # Each position is unique (no cycling)
        assert len(poly.sequence_set[0]) == 2
        assert "[*]CC[*]" in poly.sequence_set[0][0] or "_T1" in poly.sequence_set[0][0]
        assert "[*]C=C[*]" in poly.sequence_set[0][1] or "_T1" in poly.sequence_set[0][1]

    def test_multiple_unique_monomers(self):
        """Test polymer with many unique monomer types."""
        # Use monomer names (strings that look like SMILES but won't fail validation)
        sequence = ["PE", "PS", "PMMA", "PB"]
        poly = Polymer(
            chain_num=2,
            sequence=sequence,
            tacticity="atactic"
        )
        assert poly.dop == 4
        # All four monomers should be recognized
        assert len(poly.mer_set) == 4

    def test_get_chain_info_pythonic(self):
        """Test that get_chain_info() returns Pythonic keys."""
        poly = Polymer(
            chain_num=3,
            sequence=["[*]CC[*]", "[*]C=C[*]", "[*]CC[*]"],
            tacticity="isotactic"
        )
        info = poly.get_chain_info()

        # Check Pythonic naming
        assert 'chain_num' in info
        assert 'dop' in info
        assert 'sequence' in info
        assert 'mer_set' in info
        assert 'sequence_set' in info
        assert 'tacticity_set' in info

        # Check old names are gone
        assert 'ChainNum' not in info
        assert 'DOP' not in info
        assert 'SequenceLen' not in info
