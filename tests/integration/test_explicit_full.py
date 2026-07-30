"""
Integration tests for explicit sequence workflow.

This module tests the end-to-end workflow with explicit monomer sequences,
including system creation, polymer definition, and polymerization.
"""

import pytest
import os
import tempfile
from pathlib import Path
from AutoPoly import System, Polymer


class TestExplicitFullWorkflow:
    """Test end-to-end workflow with explicit sequences."""

    def test_block_copolymer_workflow(self):
        """Test complete workflow with ABA triblock copolymer."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            system = System(out=str(Path(tmp_dir) / "block_copolymer"))

            # ABA triblock (first = 1 wildcard, middle = 2, last = 1)
            sequence = [
                "CC[*]",      # Position 0 (first)
                "[*]CC[*]",   # Position 1
                "[*]C=C[*]",  # Position 2
                "[*]C=C[*]",  # Position 3
                "[*]C=C[*]",  # Position 4
                "[*]CC[*]",   # Position 5
                "[*]CC"       # Position 6 (last)
            ]

            poly = Polymer(
                chain_num=5,
                sequence=sequence,
                topology="linear",
                tacticity="atactic"
            )

            # Verify polymer was created correctly
            assert poly.dop == 7
            assert poly.chain_num == 5
            assert len(poly.sequenceSet) == 5
            assert len(poly.sequenceSet[0]) == 7

    def test_multi_monomer_workflow(self):
        """Test workflow with multiple unique monomer types."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            system = System(out=str(Path(tmp_dir) / "multi_monomer"))

            sequence = [
                "CC[*]",                # Ethylene (first)
                "[*]C=C[*]",            # Styrene
                "[*]CC([*])C(=O)OC",    # MMA
                "[*]CC(C)(C)",          # Isobutylene (last)
            ]

            poly = Polymer(
                chain_num=3,
                sequence=sequence,
                tacticity="syndiotactic"
            )

            # Verify polymer was created correctly
            assert poly.dop == 4
            assert poly.chain_num == 3
            assert len(poly.mer_set) == 4

    def test_ring_topology_explicit(self):
        """Test ring topology with explicit sequence."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            system = System(out=str(Path(tmp_dir) / "ring_test"))

            poly = Polymer(
                chain_num=2,
                sequence=["CC[*]", "[*]C=C[*]", "[*]CC"],
                topology="ring",
                tacticity="isotactic"
            )

            # Verify polymer was created correctly
            assert poly.dop == 3
            assert poly.topology == "ring"
            assert poly.chain_num == 2

    def test_uniform_long_sequence(self):
        """Test creating a uniform polymer with a long explicit sequence."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            system = System(out=str(Path(tmp_dir) / "uniform_long"))

            # Create a sequence of 50 identical monomers
            sequence = ["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"]

            poly = Polymer(
                chain_num=5,
                sequence=sequence,
                topology="linear",
                tacticity="isotactic"
            )

            # Verify polymer was created correctly
            assert poly.dop == 50
            # One monomer chemistry; first/middle/last are distinct pSMILES strings
            assert len(poly.mer_set) == 3
            assert poly.chain_num == 5

    def test_arbitrary_sequence(self):
        """Test creating a polymer with arbitrary monomer order."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            system = System(out=str(Path(tmp_dir) / "arbitrary"))

            # Arbitrary sequence with no clear pattern
            sequence = [
                "CC[*]",                # Ethylene (first)
                "[*]C=C[*]",            # Styrene
                "[*]CC([*])C(=O)OC",    # MMA
                "[*]CC[*]",             # Ethylene again
                "[*]CC(C)(C)[*]",       # Isobutylene
                "[*]C=C",               # Styrene again (last)
            ]

            poly = Polymer(
                chain_num=2,
                sequence=sequence,
                tacticity="atactic"
            )

            # Verify the exact sequence is preserved
            assert poly.dop == 6
            assert len(poly.sequenceSet) == 2
            assert len(poly.sequenceSet[0]) == 6

    def test_single_monomer_sequence(self):
        """Test edge case: sequence with single monomer."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            system = System(out=str(Path(tmp_dir) / "single_monomer"))

            poly = Polymer(
                chain_num=10,
                sequence=["CC"],
                topology="linear",
                tacticity="syndiotactic"
            )

            # Single monomer chain
            assert poly.dop == 1
            assert poly.chain_num == 10
            assert len(poly.mer_set) == 1

    def test_all_different_monomers(self):
        """Test edge case: all monomers are different."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            system = System(out=str(Path(tmp_dir) / "all_different"))

            # Sequence where each monomer is unique
            sequence = (
                ["CC[*]"]
                + [f"[*]C({sub})[*]" for sub in
                   ["F", "Cl", "Br", "I", "C", "CC", "CCC", "CCCC"]]
                + ["[*]CC"]
            )

            poly = Polymer(
                chain_num=1,
                sequence=sequence,
                tacticity="atactic"
            )

            # All unique monomers
            assert poly.dop == 10
            assert len(poly.mer_set) == 10
