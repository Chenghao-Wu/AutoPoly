"""Integration tests for full AutoPoly workflow."""

import pytest
from pathlib import Path
from AutoPoly.system import System
from AutoPoly.polymer import Polymer
from AutoPoly.molecule import Molecule


@pytest.mark.integration
class TestFullWorkflow:
    """Test complete polymer generation workflow."""

    def test_pe_linear_polymer_workflow_structure(self, tmp_path):
        """Test PE linear polymer generation structure."""
        # Setup
        output_dir = tmp_path / "pe_test"
        system = System(out=str(output_dir))

        # Create polymer
        poly = Polymer(chain_num=1, sequence=["PE"] * 5, tacticity="isotactic")
        poly.set_Sequence()

        # Verify polymer structure
        assert len(poly.sequenceSet) == 1
        assert len(poly.sequenceSet[0]) == 5
        # Isotactic: all chirality values should be the same
        assert all(poly.tacticity_set[0]) or all(not t for t in poly.tacticity_set[0])

    def test_ps_atactic_polymer_workflow_structure(self, tmp_path):
        """Test PS atactic polymer with random chirality."""
        output_dir = tmp_path / "ps_test"
        system = System(out=str(output_dir))

        poly = Polymer(chain_num=1, sequence=["PS"] * 3, tacticity="atactic")
        poly.set_Sequence()

        # Verify chirality is assigned (random but deterministic)
        assert len(poly.tacticity_set[0]) == 3
        assert all(isinstance(t, bool) for t in poly.tacticity_set[0])

    def test_water_molecule_workflow(self, tmp_path):
        """Test water molecule generation."""
        output_dir = tmp_path / "water_test"
        system = System(out=str(output_dir))

        mol = Molecule(Count=10, Smiles="O")

        # Verify structure
        assert mol.Count == 10
        assert len(mol.sequenceSet) == 10
        assert mol.DOP == 1
        assert all(item == ["molecule_O"] for item in mol.sequenceSet)

    def test_copolymer_workflow_structure(self, tmp_path):
        """Test PE-PS copolymer generation."""
        output_dir = tmp_path / "copolymer_test"
        system = System(out=str(output_dir))

        poly = Polymer(chain_num=1, sequence=["PE", "PS", "PE", "PS"])
        poly.set_Sequence()

        # Verify copolymer structure
        assert len(poly.sequenceSet[0]) == 4
        mer_set = poly.get_mer_set()
        assert "PE" in mer_set
        assert "PS" in mer_set

    def test_ring_polymer_workflow(self, tmp_path):
        """Test ring polymer topology."""
        output_dir = tmp_path / "ring_test"
        system = System(out=str(output_dir))

        poly = Polymer(chain_num=1, sequence=["PE"] * 5, topology="ring")
        poly.set_Sequence()

        # Verify ring topology
        assert poly.topology == "ring"
        assert len(poly.sequenceSet[0]) == 5

    def test_multiple_chains_workflow(self, tmp_path):
        """Test multiple chains generation."""
        output_dir = tmp_path / "multi_chain_test"
        system = System(out=str(output_dir))

        poly = Polymer(chain_num=3, sequence=["PE"] * 4, tacticity="isotactic")
        poly.set_Sequence()

        # Verify multiple chains
        assert len(poly.sequenceSet) == 3
        assert all(len(chain) == 4 for chain in poly.sequenceSet)
        assert len(poly.tacticity_set) == 3

    def test_syndiotactic_polymer_workflow(self, tmp_path):
        """Test syndiotactic polymer with alternating chirality."""
        output_dir = tmp_path / "syndiotactic_test"
        system = System(out=str(output_dir))

        poly = Polymer(chain_num=1, sequence=["PE"] * 6, tacticity="syndiotactic")
        poly.set_Sequence()

        # Verify alternating chirality
        chirality = poly.tacticity_set[0]
        for i in range(len(chirality) - 1):
            assert chirality[i] != chirality[i + 1]

    def test_benzene_molecule_workflow(self, tmp_path):
        """Test benzene molecule generation."""
        output_dir = tmp_path / "benzene_test"
        system = System(out=str(output_dir))

        mol = Molecule(Count=5, Smiles="c1ccccc1")

        # Verify structure
        assert mol.Count == 5
        assert len(mol.sequenceSet) == 5
        assert mol.DOP == 1
        # All sequence items should be the same
        assert all(item == ["molecule_c1ccccc1"] for item in mol.sequenceSet)
