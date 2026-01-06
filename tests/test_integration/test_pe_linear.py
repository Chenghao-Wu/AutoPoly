"""
Integration tests for AutoPoly package.

These tests verify the complete workflow from System creation to Polymer
definition to Polymerization, including file generation and output management.
"""

import os
import sys
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch, Mock

# Add the AutoPoly package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

import pytest
from AutoPoly.system import System
from AutoPoly.polymer import Polymer
from AutoPoly.bead_spring import BeadSpringPolymer


@pytest.mark.slow
@pytest.mark.integration
class TestSystemToBeadSpringIntegration:
    """Test integration between System and BeadSpringPolymer."""

    def test_full_bead_spring_workflow(self):
        """Test complete workflow: System → BeadSpringPolymer → file generation."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_integration_")

        try:
            # Step 1: Create System
            system = System(out=temp_dir)
            assert system.get_output_path() == temp_dir
            assert os.path.exists(temp_dir)

            # Step 2: Create BeadSpringPolymer
            polymer = BeadSpringPolymer(
                name="integration_test",
                system=system,
                n_chains=2,
                n_beads=5,
                topology="linear"
            )
            assert polymer.system == system

            # Step 3: Generate files
            polymer.generate_data_file()

            # Step 4: Verify output
            expected_path = Path(temp_dir) / "integration_test"
            assert expected_path.exists()
            assert (expected_path / "polymer.data").exists()
            assert (expected_path / "in.polymer").exists()

        finally:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)


@pytest.mark.slow
@pytest.mark.integration
class TestPolymerCreationIntegration:
    """Test Polymer object creation and configuration."""

    def test_linear_pe_polymer_creation(self):
        """Test creation of linear PE polymer."""
        polymer = Polymer(
            ChainNum=2,
            Sequence=["PE", "PE", "PE"],
            DOP=3,
            topology="linear",
            tacticity="atactic"
        )

        assert polymer.ChainNum == 2
        assert polymer.topology == "linear"
        assert polymer.tacticity == "atactic"
        assert polymer.SequenceLen == 3

        # Verify sequences are generated
        sequence_set = polymer.get_sequence_set()
        sequence_names = polymer.get_sequence_names()

        assert len(sequence_set) == 2
        assert len(sequence_names) == 2

    def test_ring_polymer_creation(self):
        """Test creation of ring polymer."""
        polymer = Polymer(
            ChainNum=1,
            Sequence=["PE", "PE", "PE", "PE"],
            DOP=4,
            topology="ring",
            tacticity="isotactic"
        )

        assert polymer.topology == "ring"
        assert polymer.tacticity == "isotactic"

        # Ring polymers should have all internal monomers
        sequence_set = polymer.get_sequence_set()
        assert all("i.lt" in monomer for monomer in sequence_set[0])

    def test_copolymer_creation(self):
        """Test creation of copolymer."""
        polymer = Polymer(
            ChainNum=1,
            Sequence=["PE", "PS", "PE", "PS"],
            DOP=4,
            topology="linear",
            tacticity="atactic"
        )

        mer_set = polymer.get_mer_set()
        assert len(mer_set) == 2  # PE and PS
        assert "PE" in mer_set
        assert "PS" in mer_set


@pytest.mark.slow
@pytest.mark.integration
class TestMultiplePolymerTypes:
    """Test multiple polymer types in a single system."""

    def test_multiple_bead_spring_polymers(self):
        """Test creating multiple bead-spring polymers in one system."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_multi_")

        try:
            system = System(out=temp_dir)

            # Create first polymer
            polymer1 = BeadSpringPolymer(
                name="polymer1",
                system=system,
                n_chains=1,
                n_beads=3,
                topology="linear"
            )
            polymer1.generate_data_file()

            # Create second polymer
            polymer2 = BeadSpringPolymer(
                name="polymer2",
                system=system,
                n_chains=1,
                n_beads=5,
                topology="ring"
            )
            polymer2.generate_data_file()

            # Verify both outputs exist
            path1 = Path(temp_dir) / "polymer1"
            path2 = Path(temp_dir) / "polymer2"

            assert path1.exists()
            assert path2.exists()
            assert (path1 / "polymer.data").exists()
            assert (path2 / "polymer.data").exists()

        finally:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)


@pytest.mark.slow
@pytest.mark.integration
class TestOutputDirectoryStructure:
    """Test output directory structure and file organization."""

    def test_bead_spring_output_structure(self):
        """Test that BeadSpringPolymer creates proper output structure."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_structure_")

        try:
            system = System(out=temp_dir)
            polymer = BeadSpringPolymer(
                name="structure_test",
                system=system,
                n_chains=1,
                n_beads=3
            )
            polymer.generate_data_file()

            output_path = Path(temp_dir) / "structure_test"

            # Check directory exists
            assert output_path.exists()
            assert output_path.is_dir()

            # Check required files
            required_files = ["polymer.data", "in.polymer"]
            for file_name in required_files:
                file_path = output_path / file_name
                assert file_path.exists(), f"Missing file: {file_name}"
                assert file_path.is_file(), f"Not a file: {file_path}"

            # Check data file content
            data_file = output_path / "polymer.data"
            with open(data_file, 'r') as f:
                content = f.read()
                # Verify essential sections
                assert "atoms" in content
                assert "bonds" in content
                assert "Masses" in content
                assert "Atoms" in content
                assert "Bonds" in content

        finally:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)


@pytest.mark.slow
@pytest.mark.integration
class TestPolymerInfoRetrieval:
    """Test polymer information retrieval across components."""

    def test_system_chain_info(self):
        """Test System path information."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_info_")

        try:
            system = System(out="test_output")
            output_path = system.get_output_path()

            assert isinstance(output_path, str)
            assert os.path.isabs(output_path)

        finally:
            # Cleanup
            if os.path.exists("test_output"):
                shutil.rmtree("test_output")
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)

    def test_bead_spring_system_info(self):
        """Test BeadSpringPolymer system info retrieval."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_bs_info_")

        try:
            system = System(out=temp_dir)
            polymer = BeadSpringPolymer(
                name="info_test",
                system=system,
                n_chains=2,
                n_beads=4,
                topology="linear"
            )

            info = polymer.get_system_info()

            # Verify all expected keys are present
            expected_keys = [
                'name', 'n_chains', 'n_beads_per_chain',
                'topology', 'total_beads', 'total_bonds',
                'bond_length', 'mass', 'epsilon', 'sigma', 'output_path'
            ]
            for key in expected_keys:
                assert key in info, f"Missing key: {key}"

            # Verify values
            assert info['n_chains'] == 2
            assert info['n_beads_per_chain'] == 4
            assert info['total_beads'] == 8
            assert info['total_bonds'] == 6  # 2 chains * (4-1) bonds

        finally:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
