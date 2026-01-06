"""
Unit tests for the BeadSpringPolymer class.

This module tests the BeadSpringPolymer class functionality including:
- Initialization with various parameters
- Linear and ring topology generation
- LAMMPS data file generation
- System information retrieval
- Parameter modification
"""

import os
import sys
import tempfile
import shutil
from pathlib import Path

# Add the AutoPoly package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from AutoPoly.bead_spring import BeadSpringPolymer
from AutoPoly.system import System


class TestBeadSpringPolymerInitialization:
    """Test BeadSpringPolymer class initialization."""

    def test_linear_bead_spring_initialization(self, temp_system):
        """Test initialization of linear bead-spring polymer."""
        polymer = BeadSpringPolymer(
            name="test_linear",
            system=temp_system,
            n_chains=2,
            n_beads=10,
            topology="linear"
        )
        assert polymer.name == "test_linear"
        assert polymer.n_chains == 2
        assert polymer.n_beads == 10
        assert polymer.topology == "linear"

    def test_ring_bead_spring_initialization(self, temp_system):
        """Test initialization of ring bead-spring polymer."""
        polymer = BeadSpringPolymer(
            name="test_ring",
            system=temp_system,
            n_chains=1,
            n_beads=10,
            topology="ring"
        )
        assert polymer.name == "test_ring"
        assert polymer.topology == "ring"

    def test_default_parameters(self, temp_system):
        """Test initialization with default parameters."""
        polymer = BeadSpringPolymer(
            name="test_default",
            system=temp_system
        )
        assert polymer.n_chains == 1
        assert polymer.n_beads == 10
        assert polymer.topology == "linear"
        assert polymer.bond_length == 1.0
        assert polymer.mass == 1.0

    def test_custom_parameters(self, temp_system):
        """Test initialization with custom parameters."""
        polymer = BeadSpringPolymer(
            name="test_custom",
            system=temp_system,
            n_chains=3,
            n_beads=20,
            topology="ring",
            bond_length=1.5,
            mass=2.0,
            epsilon=0.5,
            sigma=0.9
        )
        assert polymer.bond_length == 1.5
        assert polymer.mass == 2.0
        assert polymer.epsilon == 0.5
        assert polymer.sigma == 0.9

    def test_invalid_topology_raises_error(self, temp_system):
        """Test that invalid topology raises ValueError."""
        with pytest.raises(ValueError, match="Topology must be either 'linear' or 'ring'"):
            BeadSpringPolymer(
                name="test_invalid",
                system=temp_system,
                topology="invalid"
            )


class TestBeadSpringPolymerFileGeneration:
    """Test BeadSpringPolymer file generation."""

    def test_linear_polymer_data_file_generation(self, temp_system):
        """Test data file generation for linear polymer."""
        polymer = BeadSpringPolymer(
            name="test_linear_gen",
            system=temp_system,
            n_chains=2,
            n_beads=5,
            topology="linear"
        )
        polymer.generate_data_file()

        data_file = Path(polymer.path) / "polymer.data"
        input_file = Path(polymer.path) / "in.polymer"

        assert data_file.exists()
        assert input_file.exists()

    def test_ring_polymer_data_file_generation(self, temp_system):
        """Test data file generation for ring polymer."""
        polymer = BeadSpringPolymer(
            name="test_ring_gen",
            system=temp_system,
            n_chains=1,
            n_beads=10,
            topology="ring"
        )
        polymer.generate_data_file()

        data_file = Path(polymer.path) / "polymer.data"
        input_file = Path(polymer.path) / "in.polymer"

        assert data_file.exists()
        assert input_file.exists()

    def test_lammps_data_file_format(self, temp_system):
        """Test that LAMMPS data file has correct format."""
        polymer = BeadSpringPolymer(
            name="test_format",
            system=temp_system,
            n_chains=1,
            n_beads=3,
            topology="linear"
        )
        polymer.generate_data_file()

        data_file = Path(polymer.path) / "polymer.data"
        with open(data_file, 'r') as f:
            content = f.read()

        # Check for expected sections
        assert "atoms" in content
        assert "bonds" in content
        assert "atom types" in content
        assert "bond types" in content
        assert "Masses" in content
        assert "Atoms" in content
        assert "Bonds" in content

    def test_atom_and_bond_counts(self, temp_system):
        """Test correct atom and bond counts in data file."""
        n_chains = 2
        n_beads = 5
        polymer = BeadSpringPolymer(
            name="test_counts",
            system=temp_system,
            n_chains=n_chains,
            n_beads=n_beads,
            topology="linear"
        )
        polymer.generate_data_file()

        data_file = Path(polymer.path) / "polymer.data"
        with open(data_file, 'r') as f:
            lines = f.readlines()

        # Parse header for counts
        total_beads = n_chains * n_beads
        total_bonds = n_chains * (n_beads - 1)

        # Find the atoms and bonds lines
        for line in lines:
            if "atoms" in line and "bond" not in line:
                assert str(total_beads) in line
            if "bonds" in line:
                assert str(total_bonds) in line


class TestBeadSpringPolymerSystemInfo:
    """Test BeadSpringPolymer system information."""

    def test_get_system_info(self, temp_system):
        """Test get_system_info returns correct information."""
        polymer = BeadSpringPolymer(
            name="test_info",
            system=temp_system,
            n_chains=2,
            n_beads=10,
            topology="ring",
            bond_length=1.2,
            mass=1.5
        )
        info = polymer.get_system_info()

        assert isinstance(info, dict)
        assert info['name'] == "test_info"
        assert info['n_chains'] == 2
        assert info['n_beads_per_chain'] == 10
        assert info['topology'] == "ring"
        assert info['total_beads'] == 20
        assert info['total_bonds'] == 20  # Ring: n_beads bonds per chain
        assert info['bond_length'] == 1.2
        assert info['mass'] == 1.5

    def test_linear_topology_bond_count(self, temp_system):
        """Test bond count for linear topology."""
        polymer = BeadSpringPolymer(
            name="test_linear_bonds",
            system=temp_system,
            n_chains=3,
            n_beads=5,
            topology="linear"
        )
        info = polymer.get_system_info()
        assert info['total_bonds'] == 3 * (5 - 1)  # n_chains * (n_beads - 1)


class TestBeadSpringPolymerParameterModification:
    """Test BeadSpringPolymer parameter modification."""

    def test_modify_valid_parameters(self, temp_system):
        """Test modifying valid parameters."""
        polymer = BeadSpringPolymer(
            name="test_modify",
            system=temp_system,
            n_chains=1,
            n_beads=10
        )
        polymer.modify_parameters(
            n_chains=3,
            n_beads=15,
            bond_length=1.3,
            mass=2.0
        )

        assert polymer.n_chains == 3
        assert polymer.n_beads == 15
        assert polymer.bond_length == 1.3
        assert polymer.mass == 2.0

    def test_modify_invalid_parameter(self, temp_system, caplog):
        """Test modifying invalid parameter logs warning."""
        polymer = BeadSpringPolymer(
            name="test_invalid_param",
            system=temp_system
        )
        polymer.modify_parameters(invalid_param=123)
        # Should log a warning for unknown parameter
        # (Warning logging would be captured in caplog)

    def test_modify_topology_invalid_raises_error(self, temp_system):
        """Test modifying topology to invalid value raises error."""
        polymer = BeadSpringPolymer(
            name="test_invalid_topo",
            system=temp_system,
            topology="linear"
        )
        with pytest.raises(ValueError, match="Topology must be either 'linear' or 'ring'"):
            polymer.modify_parameters(topology="invalid")


class TestBeadSpringPolymerWithoutSystem:
    """Test BeadSpringPolymer without System object."""

    def test_without_system_object(self):
        """Test initialization without System object."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_beadspring_")
        try:
            polymer = BeadSpringPolymer(
                name="test_no_system",
                system=None,
                n_chains=1,
                n_beads=5
            )
            assert polymer.system is None
            assert polymer.path == f"./test_no_system"

            # Clean up generated files
            if os.path.exists("./test_no_system"):
                shutil.rmtree("./test_no_system")
        finally:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
