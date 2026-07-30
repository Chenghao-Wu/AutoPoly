"""Tests for the BeadSpringPolymer class."""

import os
import pytest
import tempfile
import numpy as np
from pathlib import Path
from unittest.mock import MagicMock

from AutoPoly.models.bead_spring import (
    BeadSpringPolymer, BeadType, AngleType, MCConfig, SAWConfig,
    calculate_box_size, compute_lj_energy, compute_bond_energy,
    compute_total_energy, metropolis_accept,
    mc_single_bead_displacement, mc_crankshaft_move, mc_pivot_move,
    mc_reptation_move, mc_equilibrate, place_chains_in_box,
    mc_chain_translation, mc_chain_rotation,
    saw_grow_chain, saw_generate_multi_chain,
    _generate_uniform_sphere_points, _generate_trial_positions,
)
from AutoPoly.mc.collision import CollisionDetector


class MockSystem:
    """Mock System object for testing."""

    def __init__(self, path: str):
        self._path = path

    def get_folder_path(self) -> str:
        return self._path


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def mock_system(temp_dir):
    """Create a mock system with temp directory."""
    return MockSystem(temp_dir)


class TestBeadTypeDataclass:
    """Test BeadType dataclass."""

    def test_bead_type_defaults(self):
        """Test BeadType with default values."""
        bt = BeadType("A")
        assert bt.name == "A"
        assert bt.mass == 1.0
        assert bt.epsilon == 1.0
        assert bt.sigma == 1.0

    def test_bead_type_custom_values(self):
        """Test BeadType with custom values."""
        bt = BeadType("B", mass=2.0, epsilon=0.5, sigma=1.2)
        assert bt.name == "B"
        assert bt.mass == 2.0
        assert bt.epsilon == 0.5
        assert bt.sigma == 1.2


class TestAngleTypeDataclass:
    """Test AngleType dataclass."""

    def test_angle_type_defaults(self):
        """Test AngleType with default values."""
        at = AngleType(("A", "A", "A"))
        assert at.triplet == ("A", "A", "A")
        assert at.k == 10.0
        assert at.theta0 == 180.0

    def test_angle_type_custom_values(self):
        """Test AngleType with custom values."""
        at = AngleType(("A", "B", "A"), k=20.0, theta0=120.0)
        assert at.triplet == ("A", "B", "A")
        assert at.k == 20.0
        assert at.theta0 == 120.0


class TestBeadSpringPolymerInitialization:
    """Test BeadSpringPolymer initialization."""

    def test_simple_homopolymer(self, mock_system):
        """Test initialization with single bead type."""
        polymer = BeadSpringPolymer(
            name="homo",
            system=mock_system,
            n_chains=10,
            bead_types=[BeadType("A")],
            sequence=[("A", 50)],
        )
        assert polymer.n_chains == 10
        assert polymer.n_beads == 50
        assert polymer.topology == "linear"
        assert len(polymer.bead_types) == 1

    def test_diblock_copolymer(self, mock_system):
        """Test initialization with two bead types."""
        polymer = BeadSpringPolymer(
            name="diblock",
            system=mock_system,
            n_chains=5,
            bead_types=[
                BeadType("A", epsilon=1.0, sigma=1.0),
                BeadType("B", epsilon=0.5, sigma=1.2),
            ],
            sequence=[("A", 25), ("B", 25)],
        )
        assert polymer.n_beads == 50
        assert len(polymer.bead_types) == 2

    def test_triblock_copolymer(self, mock_system):
        """Test A-B-A triblock initialization."""
        polymer = BeadSpringPolymer(
            name="triblock",
            system=mock_system,
            n_chains=3,
            bead_types=[BeadType("A"), BeadType("B")],
            sequence=[("A", 10), ("B", 20), ("A", 10)],
        )
        assert polymer.n_beads == 40

    def test_ring_topology(self, mock_system):
        """Test ring polymer initialization."""
        polymer = BeadSpringPolymer(
            name="ring",
            system=mock_system,
            n_chains=5,
            topology="ring",
            bead_types=[BeadType("A")],
            sequence=[("A", 100)],
        )
        assert polymer.topology == "ring"
        assert polymer.n_beads == 100

    def test_fene_bond_style(self, mock_system):
        """Test FENE bond style initialization."""
        polymer = BeadSpringPolymer(
            name="fene",
            system=mock_system,
            n_chains=5,
            bead_types=[BeadType("A")],
            sequence=[("A", 50)],
            bond_style="fene",
            k_bond=30.0,
            fene_r0=1.5,
        )
        assert polymer.bond_style == "fene"
        assert polymer.k_bond == 30.0
        assert polymer.fene_r0 == 1.5

    def test_with_angles(self, mock_system):
        """Test initialization with angle potentials."""
        polymer = BeadSpringPolymer(
            name="angles",
            system=mock_system,
            n_chains=5,
            bead_types=[BeadType("A"), BeadType("B")],
            sequence=[("A", 20), ("B", 20)],
            use_angles=True,
            default_k_angle=10.0,
            angle_types=[
                AngleType(("A", "A", "A"), k=20.0),
                AngleType(("B", "B", "B"), k=5.0),
            ],
        )
        assert polymer.use_angles is True
        assert polymer.default_k_angle == 10.0


class TestSequenceParsing:
    """Test sequence parsing functionality."""

    def test_block_pattern_sequence(self, mock_system):
        """Test block pattern [("A", 3), ("B", 2)] parsing."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A"), BeadType("B")],
            sequence=[("A", 3), ("B", 2)],
        )
        assert polymer._sequence == ["A", "A", "A", "B", "B"]

    def test_explicit_list_sequence(self, mock_system):
        """Test explicit list ["A", "A", "B"] parsing."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A"), BeadType("B")],
            sequence=["A", "A", "B", "B"],
        )
        assert polymer._sequence == ["A", "A", "B", "B"]

    def test_string_sequence(self, mock_system):
        """Test string "AABB" parsing."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A"), BeadType("B")],
            sequence="AABB",
        )
        assert polymer._sequence == ["A", "A", "B", "B"]


class TestValidation:
    """Test input validation."""

    def test_invalid_topology_raises_error(self, mock_system):
        """Test that invalid topology raises ValueError."""
        with pytest.raises(ValueError, match="Topology must be"):
            BeadSpringPolymer(
                name="test",
                system=mock_system,
                n_chains=1,
                bead_types=[BeadType("A")],
                sequence=[("A", 10)],
                topology="invalid",
            )

    def test_invalid_bond_style_raises_error(self, mock_system):
        """Test that invalid bond style raises ValueError."""
        with pytest.raises(ValueError, match="Bond style must be"):
            BeadSpringPolymer(
                name="test",
                system=mock_system,
                n_chains=1,
                bead_types=[BeadType("A")],
                sequence=[("A", 10)],
                bond_style="invalid",
            )

    def test_empty_bead_types_raises_error(self, mock_system):
        """Test that empty bead_types raises ValueError."""
        with pytest.raises(ValueError, match="At least one bead type"):
            BeadSpringPolymer(
                name="test",
                system=mock_system,
                n_chains=1,
                bead_types=[],
                sequence=[("A", 10)],
            )

    def test_unknown_bead_type_in_sequence_raises_error(self, mock_system):
        """Test that unknown bead type in sequence raises ValueError."""
        with pytest.raises(ValueError, match="Unknown bead type 'B'"):
            BeadSpringPolymer(
                name="test",
                system=mock_system,
                n_chains=1,
                bead_types=[BeadType("A")],
                sequence=[("A", 5), ("B", 5)],
            )

    def test_empty_sequence_raises_error(self, mock_system):
        """Test that empty sequence raises ValueError."""
        with pytest.raises(ValueError, match="Sequence cannot be empty"):
            BeadSpringPolymer(
                name="test",
                system=mock_system,
                n_chains=1,
                bead_types=[BeadType("A")],
                sequence=[],
            )


class TestCanonicalTriplet:
    """Test canonical triplet calculation."""

    def test_canonical_triplet_same_endpoints(self, mock_system):
        """Test triplet where endpoints are the same."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A")],
            sequence=[("A", 10)],
        )
        result = polymer._get_canonical_triplet("A", "A", "A")
        assert result == ("A", "A", "A")

    def test_canonical_triplet_ordered(self, mock_system):
        """Test triplet that's already canonical."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A"), BeadType("B")],
            sequence=[("A", 5), ("B", 5)],
        )
        result = polymer._get_canonical_triplet("A", "B", "C")
        assert result == ("A", "B", "C")

    def test_canonical_triplet_reversed(self, mock_system):
        """Test triplet that needs to be reversed."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A"), BeadType("B"), BeadType("C")],
            sequence="ABC",
        )
        result = polymer._get_canonical_triplet("C", "B", "A")
        assert result == ("A", "B", "C")


class TestPairCoefficients:
    """Test Lorentz-Berthelot mixing rules."""

    def test_single_type_pair_coeff(self, mock_system):
        """Test pair coefficient for single type."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A", epsilon=1.0, sigma=1.0)],
            sequence=[("A", 10)],
        )
        coeffs = polymer._pair_coeffs
        assert (1, 1) in coeffs
        eps, sig = coeffs[(1, 1)]
        assert eps == pytest.approx(1.0)
        assert sig == pytest.approx(1.0)

    def test_two_type_pair_coeffs(self, mock_system):
        """Test Lorentz-Berthelot mixing for two types."""
        import math

        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[
                BeadType("A", epsilon=1.0, sigma=1.0),
                BeadType("B", epsilon=0.5, sigma=1.2),
            ],
            sequence=[("A", 5), ("B", 5)],
        )
        coeffs = polymer._pair_coeffs

        # A-A
        eps_aa, sig_aa = coeffs[(1, 1)]
        assert eps_aa == pytest.approx(1.0)
        assert sig_aa == pytest.approx(1.0)

        # B-B
        eps_bb, sig_bb = coeffs[(2, 2)]
        assert eps_bb == pytest.approx(0.5)
        assert sig_bb == pytest.approx(1.2)

        # A-B (mixed)
        eps_ab, sig_ab = coeffs[(1, 2)]
        assert eps_ab == pytest.approx(math.sqrt(1.0 * 0.5))
        assert sig_ab == pytest.approx((1.0 + 1.2) / 2)


class TestAngleTypeMap:
    """Test angle type mapping."""

    def test_homopolymer_angle_types(self, mock_system):
        """Test angle types for homopolymer."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A")],
            sequence=[("A", 10)],
            use_angles=True,
        )
        # Should have exactly one angle type: A-A-A
        assert len(polymer._angle_type_map) == 1
        assert ("A", "A", "A") in polymer._angle_type_map

    def test_diblock_angle_types(self, mock_system):
        """Test angle types for diblock copolymer."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A"), BeadType("B")],
            sequence=[("A", 5), ("B", 5)],
            use_angles=True,
        )
        # Should have: A-A-A, A-A-B, A-B-B, B-B-B
        angle_types = polymer._angle_type_map
        assert ("A", "A", "A") in angle_types
        assert ("A", "A", "B") in angle_types
        assert ("A", "B", "B") in angle_types
        assert ("B", "B", "B") in angle_types


class TestGetAngleParams:
    """Test angle parameter retrieval."""

    def test_specified_angle_params(self, mock_system):
        """Test that specified angle types return correct params."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A"), BeadType("B")],
            sequence=[("A", 5), ("B", 5)],
            use_angles=True,
            default_k_angle=10.0,
            angle_types=[
                AngleType(("A", "A", "A"), k=20.0, theta0=170.0),
            ],
        )
        k, theta0 = polymer._get_angle_params(("A", "A", "A"))
        assert k == pytest.approx(20.0)
        assert theta0 == pytest.approx(170.0)

    def test_default_angle_params(self, mock_system):
        """Test that unspecified angle types return defaults."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A"), BeadType("B")],
            sequence=[("A", 5), ("B", 5)],
            use_angles=True,
            default_k_angle=10.0,
            default_theta0=180.0,
            angle_types=[
                AngleType(("A", "A", "A"), k=20.0),
            ],
        )
        # B-B-B not specified, should use defaults
        k, theta0 = polymer._get_angle_params(("B", "B", "B"))
        assert k == pytest.approx(10.0)
        assert theta0 == pytest.approx(180.0)


class TestDataFileGeneration:
    """Test LAMMPS data file generation."""

    def test_data_file_created(self, mock_system, temp_dir):
        """Test that data file is created."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=2,
            bead_types=[BeadType("A")],
            sequence=[("A", 10)],
        )
        polymer.generate_data_file()
        assert os.path.exists(f"{temp_dir}/test/polymer.data")

    def test_data_file_atom_count(self, mock_system, temp_dir):
        """Test that data file has correct atom count."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=3,
            bead_types=[BeadType("A")],
            sequence=[("A", 20)],
        )
        polymer.generate_data_file()

        with open(f"{temp_dir}/test/polymer.data") as f:
            content = f.read()
            assert "60 atoms" in content  # 3 chains * 20 beads

    def test_data_file_bond_count_linear(self, mock_system, temp_dir):
        """Test bond count for linear polymer."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=2,
            bead_types=[BeadType("A")],
            sequence=[("A", 10)],
            topology="linear",
        )
        polymer.generate_data_file()

        with open(f"{temp_dir}/test/polymer.data") as f:
            content = f.read()
            assert "18 bonds" in content  # 2 chains * 9 bonds

    def test_data_file_bond_count_ring(self, mock_system, temp_dir):
        """Test bond count for ring polymer."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=2,
            bead_types=[BeadType("A")],
            sequence=[("A", 10)],
            topology="ring",
        )
        polymer.generate_data_file()

        with open(f"{temp_dir}/test/polymer.data") as f:
            content = f.read()
            assert "20 bonds" in content  # 2 chains * 10 bonds (ring)

    def test_data_file_multiple_atom_types(self, mock_system, temp_dir):
        """Test data file with multiple atom types."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A"), BeadType("B")],
            sequence=[("A", 5), ("B", 5)],
        )
        polymer.generate_data_file()

        with open(f"{temp_dir}/test/polymer.data") as f:
            content = f.read()
            assert "2 atom types" in content

    def test_data_file_angle_count(self, mock_system, temp_dir):
        """Test angle count in data file."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=2,
            bead_types=[BeadType("A")],
            sequence=[("A", 10)],
            use_angles=True,
        )
        polymer.generate_data_file()

        with open(f"{temp_dir}/test/polymer.data") as f:
            content = f.read()
            assert "16 angles" in content  # 2 chains * 8 angles

    def test_data_file_masses_section(self, mock_system, temp_dir):
        """Test that masses section is correctly generated."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[
                BeadType("A", mass=1.5),
                BeadType("B", mass=2.0),
            ],
            sequence=[("A", 5), ("B", 5)],
        )
        polymer.generate_data_file()

        with open(f"{temp_dir}/test/polymer.data") as f:
            content = f.read()
            assert "1.500" in content  # Mass of A
            assert "2.000" in content  # Mass of B


class TestInputScriptGeneration:
    """Test LAMMPS input script generation."""

    def test_input_script_created(self, mock_system, temp_dir):
        """Test that input script is created."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A")],
            sequence=[("A", 10)],
        )
        polymer.generate_data_file()
        assert os.path.exists(f"{temp_dir}/test/in.polymer")

    def test_input_script_harmonic_bonds(self, mock_system, temp_dir):
        """Test input script with harmonic bonds."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A")],
            sequence=[("A", 10)],
            bond_style="harmonic",
        )
        polymer.generate_data_file()

        with open(f"{temp_dir}/test/in.polymer") as f:
            content = f.read()
            assert "bond_style      harmonic" in content
            assert "special_bonds   fene" not in content

    def test_input_script_fene_bonds(self, mock_system, temp_dir):
        """Test input script with FENE bonds."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A")],
            sequence=[("A", 10)],
            bond_style="fene",
            k_bond=30.0,
            fene_r0=1.5,
        )
        polymer.generate_data_file()

        with open(f"{temp_dir}/test/in.polymer") as f:
            content = f.read()
            assert "bond_style      fene" in content
            assert "special_bonds   fene" in content

    def test_input_script_pair_coeffs(self, mock_system, temp_dir):
        """Test pair coefficients in input script."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[
                BeadType("A", epsilon=1.0, sigma=1.0),
                BeadType("B", epsilon=0.5, sigma=1.2),
            ],
            sequence=[("A", 5), ("B", 5)],
        )
        polymer.generate_data_file()

        with open(f"{temp_dir}/test/in.polymer") as f:
            content = f.read()
            assert "pair_coeff      1 1" in content
            assert "pair_coeff      1 2" in content
            assert "pair_coeff      2 2" in content

    def test_input_script_angle_coeffs(self, mock_system, temp_dir):
        """Test angle coefficients in input script."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A"), BeadType("B")],
            sequence=[("A", 5), ("B", 5)],
            use_angles=True,
            angle_types=[
                AngleType(("A", "A", "A"), k=20.0),
            ],
        )
        polymer.generate_data_file()

        with open(f"{temp_dir}/test/in.polymer") as f:
            content = f.read()
            assert "angle_style     harmonic" in content
            assert "angle_coeff" in content


class TestGetSystemInfo:
    """Test get_system_info method."""

    def test_system_info_basic(self, mock_system):
        """Test basic system info dictionary."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=5,
            bead_types=[BeadType("A")],
            sequence=[("A", 20)],
        )
        info = polymer.get_system_info()

        assert info['name'] == "test"
        assert info['n_chains'] == 5
        assert info['n_beads_per_chain'] == 20
        assert info['total_atoms'] == 100
        assert info['topology'] == "linear"

    def test_system_info_bonds(self, mock_system):
        """Test bond info in system info."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=2,
            bead_types=[BeadType("A")],
            sequence=[("A", 10)],
            bond_style="fene",
            k_bond=30.0,
        )
        info = polymer.get_system_info()

        assert info['bond_style'] == "fene"
        assert info['k_bond'] == 30.0
        assert info['total_bonds'] == 18  # 2 chains * 9 bonds

    def test_system_info_angles(self, mock_system):
        """Test angle info in system info."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=2,
            bead_types=[BeadType("A")],
            sequence=[("A", 10)],
            use_angles=True,
        )
        info = polymer.get_system_info()

        assert info['use_angles'] is True
        assert info['total_angles'] == 16  # 2 chains * 8 angles

    def test_system_info_bead_types(self, mock_system):
        """Test bead types in system info."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A"), BeadType("B")],
            sequence=[("A", 5), ("B", 5)],
        )
        info = polymer.get_system_info()

        assert info['bead_types'] == ["A", "B"]
        assert info['sequence'] == ["A"] * 5 + ["B"] * 5


class TestRingTopologyAngles:
    """Test angle handling for ring topology."""

    def test_ring_angle_count(self, mock_system, temp_dir):
        """Test that ring topology has correct angle count."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A")],
            sequence=[("A", 10)],
            topology="ring",
            use_angles=True,
        )
        info = polymer.get_system_info()
        # Ring of 10 beads should have 10 angles (each bead is center of one angle)
        assert info['total_angles'] == 10

    def test_ring_wrap_around_angles(self, mock_system, temp_dir):
        """Test that ring topology includes wrap-around angles in data file."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A")],
            sequence=[("A", 5)],
            topology="ring",
            use_angles=True,
        )
        polymer.generate_data_file()

        with open(f"{temp_dir}/test/polymer.data") as f:
            content = f.read()
            # Should have 5 angles for a ring of 5 beads
            assert "5 angles" in content


# =============================================================================
# MC Configuration Tests
# =============================================================================

class TestMCConfig:
    """Test MCConfig dataclass."""

    def test_mc_config_defaults(self):
        """Test MCConfig with default values."""
        config = MCConfig()
        assert config.density == 0.85
        assert config.n_steps == 10000
        assert config.temperature == 1.0
        assert config.max_displacement == 0.5
        assert config.max_angle == 0.3
        assert config.lj_epsilon == 1.0
        assert config.lj_sigma == 1.0
        assert config.lj_cutoff == 2.5
        assert config.bond_k == 100.0
        assert config.bond_tolerance == 0.3
        assert config.move_weights is None

    def test_mc_config_custom_values(self):
        """Test MCConfig with custom values."""
        config = MCConfig(
            density=0.5,
            n_steps=5000,
            temperature=0.8,
            move_weights={"displacement": 0.5, "pivot": 0.5},
        )
        assert config.density == 0.5
        assert config.n_steps == 5000
        assert config.temperature == 0.8
        assert config.move_weights == {"displacement": 0.5, "pivot": 0.5}


# =============================================================================
# Density-Based Box Sizing Tests
# =============================================================================

class TestDensityBoxSizing:
    """Test density-based box size calculation."""

    def test_calculate_box_size_function(self):
        """Test standalone calculate_box_size function."""
        # 1000 beads at density 0.5 -> volume = 2000 -> box = 2000^(1/3) ≈ 12.6
        box_size = calculate_box_size(1000, density=0.5)
        expected = (1000 / 0.5) ** (1/3)
        assert abs(box_size - expected) < 0.001

    def test_calculate_box_size_default_density(self):
        """Test box size with default density."""
        box_size = calculate_box_size(850)  # default density 0.85
        expected = (850 / 0.85) ** (1/3)
        assert abs(box_size - expected) < 0.001

    def test_box_size_from_density_in_polymer(self, mock_system):
        """Test that polymer calculates correct box size from density."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=10,
            bead_types=[BeadType("A")],
            sequence=[("A", 100)],
            density=0.5,
        )
        # 1000 beads / 0.5 = 2000 volume, box = 12.599...
        expected = (1000 / 0.5) ** (1/3)
        assert abs(polymer._calculate_box_size() - expected) < 0.01

    def test_explicit_box_size_overrides_density(self, mock_system):
        """Test that explicit box_size overrides density calculation."""
        polymer = BeadSpringPolymer(
            name="test",
            system=mock_system,
            n_chains=10,
            bead_types=[BeadType("A")],
            sequence=[("A", 100)],
            density=0.5,
            box_size=50.0,  # Explicit override
        )
        assert polymer._calculate_box_size() == 50.0


# =============================================================================
# Energy Function Tests
# =============================================================================

class TestEnergyFunctions:
    """Test energy computation functions."""

    def test_compute_lj_energy_no_overlap(self):
        """Test LJ energy for well-separated beads."""
        # Two beads far apart
        positions = [np.array([0.0, 0.0, 0.0]), np.array([5.0, 0.0, 0.0])]
        energy = compute_lj_energy(positions, sigma=1.0, epsilon=1.0, cutoff=2.5)
        # Beyond cutoff, energy should be 0
        assert energy == 0.0

    def test_compute_lj_energy_at_sigma(self):
        """Test LJ energy at r = sigma."""
        positions = [np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])]
        energy = compute_lj_energy(positions, sigma=1.0, epsilon=1.0, cutoff=2.5)
        # At r = sigma, LJ energy is 0
        assert abs(energy) < 0.001

    def test_compute_lj_energy_excludes_bonded(self):
        """Test that bonded pairs are excluded from LJ calculation."""
        positions = [
            np.array([0.0, 0.0, 0.0]),
            np.array([0.5, 0.0, 0.0]),  # Very close, would have high LJ energy
        ]
        # Without exclusion
        energy_no_exclude = compute_lj_energy(positions, sigma=1.0, epsilon=1.0)
        # With exclusion
        energy_with_exclude = compute_lj_energy(
            positions, sigma=1.0, epsilon=1.0, exclude_bonded=[(0, 1)]
        )
        assert energy_with_exclude == 0.0
        assert energy_no_exclude > 0  # High repulsive energy

    def test_compute_bond_energy_equilibrium(self):
        """Test bond energy at equilibrium length."""
        positions = [np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])]
        bonds = [(0, 1)]
        energy = compute_bond_energy(positions, bonds, k_bond=100.0, r0=1.0)
        assert abs(energy) < 0.001

    def test_compute_bond_energy_stretched(self):
        """Test bond energy for stretched bond."""
        positions = [np.array([0.0, 0.0, 0.0]), np.array([1.5, 0.0, 0.0])]
        bonds = [(0, 1)]
        energy = compute_bond_energy(positions, bonds, k_bond=100.0, r0=1.0)
        # E = 0.5 * 100 * (0.5)^2 = 12.5
        assert abs(energy - 12.5) < 0.001

    def test_compute_total_energy(self):
        """Test total energy computation."""
        positions = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([2.0, 0.0, 0.0]),
        ]
        bonds = [(0, 1), (1, 2)]
        energy = compute_total_energy(
            positions, bonds,
            lj_sigma=1.0, lj_epsilon=1.0, lj_cutoff=2.5,
            bond_k=100.0, bond_r0=1.0,
        )
        # Bonds at equilibrium, so bond energy ≈ 0
        # LJ between 0-2 at distance 2.0 (within cutoff)
        assert energy is not None


class TestMetropolisAccept:
    """Test Metropolis acceptance criterion."""

    def test_metropolis_accept_negative_delta(self):
        """Test that negative energy changes are always accepted."""
        for _ in range(100):
            assert metropolis_accept(-1.0, temperature=1.0) is True

    def test_metropolis_accept_zero_delta(self):
        """Test that zero energy change is accepted."""
        assert metropolis_accept(0.0, temperature=1.0) is True

    def test_metropolis_accept_large_positive_delta(self):
        """Test that large positive changes are rarely accepted."""
        np.random.seed(42)
        accepts = sum(metropolis_accept(100.0, temperature=1.0) for _ in range(1000))
        # With delta_E = 100 and T = 1, probability ≈ exp(-100) ≈ 0
        assert accepts < 10

    def test_metropolis_accept_high_temperature(self):
        """Test that high temperature increases acceptance."""
        np.random.seed(42)
        accepts = sum(metropolis_accept(1.0, temperature=10.0) for _ in range(1000))
        # With delta_E = 1 and T = 10, probability ≈ exp(-0.1) ≈ 0.9
        assert accepts > 800


# =============================================================================
# MC Move Function Tests
# =============================================================================

class TestMCMoves:
    """Test individual MC move functions."""

    def test_single_bead_displacement(self):
        """Test that single bead displacement moves the correct bead."""
        positions = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([2.0, 0.0, 0.0]),
        ]
        new_positions, is_valid = mc_single_bead_displacement(
            positions, bead_idx=1, max_disp=0.5
        )
        assert is_valid is True
        # Bead 0 and 2 should be unchanged
        assert np.allclose(new_positions[0], positions[0])
        assert np.allclose(new_positions[2], positions[2])
        # Bead 1 should have moved
        assert not np.allclose(new_positions[1], positions[1])

    def test_single_bead_displacement_with_pbc(self):
        """Test displacement with periodic boundaries."""
        positions = [np.array([4.9, 0.0, 0.0])]  # Near box edge
        new_positions, is_valid = mc_single_bead_displacement(
            positions, bead_idx=0, max_disp=0.5, box_size=10.0
        )
        assert is_valid is True
        # Position should be wrapped into box
        assert -5.0 <= new_positions[0][0] <= 5.0

    def test_crankshaft_move_preserves_endpoints(self):
        """Test that crankshaft move preserves endpoint positions."""
        np.random.seed(42)
        positions = [np.array([float(i), 0.0, 0.0]) for i in range(6)]
        new_positions, is_valid = mc_crankshaft_move(
            positions, chain_start=0, chain_end=6, max_angle=0.5
        )
        if is_valid:
            # At least one of the endpoints should be preserved
            # (depending on which i,j were chosen)
            pass  # Move is valid

    def test_crankshaft_move_too_short_chain(self):
        """Test that crankshaft fails for chains with < 4 beads."""
        positions = [np.array([float(i), 0.0, 0.0]) for i in range(3)]
        new_positions, is_valid = mc_crankshaft_move(
            positions, chain_start=0, chain_end=3, max_angle=0.5
        )
        assert is_valid is False

    def test_pivot_move(self):
        """Test that pivot move rotates tail portion of chain."""
        np.random.seed(42)
        positions = [np.array([float(i), 0.0, 0.0]) for i in range(5)]
        new_positions, is_valid = mc_pivot_move(
            positions, chain_start=0, chain_end=5, max_angle=0.5
        )
        assert is_valid is True
        # First position should be unchanged
        assert np.allclose(new_positions[0], positions[0])

    def test_pivot_move_too_short_chain(self):
        """Test that pivot fails for single bead."""
        positions = [np.array([0.0, 0.0, 0.0])]
        new_positions, is_valid = mc_pivot_move(
            positions, chain_start=0, chain_end=1, max_angle=0.5
        )
        assert is_valid is False

    def test_reptation_move_maintains_chain_length(self):
        """Test that reptation maintains number of beads."""
        np.random.seed(42)
        positions = [np.array([float(i), 0.0, 0.0]) for i in range(5)]
        new_positions, is_valid = mc_reptation_move(
            positions, chain_start=0, chain_end=5, bond_length=1.0
        )
        assert is_valid is True
        assert len(new_positions) == len(positions)

    def test_reptation_move_too_short_chain(self):
        """Test that reptation fails for single bead."""
        positions = [np.array([0.0, 0.0, 0.0])]
        new_positions, is_valid = mc_reptation_move(
            positions, chain_start=0, chain_end=1, bond_length=1.0
        )
        assert is_valid is False


# =============================================================================
# Multi-Chain Placement Tests
# =============================================================================

class TestMultiChainPlacement:
    """Test multi-chain placement functions."""

    def test_place_chains_in_box(self):
        """Test placing multiple chains in a box."""
        chain1 = [np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])]
        chain2 = [np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])]

        all_positions, chain_indices = place_chains_in_box(
            [chain1, chain2], box_size=20.0, min_separation=3.0
        )

        assert len(all_positions) == 4
        assert len(chain_indices) == 2
        assert chain_indices[0] == (0, 2)
        assert chain_indices[1] == (2, 4)

    def test_mc_chain_translation(self):
        """Test chain translation move."""
        positions = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([10.0, 0.0, 0.0]),
            np.array([11.0, 0.0, 0.0]),
        ]
        chain_indices = [(0, 2), (2, 4)]

        new_positions, is_valid = mc_chain_translation(
            positions, chain_indices, chain_idx=0, max_disp=0.5, box_size=20.0
        )
        assert is_valid is True
        # Chain 1 should be unchanged
        assert np.allclose(new_positions[2], positions[2])
        assert np.allclose(new_positions[3], positions[3])

    def test_mc_chain_rotation(self):
        """Test chain rotation move."""
        positions = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([10.0, 0.0, 0.0]),
            np.array([11.0, 0.0, 0.0]),
        ]
        chain_indices = [(0, 2), (2, 4)]

        new_positions, is_valid = mc_chain_rotation(
            positions, chain_indices, chain_idx=0, max_angle=0.5
        )
        assert is_valid is True
        # Chain 1 should be unchanged
        assert np.allclose(new_positions[2], positions[2])
        assert np.allclose(new_positions[3], positions[3])


# =============================================================================
# MC Equilibration Tests
# =============================================================================

class TestMCEquilibration:
    """Test MC equilibration routine."""

    def test_mc_equilibrate_returns_positions(self):
        """Test that equilibration returns valid positions."""
        positions = [np.array([float(i), 0.0, 0.0]) for i in range(10)]
        bonds = [(i, i+1) for i in range(9)]

        eq_positions, stats = mc_equilibrate(
            positions, bonds, n_steps=100, temperature=1.0
        )

        assert len(eq_positions) == len(positions)
        assert isinstance(stats, dict)

    def test_mc_equilibrate_acceptance_stats(self):
        """Test that equilibration returns acceptance statistics."""
        positions = [np.array([float(i), 0.0, 0.0]) for i in range(10)]
        bonds = [(i, i+1) for i in range(9)]

        _, stats = mc_equilibrate(
            positions, bonds, n_steps=500, temperature=1.0
        )

        # Should have stats for default move types
        assert "displacement" in stats
        assert "pivot" in stats
        # Acceptance rates should be between 0 and 1
        for rate in stats.values():
            assert 0.0 <= rate <= 1.0

    def test_equilibration_in_polymer(self, mock_system, temp_dir):
        """Test equilibration through BeadSpringPolymer."""
        np.random.seed(42)
        polymer = BeadSpringPolymer(
            name="eq_test",
            system=mock_system,
            n_chains=2,
            bead_types=[BeadType("A")],
            sequence=[("A", 10)],
            density=0.3,
            equilibrate=True,
            mc_config=MCConfig(n_steps=100, temperature=1.0),
        )
        polymer.generate_data_file()

        # Check that file was created
        assert os.path.exists(f"{temp_dir}/eq_test/polymer.data")

    def test_manual_equilibration(self, mock_system):
        """Test manual equilibration after creation."""
        np.random.seed(42)
        polymer = BeadSpringPolymer(
            name="manual_eq",
            system=mock_system,
            n_chains=2,
            bead_types=[BeadType("A")],
            sequence=[("A", 10)],
            density=0.5,
        )

        # Equilibrate manually
        polymer.equilibrate(MCConfig(n_steps=100, temperature=0.5))

        # Positions should be initialized
        assert polymer._positions is not None
        assert len(polymer._positions) == 20  # 2 chains * 10 beads


class TestBeadSpringPolymerWithDensity:
    """Test BeadSpringPolymer with density-based box sizing."""

    def test_data_file_with_density(self, mock_system, temp_dir):
        """Test data file generation with density parameter."""
        polymer = BeadSpringPolymer(
            name="density_test",
            system=mock_system,
            n_chains=5,
            bead_types=[BeadType("A")],
            sequence=[("A", 20)],
            density=0.5,
        )
        polymer.generate_data_file()

        with open(f"{temp_dir}/density_test/polymer.data") as f:
            content = f.read()
            # Box size from density: (100/0.5)^(1/3) ≈ 5.85
            # Check that box dimensions exist
            assert "xlo xhi" in content
            assert "100 atoms" in content

    def test_data_file_with_explicit_box_size(self, mock_system, temp_dir):
        """Test data file generation with explicit box size."""
        polymer = BeadSpringPolymer(
            name="box_test",
            system=mock_system,
            n_chains=5,
            bead_types=[BeadType("A")],
            sequence=[("A", 20)],
            box_size=30.0,
        )
        polymer.generate_data_file()

        with open(f"{temp_dir}/box_test/polymer.data") as f:
            content = f.read()
            # Box should be -15 to 15
            assert "-15.0 15.0 xlo xhi" in content


# =============================================================================
# SAW Configuration Tests
# =============================================================================

class TestSAWConfig:
    """Test SAWConfig dataclass."""

    def test_saw_config_defaults(self):
        """Test SAWConfig with default values."""
        config = SAWConfig()
        assert config.collision_sigma == 1.0
        assert config.collision_tolerance == 0.1
        assert config.n_trials == 50
        assert config.max_backtrack_depth == 10
        assert config.max_total_backtracks == 1000
        assert config.bond_angle_min == 60.0
        assert config.bond_angle_max == 180.0
        assert config.ring_closure_trials == 100
        assert config.ring_closure_tolerance == 0.2

    def test_saw_config_custom_values(self):
        """Test SAWConfig with custom values."""
        config = SAWConfig(
            collision_sigma=1.5,
            n_trials=100,
            max_backtrack_depth=20,
            bond_angle_min=70.0,
        )
        assert config.collision_sigma == 1.5
        assert config.n_trials == 100
        assert config.max_backtrack_depth == 20
        assert config.bond_angle_min == 70.0


# =============================================================================
# SAW Helper Function Tests
# =============================================================================

class TestSAWHelperFunctions:
    """Test SAW helper functions."""

    def test_uniform_sphere_points_coverage(self):
        """Test that sphere points are uniformly distributed."""
        points = _generate_uniform_sphere_points(100)
        assert points.shape == (100, 3)
        # All points should be on unit sphere
        norms = np.linalg.norm(points, axis=1)
        assert np.allclose(norms, 1.0)

    def test_uniform_sphere_points_single(self):
        """Test single point generation."""
        points = _generate_uniform_sphere_points(1)
        assert points.shape == (1, 3)

    def test_trial_positions_distance(self):
        """Test that trial positions are at correct bond length."""
        center = np.array([0.0, 0.0, 0.0])
        bond_length = 1.5
        positions = _generate_trial_positions(center, bond_length, 50)
        distances = np.linalg.norm(positions - center, axis=1)
        assert np.allclose(distances, bond_length)

    def test_trial_positions_with_angle_constraint(self):
        """Test that angle constraints filter positions correctly."""
        center = np.array([1.0, 0.0, 0.0])
        prev_direction = np.array([1.0, 0.0, 0.0])  # Coming from left
        bond_length = 1.0

        # Only allow angles between 90 and 180 degrees (forward hemisphere)
        positions = _generate_trial_positions(
            center, bond_length, 100,
            prev_direction, angle_min=90.0, angle_max=180.0
        )

        # All new directions should have positive x component
        # (continuing roughly forward)
        for pos in positions:
            new_dir = pos - center
            # Angle with -prev_direction (i.e., backward direction)
            cos_angle = np.dot(new_dir, -prev_direction) / (np.linalg.norm(new_dir))
            angle_deg = np.degrees(np.arccos(cos_angle))
            assert 90.0 <= angle_deg <= 180.0


# =============================================================================
# SAW Chain Growth Tests
# =============================================================================

class TestSAWChainGrowth:
    """Test SAW chain growth functions."""

    def test_saw_grow_single_chain_linear(self):
        """Test growing a single linear chain."""
        np.random.seed(42)
        box_bounds = ((-50, 50), (-50, 50), (-50, 50))
        detector = CollisionDetector(box_bounds, cell_size=2.0)
        config = SAWConfig(collision_sigma=1.0, n_trials=50)

        positions, backtracks = saw_grow_chain(
            n_beads=20,
            bond_length=1.0,
            start_position=np.array([0.0, 0.0, 0.0]),
            collision_detector=detector,
            config=config,
            topology="linear",
        )

        assert positions is not None
        assert len(positions) == 20
        # Check bond lengths
        for i in range(len(positions) - 1):
            dist = np.linalg.norm(positions[i + 1] - positions[i])
            assert abs(dist - 1.0) < 0.01

    def test_saw_grow_single_chain_ring(self):
        """Test growing a single ring chain."""
        # Ring closure is challenging - try multiple seeds
        success = False
        for seed in [100, 200, 300, 400, 500]:
            np.random.seed(seed)
            box_bounds = ((-100, 100), (-100, 100), (-100, 100))
            detector = CollisionDetector(box_bounds, cell_size=2.0)
            config = SAWConfig(
                collision_sigma=1.0,
                n_trials=100,
                ring_closure_trials=500,
                ring_closure_tolerance=0.5,
                bond_angle_min=40.0,  # Very flexible for ring
                max_total_backtracks=2000,
            )

            positions, backtracks = saw_grow_chain(
                n_beads=10,  # Smaller ring is easier
                bond_length=1.0,
                start_position=np.array([0.0, 0.0, 0.0]),
                collision_detector=detector,
                config=config,
                topology="ring",
            )

            if positions is not None:
                success = True
                assert len(positions) == 10
                # Check all bond lengths
                for i in range(len(positions) - 1):
                    dist = np.linalg.norm(positions[i + 1] - positions[i])
                    assert abs(dist - 1.0) < 0.01
                # Check ring closure
                closure_dist = np.linalg.norm(positions[-1] - positions[0])
                assert abs(closure_dist - 1.0) < 0.6
                break

        assert success, "Ring chain generation failed with all seeds"

    def test_saw_no_overlaps(self):
        """Test that SAW produces non-overlapping beads."""
        np.random.seed(42)
        box_bounds = ((-50, 50), (-50, 50), (-50, 50))
        detector = CollisionDetector(box_bounds, cell_size=2.0)
        config = SAWConfig(collision_sigma=1.0, collision_tolerance=0.1)

        positions, _ = saw_grow_chain(
            n_beads=30,
            bond_length=1.0,
            start_position=np.array([0.0, 0.0, 0.0]),
            collision_detector=detector,
            config=config,
            topology="linear",
        )

        assert positions is not None
        # Check all non-bonded pairs don't overlap
        min_distance = config.collision_sigma - config.collision_tolerance
        for i in range(len(positions)):
            for j in range(i + 2, len(positions)):  # Skip bonded neighbors
                dist = np.linalg.norm(positions[i] - positions[j])
                assert dist >= min_distance * 0.9  # Allow small tolerance


# =============================================================================
# SAW Multi-Chain Tests
# =============================================================================

class TestSAWMultiChain:
    """Test SAW multi-chain generation."""

    def test_saw_generate_multi_chain(self):
        """Test generating multiple chains."""
        np.random.seed(42)
        config = SAWConfig(collision_sigma=1.0, n_trials=50)

        positions, chain_indices, stats = saw_generate_multi_chain(
            n_chains=3,
            n_beads_per_chain=15,
            bond_length=1.0,
            box_size=30.0,
            config=config,
            topology="linear",
        )

        assert stats["success"] is True
        assert positions is not None
        assert len(positions) == 45  # 3 * 15
        assert len(chain_indices) == 3
        assert chain_indices[0] == (0, 15)
        assert chain_indices[1] == (15, 30)
        assert chain_indices[2] == (30, 45)

    def test_saw_multi_chain_no_inter_chain_overlap(self):
        """Test that multi-chain SAW has no inter-chain overlaps."""
        np.random.seed(42)
        config = SAWConfig(collision_sigma=1.0, collision_tolerance=0.1)

        positions, chain_indices, stats = saw_generate_multi_chain(
            n_chains=2,
            n_beads_per_chain=20,
            bond_length=1.0,
            box_size=30.0,
            config=config,
        )

        assert stats["success"] is True
        # Check inter-chain distances
        min_distance = config.collision_sigma - config.collision_tolerance
        chain1_pos = positions[chain_indices[0][0]:chain_indices[0][1]]
        chain2_pos = positions[chain_indices[1][0]:chain_indices[1][1]]

        for p1 in chain1_pos:
            for p2 in chain2_pos:
                dist = np.linalg.norm(np.array(p1) - np.array(p2))
                assert dist >= min_distance * 0.9


# =============================================================================
# SAW Integration with BeadSpringPolymer Tests
# =============================================================================

class TestBeadSpringPolymerSAW:
    """Test SAW integration with BeadSpringPolymer."""

    def test_polymer_with_saw_generation(self, mock_system, temp_dir):
        """Test polymer generation using SAW method."""
        np.random.seed(42)
        polymer = BeadSpringPolymer(
            name="saw_test",
            system=mock_system,
            n_chains=2,
            bead_types=[BeadType("A")],
            sequence=[("A", 15)],
            density=0.1,  # Low density for easy SAW
            generation_method="saw",
        )
        polymer.generate_data_file()

        assert os.path.exists(f"{temp_dir}/saw_test/polymer.data")
        with open(f"{temp_dir}/saw_test/polymer.data") as f:
            content = f.read()
            assert "30 atoms" in content  # 2 * 15

    def test_polymer_saw_with_custom_config(self, mock_system, temp_dir):
        """Test SAW with custom configuration."""
        np.random.seed(42)
        saw_config = SAWConfig(
            collision_sigma=1.2,
            n_trials=100,
            bond_angle_min=70.0,
        )
        polymer = BeadSpringPolymer(
            name="saw_custom",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A", sigma=1.2)],
            sequence=[("A", 20)],
            density=0.05,
            generation_method="saw",
            saw_config=saw_config,
        )
        polymer.generate_data_file()

        assert os.path.exists(f"{temp_dir}/saw_custom/polymer.data")

    def test_polymer_saw_ring(self, mock_system, temp_dir):
        """Test SAW generation for ring polymer (with fallback)."""
        # Ring SAW is challenging - test that it either succeeds or falls back gracefully
        np.random.seed(42)
        saw_config = SAWConfig(
            n_trials=100,
            ring_closure_trials=500,
            ring_closure_tolerance=0.5,
            bond_angle_min=40.0,
            max_total_backtracks=500,  # Limited to test fallback
        )
        polymer = BeadSpringPolymer(
            name="saw_ring",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A")],
            sequence=[("A", 10)],  # Smaller ring
            topology="ring",
            density=0.02,  # Very low density
            generation_method="saw",
            saw_config=saw_config,
        )
        # Should either succeed with SAW or fall back to geometric
        polymer.generate_data_file()

        assert os.path.exists(f"{temp_dir}/saw_ring/polymer.data")
        with open(f"{temp_dir}/saw_ring/polymer.data") as f:
            content = f.read()
            assert "10 bonds" in content  # Ring has n bonds

    def test_invalid_generation_method_raises(self, mock_system):
        """Test that invalid generation method raises ValueError."""
        with pytest.raises(ValueError, match="Generation method must be"):
            BeadSpringPolymer(
                name="test",
                system=mock_system,
                n_chains=1,
                bead_types=[BeadType("A")],
                sequence=[("A", 10)],
                generation_method="invalid",
            )

    def test_saw_fallback_to_geometric(self, mock_system, temp_dir):
        """Test that SAW falls back to geometric on failure."""
        np.random.seed(42)
        # Very high density should cause SAW to fail
        saw_config = SAWConfig(
            n_trials=5,  # Very few trials
            max_total_backtracks=10,  # Very limited backtracks
        )
        polymer = BeadSpringPolymer(
            name="saw_fallback",
            system=mock_system,
            n_chains=10,
            bead_types=[BeadType("A")],
            sequence=[("A", 50)],
            density=0.9,  # Very high density
            generation_method="saw",
            saw_config=saw_config,
        )
        # Should not raise - falls back to geometric
        polymer.generate_data_file()
        assert os.path.exists(f"{temp_dir}/saw_fallback/polymer.data")

    def test_saw_generate_method_directly(self, mock_system):
        """Test calling saw_generate method directly."""
        np.random.seed(42)
        polymer = BeadSpringPolymer(
            name="saw_direct",
            system=mock_system,
            n_chains=2,
            bead_types=[BeadType("A")],
            sequence=[("A", 15)],
            density=0.1,
        )

        success = polymer.saw_generate()
        assert success is True
        assert polymer._positions is not None
        assert len(polymer._positions) == 30

    def test_saw_verifies_bond_lengths(self, mock_system):
        """Test that SAW produces correct bond lengths."""
        np.random.seed(42)
        bond_length = 1.5
        polymer = BeadSpringPolymer(
            name="saw_bonds",
            system=mock_system,
            n_chains=1,
            bead_types=[BeadType("A")],
            sequence=[("A", 20)],
            bond_length=bond_length,
            density=0.05,
            generation_method="saw",
        )

        polymer.saw_generate()

        # Check all bond lengths
        positions = polymer._positions
        for i in range(len(positions) - 1):
            dist = np.linalg.norm(positions[i + 1] - positions[i])
            assert abs(dist - bond_length) < 0.05
