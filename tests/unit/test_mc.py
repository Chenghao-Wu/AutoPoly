#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit Tests for Monte Carlo Placement Module

Tests for collision detection, chain growth, and molecular placement.
"""
import pytest
import numpy as np
from pathlib import Path
import tempfile
import os

# Import MC module components
from AutoPoly.mc.collision import (
    CollisionDetector,
    MonomerSphere,
    calculate_box_size,
)

from AutoPoly.mc.chain_growth import (
    AtomData,
    MonomerTemplate,
    MonomerPlacement,
    ChainGrowthMC,
    parse_lt_file,
    rotation_matrix_from_axis_angle,
    rotation_matrix_align_vectors,
    rotation_matrix_to_axis_angle,
    random_rotation_matrix,
    _parse_atom_line,
    _determine_monomer_type,
)

from AutoPoly.mc.placement import (
    PolymerPlacement,
    MoleculePlacement,
    MolecularPlacementMC,
)


class TestCollisionDetector:
    """Tests for CollisionDetector class."""

    def test_initialization(self):
        """Test CollisionDetector initialization."""
        bounds = ((-50, 50), (-50, 50), (-50, 50))
        detector = CollisionDetector(bounds, cell_size=5.0)

        assert detector.cell_size == 5.0
        assert detector.get_sphere_count() == 0
        assert detector._xmin == -50
        assert detector._xmax == 50

    def test_add_monomer(self):
        """Test adding monomers to the detector."""
        bounds = ((-50, 50), (-50, 50), (-50, 50))
        detector = CollisionDetector(bounds, cell_size=5.0)

        detector.add_monomer(0, np.array([0, 0, 0]), 2.0)
        detector.add_monomer(1, np.array([10, 0, 0]), 2.0)

        assert detector.get_sphere_count() == 2

    def test_remove_monomer(self):
        """Test removing monomers from the detector."""
        bounds = ((-50, 50), (-50, 50), (-50, 50))
        detector = CollisionDetector(bounds, cell_size=5.0)

        detector.add_monomer(0, np.array([0, 0, 0]), 2.0)
        detector.add_monomer(1, np.array([10, 0, 0]), 2.0)

        detector.remove_monomer(0)
        assert detector.get_sphere_count() == 1

    def test_check_collision_true(self):
        """Test collision detection when spheres overlap."""
        bounds = ((-50, 50), (-50, 50), (-50, 50))
        detector = CollisionDetector(bounds, cell_size=5.0)

        detector.add_monomer(0, np.array([0, 0, 0]), 2.0)

        # This sphere overlaps with the first
        result = detector.check_collision(np.array([3, 0, 0]), 2.0)
        assert result is True

    def test_check_collision_false(self):
        """Test collision detection when spheres don't overlap."""
        bounds = ((-50, 50), (-50, 50), (-50, 50))
        detector = CollisionDetector(bounds, cell_size=5.0)

        detector.add_monomer(0, np.array([0, 0, 0]), 2.0)

        # This sphere doesn't overlap
        result = detector.check_collision(np.array([10, 0, 0]), 2.0)
        assert result is False

    def test_check_collision_with_exclusion(self):
        """Test collision detection with excluded IDs."""
        bounds = ((-50, 50), (-50, 50), (-50, 50))
        detector = CollisionDetector(bounds, cell_size=5.0)

        detector.add_monomer(0, np.array([0, 0, 0]), 2.0)

        # This would collide, but we exclude ID 0
        result = detector.check_collision(
            np.array([3, 0, 0]), 2.0, exclude_ids={0}
        )
        assert result is False

    def test_check_bounds(self):
        """Test boundary checking."""
        bounds = ((-10, 10), (-10, 10), (-10, 10))
        detector = CollisionDetector(bounds, cell_size=5.0)

        # Inside bounds
        assert detector.check_bounds(np.array([0, 0, 0]), 2.0) is True

        # Outside bounds (radius extends beyond)
        assert detector.check_bounds(np.array([9, 0, 0]), 2.0) is False

    def test_estimate_radius(self):
        """Test radius estimation from connection atoms."""
        left = np.array([0, 0, 0])
        right = np.array([3, 0, 0])  # Distance = 3

        radius = CollisionDetector.estimate_radius(left, right, buffer=1.0)

        # Expected: 3/2 + 1 = 2.5
        assert np.isclose(radius, 2.5)

    def test_estimate_radius_from_coords(self):
        """Test radius estimation from atom coordinates."""
        # Create a simple cube of atoms
        coords = np.array([
            [-1, -1, -1],
            [1, -1, -1],
            [-1, 1, -1],
            [1, 1, -1],
            [-1, -1, 1],
            [1, -1, 1],
            [-1, 1, 1],
            [1, 1, 1],
        ], dtype=float)

        radius = CollisionDetector.estimate_radius_from_coords(coords, buffer=0.5)

        # Max distance from center (0,0,0) to corner is sqrt(3) ~ 1.732
        expected = np.sqrt(3) + 0.5
        assert np.isclose(radius, expected, rtol=0.01)

    def test_clear(self):
        """Test clearing all spheres."""
        bounds = ((-50, 50), (-50, 50), (-50, 50))
        detector = CollisionDetector(bounds, cell_size=5.0)

        detector.add_monomer(0, np.array([0, 0, 0]), 2.0)
        detector.add_monomer(1, np.array([10, 0, 0]), 2.0)

        detector.clear()
        assert detector.get_sphere_count() == 0


class TestCalculateBoxSize:
    """Tests for calculate_box_size function."""

    def test_basic_calculation(self):
        """Test basic box size calculation."""
        # 1000 monomers at density 1.0 -> volume = 1000 -> box = 10
        box_size = calculate_box_size(1000, monomer_density=1.0)
        assert np.isclose(box_size, 10.0)

    def test_low_density(self):
        """Test with low density."""
        # 100 monomers at density 0.1 -> volume = 1000 -> box = 10
        box_size = calculate_box_size(100, monomer_density=0.1)
        assert np.isclose(box_size, 10.0)

    def test_invalid_density(self):
        """Test with invalid density."""
        with pytest.raises(ValueError):
            calculate_box_size(100, monomer_density=0)

    def test_invalid_monomers(self):
        """Test with invalid monomer count."""
        with pytest.raises(ValueError):
            calculate_box_size(0, monomer_density=0.1)


class TestRotationFunctions:
    """Tests for rotation matrix functions."""

    def test_rotation_matrix_from_axis_angle(self):
        """Test creating rotation matrix from axis-angle."""
        # 90 degree rotation around z-axis
        axis = np.array([0, 0, 1])
        angle = np.pi / 2

        R = rotation_matrix_from_axis_angle(axis, angle)

        # Should rotate [1,0,0] to [0,1,0]
        v = np.array([1, 0, 0])
        v_rotated = R @ v

        assert np.allclose(v_rotated, [0, 1, 0], atol=1e-10)

    def test_rotation_matrix_align_vectors(self):
        """Test aligning two vectors."""
        v1 = np.array([1, 0, 0])
        v2 = np.array([0, 1, 0])

        R = rotation_matrix_align_vectors(v1, v2)
        v1_rotated = R @ v1

        # Result should be parallel to v2
        assert np.allclose(v1_rotated, v2, atol=1e-10)

    def test_rotation_matrix_align_same_vectors(self):
        """Test aligning a vector to itself."""
        v1 = np.array([1, 0, 0])
        v2 = np.array([1, 0, 0])

        R = rotation_matrix_align_vectors(v1, v2)

        # Should be identity
        assert np.allclose(R, np.eye(3))

    def test_rotation_matrix_align_opposite_vectors(self):
        """Test aligning opposite vectors."""
        v1 = np.array([1, 0, 0])
        v2 = np.array([-1, 0, 0])

        R = rotation_matrix_align_vectors(v1, v2)
        v1_rotated = R @ v1

        # Result should be parallel to v2
        assert np.allclose(v1_rotated, v2, atol=1e-10)

    def test_rotation_matrix_to_axis_angle(self):
        """Test converting rotation matrix to axis-angle."""
        # Create a known rotation
        axis = np.array([0, 0, 1])
        angle_rad = np.pi / 3  # 60 degrees
        R = rotation_matrix_from_axis_angle(axis, angle_rad)

        angle_deg, ax, ay, az = rotation_matrix_to_axis_angle(R)

        assert np.isclose(angle_deg, 60.0, atol=0.1)
        assert np.isclose(abs(az), 1.0, atol=0.01)  # Axis should be [0,0,1] or [0,0,-1]

    def test_random_rotation_matrix(self):
        """Test random rotation matrix generation."""
        R = random_rotation_matrix()

        # Check orthogonality: R @ R.T = I
        assert np.allclose(R @ R.T, np.eye(3), atol=1e-10)

        # Check determinant = 1 (proper rotation)
        assert np.isclose(np.linalg.det(R), 1.0, atol=1e-10)


class TestParseAtomLine:
    """Tests for _parse_atom_line function."""

    def test_parse_valid_line(self):
        """Test parsing a valid atom line."""
        line = "$atom:C1 $mol:... @atom:81 -0.1200    -0.754   0.000   -0.000"
        atom = _parse_atom_line(line)

        assert atom is not None
        assert atom.atom_id == "C1"
        assert atom.element == "C"
        assert atom.atom_type == "81"
        assert np.isclose(atom.charge, -0.12)
        assert np.allclose(atom.coords, [-0.754, 0.0, 0.0])

    def test_parse_hydrogen_line(self):
        """Test parsing a hydrogen atom line."""
        line = "$atom:H3 $mol:... @atom:85 0.0600    -1.093   -1.061   -0.070"
        atom = _parse_atom_line(line)

        assert atom is not None
        assert atom.atom_id == "H3"
        assert atom.element == "H"

    def test_parse_invalid_line(self):
        """Test parsing an invalid line."""
        line = "invalid line"
        atom = _parse_atom_line(line)

        assert atom is None


class TestDetermineMonomerType:
    """Tests for _determine_monomer_type function."""

    def test_first_monomer(self):
        """Test detection of first (left-end) monomer."""
        assert _determine_monomer_type("monomer_0_0le") == "first"
        assert _determine_monomer_type("monomer_0_0le_T1") == "first"

    def test_last_monomer(self):
        """Test detection of last (right-end) monomer."""
        assert _determine_monomer_type("monomer_0_49re") == "last"
        assert _determine_monomer_type("monomer_0_49re_T1") == "last"

    def test_middle_monomer(self):
        """Test detection of middle (internal) monomer."""
        assert _determine_monomer_type("monomer_0_1i") == "middle"
        assert _determine_monomer_type("monomer_0_1i_T1") == "middle"

    def test_default_middle(self):
        """Test default to middle for unknown patterns."""
        assert _determine_monomer_type("monomer_unknown") == "middle"


class TestMolecularPlacementMC:
    """Tests for MolecularPlacementMC class."""

    def test_initialization(self):
        """Test MolecularPlacementMC initialization."""
        bounds = ((-50, 50), (-50, 50), (-50, 50))
        placer = MolecularPlacementMC(bounds, max_attempts=1000)

        assert placer.max_attempts == 1000
        assert placer.collision_detector is not None

    def test_random_position(self):
        """Test random position generation."""
        bounds = ((-50, 50), (-50, 50), (-50, 50))
        placer = MolecularPlacementMC(bounds)

        pos = placer.random_position()

        assert -50 <= pos[0] <= 50
        assert -50 <= pos[1] <= 50
        assert -50 <= pos[2] <= 50

    def test_random_position_with_margin(self):
        """Test random position generation with margin."""
        bounds = ((-50, 50), (-50, 50), (-50, 50))
        placer = MolecularPlacementMC(bounds)

        pos = placer.random_position(margin=10.0)

        assert -40 <= pos[0] <= 40
        assert -40 <= pos[1] <= 40
        assert -40 <= pos[2] <= 40

    def test_random_orientation(self):
        """Test random orientation generation."""
        bounds = ((-50, 50), (-50, 50), (-50, 50))
        placer = MolecularPlacementMC(bounds)

        R, axis_angle = placer.random_orientation()

        # Check it's a valid rotation matrix
        assert np.allclose(R @ R.T, np.eye(3), atol=1e-10)
        assert np.isclose(np.linalg.det(R), 1.0, atol=1e-10)

    def test_place_molecule(self):
        """Test placing a single molecule."""
        bounds = ((-100, 100), (-100, 100), (-100, 100))
        placer = MolecularPlacementMC(bounds, max_attempts=1000)

        placement = placer.place_molecule("water", radius=2.0)

        assert placement is not None
        assert placement.molecule_name == "water"
        assert placement.radius == 2.0

    def test_place_polymer(self):
        """Test placing a single polymer."""
        bounds = ((-100, 100), (-100, 100), (-100, 100))
        placer = MolecularPlacementMC(bounds, max_attempts=1000)

        placement = placer.place_polymer("poly_1", radius=10.0)

        assert placement is not None
        assert placement.poly_name == "poly_1"
        assert placement.radius == 10.0

    def test_place_multiple_molecules(self):
        """Test placing multiple molecules."""
        bounds = ((-100, 100), (-100, 100), (-100, 100))
        placer = MolecularPlacementMC(bounds, max_attempts=1000)

        specs = [
            {"molecule_name": "water", "radius": 2.0},
            {"molecule_name": "water", "radius": 2.0},
            {"molecule_name": "water", "radius": 2.0},
        ]

        placements = placer.place_all_molecules(specs)

        assert len(placements) == 3

    def test_generate_molecule_lt_commands(self):
        """Test generating moltemplate commands."""
        bounds = ((-100, 100), (-100, 100), (-100, 100))
        placer = MolecularPlacementMC(bounds, max_attempts=1000)

        placement = placer.place_molecule("water", radius=2.0)
        commands = placer.generate_molecule_lt_commands([placement])

        assert len(commands) == 1
        assert "water" in commands[0]
        assert ".move(" in commands[0]

    def test_placement_stats(self):
        """Test getting placement statistics."""
        bounds = ((-100, 100), (-100, 100), (-100, 100))
        placer = MolecularPlacementMC(bounds, max_attempts=1000)

        placer.place_molecule("water", radius=2.0)
        placer.place_polymer("poly_1", radius=10.0)

        stats = placer.get_placement_stats()

        assert stats['molecules'] == 1
        assert stats['polymers'] == 1

    def test_clear(self):
        """Test clearing placements."""
        bounds = ((-100, 100), (-100, 100), (-100, 100))
        placer = MolecularPlacementMC(bounds, max_attempts=1000)

        placer.place_molecule("water", radius=2.0)
        placer.clear()

        stats = placer.get_placement_stats()
        assert stats['total_spheres'] == 0


class TestMonomerTemplate:
    """Tests for MonomerTemplate dataclass."""

    def test_get_center(self):
        """Test calculating geometric center."""
        atoms = [
            AtomData("C1", "C", np.array([0, 0, 0]), "81", -0.12),
            AtomData("C2", "C", np.array([2, 0, 0]), "81", -0.12),
        ]

        template = MonomerTemplate(
            lt_file="test.lt",
            monomer_name="test",
            monomer_type="middle",
            atoms=atoms,
            left_conn_coords=np.array([0, 0, 0]),
            right_conn_coords=np.array([2, 0, 0]),
            left_conn_id="C1",
            right_conn_id="C2"
        )

        center = template.get_center()
        assert np.allclose(center, [1, 0, 0])

    def test_bond_vector_calculation(self):
        """Test automatic bond vector calculation."""
        template = MonomerTemplate(
            lt_file="test.lt",
            monomer_name="test",
            monomer_type="middle",
            atoms=[],
            left_conn_coords=np.array([0, 0, 0]),
            right_conn_coords=np.array([1.5, 0, 0]),
            left_conn_id="C1",
            right_conn_id="C2"
        )

        assert np.allclose(template.bond_vector, [1.5, 0, 0])


class TestParseLtFile:
    """Tests for parse_lt_file function."""

    def test_parse_valid_file(self, tmp_path):
        """Test parsing a valid .lt file."""
        lt_content = '''import "oplsaa.lt"
monomer_0_1i inherits OPLSAA {
  write("Data Atoms") {
    $atom:C1 $mol:... @atom:81 -0.1200    -0.754   0.000   -0.000
    $atom:C2 $mol:... @atom:81 -0.1200    0.754   -0.000   0.000
    $atom:H3 $mol:... @atom:85 0.0600    -1.093   -1.061   -0.070
  }
}
'''
        lt_file = tmp_path / "monomer_0_1i.lt"
        lt_file.write_text(lt_content)

        template = parse_lt_file(str(lt_file))

        assert template.monomer_name == "monomer_0_1i"
        assert template.monomer_type == "middle"
        assert len(template.atoms) == 3
        assert np.allclose(template.left_conn_coords, [-0.754, 0, 0])
        assert np.allclose(template.right_conn_coords, [0.754, 0, 0])

    def test_parse_file_not_found(self):
        """Test parsing a non-existent file."""
        with pytest.raises(FileNotFoundError):
            parse_lt_file("/nonexistent/path/file.lt")

    def test_parse_file_too_few_atoms(self, tmp_path):
        """Test parsing a file with too few atoms."""
        lt_content = '''import "oplsaa.lt"
monomer_test inherits OPLSAA {
  write("Data Atoms") {
    $atom:C1 $mol:... @atom:81 -0.1200    0.0   0.000   0.000
  }
}
'''
        lt_file = tmp_path / "monomer_test.lt"
        lt_file.write_text(lt_content)

        with pytest.raises(ValueError):
            parse_lt_file(str(lt_file))


class TestChainGrowthMC:
    """Tests for ChainGrowthMC class."""

    def test_initialization(self):
        """Test ChainGrowthMC initialization."""
        bounds = ((-100, 100), (-100, 100), (-100, 100))
        detector = CollisionDetector(bounds, cell_size=5.0)
        chain_mc = ChainGrowthMC(detector, max_attempts=1000)

        assert chain_mc.max_attempts == 1000
        assert chain_mc.collision_detector is detector

    def test_load_monomer_template(self, tmp_path):
        """Test loading and caching monomer templates."""
        lt_content = '''import "oplsaa.lt"
monomer_test inherits OPLSAA {
  write("Data Atoms") {
    $atom:C1 $mol:... @atom:81 -0.1200    0.0   0.000   0.000
    $atom:C2 $mol:... @atom:81 -0.1200    1.5   0.000   0.000
  }
}
'''
        lt_file = tmp_path / "monomer_test.lt"
        lt_file.write_text(lt_content)

        bounds = ((-100, 100), (-100, 100), (-100, 100))
        detector = CollisionDetector(bounds, cell_size=5.0)
        chain_mc = ChainGrowthMC(detector, max_attempts=1000)

        template = chain_mc.load_monomer_template(str(lt_file))
        assert template.monomer_name == "monomer_test"

        # Test caching
        template2 = chain_mc.load_monomer_template(str(lt_file))
        assert template is template2  # Same object (cached)

    def test_generate_lt_commands(self, tmp_path):
        """Test generating moltemplate commands from placements."""
        bounds = ((-100, 100), (-100, 100), (-100, 100))
        detector = CollisionDetector(bounds, cell_size=5.0)
        chain_mc = ChainGrowthMC(detector, max_attempts=1000)

        # Create a mock template
        template = MonomerTemplate(
            lt_file="test.lt",
            monomer_name="monomer_test",
            monomer_type="middle",
            atoms=[],
            left_conn_coords=np.array([0, 0, 0]),
            right_conn_coords=np.array([1.5, 0, 0]),
            left_conn_id="C1",
            right_conn_id="C2"
        )

        # Create a mock placement
        placement = MonomerPlacement(
            template=template,
            position=np.array([5.0, 0.0, 0.0]),
            rotation_matrix=np.eye(3),
            rotation_axis_angle=(0.0, 1.0, 0.0, 0.0),
            monomer_index=0
        )

        commands = chain_mc.generate_lt_commands([placement])

        assert len(commands) == 1
        assert "monomer_test" in commands[0]
        assert ".move(5.0000,0.0000,0.0000)" in commands[0]
