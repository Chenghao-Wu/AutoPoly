#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Molecular Placement Module for Monte Carlo Random Placement

This module provides Monte Carlo random placement of polymers and molecules
in a simulation box with collision detection.

Key Features:
- Random position generation within box bounds
- Uniform random orientation using quaternion sampling
- Collision-aware sequential placement
- Generation of moltemplate .rot().move() commands

Created on 2026-01-29
@author: zwu
"""
import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Any
from pathlib import Path

from .collision import CollisionDetector
from .chain_growth import (
    MonomerTemplate,
    MonomerPlacement,
    parse_lt_file,
    random_rotation_matrix,
    rotation_matrix_to_axis_angle,
)


@dataclass
class PolymerPlacement:
    """
    Represents a placed polymer chain with its transformation.

    Attributes:
        polymer_id: Unique identifier for this polymer
        poly_name: Name of the polymer (e.g., "poly_1")
        position: Translation vector for the entire polymer
        rotation_matrix: 3x3 rotation matrix for the entire polymer
        rotation_axis_angle: (angle_deg, ax, ay, az) for moltemplate
        monomer_placements: List of MonomerPlacement for chain growth (optional)
        center: Center of mass of the polymer
        radius: Bounding sphere radius
    """
    polymer_id: int
    poly_name: str
    position: np.ndarray
    rotation_matrix: np.ndarray
    rotation_axis_angle: Tuple[float, float, float, float]
    monomer_placements: Optional[List[MonomerPlacement]] = None
    center: Optional[np.ndarray] = None
    radius: float = 0.0


@dataclass
class MoleculePlacement:
    """
    Represents a placed molecule with its transformation.

    Attributes:
        molecule_id: Unique identifier for this molecule
        molecule_name: Name of the molecule type (e.g., "water")
        instance_name: Instance name in moltemplate (e.g., "molecule_1")
        position: Translation vector
        rotation_matrix: 3x3 rotation matrix
        rotation_axis_angle: (angle_deg, ax, ay, az) for moltemplate
        center: Center of the molecule after placement
        radius: Collision radius
    """
    molecule_id: int
    molecule_name: str
    instance_name: str
    position: np.ndarray
    rotation_matrix: np.ndarray
    rotation_axis_angle: Tuple[float, float, float, float]
    center: Optional[np.ndarray] = None
    radius: float = 0.0


class MolecularPlacementMC:
    """
    Monte Carlo random placement of polymers and molecules.

    This class handles placing multiple polymers and molecules in a simulation
    box using random positions and orientations while avoiding collisions.

    Attributes:
        box_bounds: Simulation box boundaries
        collision_detector: CollisionDetector for overlap checking
        max_attempts: Maximum placement attempts per molecule
    """

    def __init__(
        self,
        box_bounds: Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]],
        collision_detector: Optional[CollisionDetector] = None,
        max_attempts: int = 10000
    ):
        """
        Initialize the molecular placement MC sampler.

        Args:
            box_bounds: Simulation box as ((xmin, xmax), (ymin, ymax), (zmin, zmax))
            collision_detector: Optional CollisionDetector. If None, creates one.
            max_attempts: Maximum placement attempts per molecule/polymer
        """
        self.box_bounds = box_bounds
        self.max_attempts = max_attempts

        if collision_detector is None:
            # Estimate cell size from box dimensions
            box_size = min(
                box_bounds[0][1] - box_bounds[0][0],
                box_bounds[1][1] - box_bounds[1][0],
                box_bounds[2][1] - box_bounds[2][0]
            )
            cell_size = max(5.0, box_size / 20)
            self.collision_detector = CollisionDetector(box_bounds, cell_size)
        else:
            self.collision_detector = collision_detector

        # Counters for unique IDs
        self._polymer_counter = 0
        self._molecule_counter = 0

    def random_position(self, margin: float = 0.0) -> np.ndarray:
        """
        Generate a uniform random position within the box bounds.

        Args:
            margin: Minimum distance from box walls

        Returns:
            3D position as numpy array
        """
        xmin, xmax = self.box_bounds[0]
        ymin, ymax = self.box_bounds[1]
        zmin, zmax = self.box_bounds[2]

        x = np.random.uniform(xmin + margin, xmax - margin)
        y = np.random.uniform(ymin + margin, ymax - margin)
        z = np.random.uniform(zmin + margin, zmax - margin)

        return np.array([x, y, z])

    def random_orientation(self) -> Tuple[np.ndarray, Tuple[float, float, float, float]]:
        """
        Generate a uniform random orientation.

        Uses quaternion sampling for uniform distribution on SO(3).

        Returns:
            Tuple of (3x3 rotation matrix, (angle_deg, ax, ay, az))
        """
        R = random_rotation_matrix()
        axis_angle = rotation_matrix_to_axis_angle(R)
        return R, axis_angle

    @staticmethod
    def estimate_polymer_radius(
        monomer_placements: List[MonomerPlacement]
    ) -> Tuple[np.ndarray, float]:
        """
        Estimate bounding sphere for a polymer chain.

        Computes the center of mass and maximum distance from center
        to any atom in the chain.

        Args:
            monomer_placements: List of MonomerPlacement from chain growth

        Returns:
            Tuple of (center position, bounding sphere radius)
        """
        if not monomer_placements:
            return np.zeros(3), 0.0

        # Collect all atom coordinates
        all_coords = []
        for placement in monomer_placements:
            all_coords.append(placement.world_coords)

        all_coords = np.vstack(all_coords)
        center = np.mean(all_coords, axis=0)

        # Find maximum distance from center
        distances = np.linalg.norm(all_coords - center, axis=1)
        radius = np.max(distances)

        return center, radius

    @staticmethod
    def estimate_polymer_radius_from_templates(
        lt_files: List[str],
        spacing: float = 3.5
    ) -> float:
        """
        Estimate polymer radius from template files and spacing.

        Uses a simple linear model: radius = N * spacing / 2

        Args:
            lt_files: List of .lt file paths
            spacing: Expected monomer spacing in Angstroms

        Returns:
            Estimated radius in Angstroms
        """
        n_monomers = len(lt_files)
        # Approximate as extended chain divided by 2
        return n_monomers * spacing / 2 + 2.0  # Extra buffer

    def place_polymer(
        self,
        poly_name: str,
        radius: float,
        center: Optional[np.ndarray] = None
    ) -> Optional[PolymerPlacement]:
        """
        Place a polymer with random position and orientation.

        Args:
            poly_name: Name of the polymer (e.g., "poly_1")
            radius: Bounding sphere radius for collision detection
            center: Optional center offset within polymer coordinates

        Returns:
            PolymerPlacement if successful, None if failed after max_attempts
        """
        if center is None:
            center = np.zeros(3)

        for attempt in range(self.max_attempts):
            # Random position with margin for polymer radius
            position = self.random_position(margin=radius)

            # Random orientation
            R, axis_angle = self.random_orientation()

            # World center after transformation
            world_center = R @ center + position

            # Check collision
            if not self.collision_detector.check_collision(world_center, radius):
                # Add to collision detector
                polymer_id = self._polymer_counter
                self._polymer_counter += 1

                self.collision_detector.add_monomer(
                    polymer_id + 100000,  # Offset to avoid collision with monomer IDs
                    world_center,
                    radius
                )

                return PolymerPlacement(
                    polymer_id=polymer_id,
                    poly_name=poly_name,
                    position=position,
                    rotation_matrix=R,
                    rotation_axis_angle=axis_angle,
                    center=world_center,
                    radius=radius
                )

        return None

    def place_molecule(
        self,
        molecule_name: str,
        radius: float,
        instance_name: Optional[str] = None
    ) -> Optional[MoleculePlacement]:
        """
        Place a molecule with random position and orientation.

        Args:
            molecule_name: Type name of the molecule (e.g., "water")
            radius: Collision radius for the molecule
            instance_name: Optional instance name. If None, auto-generated.

        Returns:
            MoleculePlacement if successful, None if failed after max_attempts
        """
        for attempt in range(self.max_attempts):
            # Random position with margin
            position = self.random_position(margin=radius)

            # Random orientation
            R, axis_angle = self.random_orientation()

            # Check collision
            if not self.collision_detector.check_collision(position, radius):
                # Add to collision detector
                molecule_id = self._molecule_counter
                self._molecule_counter += 1

                if instance_name is None:
                    instance_name = f"molecule_{molecule_id + 1}"

                self.collision_detector.add_monomer(
                    molecule_id + 200000,  # Offset to avoid collision with other IDs
                    position,
                    radius
                )

                return MoleculePlacement(
                    molecule_id=molecule_id,
                    molecule_name=molecule_name,
                    instance_name=instance_name,
                    position=position,
                    rotation_matrix=R,
                    rotation_axis_angle=axis_angle,
                    center=position,
                    radius=radius
                )

        return None

    def place_all_polymers(
        self,
        polymer_specs: List[Dict[str, Any]]
    ) -> List[PolymerPlacement]:
        """
        Place multiple polymers sequentially with collision avoidance.

        Args:
            polymer_specs: List of dicts with keys:
                - poly_name: Name of the polymer (e.g., "poly_1")
                - radius: Bounding sphere radius
                - center: Optional center offset (default: origin)

        Returns:
            List of PolymerPlacement objects

        Raises:
            RuntimeError: If any polymer fails to place
        """
        placements = []

        for spec in polymer_specs:
            placement = self.place_polymer(
                poly_name=spec['poly_name'],
                radius=spec['radius'],
                center=spec.get('center')
            )

            if placement is None:
                raise RuntimeError(
                    f"Failed to place polymer {spec['poly_name']} "
                    f"after {self.max_attempts} attempts"
                )

            placements.append(placement)

        return placements

    def place_all_molecules(
        self,
        molecule_specs: List[Dict[str, Any]]
    ) -> List[MoleculePlacement]:
        """
        Place multiple molecules sequentially with collision avoidance.

        Args:
            molecule_specs: List of dicts with keys:
                - molecule_name: Type name of the molecule
                - radius: Collision radius
                - instance_name: Optional instance name

        Returns:
            List of MoleculePlacement objects

        Raises:
            RuntimeError: If any molecule fails to place
        """
        placements = []

        for spec in molecule_specs:
            placement = self.place_molecule(
                molecule_name=spec['molecule_name'],
                radius=spec['radius'],
                instance_name=spec.get('instance_name')
            )

            if placement is None:
                raise RuntimeError(
                    f"Failed to place molecule {spec['molecule_name']} "
                    f"after {self.max_attempts} attempts"
                )

            placements.append(placement)

        return placements

    def generate_polymer_lt_commands(
        self,
        placements: List[PolymerPlacement],
        use_rotation: bool = True
    ) -> List[str]:
        """
        Generate moltemplate instantiation commands for polymers.

        Output format:
        polymer_1 = new poly_1.rot(angle,ax,ay,az).move(x,y,z)

        Args:
            placements: List of PolymerPlacement
            use_rotation: Whether to include .rot() commands

        Returns:
            List of moltemplate command strings
        """
        commands = []

        for placement in placements:
            poly_name = placement.poly_name
            pos = placement.position

            # Derive instance name from poly_name (e.g., "poly_1" -> "polymer_1")
            # This ensures correct naming even when some placements fail
            if poly_name.startswith("poly_"):
                poly_num = poly_name[5:]  # Extract number after "poly_"
                instance_name = f"polymer_{poly_num}"
            else:
                instance_name = f"polymer_{placement.polymer_id + 1}"

            if use_rotation and not np.allclose(placement.rotation_matrix, np.eye(3)):
                angle, ax, ay, az = placement.rotation_axis_angle
                cmd = (
                    f"{instance_name} = new {poly_name}"
                    f".rot({angle:.4f},{ax:.4f},{ay:.4f},{az:.4f})"
                    f".move({pos[0]:.4f},{pos[1]:.4f},{pos[2]:.4f})"
                )
            else:
                cmd = (
                    f"{instance_name} = new {poly_name}"
                    f".move({pos[0]:.4f},{pos[1]:.4f},{pos[2]:.4f})"
                )

            commands.append(cmd)

        return commands

    def generate_molecule_lt_commands(
        self,
        placements: List[MoleculePlacement],
        use_rotation: bool = True
    ) -> List[str]:
        """
        Generate moltemplate instantiation commands for molecules.

        Output format:
        molecule_1 = new water.rot(angle,ax,ay,az).move(x,y,z)

        Args:
            placements: List of MoleculePlacement
            use_rotation: Whether to include .rot() commands

        Returns:
            List of moltemplate command strings
        """
        commands = []

        for placement in placements:
            mol_name = placement.molecule_name
            instance = placement.instance_name
            pos = placement.position

            if use_rotation and not np.allclose(placement.rotation_matrix, np.eye(3)):
                angle, ax, ay, az = placement.rotation_axis_angle
                cmd = (
                    f"{instance} = new {mol_name}"
                    f".rot({angle:.4f},{ax:.4f},{ay:.4f},{az:.4f})"
                    f".move({pos[0]:.4f},{pos[1]:.4f},{pos[2]:.4f})"
                )
            else:
                cmd = (
                    f"{instance} = new {mol_name}"
                    f".move({pos[0]:.4f},{pos[1]:.4f},{pos[2]:.4f})"
                )

            commands.append(cmd)

        return commands

    def get_placement_stats(self) -> Dict[str, int]:
        """
        Get statistics about placed entities.

        Returns:
            Dict with counts of placed polymers and molecules
        """
        return {
            'polymers': self._polymer_counter,
            'molecules': self._molecule_counter,
            'total_spheres': self.collision_detector.get_sphere_count()
        }

    def clear(self) -> None:
        """Reset the placement state."""
        self.collision_detector.clear()
        self._polymer_counter = 0
        self._molecule_counter = 0
