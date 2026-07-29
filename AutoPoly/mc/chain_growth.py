#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chain Growth Monte Carlo Module for Self-Avoiding Random Walk

This module implements self-avoiding random walk (SAW) for building polymer chains
monomer-by-monomer with collision detection. The chain grows by aligning each new
monomer's left connection point to the previous monomer's right connection point.

Key Features:
- Parse monomer templates from .lt files
- 3D transformation (rotation + translation) for monomer alignment
- Self-avoiding random walk with dihedral angle sampling
- Generation of moltemplate .rot().move() commands

Connection Point Convention (from .lt files):
- First atom in "Data Atoms": Left connection point (C1)
- Second atom in "Data Atoms": Right connection point (C2)

Created on 2026-01-29
@author: zwu
"""
import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Any
from pathlib import Path
import re
import logging

from .collision import CollisionDetector

# Set up module logger
logger = logging.getLogger(__name__)


@dataclass
class AtomData:
    """
    Data for a single atom parsed from .lt file.

    Attributes:
        atom_id: Atom identifier (e.g., "C1", "H3")
        element: Element symbol (e.g., "C", "H")
        coords: 3D coordinates as numpy array
        atom_type: Force field atom type
        charge: Partial charge
    """
    atom_id: str
    element: str
    coords: np.ndarray
    atom_type: str
    charge: float


@dataclass
class MonomerTemplate:
    """
    Template for a monomer parsed from .lt file.

    Attributes:
        lt_file: Path to the .lt file
        monomer_name: Name of the monomer (e.g., "monomer_0_1i")
        monomer_type: Type of monomer ("first", "middle", "last", or "ring")
        atoms: List of all atoms with local coordinates
        left_conn_coords: Coordinates of left connection atom (first atom, C1)
        right_conn_coords: Coordinates of right connection atom (second atom, C2)
        left_conn_id: Atom ID of left connection (e.g., "C1") or None for first monomer
        right_conn_id: Atom ID of right connection (e.g., "C2") or None for last monomer
        bond_vector: Vector from left to right connection (for middle monomers)
    """
    lt_file: str
    monomer_name: str
    monomer_type: str  # "first", "middle", "last", "ring"
    atoms: List[AtomData]
    left_conn_coords: Optional[np.ndarray]
    right_conn_coords: Optional[np.ndarray]
    left_conn_id: Optional[str]
    right_conn_id: Optional[str]
    bond_vector: Optional[np.ndarray] = None

    def __post_init__(self):
        """Calculate bond vector if both connection points exist."""
        if self.left_conn_coords is not None and self.right_conn_coords is not None:
            self.bond_vector = self.right_conn_coords - self.left_conn_coords

    def get_center(self) -> np.ndarray:
        """Calculate the geometric center of all atoms."""
        if not self.atoms:
            return np.zeros(3)
        coords = np.array([a.coords for a in self.atoms])
        return np.mean(coords, axis=0)

    def get_all_coords(self) -> np.ndarray:
        """Get array of all atom coordinates."""
        return np.array([a.coords for a in self.atoms])


@dataclass
class MonomerPlacement:
    """
    Represents a placed monomer with its transformation.

    Attributes:
        template: The monomer template used
        position: Translation vector applied
        rotation_matrix: 3x3 rotation matrix applied
        rotation_axis_angle: (angle_deg, ax, ay, az) for moltemplate .rot() command
        monomer_index: Index in the chain
        world_coords: Atom coordinates after transformation
        world_left_conn: Left connection point after transformation
        world_right_conn: Right connection point after transformation
    """
    template: MonomerTemplate
    position: np.ndarray
    rotation_matrix: np.ndarray
    rotation_axis_angle: Tuple[float, float, float, float]
    monomer_index: int
    world_coords: np.ndarray = field(default_factory=lambda: np.array([]))
    world_left_conn: Optional[np.ndarray] = None
    world_right_conn: Optional[np.ndarray] = None


def parse_lt_file(lt_file: str) -> MonomerTemplate:
    """
    Parse a monomer .lt file to extract atom coordinates and connection points.

    The .lt file format has atoms in the "Data Atoms" block with the format:
    $atom:ID $mol:... @atom:TYPE CHARGE X Y Z

    Connection point convention:
    - First atom (C1): Left connection point
    - Second atom (C2): Right connection point

    Args:
        lt_file: Path to the .lt file

    Returns:
        MonomerTemplate with parsed data

    Raises:
        FileNotFoundError: If the .lt file doesn't exist
        ValueError: If the file format is invalid
    """
    lt_path = Path(lt_file)
    if not lt_path.exists():
        raise FileNotFoundError(f"Monomer .lt file not found: {lt_file}")

    # Extract monomer name from file
    monomer_name = lt_path.stem

    atoms = []
    in_atoms_block = False

    with open(lt_path, 'r') as f:
        for line in f:
            line_stripped = line.strip()

            # Detect start of Data Atoms block
            if line_stripped == 'write("Data Atoms") {':
                in_atoms_block = True
                continue
            elif line_stripped == '}':
                if in_atoms_block:
                    in_atoms_block = False
                continue

            if in_atoms_block and line_stripped:
                atom = _parse_atom_line(line_stripped)
                if atom:
                    atoms.append(atom)

    if len(atoms) < 2:
        raise ValueError(f"Expected at least 2 atoms in {lt_file}, found {len(atoms)}")

    # Determine monomer type from filename
    monomer_type = _determine_monomer_type(monomer_name)

    # Connection points are first two atoms
    left_conn_coords = atoms[0].coords if monomer_type != "first" else atoms[0].coords
    right_conn_coords = atoms[1].coords if monomer_type != "last" else atoms[1].coords

    left_conn_id = atoms[0].atom_id
    right_conn_id = atoms[1].atom_id

    # For first/last monomers, one connection is "virtual" (used only for chain bonding)
    if monomer_type == "first":
        # First monomer has no left connection to previous
        left_conn_id = None
        left_conn_coords = None
    elif monomer_type == "last":
        # Last monomer has no right connection to next
        right_conn_id = None
        right_conn_coords = None

    return MonomerTemplate(
        lt_file=str(lt_path),
        monomer_name=monomer_name,
        monomer_type=monomer_type,
        atoms=atoms,
        left_conn_coords=left_conn_coords,
        right_conn_coords=right_conn_coords,
        left_conn_id=left_conn_id,
        right_conn_id=right_conn_id
    )


def _parse_atom_line(line: str) -> Optional[AtomData]:
    """
    Parse a single atom line from the Data Atoms block.

    Format: $atom:ID $mol:... @atom:TYPE CHARGE X Y Z

    Args:
        line: The line to parse

    Returns:
        AtomData if successfully parsed, None otherwise
    """
    parts = line.split()
    if len(parts) < 7:
        return None

    # Extract atom ID (e.g., "$atom:C1" -> "C1")
    atom_id_match = re.match(r'\$atom:(\w+)', parts[0])
    if not atom_id_match:
        return None
    atom_id = atom_id_match.group(1)

    # Extract element from atom ID (e.g., "C1" -> "C", "H16" -> "H")
    element_match = re.match(r'([A-Z][a-z]?)', atom_id)
    element = element_match.group(1) if element_match else "X"

    # Extract atom type (e.g., "@atom:81" -> "81")
    atom_type_match = re.match(r'@atom:(\w+)', parts[2])
    atom_type = atom_type_match.group(1) if atom_type_match else ""

    try:
        charge = float(parts[3])
        x = float(parts[4])
        y = float(parts[5])
        z = float(parts[6])
    except (ValueError, IndexError):
        return None

    return AtomData(
        atom_id=atom_id,
        element=element,
        coords=np.array([x, y, z]),
        atom_type=atom_type,
        charge=charge
    )


def _determine_monomer_type(monomer_name: str) -> str:
    """
    Determine monomer type from filename convention.

    Conventions:
    - "*le*" or "*_0le*": Left-end (first monomer)
    - "*re*" or "*_re*": Right-end (last monomer)
    - "*i*" or "*_i*": Internal (middle monomer)
    - Otherwise: middle (default)

    Args:
        monomer_name: Name of the monomer file (without extension)

    Returns:
        Monomer type: "first", "last", "middle", or "ring"
    """
    name_lower = monomer_name.lower()

    if "le" in name_lower and ("0le" in name_lower or "_le" in name_lower or name_lower.endswith("le")):
        return "first"
    elif "re" in name_lower and ("_re" in name_lower or name_lower.endswith("re")):
        return "last"
    elif "_i" in name_lower or name_lower.endswith("i"):
        return "middle"
    else:
        return "middle"


def rotation_matrix_from_axis_angle(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    """
    Create a 3x3 rotation matrix from axis-angle representation.

    Uses Rodrigues' rotation formula.

    Args:
        axis: Unit vector defining rotation axis
        angle_rad: Rotation angle in radians

    Returns:
        3x3 rotation matrix
    """
    axis = axis / np.linalg.norm(axis)
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    t = 1 - c

    x, y, z = axis

    return np.array([
        [t*x*x + c,    t*x*y - s*z,  t*x*z + s*y],
        [t*x*y + s*z,  t*y*y + c,    t*y*z - s*x],
        [t*x*z - s*y,  t*y*z + s*x,  t*z*z + c]
    ])


def rotation_matrix_align_vectors(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """
    Create rotation matrix that aligns vector v1 to vector v2.

    Args:
        v1: Source vector
        v2: Target vector

    Returns:
        3x3 rotation matrix R such that R @ v1 is parallel to v2
    """
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)

    # Check if vectors are already aligned or opposite
    dot = np.dot(v1, v2)

    # For nearly parallel vectors, add small perturbation for numerical stability
    if dot > 0.9999:
        # Add tiny random perturbation to break symmetry and avoid numerical issues
        perturbation = np.random.randn(3) * 1e-6
        v2_perturbed = v2 + perturbation
        v2_perturbed = v2_perturbed / np.linalg.norm(v2_perturbed)
        dot = np.dot(v1, v2_perturbed)
        if dot > 0.9999:
            return np.eye(3)
        v2 = v2_perturbed
    elif dot < -0.9999:
        # Find a perpendicular axis
        perp = np.array([1, 0, 0]) if abs(v1[0]) < 0.9 else np.array([0, 1, 0])
        axis = np.cross(v1, perp)
        axis = axis / np.linalg.norm(axis)
        return rotation_matrix_from_axis_angle(axis, np.pi)

    # Rotation axis is cross product
    axis = np.cross(v1, v2)
    axis_norm = np.linalg.norm(axis)

    # Handle edge case where cross product is very small
    if axis_norm < 1e-10:
        return np.eye(3)

    axis = axis / axis_norm

    # Rotation angle from dot product
    angle = np.arccos(np.clip(dot, -1, 1))

    return rotation_matrix_from_axis_angle(axis, angle)


def rotation_matrix_to_axis_angle(R: np.ndarray) -> Tuple[float, float, float, float]:
    """
    Convert a 3x3 rotation matrix to axis-angle representation.

    Returns format suitable for moltemplate: (angle_degrees, ax, ay, az)

    Args:
        R: 3x3 rotation matrix

    Returns:
        Tuple of (angle_degrees, axis_x, axis_y, axis_z)
    """
    # Handle identity matrix
    trace = np.trace(R)
    if np.isclose(trace, 3.0):
        return (0.0, 1.0, 0.0, 0.0)

    # Handle 180 degree rotation
    if np.isclose(trace, -1.0):
        # Find the column with largest diagonal element
        diag = np.diag(R)
        i = np.argmax(diag)
        axis = np.zeros(3)
        axis[i] = 1.0
        return (180.0, axis[0], axis[1], axis[2])

    # General case
    angle = np.arccos(np.clip((trace - 1) / 2, -1, 1))

    # Axis from skew-symmetric part
    axis = np.array([
        R[2, 1] - R[1, 2],
        R[0, 2] - R[2, 0],
        R[1, 0] - R[0, 1]
    ])

    axis_norm = np.linalg.norm(axis)
    if axis_norm < 1e-10:
        return (0.0, 1.0, 0.0, 0.0)

    axis = axis / axis_norm
    angle_deg = np.degrees(angle)

    return (angle_deg, axis[0], axis[1], axis[2])


def random_rotation_matrix() -> np.ndarray:
    """
    Generate a uniformly distributed random rotation matrix.

    Uses the algorithm from Graphics Gems III for uniform random rotations
    based on quaternion sampling.

    Returns:
        3x3 rotation matrix
    """
    # Generate uniform random quaternion
    u1, u2, u3 = np.random.random(3)

    q = np.array([
        np.sqrt(1 - u1) * np.sin(2 * np.pi * u2),
        np.sqrt(1 - u1) * np.cos(2 * np.pi * u2),
        np.sqrt(u1) * np.sin(2 * np.pi * u3),
        np.sqrt(u1) * np.cos(2 * np.pi * u3)
    ])

    # Convert quaternion to rotation matrix
    q0, q1, q2, q3 = q
    return np.array([
        [1 - 2*(q2**2 + q3**2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)],
        [2*(q1*q2 + q0*q3), 1 - 2*(q1**2 + q3**2), 2*(q2*q3 - q0*q1)],
        [2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1 - 2*(q1**2 + q2**2)]
    ])


class ChainGrowthMC:
    """
    Self-Avoiding Random Walk chain growth using Monte Carlo.

    This class builds polymer chains by placing monomers one at a time,
    using collision detection to ensure self-avoiding conformations.

    The algorithm:
    1. Place first monomer at origin with random orientation
    2. For each subsequent monomer:
       a. Align left connection to previous monomer's right connection
       b. Sample random dihedral angles
       c. Check for collisions
       d. Accept if no collision, retry otherwise

    Attributes:
        collision_detector: CollisionDetector for checking overlaps
        max_attempts: Maximum placement attempts per monomer
        templates: Cache of loaded monomer templates
    """

    def __init__(
        self,
        collision_detector: CollisionDetector,
        max_attempts: int = 1000,
        bond_angle_min: float = 95.0,
        bond_angle_max: float = 150.0,
        intrachain_exclude_neighbors: int = 2,
        junction_bond_length: float = 1.54
    ):
        """
        Initialize the chain growth Monte Carlo sampler.

        Args:
            collision_detector: CollisionDetector instance for checking overlaps
            max_attempts: Maximum number of placement attempts per monomer
            bond_angle_min: Minimum bond angle in degrees (C-C-C angle between consecutive monomers)
            bond_angle_max: Maximum bond angle in degrees (C-C-C angle between consecutive monomers)
            intrachain_exclude_neighbors: Number of neighbors to exclude from intra-chain collision detection.
                                         0 = exclude only current monomer (maximum collision checking)
                                         1 = exclude i-1, i, i+1 (immediate neighbors)
                                         2 = exclude i-2, i-1, i, i+1, i+2 (default, recommended)
                                         3 = exclude i-3 through i+3
                                         Higher values allow tighter packing but must maintain zero self-intersections.
            junction_bond_length: Equilibrium bond length (Angstrom) of the
                                         inter-monomer bond formed at each junction.
                                         The incoming monomer's left connection is
                                         placed this far beyond the previous right
                                         connection (previously it was placed exactly
                                         ON it, giving zero-length junction bonds).
        """
        self.collision_detector = collision_detector
        self.max_attempts = max_attempts
        self.bond_angle_min = bond_angle_min
        self.bond_angle_max = bond_angle_max
        self.intrachain_exclude_neighbors = intrachain_exclude_neighbors
        self.junction_bond_length = junction_bond_length
        self.templates: Dict[str, MonomerTemplate] = {}

    def load_monomer_template(self, lt_file: str) -> MonomerTemplate:
        """
        Load and cache a monomer template from .lt file.

        Args:
            lt_file: Path to the .lt file

        Returns:
            MonomerTemplate with parsed data
        """
        if lt_file not in self.templates:
            self.templates[lt_file] = parse_lt_file(lt_file)
        return self.templates[lt_file]

    def align_monomer_to_connection(
        self,
        template: MonomerTemplate,
        target_position: np.ndarray,
        incoming_direction: np.ndarray,
        dihedral_angle: float = 0.0
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Align monomer so its left connection is at target position.

        Note: In the new bond angle implementation, incoming_direction already
        includes both bond angle bend and dihedral rotation. The dihedral_angle
        parameter is kept for backward compatibility with legacy mode.

        Steps:
        1. Translate monomer so left_conn is at origin
        2. Rotate to align bond_vector with incoming_direction
        3. (Legacy) Apply dihedral rotation around backbone axis
        4. Translate so left_conn is at target_position

        Args:
            template: MonomerTemplate to align
            target_position: Where left_conn should be placed
            incoming_direction: Direction the backbone is coming from
                              (already includes bond angle and dihedral)
            dihedral_angle: (Deprecated) Dihedral angle in radians

        Returns:
            Tuple of (3x3 rotation matrix, translation vector)
        """
        if template.left_conn_coords is None:
            raise ValueError("Cannot align monomer without left connection point")

        # Get the bond vector (direction from left to right connection)
        if template.bond_vector is not None:
            bond_vec = template.bond_vector
        else:
            # Use a default direction if no bond vector
            bond_vec = np.array([1.0, 0.0, 0.0])

        # Normalize incoming direction
        incoming_dir = incoming_direction / np.linalg.norm(incoming_direction)

        # Align bond_vector with incoming_direction
        # (incoming_dir already includes bond angle and dihedral from _apply_bond_angle_bend)
        R_align = rotation_matrix_align_vectors(bond_vec, incoming_dir)

        # Apply additional dihedral rotation only if provided (legacy compatibility)
        if dihedral_angle != 0.0:
            R_dihedral = rotation_matrix_from_axis_angle(incoming_dir, dihedral_angle)
            R_total = R_dihedral @ R_align
        else:
            R_total = R_align

        # Calculate translation
        # After rotation, left_conn moves to: R_total @ left_conn_coords
        rotated_left_conn = R_total @ template.left_conn_coords
        translation = target_position - rotated_left_conn

        return R_total, translation

    def sample_bond_angle(self) -> float:
        """
        Sample a bond angle randomly within the configured range.

        Returns:
            Bond angle in radians (C-C-C angle between consecutive monomers)
        """
        # Uniform random sampling within range
        angle_deg = np.random.uniform(self.bond_angle_min, self.bond_angle_max)
        return np.radians(angle_deg)

    def _apply_bond_angle_bend(
        self,
        prev_bond_dir: np.ndarray,
        bond_angle: float,
        dihedral: float
    ) -> np.ndarray:
        """
        Apply bond angle bend to previous bond direction using spherical coordinates.

        Constructs direction on a cone around prev_bond_dir:
        - Cone half-angle = bond_angle (angle between consecutive bonds)
        - Azimuthal position = dihedral

        Args:
            prev_bond_dir: Normalized direction of previous bond
            bond_angle: C-C-C bond angle in radians (e.g., 109.5° for tetrahedral)
            dihedral: Dihedral rotation angle in radians

        Returns:
            Normalized incoming direction vector for new monomer
        """
        # The angle between consecutive bond vectors equals the bond angle
        # (not the supplement)

        # Find two orthogonal vectors perpendicular to prev_bond_dir
        if abs(prev_bond_dir[2]) < 0.9:
            u = np.cross(prev_bond_dir, np.array([0, 0, 1]))
        else:
            u = np.cross(prev_bond_dir, np.array([1, 0, 0]))
        u = u / np.linalg.norm(u)

        v = np.cross(prev_bond_dir, u)
        v = v / np.linalg.norm(v)

        # Construct direction on cone using spherical coordinates
        incoming_dir = (
            np.cos(bond_angle) * prev_bond_dir +
            np.sin(bond_angle) * (np.cos(dihedral) * u + np.sin(dihedral) * v)
        )

        return incoming_dir / np.linalg.norm(incoming_dir)

    def _transform_coords(
        self,
        coords: np.ndarray,
        rotation: np.ndarray,
        translation: np.ndarray
    ) -> np.ndarray:
        """Apply rotation and translation to coordinates."""
        return (rotation @ coords.T).T + translation

    def _estimate_monomer_radius(self, template: MonomerTemplate) -> float:
        """Estimate collision radius using backbone atoms only (C, O, N, S).

        Using backbone-only atoms produces a smaller radius (~1.0-1.5 Å) than
        the full bounding sphere (~3-5 Å). This keeps the radius below the
        monomer connection distance (~2.35 Å for PEO), so the auto-calibration
        condition ``collision_diameter > conn_dist`` stays false and intra-chain
        collision detection remains active.
        """
        backbone_elements = {'C', 'O', 'N', 'S'}
        backbone_coords = [a.coords for a in template.atoms if a.element in backbone_elements]

        if not backbone_coords:
            # Fallback to all heavy atoms
            backbone_coords = [a.coords for a in template.atoms if a.element != 'H']

        if not backbone_coords:
            return 1.5  # Minimum default

        coords = np.array(backbone_coords)
        center = np.mean(coords, axis=0)
        distances = np.linalg.norm(coords - center, axis=1)
        max_distance = np.max(distances) if len(distances) > 0 else 0

        # Small buffer (0.3 Å) — backbone-only radius should be smaller than
        # the connection distance to keep intra-chain collision active
        radius = max_distance + 0.3
        logger.debug(f"Monomer {template.monomer_name}: backbone-only radius = {radius:.2f} Å")
        return radius

    def grow_chain(
        self,
        monomer_lt_files: List[str],
        chain_id: int = 0,
        start_position: Optional[np.ndarray] = None
    ) -> List[MonomerPlacement]:
        """
        Build a polymer chain using self-avoiding random walk.

        Algorithm:
        1. Place first monomer at start_position with random 3D orientation
        2. For each subsequent monomer:
           - Get right_conn of previous monomer as target for left_conn
           - Sample random dihedral angles
           - Check collision with all placed atoms
           - Accept if no collision, retry otherwise

        Args:
            monomer_lt_files: List of .lt file paths for each monomer in sequence
            chain_id: Unique identifier for this chain
            start_position: Starting position (default: origin)

        Returns:
            List of MonomerPlacement objects for each placed monomer

        Raises:
            RuntimeError: If chain growth fails after max_attempts
        """
        if start_position is None:
            start_position = np.zeros(3)

        placements = []

        for i, lt_file in enumerate(monomer_lt_files):
            template = self.load_monomer_template(lt_file)
            monomer_radius = self._estimate_monomer_radius(template)

            placement = None

            if i == 0:
                # First monomer: random orientation at start position
                placement = self._place_first_monomer(
                    template, start_position, chain_id, monomer_radius
                )
            else:
                # Subsequent monomers: align to previous
                prev_placement = placements[-1]
                placement = self._place_subsequent_monomer(
                    template, prev_placement, i, chain_id, monomer_radius
                )

            if placement is None:
                raise RuntimeError(
                    f"Failed to place monomer {i} ({template.monomer_name}) "
                    f"after {self.max_attempts} attempts"
                )

            placements.append(placement)

        return placements

    def _place_first_monomer(
        self,
        template: MonomerTemplate,
        position: np.ndarray,
        chain_id: int,
        radius: float
    ) -> Optional[MonomerPlacement]:
        """Place the first monomer with random orientation."""
        for _ in range(self.max_attempts):
            # Random orientation
            R = random_rotation_matrix()

            # Transform coordinates
            coords = template.get_all_coords()
            center = template.get_center()

            # Rotate around center, then translate to position
            rotated_coords = (R @ (coords - center).T).T + position

            # Calculate world connection points
            world_left_conn = None
            world_right_conn = None
            if template.left_conn_coords is not None:
                world_left_conn = R @ (template.left_conn_coords - center) + position
            if template.right_conn_coords is not None:
                world_right_conn = R @ (template.right_conn_coords - center) + position

            # Check collision
            monomer_center = np.mean(rotated_coords, axis=0)
            if not self.collision_detector.check_collision(monomer_center, radius):
                # Add to collision detector
                mon_id = chain_id * 10000  # Unique ID scheme
                self.collision_detector.add_monomer(mon_id, monomer_center, radius)

                axis_angle = rotation_matrix_to_axis_angle(R)

                return MonomerPlacement(
                    template=template,
                    position=position - R @ center,  # Offset for moltemplate
                    rotation_matrix=R,
                    rotation_axis_angle=axis_angle,
                    monomer_index=0,
                    world_coords=rotated_coords,
                    world_left_conn=world_left_conn,
                    world_right_conn=world_right_conn
                )

        return None

    def _place_subsequent_monomer(
        self,
        template: MonomerTemplate,
        prev_placement: MonomerPlacement,
        monomer_index: int,
        chain_id: int,
        radius: float
    ) -> Optional[MonomerPlacement]:
        """Place a subsequent monomer aligned to the previous one."""
        if prev_placement.world_right_conn is None:
            raise ValueError("Previous monomer has no right connection point")

        # Target position for this monomer's left connection
        target_pos = prev_placement.world_right_conn

        # Previous bond direction (to apply bend angle to)
        if prev_placement.world_left_conn is not None:
            prev_bond_dir = prev_placement.world_right_conn - prev_placement.world_left_conn
        else:
            prev_bond_dir = prev_placement.world_right_conn - np.mean(prev_placement.world_coords, axis=0)

        prev_bond_dir = prev_bond_dir / np.linalg.norm(prev_bond_dir)

        # Build exclusion set for selective intra-chain collision detection
        # Exclude only nearby neighbors (±intrachain_exclude_neighbors) - this prevents
        # the chain from folding back on itself while allowing realistic bond angles
        exclude_ids = set()
        for offset in range(-self.intrachain_exclude_neighbors,
                           self.intrachain_exclude_neighbors + 1):
            neighbor_idx = monomer_index + offset
            if 0 <= neighbor_idx < monomer_index:  # Only exclude already-placed monomers
                exclude_ids.add(chain_id * 10000 + neighbor_idx)

        # Always exclude current monomer
        exclude_ids.add(chain_id * 10000 + monomer_index)

        for attempt in range(self.max_attempts):
            # 1. Sample bond angle
            bond_angle = self.sample_bond_angle()

            # 2. Sample dihedral angle — uniform sampling for generic SAW.
            # Uniform dihedral + excluded volume naturally produces the correct
            # Flory exponent (ν ≈ 0.588) for self-avoiding walks.
            dihedral = np.random.uniform(0, 2 * np.pi)

            # 3. Construct incoming direction with bond angle bend
            incoming_dir = self._apply_bond_angle_bend(prev_bond_dir, bond_angle, dihedral)

            # 4. Align monomer (dihedral=0 since already applied in incoming_dir)
            #    The left connection is placed one equilibrium bond length
            #    BEYOND the previous right connection, along the growth
            #    direction — otherwise the junction atoms coincide exactly
            #    and the inter-monomer bond written into the .lt has zero length.
            junction_target = target_pos + self.junction_bond_length * incoming_dir
            R, translation = self.align_monomer_to_connection(
                template, junction_target, incoming_dir, 0.0
            )

            # Transform all coordinates
            coords = template.get_all_coords()
            world_coords = self._transform_coords(coords, R, translation)

            # Calculate world connection points
            world_left_conn = None
            world_right_conn = None
            if template.left_conn_coords is not None:
                world_left_conn = R @ template.left_conn_coords + translation
            if template.right_conn_coords is not None:
                world_right_conn = R @ template.right_conn_coords + translation

            # Check collision
            monomer_center = np.mean(world_coords, axis=0)

            # Check bounds
            if not self.collision_detector.check_bounds(monomer_center, radius):
                continue

            if not self.collision_detector.check_collision(
                monomer_center, radius, exclude_ids=exclude_ids
            ):
                # Add to collision detector
                mon_id = chain_id * 10000 + monomer_index
                self.collision_detector.add_monomer(mon_id, monomer_center, radius)

                axis_angle = rotation_matrix_to_axis_angle(R)

                # Enhanced logging with collision info
                logger.debug(
                    f"Chain {chain_id} monomer {monomer_index}: placed after {attempt + 1} attempts, "
                    f"bond_angle = {np.degrees(bond_angle):.1f}°, dihedral = {np.degrees(dihedral):.1f}°, "
                    f"excluded {len(exclude_ids)} neighbors"
                )

                return MonomerPlacement(
                    template=template,
                    position=translation,
                    rotation_matrix=R,
                    rotation_axis_angle=axis_angle,
                    monomer_index=monomer_index,
                    world_coords=world_coords,
                    world_left_conn=world_left_conn,
                    world_right_conn=world_right_conn
                )

        # Log failure for diagnostic purposes
        logger.debug(
            f"Chain {chain_id} monomer {monomer_index}: FAILED after {self.max_attempts} attempts"
        )
        return None

    def generate_lt_commands(
        self,
        placements: List[MonomerPlacement],
        use_rotation: bool = True
    ) -> List[str]:
        """
        Generate moltemplate instantiation commands from placements.

        Output format:
        monomer[0] = new MonomerA.rot(angle,ax,ay,az).move(x,y,z)

        Args:
            placements: List of MonomerPlacement from grow_chain()
            use_rotation: Whether to include .rot() commands

        Returns:
            List of moltemplate command strings
        """
        commands = []

        for placement in placements:
            monomer_name = placement.template.monomer_name
            idx = placement.monomer_index
            pos = placement.position

            if use_rotation and not np.allclose(placement.rotation_matrix, np.eye(3)):
                angle, ax, ay, az = placement.rotation_axis_angle
                cmd = (
                    f"    monomer[{idx}] = new {monomer_name}"
                    f".rot({angle:.4f},{ax:.4f},{ay:.4f},{az:.4f})"
                    f".move({pos[0]:.4f},{pos[1]:.4f},{pos[2]:.4f})"
                )
            else:
                cmd = (
                    f"    monomer[{idx}] = new {monomer_name}"
                    f".move({pos[0]:.4f},{pos[1]:.4f},{pos[2]:.4f})"
                )

            commands.append(cmd)

        return commands
