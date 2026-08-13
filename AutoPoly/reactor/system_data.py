#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LAMMPS system.data Parser for the AutoPoly Reactor

Reads the ``system.data`` produced by AutoPoly's stage 3 (moltemplate) and
derives *by-example* numeric type assignments:

- atom types: numeric id -> symbolic name (from ``Masses`` ``# name``
  comments that moltemplate writes),
- bond/angle/dihedral/improper types: canonical symbolic key (the atom
  types of the participating atoms) -> numeric type id.

The reactor uses these maps so that template molecule files reference the
*same* numeric types the simulation already uses whenever a matching
interaction exists in the data file; only genuinely new interactions
(e.g. a cross-monomer bond created by a reaction) receive fresh ids.

Also parses ``system.in.init`` for the LAMMPS style commands that the
reaction input script reproduces.

Created on 2026-08-04
@author: zwu
"""
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from ..core.exceptions import GenerationError

_STYLE_COMMANDS = (
    "units", "dimension", "boundary", "atom_style",
    "bond_style", "angle_style", "dihedral_style", "improper_style",
    "pair_style", "kspace_style", "pair_modify", "special_bonds",
    "neighbor", "neigh_modify", "timestep",
)

_SECTION_ATOM_COUNT = {"Bonds": 2, "Angles": 3, "Dihedrals": 4, "Impropers": 4}
_SECTION_KEY = {"Bonds": "bond", "Angles": "angle",
                "Dihedrals": "dihedral", "Impropers": "improper"}


def _canonical(kind: str, types: Tuple[str, ...]) -> Tuple[str, ...]:
    """
    Canonicalize a symbolic interaction key.

    Bonds, angles and dihedrals are reversal-symmetric. Impropers use the
    cvff convention (atom 1 central; atoms 3 and 4 swappable).
    """
    if kind == "improper":
        return (types[0], types[1]) + tuple(sorted(types[2:4]))
    rev = tuple(reversed(types))
    return min(types, rev)


class SystemData:
    """
    Parsed view of one AutoPoly ``system.data`` file.

    Attributes:
        path: The data file.
        atom_type_names: numeric atom type -> symbolic name.
        bond_types / angle_types / dihedral_types / improper_types:
            canonical symbolic key -> numeric type id.
        header_counts: ``<n> <keyword>`` pairs from the header.
        box_bounds: {"x": (xlo, xhi), "y": ..., "z": ...} in Angstrom.
    """

    def __init__(self, path: Path) -> None:
        self.path = Path(path)
        if not self.path.is_file():
            raise GenerationError(f"system.data not found: {self.path}")
        self.atom_type_names: Dict[int, str] = {}
        self.bond_types: Dict[Tuple[str, ...], int] = {}
        self.angle_types: Dict[Tuple[str, ...], int] = {}
        self.dihedral_types: Dict[Tuple[str, ...], int] = {}
        self.improper_types: Dict[Tuple[str, ...], int] = {}
        self.header_counts: List[Tuple[int, str]] = []
        self.box_bounds: Dict[str, Tuple[float, float]] = {}
        self._parse()

    # ------------------------------------------------------------------
    # Parsing
    # ------------------------------------------------------------------
    def _parse(self) -> None:
        text = self.path.read_text()
        # Box bounds live in the header ("<lo> <hi> xlo xhi" lines).
        for axis in ("x", "y", "z"):
            m = re.search(
                rf"^\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+{axis}lo\s+{axis}hi",
                text, re.MULTILINE)
            if m:
                self.box_bounds[axis] = (float(m.group(1)), float(m.group(2)))

        sections: Dict[str, List[str]] = {}
        current: Optional[str] = None
        for raw in text.splitlines():
            stripped = raw.strip()
            if not stripped:
                continue
            header_match = re.match(r"^(\d+)\s+(\S.*)$", stripped)
            if current is None and header_match and not sections:
                count, keyword = int(header_match.group(1)), header_match.group(2)
                self.header_counts.append((count, keyword))
                continue
            # Section headers are single capitalized words (optionally with
            # a trailing comment), e.g. "Atoms  # full".
            if re.match(r"^[A-Z][A-Za-z ]*(\s+#.*)?$", stripped) and \
                    stripped.split()[0] in (
                        "Masses", "Atoms", "Bonds", "Angles", "Dihedrals",
                        "Impropers", "Velocities", "Pair", "Bond",
                        "Angle", "Dihedral", "Improper",
                    ):
                current = stripped.split("#")[0].strip()
                sections.setdefault(current, [])
                continue
            if current:
                sections[current].append(stripped)

        self._parse_masses(sections.get("Masses", []))
        atoms = self._parse_atoms(sections.get("Atoms", []))
        for section_name, n_atoms in _SECTION_ATOM_COUNT.items():
            self._parse_topology(sections.get(section_name, []),
                                 section_name, n_atoms, atoms)

    def _parse_masses(self, lines: List[str]) -> None:
        for line in lines:
            match = re.match(r"^(\d+)\s+([\d.eE+-]+)\s*(?:#\s*(\S+))?", line)
            if match:
                type_id = int(match.group(1))
                name = match.group(3) or str(type_id)
                self.atom_type_names[type_id] = name

    def _parse_atoms(self, lines: List[str]) -> Dict[int, int]:
        """atom id -> numeric atom type (style-agnostic column scan)."""
        atoms: Dict[int, int] = {}
        for line in lines:
            comment_free = line.split("#")[0].split()
            if len(comment_free) < 3:
                continue
            atoms[int(comment_free[0])] = int(comment_free[2])
        return atoms

    def _parse_topology(
        self,
        lines: List[str],
        section_name: str,
        n_atoms: int,
        atoms: Dict[int, int],
    ) -> None:
        target = getattr(self, f"{_SECTION_KEY[section_name]}_types")
        for line in lines:
            comment_free = line.split("#")[0].split()
            if len(comment_free) < 2 + n_atoms:
                continue
            type_id = int(comment_free[1])
            atom_ids = [int(x) for x in comment_free[2:2 + n_atoms]]
            try:
                sym = tuple(self.atom_type_names[atoms[a]] for a in atom_ids)
            except KeyError:
                continue  # atom id missing from Atoms (shouldn't happen)
            key = _canonical(_SECTION_KEY[section_name], sym)
            target.setdefault(key, type_id)

    # ------------------------------------------------------------------
    # By-example lookups
    # ------------------------------------------------------------------
    def atom_type_id(self, name: str) -> Optional[int]:
        """Numeric id for a symbolic atom type name."""
        for type_id, type_name in self.atom_type_names.items():
            if type_name == name:
                return type_id
        return None

    def max_atom_type(self) -> int:
        return max(self.atom_type_names, default=0)

    def max_type(self, kind: str) -> int:
        """Largest numeric type id for a bonded kind."""
        table = getattr(self, f"{kind}_types")
        return max(table.values(), default=0)

    def lookup(self, kind: str, types: Tuple[str, ...]) -> Optional[int]:
        """Numeric type id for a bonded interaction (by example)."""
        key = _canonical(kind, types)
        return getattr(self, f"{kind}_types").get(key)

    def header_types_count(self, keyword: str) -> int:
        """The ``<n> <keyword>`` header count (e.g. "bond types")."""
        for count, kw in self.header_counts:
            if kw == keyword:
                return count
        return 0


def parse_in_init(path: Path) -> Dict[str, str]:
    """
    Parse style commands from ``system.in.init``.

    Returns:
        Dict of command name -> full argument text (comment stripped),
        e.g. {"bond_style": "hybrid harmonic", ...}.
    """
    path = Path(path)
    styles: Dict[str, str] = {}
    if not path.is_file():
        return styles
    for raw in path.read_text().splitlines():
        line = raw.split("#")[0].strip()
        if not line:
            continue
        parts = line.split(None, 1)
        if parts and parts[0] in _STYLE_COMMANDS:
            styles[parts[0]] = parts[1] if len(parts) > 1 else ""
    return styles
