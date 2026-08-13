#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Force-Field Parameter Tables for the AutoPoly Reactor

Parses a force-field ``.lt`` file (gaff.lt, gaff2.lt, dreiding.lt,
compass_published.lt, oplsaa.lt, loplsaa.lt) into symbolically-keyed
parameter tables so the reactor can:

- assign the *same* bonded parameter the moltemplate "By Type" machinery
  would assign, by looking up ``bond_defs``/``angle_defs``/... entries whose
  ``@atom:`` patterns match the typed atoms of a reaction template, and
- emit supplementary ``bond_coeff``/``angle_coeff``/... lines for parameter
  types that do not yet exist in the generated ``system.data`` (e.g. the
  new cross-monomer bond created by a reaction).

Matching is wildcard-aware: a definition pattern such as ``@atom:*`` or
``@atom:*_bCT*_a*_d*_i*`` matches any concrete type consistent with the
non-wildcard parts. Definitions are tried in file order and the first
match wins (matching moltemplate's By-Type semantics).

Created on 2026-08-04
@author: zwu
"""
import re
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

from ..core.exceptions import GenerationError

_KIND_PREFIX = {
    "bond": "@bond:",
    "angle": "@angle:",
    "dihedral": "@dihedral:",
    "improper": "@improper:",
}
_N_ATOMS = {"bond": 2, "angle": 3, "dihedral": 4, "improper": 4}


class ForceFieldTables:
    """
    Symbolic parameter tables parsed from one force-field ``.lt`` file.

    Attributes:
        path: Source .lt file.
        masses: ``@atom:<type>`` -> mass.
        pair_coeffs: ``@atom:<t1> @atom:<t2>`` (canonical order) -> coeff text.
        <kind>_defs: list of (def_name, [atom patterns]) per bonded kind.
        <kind>_coeffs: def_name -> coeff text (after the ``@kind:name`` token).
    """

    def __init__(self, path: Path) -> None:
        self.path = Path(path)
        if not self.path.is_file():
            raise GenerationError(f"Force field file not found: {self.path}")
        self.masses: Dict[str, float] = {}
        self.pair_coeffs: Dict[Tuple[str, str], str] = {}
        self.bond_defs: List[Tuple[str, List[str]]] = []
        self.bond_coeffs: Dict[str, str] = {}
        self.angle_defs: List[Tuple[str, List[str]]] = []
        self.angle_coeffs: Dict[str, str] = {}
        self.dihedral_defs: List[Tuple[str, List[str]]] = []
        self.dihedral_coeffs: Dict[str, str] = {}
        self.improper_defs: List[Tuple[str, List[str]]] = []
        self.improper_coeffs: Dict[str, str] = {}
        self._parse()

    # ------------------------------------------------------------------
    # Parsing
    # ------------------------------------------------------------------
    def _parse(self) -> None:
        section: Optional[str] = None
        with open(self.path) as f:
            for raw in f:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                if 'write_once("Data Masses")' in line:
                    section = "masses"; continue
                if 'write_once("Data Bonds By Type")' in line:
                    section = "bond_defs"; continue
                if 'write_once("Data Angles By Type")' in line:
                    section = "angle_defs"; continue
                if 'write_once("Data Dihedrals By Type")' in line:
                    section = "dihedral_defs"; continue
                if 'write_once("Data Impropers By Type' in line:
                    section = "improper_defs"; continue
                if 'write_once("In Settings")' in line:
                    section = "settings"; continue
                if line.startswith("}"):
                    section = None; continue

                if section == "masses" and "@atom:" in line:
                    parts = line.split()
                    self.masses[parts[0].split("@atom:")[1]] = float(parts[1])
                elif section == "bond_defs" and "@bond:" in line:
                    self.bond_defs.append(self._parse_def(line, "bond"))
                elif section == "angle_defs" and "@angle:" in line:
                    self.angle_defs.append(self._parse_def(line, "angle"))
                elif section == "dihedral_defs" and "@dihedral:" in line:
                    self.dihedral_defs.append(self._parse_def(line, "dihedral"))
                elif section == "improper_defs" and "@improper:" in line:
                    self.improper_defs.append(self._parse_def(line, "improper"))
                elif section == "settings":
                    self._parse_coeff(line)

    @staticmethod
    def _parse_def(line: str, kind: str) -> Tuple[str, List[str]]:
        """Parse ``@kind:name @atom:p1 @atom:p2 ...`` into (name, patterns)."""
        tokens = line.split()
        name = tokens[0].split(_KIND_PREFIX[kind])[1]
        patterns = [t.split("@atom:")[1] for t in tokens[1:] if "@atom:" in t]
        return name, patterns

    def _parse_coeff(self, line: str) -> None:
        if "pair_coeff" in line and "@atom:" in line:
            types = re.findall(r"@atom:([^\s]+)", line)[:2]
            if len(types) == 2:
                key = tuple(sorted(types))
                # Coeff text = everything after the second "@atom:<type>" token.
                after_first = line.split(f"@atom:{types[0]}", 1)[1]
                coeff = after_first.split(f"@atom:{types[1]}", 1)[1].strip()
                self.pair_coeffs[key] = coeff
            return
        for kind, prefix in _KIND_PREFIX.items():
            coeff_kw = f"{kind}_coeff"
            if line.startswith(coeff_kw) and prefix in line:
                name = line.split(prefix)[1].split()[0]
                # Split on the full "@kind:<name>" token, not the bare name
                # (bare names can appear inside longer coeff text).
                marker = f"{prefix}{name}"
                coeff = line.split(marker, 1)[1].strip()
                getattr(self, f"{kind}_coeffs")[name] = coeff
                return

    # ------------------------------------------------------------------
    # Symbolic lookup
    # ------------------------------------------------------------------
    @staticmethod
    def _pattern_matches(pattern: str, atom_type: str) -> bool:
        """True if a By-Type ``@atom:`` pattern matches a concrete type.

        Exact matches always succeed; patterns containing ``*`` match when
        every literal (non-wildcard) fragment appears in order. This covers
        GAFF ``@atom:*`` wildcards, DREIDING prefixes (``C_3*``), and
        OPLSAA/COMPASS extended patterns (``*_bCT*_a*_d*_i*``,
        ``*~pc4~b*~a*~d*~i*``).
        """
        if pattern == atom_type or pattern == "*":
            return True
        if "*" not in pattern:
            return False
        pos = 0
        for frag in pattern.split("*"):
            if not frag:
                continue
            idx = atom_type.find(frag, pos)
            if idx == -1:
                return False
            pos = idx + len(frag)
        return True

    @staticmethod
    def _specificity(pattern: str) -> Tuple[int, int]:
        """
        Rank how specific a pattern is: (literal chars, -wildcards).
        Higher beats lower when comparing candidate definitions.
        """
        wildcards = pattern.count("*")
        literal = len(pattern) - wildcards
        return (literal, -wildcards)

    @staticmethod
    def _def_matches(patterns: List[str], types: Tuple[str, ...]) -> bool:
        """All patterns match the corresponding concrete types (in order)."""
        if len(patterns) != len(types):
            return False
        return all(
            ForceFieldTables._pattern_matches(p, t)
            for p, t in zip(patterns, types)
        )

    def _lookup(
        self, kind: str, types: Tuple[str, ...]
    ) -> Optional[Tuple[str, str]]:
        """
        Find the (def_name, coeff_text) moltemplate's By-Type machinery
        would assign for a bonded interaction over ``types``.

        Both the forward and reversed type sequences are considered. Among
        all matching definitions the *most specific* one wins (measured by
        the total literal content of its atom patterns), preferring exact
        matches over wildcard entries; ties keep file order.
        """
        defs = getattr(self, f"{kind}_defs")
        coeffs = getattr(self, f"{kind}_coeffs")
        # Specificity must beat file order: exact/concrete definitions take
        # priority over generic wildcard entries (e.g. "c3-c3-c3-c3" over
        # "X-c3-c3-X"). Score as (literal chars, -wildcards, -file order).
        best: Optional[Tuple[Tuple[int, int, int], str]] = None
        for order, (name, patterns) in enumerate(defs):
            if not (self._def_matches(patterns, types) or
                    self._def_matches(patterns, tuple(reversed(types)))):
                continue
            lit = sum(self._specificity(p)[0] for p in patterns)
            wild = sum(self._specificity(p)[1] for p in patterns)
            score = (lit, wild, -order)
            if best is None or score > best[0]:
                best = (score, name)
        if best is None:
            return None
        name = best[1]
        return name, coeffs.get(name, "")

    def bond_for(self, t1: str, t2: str) -> Optional[Tuple[str, str]]:
        """(def_name, coeff) for a bond between the two atom types."""
        return self._lookup("bond", (t1, t2))

    def angle_for(self, t1: str, t2: str, t3: str) -> Optional[Tuple[str, str]]:
        """(def_name, coeff) for an angle over the three atom types."""
        return self._lookup("angle", (t1, t2, t3))

    def dihedral_for(self, t1: str, t2: str, t3: str, t4: str) -> Optional[Tuple[str, str]]:
        """(def_name, coeff) for a dihedral over the four atom types."""
        return self._lookup("dihedral", (t1, t2, t3, t4))

    def improper_for(self, t1: str, t2: str, t3: str, t4: str) -> Optional[Tuple[str, str]]:
        """
        (def_name, coeff) for an improper over the four atom types.

        The reactor enumerates impropers in LAMMPS cvff order (central atom
        first, positions 3 and 4 swappable). GAFF/CVFF def names list the
        central atom at position 2 (e.g. ``X-o-c-os``, central ``c``); the
        def lookup tries both the (t2, central, t3, t4) GAFF layout and the
        direct central-first layout.
        """
        # GAFF def layout: (t2, central, t3, t4) and the t3<->t4 swap.
        for cand in ((t2, t1, t3, t4), (t2, t1, t4, t3),
                     (t1, t2, t3, t4), (t1, t2, t4, t3)):
            result = self._lookup("improper", cand)
            if result is not None:
                return result
        return None

    def pair_for(self, t1: str, t2: str) -> Optional[str]:
        """Pair coeff text for two atom types (exact or wildcard match)."""
        key = tuple(sorted((t1, t2)))
        if key in self.pair_coeffs:
            return self.pair_coeffs[key]
        for (a, b), coeff in self.pair_coeffs.items():
            if (self._pattern_matches(a, t1) and self._pattern_matches(b, t2)) or \
               (self._pattern_matches(a, t2) and self._pattern_matches(b, t1)):
                return coeff
        return None

    def pair_coeff_line(self, t1: str, t2: str, id1: int, id2: int) -> Optional[str]:
        """
        Formatted ``pair_coeff`` line for two numeric ids.

        The stored coeff text starts with the style name; any trailing
        comment from the source .lt is stripped.
        """
        coeff = self.pair_for(t1, t2)
        if coeff is None:
            return None
        return f"pair_coeff {id1} {id2} {coeff.split('#')[0].strip()}"

    def mass_for(self, atom_type: str) -> Optional[float]:
        """Mass for an atom type (exact or wildcard pattern match)."""
        if atom_type in self.masses:
            return self.masses[atom_type]
        for pattern, mass in self.masses.items():
            if self._pattern_matches(pattern, atom_type):
                return mass
        return None


def load_force_field_tables(force_field: str, extern_dir: Optional[Path] = None) -> ForceFieldTables:
    """
    Load the parameter tables for a registered force field.

    Args:
        force_field: One of FORCE_FIELD_REGISTRY keys.
        extern_dir: Override for the extern/ directory (defaults to the
                    package's bundled moltemplate force fields).

    Returns:
        Parsed ForceFieldTables.
    """
    from ..core.conf import FORCE_FIELD_REGISTRY

    if extern_dir is None:
        extern_dir = Path(__file__).parent.parent.resolve() / "extern"
    lt_file = FORCE_FIELD_REGISTRY[force_field]["lt_file"]
    return ForceFieldTables(Path(extern_dir) / "moltemplate" / "force_fields" / lt_file)
