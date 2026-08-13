#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Base machinery for built-in crystalline substrate slab builders.

SlabBuilder implements the carve-and-functionalize pipeline shared by all
surface builders (following mBuild's SilicaInterface recipe, mosdef-hub/
mbuild, BSD-3; adapted to crystalline bulks generated in code):

1. Replicate a periodic surface cell to cover the requested lateral area
   and thickness.
2. Cleave the slab: the z-shift is chosen automatically so that the cut
   runs along the natural cleavage plane (fewest broken bonds) while
   keeping every Si tetrahedrally coordinated.
3. Bond perception by Si-O distance cutoff (lateral minimum image).
4. Strip disconnected fragments (largest connected component).
5. Optionally form Si-O-Si bridges between dangling top-face oxygens to
   reduce the silanol density toward a target (only where the surface
   geometry allows it).
6. Cap remaining dangling oxygens with H (top and, by default, bottom).
7. Type atoms by role and emit a self-contained moltemplate .lt class.

Subclasses provide the surface cell (fractional coordinates + lattice
constants) via class attributes.

Created on 2026-08-12
@author: zwu
"""
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from ..core.exceptions import ValidationError
from ..core.system import logger
from .ff_tables import OH_BOND_LENGTH, SLAB_FF_REGISTRY, SLAB_ROLES

#: Si-O bond perception cutoff (A); Si-O bonds are ~1.6 A, next shell
#: starts above 2.4 A.
SI_O_BOND_CUTOFF = 1.9
#: Max Si-Si distance (lateral minimum image) for forming a Si-O-Si bridge
#: between dangling oxygens (mbuild: 0.45 nm).
BRIDGE_SI_SI_MAX = 4.5
#: z margin reserved on each slab face for the oxygen coating plus the
#: hydroxyl caps (A).  Si are kept inside [Z_MARGIN, thickness-Z_MARGIN],
#: which guarantees the final physical slab fits inside [0, thickness].
Z_MARGIN = 2.7


@dataclass
class SlabAtom:
    """One atom of a built slab."""
    name: str          # unique atom name (e.g. "si_3", "oh_7")
    role: str          # "Si", "OB", "OH" or "HO"
    type_name: str     # moltemplate @atom type
    charge: float
    position: np.ndarray


@dataclass
class SurfaceSlab:
    """
    A built crystalline substrate slab.

    Attributes:
        class_name: Moltemplate class name (single instance per system).
        lt_filename: File name the slab .lt is written to.
        lx, ly: Lateral slab dimensions (A); integer multiples of the
                surface cell.  The simulation box lateral sides must equal
                these values (seamless periodic slab).
        thickness: Requested slab envelope thickness (A).  The physical
                   slab (including hydroxyl caps) is guaranteed to fit
                   inside [0, thickness].
        z_extent: Actual physical z-extent of the slab atoms (A).
        atoms: Slab atoms (roles/types/charges/positions), coordinates
               centered laterally and at the slab mid-plane in z.
        n_silanol_top / n_silanol_bottom: silanol counts per face.
        slab_ff: Force-field table name used ("interface"/"clayff"/"custom").
        pair_substyle: LAMMPS pair sub-style prefix for the pair_coeff
                lines (e.g. "lj/charmm/coul/long"); required when the
                film force field uses a hybrid pair style, must be None
                otherwise.  Set automatically by BoxPacker from the film
                force field.
        bond_substyle: same for the bond_coeff lines ("harmonic" when
                the film force field uses a hybrid bond style).
        bonds: Si-O / O-H bond topology (atom index pairs).  Emitted
                with zero-force-coefficient harmonic terms: the slab is
                meant to be frozen, so the bonds add no forces but make
                special_bonds exclude the intra-slab 1-2/1-3/1-4
                nonbonded interactions (otherwise bonded O-H pairs at
                0.945 A would contribute enormous spurious LJ/Coulomb
                terms).
    """
    class_name: str
    lt_filename: str
    lx: float
    ly: float
    thickness: float
    z_extent: float
    atoms: List[SlabAtom]
    n_silanol_top: int
    n_silanol_bottom: int
    slab_ff: str
    pair_substyle: Optional[str] = None
    bond_substyle: Optional[str] = None
    bonds: Optional[List[Tuple[int, int]]] = None

    # ------------------------------------------------------------------
    def write_lt(self, directory) -> Path:
        """Write the slab as a self-contained moltemplate .lt file."""
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / self.lt_filename
        lines = [
            "# crystalline substrate slab built by AutoPoly",
            f"# lateral {self.lx:.4f} x {self.ly:.4f} A, envelope "
            f"thickness {self.thickness:.2f} A (physical z-extent "
            f"{self.z_extent:.2f} A)",
            f"# silanols: top {self.n_silanol_top}, "
            f"bottom {self.n_silanol_bottom}; slab_ff='{self.slab_ff}'",
            "# self-contained: masses and pair coefficients included;",
            "# no bonded terms (slab is meant to be held rigid/frozen)",
            "",
            f"{self.class_name} {{",
            "",
            '  write_once("Data Masses") {',
        ]
        seen = {}
        for atom in self.atoms:
            seen[atom.type_name] = atom
        mass_by_type = _mass_table(self.atoms)
        for tname in sorted(seen):
            lines.append(f"    @atom:{tname} {mass_by_type[tname]:.4f}")
        lines += ["  }", "", '  write_once("In Settings") {']
        lj_by_type = _lj_table(self.atoms)
        sub = f" {self.pair_substyle}" if self.pair_substyle else ""
        for tname in sorted(seen):
            eps, sig = lj_by_type[tname]
            lines.append(
                f"    pair_coeff @atom:{tname} @atom:{tname}{sub} "
                f"{eps:.6g} {sig:.6g}"
            )
        lines += ["  }", "", '  write("Data Atoms") {']
        for atom in self.atoms:
            x, y, z = atom.position
            lines.append(
                f"    $atom:{atom.name} $mol:. @atom:{atom.type_name} "
                f"{atom.charge:.4f} {x:.4f} {y:.4f} {z:.4f}"
            )
        lines += ["  }"]

        if self.bonds:
            # Bond types: Si-O ("*_si_o") and O-H ("*_o_h")
            prefix = self.class_name.lower()

            def btype_name(roles):
                return (f"{prefix}_o_h" if "HO" in roles
                        else f"{prefix}_si_o")

            roles_of = [
                frozenset((self.atoms[i].role, self.atoms[j].role))
                for i, j in self.bonds
            ]
            type_names = sorted({btype_name(r) for r in roles_of},
                                reverse=True)  # si_o before o_h

            bsub = f" {self.bond_substyle}" if self.bond_substyle else ""
            lines += ["", '  write_once("In Settings") {']
            for tname in type_names:
                # equilibrium length: mean length of bonds of this type
                # (lateral minimum image: bonds crossing the slab's
                # periodic edges would otherwise appear stretched)
                d = []
                for (i, j), r in zip(self.bonds, roles_of):
                    if btype_name(r) != tname:
                        continue
                    dv = (self.atoms[i].position - self.atoms[j].position)
                    dv[0] -= self.lx * round(dv[0] / self.lx)
                    dv[1] -= self.ly * round(dv[1] / self.ly)
                    d.append(float(np.linalg.norm(dv)))
                lines.append(
                    f"    bond_coeff @bond:{tname}{bsub} "
                    f"0.0 {float(np.mean(d)):.4f}"
                )
            lines += ["  }", "", '  write("Data Bonds") {']
            for k, ((i, j), r) in enumerate(
                    zip(self.bonds, roles_of), start=1):
                lines.append(
                    f"    $bond:b{k} @bond:{btype_name(r)} "
                    f"$atom:{self.atoms[i].name} $atom:{self.atoms[j].name}"
                )
            lines += ["  }"]

        lines += ["", "}", ""]
        path.write_text("\n".join(lines))
        logger.info(f"Slab written to {path} ({len(self.atoms)} atoms)")
        return path


def _mass_table(atoms) -> Dict[str, float]:
    masses = {"Si": 28.0855, "OB": 15.9994, "OH": 15.9994, "HO": 1.008}
    return {a.type_name: masses[a.role] for a in atoms}


def _lj_table(atoms) -> Dict[str, Tuple[float, float]]:
    return {a.type_name: (a._epsilon, a._sigma) for a in atoms}


class SlabBuilder:
    """
    Base class for crystalline substrate slab builders.

    Class attributes provided by subclasses:
        CELL_FRACTIONAL: ((elem, (fx, fy, fz)), ...) surface-cell atoms.
        CELL_A / CELL_B: lateral lattice constants (A).
        CELL_C: surface-cell z period (A).
        BUILDER_NAME: registry name (e.g. "alpha_quartz").
        DEFAULT_CLASS_NAME / DEFAULT_LT_FILENAME: moltemplate identifiers.

    Args:
        lx_target: Requested lateral x size (A); snapped to a multiple of
                   CELL_A.
        ly_target: Requested lateral y size (A); snapped to a multiple of
                   CELL_B.
        thickness: Slab envelope thickness (A).  The physical slab
                   including hydroxyl caps fits inside this extent; the Si
                   core occupies thickness - ~5.4 A.
        oh_density: Target top-surface silanol density (sites/nm^2).  The
                    builder attempts to reduce the intrinsic dangling-O
                    density to the target by forming Si-O-Si bridges; if
                    the surface geometry does not allow enough bridges
                    (or the target exceeds the intrinsic density), all
                    dangling oxygens are hydroxylated and a warning is
                    logged.
        hydroxylate_bottom: Cap bottom-face dangling oxygens with H
                    (default True; both faces are exposed in the
                    fully-periodic slab model).
        slab_ff: "interface", "clayff", or "custom".
        slab_types / slab_charges / slab_lj: per-role overrides of the
                    built-in tables (roles "Si", "OB", "OH", "HO"); all
                    three are required when slab_ff="custom".
        seed: RNG seed for the bridge pairing.
        pair_substyle: LAMMPS pair sub-style for the pair_coeff lines
                    (e.g. "lj/charmm/coul/long"); required when the film
                    force field uses a hybrid pair style.  Normally left
                    to BoxPacker, which sets it from the film force field.
        class_name / lt_filename: moltemplate identifiers for the slab.
    """

    CELL_FRACTIONAL: Tuple = ()
    CELL_A: float = 0.0
    CELL_B: float = 0.0
    CELL_C: float = 0.0
    BUILDER_NAME: str = "base"
    ATOM_PREFIX: str = "slab"
    DEFAULT_CLASS_NAME: str = "SurfaceSlab"
    DEFAULT_LT_FILENAME: str = "surface_slab.lt"

    def __init__(
        self,
        lx_target: float,
        ly_target: float,
        thickness: float,
        oh_density: float = 4.6,
        hydroxylate_bottom: bool = True,
        slab_ff: str = "interface",
        slab_types: Optional[Dict[str, str]] = None,
        slab_charges: Optional[Dict[str, float]] = None,
        slab_lj: Optional[Dict[str, Tuple[float, float]]] = None,
        seed: int = 42,
        pair_substyle: Optional[str] = None,
        bond_substyle: Optional[str] = None,
        class_name: Optional[str] = None,
        lt_filename: Optional[str] = None,
    ) -> None:
        self.lx_target = float(lx_target)
        self.ly_target = float(ly_target)
        self.thickness = float(thickness)
        self.oh_density = float(oh_density)
        self.hydroxylate_bottom = bool(hydroxylate_bottom)
        self.slab_ff = slab_ff
        self.seed = seed
        self.pair_substyle = pair_substyle
        self.bond_substyle = bond_substyle
        self.class_name = class_name or self.DEFAULT_CLASS_NAME
        self.lt_filename = lt_filename or self.DEFAULT_LT_FILENAME
        self._table = self._resolve_table(
            slab_ff, slab_types, slab_charges, slab_lj
        )
        self._validate()

    # ------------------------------------------------------------------
    # Validation / parameter resolution
    # ------------------------------------------------------------------
    @property
    def min_thickness(self) -> float:
        return 2.0 * Z_MARGIN + self.CELL_C / 2.0

    @staticmethod
    def _resolve_table(slab_ff, slab_types, slab_charges, slab_lj):
        if slab_ff == "custom":
            missing = [
                r for r in SLAB_ROLES
                if not (slab_types and slab_types.get(r)
                        and slab_charges and r in slab_charges
                        and slab_lj and r in slab_lj)
            ]
            if missing:
                raise ValidationError(
                    f"slab_ff='custom' requires slab_types, slab_charges "
                    f"and slab_lj covering all roles {SLAB_ROLES}; "
                    f"missing: {missing}"
                )
            return {
                r: {
                    "type": slab_types[r],
                    "charge": float(slab_charges[r]),
                    "epsilon": float(slab_lj[r][0]),
                    "sigma": float(slab_lj[r][1]),
                }
                for r in SLAB_ROLES
            }
        if slab_ff not in SLAB_FF_REGISTRY:
            raise ValidationError(
                f"Unknown slab_ff '{slab_ff}'; choose from "
                f"{sorted(SLAB_FF_REGISTRY)} or 'custom'"
            )
        base = {r: dict(e) for r, e in SLAB_FF_REGISTRY[slab_ff].items()}
        # per-role overrides on top of a built-in table
        for role, tname in (slab_types or {}).items():
            if role in base:
                base[role]["type"] = tname
        for role, q in (slab_charges or {}).items():
            if role in base:
                base[role]["charge"] = float(q)
        for role, lj in (slab_lj or {}).items():
            if role in base:
                base[role]["epsilon"], base[role]["sigma"] = (
                    float(lj[0]), float(lj[1])
                )
        return base

    def _validate(self) -> None:
        if self.lx_target <= 0 or self.ly_target <= 0:
            raise ValidationError(
                f"{self.BUILDER_NAME} builder requires positive lateral "
                f"targets, got ({self.lx_target}, {self.ly_target})"
            )
        if self.thickness < self.min_thickness:
            raise ValidationError(
                f"Slab thickness must be >= {self.min_thickness:.1f} A "
                f"(envelope incl. hydroxyl coatings), got {self.thickness}"
            )
        if self.oh_density < 0:
            raise ValidationError(
                f"oh_density must be >= 0, got {self.oh_density}"
            )

    # ------------------------------------------------------------------
    # Public build
    # ------------------------------------------------------------------
    def build(self) -> SurfaceSlab:
        nx = max(1, int(round(self.lx_target / self.CELL_A)))
        ny = max(1, int(round(self.ly_target / self.CELL_B)))
        lx, ly = nx * self.CELL_A, ny * self.CELL_B

        species, pos = self._tile(nx, ny)
        keep = self._cleave(species, pos, lx, ly)
        species, pos = species[keep], pos[keep]
        bonds = self._perceive_bonds(species, pos, lx, ly)
        species, pos, bonds = self._strip_strays(species, pos, bonds)

        coord = _coordination(len(species), bonds)
        si_idx = np.where(species == "Si")[0]
        if not np.all(coord[si_idx] == 4):
            bad = int(np.sum(coord[si_idx] != 4))
            raise ValidationError(
                f"{self.BUILDER_NAME} slab cleave left {bad} Si without "
                "full tetrahedral coordination; increase the slab "
                "thickness"
            )

        z_si_top = pos[si_idx][:, 2].max()
        dangling_top = [
            i for i in np.where(species == "O")[0]
            if coord[i] == 1 and pos[i][2] > z_si_top
        ]

        species, pos, bonds, dangling_top = self._bridge_dangling(
            species, pos, bonds, dangling_top, lx, ly
        )

        # indices were remapped by bridging: recompute coordination and
        # the bottom-face dangling oxygens
        coord = _coordination(len(species), bonds)
        dangling_set = set(dangling_top)
        dangling_bottom = [
            i for i in np.where(species == "O")[0]
            if coord[i] == 1 and i not in dangling_set
        ]

        # Hydroxylation: cap remaining dangling oxygens with H
        oh_len = OH_BOND_LENGTH.get(self.slab_ff, 0.98)
        species_l = list(species)
        pos_l = list(pos)
        bonds_l = list(bonds)
        n_top_capped = 0
        for i in dangling_top:
            self._cap_with_h(species_l, pos_l, bonds_l, i, oh_len)
            n_top_capped += 1
        n_bottom_capped = 0
        for i in dangling_bottom:
            if self.hydroxylate_bottom:
                self._cap_with_h(species_l, pos_l, bonds_l, i, oh_len)
                n_bottom_capped += 1
            # else: leave bare (documented; slab no longer neutral)

        atoms = self._assign_atoms(species_l, pos_l, bonds_l)
        slab_bonds = list(bonds_l)

        # Center coordinates: lateral at origin, z at the slab mid-plane.
        # Wrap lateral coordinates into [-l/2, l/2): edge hydroxyl caps
        # can point slightly outside the cell (harmless for a periodic
        # slab, but tidy in the data file).
        arr = np.array([a.position for a in atoms])
        arr[:, 0] -= lx / 2.0
        arr[:, 1] -= ly / 2.0
        arr[:, 0] -= lx * np.floor(arr[:, 0] / lx + 0.5)
        arr[:, 1] -= ly * np.floor(arr[:, 1] / ly + 0.5)
        z_mid = 0.5 * (arr[:, 2].min() + arr[:, 2].max())
        arr[:, 2] -= z_mid
        for a, p in zip(atoms, arr):
            a.position = p
        z_extent = float(arr[:, 2].max() - arr[:, 2].min())

        total_q = float(sum(a.charge for a in atoms))
        if abs(total_q) > 1e-4:
            logger.warning(
                f"{self.BUILDER_NAME} slab net charge {total_q:+.4f} e "
                "(expected ~0 for a fully hydroxylated slab)"
            )

        achieved = n_top_capped / (lx * ly / 100.0)
        logger.info(
            f"{self.BUILDER_NAME} slab: {nx}x{ny} cells "
            f"({lx:.3f} x {ly:.3f} A), {len(atoms)} atoms, "
            f"z-extent {z_extent:.2f} A, silanols top {n_top_capped} / "
            f"bottom {n_bottom_capped} (achieved top density "
            f"{achieved:.2f}/nm^2, target {self.oh_density}/nm^2)"
        )
        return SurfaceSlab(
            class_name=self.class_name,
            lt_filename=self.lt_filename,
            lx=lx, ly=ly,
            thickness=self.thickness,
            z_extent=z_extent,
            atoms=atoms,
            n_silanol_top=n_top_capped,
            n_silanol_bottom=n_bottom_capped,
            slab_ff=self.slab_ff,
            pair_substyle=self.pair_substyle,
            bond_substyle=self.bond_substyle,
            bonds=slab_bonds,
        )

    # ------------------------------------------------------------------
    # Pipeline steps
    # ------------------------------------------------------------------
    def _tile(self, nx, ny):
        """Replicate the surface cell; z as needed for the cut."""
        nz = int(math.ceil(self.thickness / self.CELL_C)) + 1
        species, positions = [], []
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    for elem, (fx, fy, fz) in self.CELL_FRACTIONAL:
                        species.append(elem)
                        positions.append((
                            (ix + fx) * self.CELL_A,
                            (iy + fy) * self.CELL_B,
                            (iz + fz) * self.CELL_C,
                        ))
        return np.array(species), np.array(positions)

    def _cleave(self, species, pos, lx, ly):
        """
        Choose the z-shift whose cut runs along the natural cleavage
        plane, then keep Si in the interior window and O in the
        buffer-coated range.

        Candidate shifts align each distinct Si layer with the window
        bottom.  The winner keeps every Si tetrahedrally coordinated
        (given the O coating) and breaks the fewest bonds (smallest
        dangling-O count).  Raises ValidationError if no candidate works.
        """
        z = pos[:, 2]
        si_mask = species == "Si"
        si_z = np.sort(z[si_mask])
        planes = []
        for zv in si_z:
            if not planes or zv - planes[-1] > 0.2:
                planes.append(zv)

        si_top = self.thickness - Z_MARGIN
        best = None  # (n_dangling, shift, keep, n_dangling_invalid)
        for p in planes:
            shift = Z_MARGIN - p
            zz = z + shift
            keep = (
                (si_mask & (zz >= Z_MARGIN - 1e-6) & (zz <= si_top))
                | ((species == "O") & (zz >= Z_MARGIN - 1.65)
                   & (zz <= si_top + 1.65))
            )
            if not np.any(keep & si_mask):
                continue
            sub_species = species[keep]
            sub_pos = pos[keep].copy()
            sub_pos[:, 2] = zz[keep]
            bonds = self._perceive_bonds(sub_species, sub_pos, lx, ly)
            coord = _coordination(len(sub_species), bonds)
            sub_si = np.where(sub_species == "Si")[0]
            n_bad = int(np.sum(coord[sub_si] != 4))
            n_dangling = int(np.sum(
                coord[sub_species == "O"] != 2
            ))
            candidate = (n_bad, n_dangling, keep)
            if best is None or (n_bad, n_dangling) < best[:2]:
                best = candidate

        if best is None or best[0] > 0:
            raise ValidationError(
                f"{self.BUILDER_NAME}: no cleavage plane keeps all Si "
                "tetrahedrally coordinated for thickness "
                f"{self.thickness} A; increase the slab thickness"
            )
        return best[2]

    @staticmethod
    def _perceive_bonds(species, pos, lx, ly):
        """Si-O bonds by distance cutoff (lateral minimum image)."""
        si = np.where(species == "Si")[0]
        ox = np.where(species == "O")[0]
        bonds = []
        for i in si:
            d = pos[ox] - pos[i]
            d[:, 0] -= lx * np.round(d[:, 0] / lx)
            d[:, 1] -= ly * np.round(d[:, 1] / ly)
            dist = np.linalg.norm(d, axis=1)
            for j, dij in zip(ox, dist):
                if dij < SI_O_BOND_CUTOFF:
                    bonds.append((int(i), int(j)))
        return bonds

    @staticmethod
    def _strip_strays(species, pos, bonds):
        """Keep the largest connected component of the bond graph."""
        neighbors = {}
        for i, j in bonds:
            neighbors.setdefault(i, set()).add(j)
            neighbors.setdefault(j, set()).add(i)
        seen, components = set(), []
        for start in list(neighbors):
            if start in seen:
                continue
            stack, comp = [start], set()
            while stack:
                node = stack.pop()
                if node in comp:
                    continue
                comp.add(node)
                stack.extend(neighbors.get(node, ()) - comp)
            seen |= comp
            components.append(comp)
        if not components:
            raise ValidationError(
                "Substrate slab: no bonded network found"
            )
        major = max(components, key=len)
        n_dropped = len(species) - len(major)
        if n_dropped:
            logger.info(f"Slab: stripped {n_dropped} stray atoms")
        keep_idx = sorted(major)
        remap = {old: new for new, old in enumerate(keep_idx)}
        return (
            species[keep_idx],
            pos[keep_idx],
            [(remap[i], remap[j]) for i, j in bonds
             if i in remap and j in remap],
        )

    def _bridge_dangling(self, species, pos, bonds, dangling_top, lx, ly):
        """
        Form Si-O-Si bridges to reduce the dangling-O count to the target
        silanol density (mbuild SilicaInterface._bridge_dangling_Os).
        """
        neighbors = {}
        for i, j in bonds:
            neighbors.setdefault(i, set()).add(j)
            neighbors.setdefault(j, set()).add(i)

        area_nm2 = lx * ly / 100.0
        target = int(round(self.oh_density * area_nm2))
        if len(dangling_top) <= target:
            if len(dangling_top) < target:
                logger.warning(
                    f"{self.BUILDER_NAME} slab has only "
                    f"{len(dangling_top)} dangling surface oxygens, "
                    f"below the target silanol count {target} "
                    f"({self.oh_density}/nm^2); keeping all"
                )
            return species, pos, bonds, dangling_top

        rng = np.random.default_rng(self.seed)
        species = species.copy()
        pos = pos.copy()
        bonds = list(bonds)
        dangling = list(dangling_top)
        n_bridges = int((len(dangling) - target) / 2)

        for _ in range(n_bridges):
            bridged = False
            attempts = 0
            while not bridged and attempts < 1000:
                attempts += 1
                o1 = dangling[rng.integers(len(dangling))]
                si1 = next(iter(neighbors[o1]))
                order = rng.permutation(len(dangling))
                for k in order:
                    o2 = dangling[k]
                    if o2 == o1:
                        continue
                    si2 = next(iter(neighbors[o2]))
                    if si1 == si2:
                        continue
                    if neighbors[si1] & neighbors[si2]:
                        continue
                    d = pos[si2] - pos[si1]
                    d[0] -= lx * round(d[0] / lx)
                    d[1] -= ly * round(d[1] / ly)
                    if np.linalg.norm(d) < BRIDGE_SI_SI_MAX:
                        bonds.append((o1, si2))
                        neighbors[o1].add(si2)
                        neighbors[si2].add(o1)
                        neighbors[si2].discard(o2)
                        neighbors.pop(o2, None)
                        # delete o2: mark and compact at the end
                        pos[o2] = np.nan
                        species[o2] = "X"
                        dangling.remove(o1)
                        dangling.remove(o2)
                        bridged = True
                        break
            if not bridged:
                logger.warning(
                    f"{self.BUILDER_NAME} slab: no geometrically valid "
                    "Si-O-Si bridge pairs found; the surface stays at "
                    "its intrinsic termination. Use a larger oh_density "
                    "or accept the intrinsic density."
                )
                break

        keep = species != "X"
        remap = {}
        new_idx = 0
        for old, k in enumerate(keep):
            if k:
                remap[old] = new_idx
                new_idx += 1
        species = species[keep]
        pos = pos[keep]
        bonds = [(remap[i], remap[j]) for i, j in bonds
                 if i in remap and j in remap]
        dangling_top = [remap[i] for i in dangling]
        return species, pos, bonds, dangling_top

    @staticmethod
    def _cap_with_h(species_l, pos_l, bonds_l, o_idx, oh_len):
        """Cap a dangling oxygen with H along the O-Si bond direction."""
        si_idx = None
        for i, j in bonds_l:
            if i == o_idx and species_l[j] == "Si":
                si_idx = j
                break
            if j == o_idx and species_l[i] == "Si":
                si_idx = i
                break
        if si_idx is None:
            raise ValidationError("Substrate slab: dangling O without Si")
        direction = pos_l[o_idx] - pos_l[si_idx]
        direction = direction / np.linalg.norm(direction)
        h_idx = len(species_l)
        species_l.append("H")
        pos_l.append(pos_l[o_idx] + oh_len * direction)
        bonds_l.append((o_idx, h_idx))

    def _assign_atoms(self, species_l, pos_l, bonds_l):
        """Assign roles/types/charges; every Si must stay 4-coordinated."""
        si_neighbors = {}
        for i, j in bonds_l:
            if species_l[j] == "Si":
                si_neighbors.setdefault(i, set()).add(j)
            if species_l[i] == "Si":
                si_neighbors.setdefault(j, set()).add(i)
        counters = {"Si": 0, "OB": 0, "OH": 0, "HO": 0}
        atoms = []
        table = self._table
        for idx, (elem, p) in enumerate(zip(species_l, pos_l)):
            if elem == "Si":
                role = "Si"
            elif elem == "H":
                role = "HO"
            elif len(si_neighbors.get(idx, ())) >= 2:
                role = "OB"
            else:
                role = "OH"
            counters[role] += 1
            entry = table[role]
            atom = SlabAtom(
                name=f"{self.ATOM_PREFIX}_{role.lower()}_{counters[role]}",
                role=role,
                type_name=entry["type"],
                charge=entry["charge"],
                position=np.asarray(p, dtype=float),
            )
            atom._epsilon = entry["epsilon"]
            atom._sigma = entry["sigma"]
            atoms.append(atom)
        return atoms


def _coordination(n, bonds):
    coord = np.zeros(n, dtype=int)
    for i, j in bonds:
        coord[i] += 1
        coord[j] += 1
    return coord
