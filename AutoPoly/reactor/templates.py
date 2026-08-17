#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reaction Template Writer for the AutoPoly Reactor

Turns :class:`~AutoPoly.reactor.mapper.ReactionMetadata` into the file
triplet LAMMPS ``fix bond/react`` consumes:

- ``template_pre_N.molecule`` / ``template_post_N.molecule`` — LAMMPS
  molecule files restricted to the template atoms, typed with the
  simulation's force field, carrying embedded 3D coordinates and the
  bonded topology (bonds/angles/dihedrals/impropers) of the reaction site.
- ``RXN_N.map`` — the superimpose map file (InitiatorIDs, EdgeIDs,
  Equivalences, and a ``RXN_N_with_delete_ids.map`` variant carrying
  DeleteIDs when the reaction splits off a byproduct).

Type resolution strategy:

- Atom types come from :class:`~AutoPoly.monomers.monomer_generator.SMARTSTyper`
  on the combined reactant/product molecules (full chemical environment).
- Numeric ids reuse the ids already present in ``system.data`` whenever a
  matching interaction exists there (by-example matching via
  :class:`~AutoPoly.reactor.system_data.SystemData`). Genuinely new
  interactions (e.g. the cross-monomer bond a reaction creates) get fresh
  ids appended after the simulation's counts, and their coefficients are
  written to ``system.in.settings.reactor`` from the force-field tables.

Coordinates are generated per template pair: RDKit embeds a 3D conformer
of the combined reactants, a constrained embedding keeps the product's
shared heavy-atom frame consistent with the reactant so the pre/post
superimposition LAMMPS performs is well-defined.

Created on 2026-08-04
@author: zwu
"""
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

from ..core.exceptions import GenerationError
from ..core.system import logger
from ..monomers.monomer_generator import SMARTSTyper
from .fftables import ForceFieldTables
from .mapper import ReactionMetadata
from .system_data import SystemData


@dataclass
class TypeAssignment:
    """
    Numeric type assignment for one template pair against one data file.

    The pre and post templates of one reaction resolve against a SINGLE
    assignment with one global counter: atom types are keyed by symbolic
    name (a changed atom type keeps one id across pre/post), while new
    bonded interactions are keyed by (coeff, canonical-type) so a
    reaction-changed interaction gets one id shared by both templates.
    Existing data-file types are reused by example.

    Attributes:
        atom_type_ids: template atom index (0-based) -> numeric atom type id.
        bond_types / angle_types / dihedral_types / improper_types:
            canonical symbolic key -> numeric id (new ids included).
        new_atom_types: symbolic name -> new numeric id.
        new_coeff_lines: supplementary ``*_coeff`` lines for new types.
        new_mass_lines: ``mass`` lines for new atom types.
        max_ids: per-kind highest numeric id used.
    """
    atom_type_ids: Dict[int, int] = field(default_factory=dict)
    bond_types: Dict[Tuple[str, ...], int] = field(default_factory=dict)
    angle_types: Dict[Tuple[str, ...], int] = field(default_factory=dict)
    dihedral_types: Dict[Tuple[str, ...], int] = field(default_factory=dict)
    improper_types: Dict[Tuple[str, ...], int] = field(default_factory=dict)
    new_atom_types: Dict[str, int] = field(default_factory=dict)
    new_coeff_lines: List[str] = field(default_factory=list)
    new_mass_lines: List[str] = field(default_factory=list)
    max_ids: Dict[str, int] = field(default_factory=dict)


@dataclass
class TemplateSet:
    """Paths and metadata for one built reaction template triplet."""
    reaction_id: int
    reaction_name: str
    pre_file: Path
    post_file: Path
    map_file: Path
    map_file_delete_ids: Optional[Path]
    delete_atom: bool


def _canonical(kind: str, types: Tuple[str, ...]) -> Tuple[str, ...]:
    if kind == "improper":
        return (types[0], types[1]) + tuple(sorted(types[2:4]))
    return min(types, tuple(reversed(types)))


def _embed_reactant_product(
    reactant: Chem.Mol,
    product: Chem.Mol,
    r_to_p: Dict[int, int],
    seed: int = 42,
) -> Chem.Mol:
    """
    Embed 3D coordinates for the reactant, then the product constrained to it.

    The product is embedded with a coordinate constraint on the atoms it
    shares with the reactant (all mapped atoms), so pre/post templates are
    spatially superimposable. Returns the embedded product.
    """
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(reactant, params) != 0:
        AllChem.EmbedMolecule(reactant, randomSeed=seed, useRandomCoords=True)

    r_conf = reactant.GetConformer()
    # Build a coordinate constraint map: product heavy-atom index -> reactant coords
    cmap: Dict[int, Tuple[float, float, float]] = {}
    for r_idx, p_idx in r_to_p.items():
        if reactant.GetAtomWithIdx(r_idx).GetAtomicNum() > 1:
            pos = r_conf.GetAtomPosition(r_idx)
            cmap[p_idx] = (pos.x, pos.y, pos.z)

    p_params = AllChem.ETKDGv3()
    p_params.randomSeed = seed
    try:
        p_params.SetCoordMap(cmap)
        if AllChem.EmbedMolecule(product, p_params) != 0:
            AllChem.EmbedMolecule(product, randomSeed=seed, useRandomCoords=True)
    except Exception:
        # Constrained embedding unsupported/failed: fall back to free embedding
        if AllChem.EmbedMolecule(product, p_params) != 0:
            AllChem.EmbedMolecule(product, randomSeed=seed, useRandomCoords=True)
    return product


class TemplateBuilder:
    """
    Builds fix bond/react template triplets from reaction metadata.

    Args:
        typer: SMARTSTyper for the simulation's force field.
        ff_tables: Parsed force-field tables (symbolic coeff lookup).
        system_data: Parsed system.data (by-example numeric type ids).
        charge_fallback: Value used when an atom type has no charge entry.
    """

    def __init__(
        self,
        typer: SMARTSTyper,
        ff_tables: ForceFieldTables,
        system_data: SystemData,
    ) -> None:
        self.typer = typer
        self.ff = ff_tables
        self.data = system_data

    # ------------------------------------------------------------------
    # Typing
    # ------------------------------------------------------------------
    def _type_molecule(self, mol: Chem.Mol) -> List[str]:
        """Assign force-field atom types; return per-atom symbolic types."""
        self.typer.assign_atom_types(mol)
        types: List[str] = []
        for atom in mol.GetAtoms():
            try:
                raw = atom.GetProp("AtomType")
            except KeyError:
                raw = f"@atom:{atom.GetSymbol().lower()}"
            # .fdefn features carry a leading "@atom:"; system.data Masses
            # comments store the bare name (e.g. "c3"), so strip the prefix.
            types.append(raw.split("@atom:", 1)[-1])
        return types

    def _charge_of(self, atom_type: str, atom) -> float:
        """Per-atom charge mirroring UnitTyper's policy."""
        if self.typer.force_field in ("gaff", "gaff2"):
            try:
                charge = float(atom.GetProp("_GasteigerCharge"))
                if charge != charge or charge in (float("inf"), float("-inf")):
                    return 0.0
                return charge
            except (KeyError, ValueError):
                return 0.0
        return self.typer.charge_dict.get(atom_type, 0.0)

    # ------------------------------------------------------------------
    # Topology enumeration
    # ------------------------------------------------------------------
    @staticmethod
    def _enumerate_topology(
        mol: Chem.Mol,
        keep: List[int],
    ) -> Dict[str, List[Tuple[int, ...]]]:
        """
        Bonds/angles/dihedrals/impropers over the kept atoms.

        Impropers follow the cvff convention used by moltemplate: the
        central atom is first, and the two outer atoms that LAMMPS treats
        as swappable are the last two (sorted for canonicalization).
        """
        keep_set = set(keep)
        topo: Dict[str, List[Tuple[int, ...]]] = {"bond": [], "angle": [],
                                                  "dihedral": [], "improper": []}
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if i in keep_set and j in keep_set:
                topo["bond"].append((i, j))

        # angles + dihedrals from paths
        for center in keep:
            neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(center).GetNeighbors()
                         if n.GetIdx() in keep_set]
            for a in neighbors:
                # angle a-center-b
                for b in neighbors:
                    if b > a:
                        topo["angle"].append((a, center, b))
                # dihedral x-a-center-b / a-center-b-y handled by path walk
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if i not in keep_set or j not in keep_set:
                continue
            i_nbrs = [n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors()
                      if n.GetIdx() in keep_set and n.GetIdx() != j]
            j_nbrs = [n.GetIdx() for n in mol.GetAtomWithIdx(j).GetNeighbors()
                      if n.GetIdx() in keep_set and n.GetIdx() != i]
            for a in i_nbrs:
                for b in j_nbrs:
                    if a != b:
                        topo["dihedral"].append((a, i, j, b))

        # impropers: sp2 (trigonal) centers with exactly 3 kept neighbors.
        # sp3 centers are tetrahedral and do not get cvff impropers.
        for center in keep:
            center_atom = mol.GetAtomWithIdx(center)
            if center_atom.GetHybridization() != Chem.HybridizationType.SP2:
                continue
            neighbors = sorted(n.GetIdx() for n in center_atom.GetNeighbors()
                               if n.GetIdx() in keep_set)
            if len(neighbors) == 3:
                n1, n2, n3 = neighbors
                # cvff order: central first; swappable pair last (sorted)
                topo["improper"].append((center, n1, min(n2, n3), max(n2, n3)))
        return topo

    # ------------------------------------------------------------------
    # Numeric type resolution
    # ------------------------------------------------------------------
    def _resolve_atom_types(
        self,
        types: List[str],
        assignment: TypeAssignment,
    ) -> None:
        """Resolve numeric atom type ids, keyed by symbolic name (shared)."""
        next_id = assignment.max_ids.get(
            "atom",
            max(self.data.max_atom_type(), self.data.header_types_count("atom types")),
        )
        for sym in types:
            if sym in assignment.new_atom_types or self.data.atom_type_id(sym) is not None:
                continue
            next_id += 1
            assignment.new_atom_types[sym] = next_id
            mass = self.ff.mass_for(sym)
            if mass is not None:
                assignment.new_mass_lines.append(f"{next_id} {mass}  # {sym}")
            pair_line = self.ff.pair_coeff_line(sym, sym, next_id, next_id)
            if pair_line is not None:
                assignment.new_coeff_lines.append(pair_line)
        assignment.max_ids["atom"] = next_id

    def _atom_type_id(self, sym: str, assignment: TypeAssignment) -> int:
        existing = self.data.atom_type_id(sym)
        if existing is not None:
            return existing
        return assignment.new_atom_types[sym]

    def _resolve_bonded_types(
        self,
        kind: str,
        interactions: List[Tuple[int, ...]],
        types: List[str],
        assignment: TypeAssignment,
    ) -> None:
        """Resolve numeric ids for one bonded kind (by example, else new).

        Wildcard definitions (e.g. GAFF ``X-c3-c3-X`` dihedrals) collapse
        many concrete interactions onto one parameter; identical
        (concrete-types -> resolved coeff) interactions share a single new id.
        """
        table = getattr(assignment, f"{kind}_types")
        next_id = assignment.max_ids.get(
            kind,
            max(self.data.max_type(kind), self.data.header_types_count(f"{kind} types")),
        )
        ff_lookup = getattr(self.ff, f"{kind}_for")
        # Track (coeff body) -> id across BOTH pre and post calls (stored on
        # the assignment) so identical parameters share one id globally.
        resolved: Dict[str, int] = assignment.__dict__.setdefault(
            f"_resolved_{kind}", {})
        for inter in interactions:
            sym_key = tuple(types[i] for i in inter)
            canon = _canonical(kind, sym_key)
            if canon in table:
                continue
            existing = self.data.lookup(kind, sym_key)
            if existing is not None:
                table[canon] = existing
                continue
            found = ff_lookup(*sym_key)
            if found is None:
                logger.warning(
                    f"No {kind} parameter for {'-'.join(sym_key)}; skipping"
                )
                continue
            def_name, coeff = found
            coeff_body = coeff.split("#")[0].strip()
            if coeff_body in resolved:
                table[canon] = resolved[coeff_body]
                continue
            next_id += 1
            table[canon] = next_id
            resolved[coeff_body] = next_id
            assignment.new_coeff_lines.append(
                f"{kind}_coeff {next_id} {coeff_body}  # {def_name}")
        assignment.max_ids[kind] = next_id

    # ------------------------------------------------------------------
    # Molecule file writing
    # ------------------------------------------------------------------
    def _write_molecule(
        self,
        path: Path,
        mol: Chem.Mol,
        keep: List[int],
        types: List[str],
        assignment: TypeAssignment,
        title: str,
    ) -> None:
        """Write one LAMMPS molecule template file for the kept atoms."""
        keep_sorted = sorted(keep)
        renumber = {old: new for new, old in enumerate(keep_sorted, start=1)}
        topo = self._enumerate_topology(mol, keep_sorted)
        conf = mol.GetConformer() if mol.GetNumConformers() else None

        def sym_types(i: int) -> str:
            return types[i]

        lines: List[str] = [f"{title}", ""]
        lines.append(f"{len(keep_sorted)} atoms")
        lines.append(f"{len(topo['bond'])} bonds")
        lines.append(f"{len(topo['angle'])} angles")
        lines.append(f"{len(topo['dihedral'])} dihedrals")
        lines.append(f"{len(topo['improper'])} impropers")
        lines.append("")

        lines.append("Types\n")
        for old in keep_sorted:
            lines.append(f"{renumber[old]} {self._atom_type_id(types[old], assignment)}")
        lines.append("")

        lines.append("Charges\n")
        for old in keep_sorted:
            charge = self._charge_of(types[old], mol.GetAtomWithIdx(old))
            lines.append(f"{renumber[old]} {charge:.4f}")
        lines.append("")

        if conf is not None:
            lines.append("Coords\n")
            for old in keep_sorted:
                pos = conf.GetAtomPosition(old)
                lines.append(f"{renumber[old]} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}")
            lines.append("")

        def write_topology(section: str, kind: str) -> None:
            lines.append(f"{section}\n")
            table = getattr(assignment, f"{kind}_types")
            for n, inter in enumerate(topo[kind], start=1):
                canon = _canonical(kind, tuple(sym_types(i) for i in inter))
                type_id = table.get(canon)
                if type_id is None:
                    continue
                atom_ids = " ".join(str(renumber[i]) for i in inter)
                lines.append(f"{n} {type_id} {atom_ids}")
            lines.append("")

        if topo["bond"]:
            write_topology("Bonds", "bond")
        if topo["angle"]:
            write_topology("Angles", "angle")
        if topo["dihedral"]:
            write_topology("Dihedrals", "dihedral")
        if topo["improper"]:
            write_topology("Impropers", "improper")

        path.write_text("\n".join(lines) + "\n")

    # ------------------------------------------------------------------
    # Map file writing
    # ------------------------------------------------------------------
    @staticmethod
    def write_map_file(
        path: Path,
        equivalences: Dict[int, int],
        initiators: List[int],
        edge_atoms: List[int],
        delete_ids: Optional[List[int]] = None,
        title: str = "reaction",
    ) -> None:
        """
        Write a fix bond/react map file (1-based template atom ids).
        """
        lines: List[str] = [f"# map file for {title} (generated by AutoPoly reactor)", ""]
        lines.append(f"{len(edge_atoms)} edgeIDs")
        lines.append(f"{len(equivalences)} equivalences")
        if delete_ids:
            lines.append(f"{len(delete_ids)} deleteIDs")
        lines.append("")
        lines.append("InitiatorIDs\n")
        for i in sorted(initiators):
            lines.append(str(i))
        lines.append("")
        lines.append("EdgeIDs\n")
        for i in sorted(edge_atoms):
            lines.append(str(i))
        lines.append("")
        lines.append("Equivalences\n")
        for r, p in sorted(equivalences.items()):
            lines.append(f"{r:<5} {p}")
        if delete_ids:
            lines.append("")
            lines.append("DeleteIDs\n")
            for i in sorted(delete_ids):
                lines.append(str(i))
        path.write_text("\n".join(lines) + "\n")

    # ------------------------------------------------------------------
    # Top-level build
    # ------------------------------------------------------------------
    def build(
        self,
        metadata: ReactionMetadata,
        out_dir: Path,
        seed: int = 42,
    ) -> Tuple[TemplateSet, TypeAssignment]:
        """
        Build the pre/post molecule templates and map file for one reaction.

        Returns:
            (TemplateSet, TypeAssignment): file paths and the type
            assignment (supplementary coeffs/masses to append to the
            simulation's settings).
        """
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        rid = metadata.reaction_id

        reactant = metadata.reactant_mol
        product = metadata.product_mol
        r_to_p = metadata.reactant_to_product
        p_to_r = metadata.product_to_reactant

        # Type both molecules with full chemical environment
        pre_types = self._type_molecule(reactant)
        post_types = self._type_molecule(product)

        # GAFF/GAFF2 have no charge table: Gasteiger charges on both sides
        # of the reaction (mirrors UnitTyper's full-molecule policy).
        if self.typer.force_field in ("gaff", "gaff2"):
            for side, m in (("reactant", reactant), ("product", product)):
                try:
                    AllChem.ComputeGasteigerCharges(m)
                except Exception as e:
                    logger.error(
                        f"Failed to compute Gasteiger charges for {side} "
                        f"of reaction {rid}: {e}"
                    )

        # Embed coordinates (product constrained onto reactant frame)
        product = _embed_reactant_product(reactant, product, r_to_p, seed=seed)
        metadata.product_mol = product

        # Template atom subsets (reactant space and their product images)
        pre_keep = sorted(metadata.template_atoms)
        post_keep = sorted(r_to_p[i] for i in pre_keep if i in r_to_p)

        # One SHARED assignment: atom types keyed by symbolic name (a changed
        # atom type keeps one id), new bonded interactions keyed by
        # (coeff, canonical types) so a reaction-changed interaction shares
        # one id across pre and post. One global counter per kind means no id
        # collisions between the two templates.
        assignment = TypeAssignment()
        self._resolve_atom_types([pre_types[i] for i in pre_keep], assignment)
        self._resolve_atom_types([post_types[i] for i in post_keep], assignment)

        # Bonded types: pre (reactant types) then post (product types).
        pre_topo = self._enumerate_topology(reactant, pre_keep)
        post_topo = self._enumerate_topology(product, post_keep)
        for kind in ("bond", "angle", "dihedral", "improper"):
            self._resolve_bonded_types(kind, pre_topo[kind], pre_types, assignment)
            self._resolve_bonded_types(kind, post_topo[kind], post_types, assignment)

        # Map template-space equivalences (1-based ids within each template)
        pre_renumber = {old: new for new, old in enumerate(pre_keep, start=1)}
        post_renumber = {old: new for new, old in enumerate(post_keep, start=1)}
        equivalences = {
            pre_renumber[r]: post_renumber[r_to_p[r]]
            for r in pre_keep if r in r_to_p and r_to_p[r] in post_renumber
        }
        initiators = [pre_renumber[i] for i in metadata.initiators if i in pre_renumber]
        edge_atoms = [pre_renumber[i] for i in metadata.edge_atoms if i in pre_renumber]
        delete_ids = [pre_renumber[i] for i in metadata.byproduct_indices
                      if i in pre_renumber] if metadata.delete_atom else []

        if len(initiators) != 2:
            raise GenerationError(
                f"Reaction {rid}: expected 2 initiators in template space, "
                f"got {initiators}"
            )

        pre_path = out_dir / f"template_pre_{rid}.molecule"
        post_path = out_dir / f"template_post_{rid}.molecule"
        map_path = out_dir / f"RXN_{rid}.map"

        # The molecule files store the template-space renumbered atoms; the
        # writer re-enumerates topology itself, so pass full types + keep.
        self._write_molecule(pre_path, reactant, pre_keep, pre_types,
                             assignment, f"pre-reaction template {rid}")
        self._write_molecule(post_path, product, post_keep, post_types,
                             assignment, f"post-reaction template {rid}")

        self.write_map_file(map_path, equivalences, initiators, edge_atoms,
                            title=f"RXN_{rid}")
        map_delete_path: Optional[Path] = None
        if delete_ids:
            map_delete_path = out_dir / f"RXN_{rid}_with_delete_ids.map"
            self.write_map_file(map_delete_path, equivalences, initiators,
                                edge_atoms, delete_ids=delete_ids,
                                title=f"RXN_{rid} (with DeleteIDs)")

        template_set = TemplateSet(
            reaction_id=rid,
            reaction_name=metadata.instance.reaction_name,
            pre_file=pre_path,
            post_file=post_path,
            map_file=map_path,
            map_file_delete_ids=map_delete_path,
            delete_atom=metadata.delete_atom,
        )
        return template_set, assignment
