#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit Typer — Stage 2 of the AutoPoly three-stage pipeline.

Force-field assignment on top of stored geometry:

    geometry/geometry.json + force_field ──▶ build/<ff>/*.lt + units.json

Typing is graph-based: the full chain molecule is rebuilt from the stored
mapped SMILES (no coordinates needed), typed with priority-based SMARTS
matching, and charges are assigned (Gasteiger on the full chain for GAFF,
.fdefn charge tables otherwise). Every stored variant atom carries its
chain atom-map number, so types/charges transfer by lookup — variant atoms
*are* chain atoms.

The same geometry can be typed under multiple force fields into separate
build/<ff>/ directories (typing is the cheap stage).

Created on 2026-07-30
@author: zwu
"""
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from rdkit import Chem
from rdkit.Chem import AllChem

from .system import logger
from .conf import FORCE_FIELD_REGISTRY
from .exceptions import ValidationError
from .geometry import GEOMETRY_DIRNAME, GeometryBuilder
from .monomer_generator import (
    AtomTypingError,
    SMARTSTyper,
    write_lt_footer,
    write_lt_header,
)
from .monomer_processing import read_lt_end_atoms
from .units import UnitLibrary, UnitSpec

BUILD_DIRNAME = "build"
DEFAULT_MOLECULE_RADIUS = 3.0


def mol_from_mapped_smiles(smiles_mapped: str) -> Chem.Mol:
    """
    Rebuild an RDKit Mol from a mapped SMILES, preserving explicit Hs.

    The default MolFromSmiles strips explicit hydrogens (and their map
    numbers), which would break the map-number join to stored geometry, so
    parsing goes through sanitize=False + explicit property cache update.
    """
    mol = Chem.MolFromSmiles(smiles_mapped, sanitize=False)
    if mol is None:
        return None
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol)
    return mol


class UnitTyper:
    """
    Stage 2: assign force-field types/charges onto stored geometry.

    Example:
        >>> UnitTyper(geom.dir, "oplsaa").type()   # build/oplsaa/
        >>> UnitTyper(geom.dir, "gaff2").type()    # build/gaff2/ — same geometry
    """

    def __init__(
        self,
        geometry_dir: Union[str, Path],
        force_field: str,
        output_dir: Optional[Union[str, Path]] = None,
    ) -> None:
        """
        Args:
            geometry_dir: Directory containing geometry.json (stage 1 output).
            force_field: One of oplsaa, lopls, gaff, gaff2, dreiding, compass.
            output_dir: Build directory (default: <project>/build/<ff>).

        Raises:
            ValidationError: Unknown force field, or unsupported/missing
                             geometry artifact.
        """
        if force_field not in FORCE_FIELD_REGISTRY:
            raise ValidationError(
                f"Invalid force_field '{force_field}'. "
                f"Must be one of: {list(FORCE_FIELD_REGISTRY)}"
            )
        self.geometry_dir = Path(geometry_dir)
        self.geometry = GeometryBuilder.load(self.geometry_dir)
        self.force_field = force_field
        if output_dir is None:
            output_dir = self.geometry_dir.parent / BUILD_DIRNAME / force_field
        self.output_dir = Path(output_dir)
        self.typer = SMARTSTyper(force_field, verbose=False)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def type(self) -> UnitLibrary:
        """
        Type all variants and molecules; write .lt files and units.json.

        Returns:
            The UnitLibrary manifest (also saved to the build directory).

        Raises:
            AtomTypingError: If chain/molecule typing fails.
            ValidationError: If variant map numbers are missing from the
                             typed chain.
        """
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 1. Type every chain graph (whole chain, for correct environments)
        typed_chains = self._type_chain_graphs()

        # 2. Write typed monomer variant .lt files
        variant_files = self._write_variant_files(typed_chains)

        # 3. Write poly_N.lt per multi-monomer chain and collect polymer units
        units, dop1_counts = self._write_chain_files()

        # 4. Type and write small-molecule .lt files; collect molecule units
        molecule_files = self._write_molecule_files(units)

        # 5. DOP=1 chains become molecule-like units grouped by variant
        self._add_single_monomer_units(units, dop1_counts)

        library = UnitLibrary(
            force_field=self.force_field,
            units=units,
            monomer_files=sorted(set(variant_files + molecule_files)),
            build_config={
                "use_mc_chain_growth": self.geometry["mc_config"].get(
                    "use_mc_chain_growth", True
                ),
                "offset": self.geometry["mc_config"].get("offset", 4.0),
                "rotate": self.geometry["mc_config"].get("rotate", 90.0),
            },
            geometry_source=str(
                Path("..") / GEOMETRY_DIRNAME / "geometry.json"
            ),
        )
        library.save(self.output_dir)
        library.validate(self.output_dir)
        library.source_dir = str(self.output_dir)
        logger.info(
            f"Typing complete ({self.force_field}): {len(units)} units, "
            f"{len(library.monomer_files)} monomer files -> {self.output_dir}"
        )
        return library

    # ------------------------------------------------------------------
    # Chain graph typing
    # ------------------------------------------------------------------
    def _type_chain_graphs(self) -> Dict[str, Chem.Mol]:
        """
        Rebuild and type the full chain molecule per polymer model.

        Returns:
            Dict of model_id -> typed RDKit Mol (map numbers intact).
        """
        typed = {}
        for model_id, graph in self.geometry["chain_graphs"].items():
            mol = mol_from_mapped_smiles(graph["smiles_mapped"])
            if mol is None:
                raise AtomTypingError(
                    f"Could not rebuild chain for {model_id} from mapped SMILES"
                )
            try:
                self.typer.assign_atom_types(mol)
            except Exception as e:
                raise AtomTypingError(
                    f"Atom typing failed for {model_id} "
                    f"(force field {self.force_field}): {e}"
                ) from e

            # Gasteiger charges on the FULL chain (before any splitting),
            # matching the legacy behavior for GAFF.
            if self.force_field == "gaff":
                try:
                    AllChem.ComputeGasteigerCharges(mol)
                except Exception as e:
                    logger.error(
                        f"Failed to compute Gasteiger charges for {model_id}: {e}"
                    )
            typed[model_id] = mol
        return typed

    def _chain_lookup(self, chain_mol: Chem.Mol, model_id: str) -> Dict[int, int]:
        """Map number -> atom index lookup for a typed chain."""
        lookup = {}
        for atom in chain_mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num > 0:
                lookup[map_num] = atom.GetIdx()
        if not lookup:
            raise ValidationError(
                f"Chain graph for {model_id} has no atom map numbers; "
                "cannot join types onto geometry"
            )
        return lookup

    # ------------------------------------------------------------------
    # Variant .lt files
    # ------------------------------------------------------------------
    def _write_variant_files(self, typed_chains: Dict[str, Chem.Mol]) -> List[str]:
        """Write one typed .lt file per geometry variant. Returns file names."""
        lookups = {
            model_id: self._chain_lookup(mol, model_id)
            for model_id, mol in typed_chains.items()
        }
        written = []
        for name, variant in self.geometry["variants"].items():
            chain_mol = typed_chains[variant["model_id"]]
            lookup = lookups[variant["model_id"]]
            rows = self._typed_atom_rows(variant, chain_mol, lookup, name)
            self._write_monomer_lt(name, variant, rows)
            written.append(f"{name}.lt")
        return written

    def _typed_atom_rows(
        self,
        variant: Dict[str, Any],
        chain_mol: Chem.Mol,
        lookup: Dict[int, int],
        name: str,
    ) -> List[Dict[str, Any]]:
        """
        Join variant atoms onto the typed chain via map numbers.

        Returns one row per atom: atom_id, atom_type, charge, coords.
        """
        rows = []
        for i, a in enumerate(variant["atoms"]):
            map_num = a["map_num"]
            idx = lookup.get(map_num)
            if idx is None:
                raise ValidationError(
                    f"Variant '{name}' atom {i} references map number "
                    f"{map_num}, which is missing from the typed chain "
                    f"({variant['model_id']})"
                )
            chain_atom = chain_mol.GetAtomWithIdx(idx)
            atom_type = self._atom_type(chain_atom, a["element"])
            charge = self._atom_charge(chain_atom, atom_type)
            rows.append({
                "atom_id": f"{a['element']}{i + 1}",
                "atom_type": atom_type,
                "charge": charge,
                "coords": a["coords"],
            })
        return rows

    def _atom_type(self, atom, element: str) -> str:
        """Atom type property with the legacy element fallback."""
        try:
            return atom.GetProp("AtomType")
        except KeyError:
            return f"@atom:{element.lower()}"

    def _atom_charge(self, atom, atom_type: str) -> float:
        """Per-atom charge: Gasteiger (gaff) or .fdefn charge table."""
        if self.force_field == "gaff":
            try:
                charge = float(atom.GetProp("_GasteigerCharge"))
                if charge != charge or charge in (float("inf"), float("-inf")):
                    return 0.0
                return charge
            except (KeyError, ValueError):
                return 0.0
        return self.typer.charge_dict.get(atom_type, 0.0)

    def _write_monomer_lt(
        self, name: str, variant: Dict[str, Any], rows: List[Dict[str, Any]]
    ) -> None:
        """Write one typed monomer variant .lt file."""
        path = self.output_dir / f"{name}.lt"
        with open(path, "w") as f:
            write_lt_header(f, name, self.force_field)
            f.write('  write("Data Atoms") {\n')
            for row in rows:
                x, y, z = row["coords"]
                f.write(
                    f"\t$atom:{row['atom_id']} $mol:... {row['atom_type']} "
                    f"{row['charge']:.4f}    {x:.3f}   {y:.3f}   {z:.3f}\n"
                )
            f.write("  }\n\n")
            f.write("  write('Data Bond List') {\n")
            for i, j in variant["bonds"]:
                id1 = rows[i]["atom_id"]
                id2 = rows[j]["atom_id"]
                f.write(f"\t$bond:{id1}{id2}\t$atom:{id1}\t$atom:{id2}\n")
            f.write("  }\n")
            write_lt_footer(f, name)

    # ------------------------------------------------------------------
    # poly_N.lt chain files
    # ------------------------------------------------------------------
    def _write_chain_files(self):
        """
        Write poly_N.lt for every multi-monomer chain.

        Returns:
            (units, dop1_counts): polymer UnitSpec list, and a dict of
            variant name -> chain count for DOP=1 chains.
        """
        entry = FORCE_FIELD_REGISTRY[self.force_field]
        units: List[UnitSpec] = []
        dop1_counts: Dict[str, int] = {}

        poly_index = 0
        for chain in self.geometry["chains"]:
            placements = chain["placements"]
            n = len(placements)
            if n <= 1:
                variant_name = placements[0]["variant"]
                dop1_counts[variant_name] = dop1_counts.get(variant_name, 0) + 1
                continue

            poly_index += 1
            poly_id = f"poly_{poly_index}"
            variant_names = [p["variant"] for p in placements]
            monomer_files = sorted({f"{v}.lt" for v in variant_names})

            self._write_poly_lt(
                self.output_dir / f"{poly_id}.lt",
                poly_id, entry, placements,
                ring=(chain["topology"] == "ring"),
            )

            anchors = self._chain_anchors(variant_names, n)
            units.append(UnitSpec(
                id=poly_id,
                kind="polymer",
                lt_file=f"{poly_id}.lt",
                count=1,
                topology=chain["topology"],
                n_monomers=chain["n_monomers"],
                radius=chain["radius"],
                anchors=anchors,
                monomer_files=monomer_files,
            ))

        return units, dop1_counts

    def _write_poly_lt(
        self,
        path: Path,
        poly_id: str,
        ff_entry: Dict[str, str],
        placements: List[Dict[str, Any]],
        ring: bool,
    ) -> None:
        """Write one poly_N.lt from stored placements (verbatim transforms)."""
        n = len(placements)
        variant_names = [p["variant"] for p in placements]

        with open(path, "w") as f:
            f.write(f'import "{ff_entry["lt_file"]}"\n')
            for name in dict.fromkeys(variant_names):
                f.write(f'import "{name}.lt"\n')
            f.write("\n")
            f.write(f"{poly_id} inherits {ff_entry['inherits']} {{\n\n")
            f.write("    create_var {$mol}\n\n")

            for i, p in enumerate(placements):
                f.write(f"    monomer[{i}] = new {p['variant']}")
                if p["rotation"] is not None:
                    angle = p["rotation"]["angle_deg"]
                    ax, ay, az = p["rotation"]["axis"]
                    f.write(f".rot({angle:.4f},{ax:.4f},{ay:.4f},{az:.4f})")
                if p["translation"] is not None:
                    x, y, z = p["translation"]
                    f.write(f".move({x:.4f},{y:.4f},{z:.4f})")
                f.write("\n")

            f.write("\n    write('Data Bond List') {\n")
            bond_index = 0
            for i in range(n if ring else n - 1):
                next_i = (i + 1) % n
                _, second_atom = read_lt_end_atoms(
                    self.output_dir / f"{variant_names[i]}.lt"
                )
                first_atom, _ = read_lt_end_atoms(
                    self.output_dir / f"{variant_names[next_i]}.lt"
                )
                bond_index += 1
                f.write(
                    f"      $bond:b{bond_index}  "
                    f"$atom:monomer[{i}]/{second_atom}  "
                    f"$atom:monomer[{next_i}]/{first_atom}\n"
                )
            f.write("    }\n")

            f.write(f"\n}} # {poly_id}\n")

    def _chain_anchors(
        self, variant_names: List[str], n: int
    ) -> Dict[str, str]:
        """
        Head/tail atom references (reserved for future grafting strategies).
        head = first atom of the first monomer; tail = second atom of the last.
        """
        first_atom, _ = read_lt_end_atoms(
            self.output_dir / f"{variant_names[0]}.lt"
        )
        _, second_atom = read_lt_end_atoms(
            self.output_dir / f"{variant_names[-1]}.lt"
        )
        return {
            "head": f"monomer[0]/{first_atom}",
            "tail": f"monomer[{n - 1}]/{second_atom}",
        }

    # ------------------------------------------------------------------
    # Small molecules
    # ------------------------------------------------------------------
    def _write_molecule_files(self, units: List[UnitSpec]) -> List[str]:
        """Type and write small-molecule .lt files; append molecule units."""
        written = []
        for entry in self.geometry["molecules"]:
            name = entry["name"]
            mol = mol_from_mapped_smiles(entry["smiles_mapped"])
            if mol is None:
                raise AtomTypingError(
                    f"Could not rebuild molecule '{name}' from mapped SMILES"
                )
            try:
                self.typer.assign_atom_types(mol)
            except Exception as e:
                raise AtomTypingError(
                    f"Atom typing failed for molecule '{name}' "
                    f"(force field {self.force_field}): {e}"
                ) from e

            lookup = {
                atom.GetAtomMapNum(): atom.GetIdx()
                for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0
            }
            rows = []
            for i, a in enumerate(entry["atoms"]):
                idx = lookup.get(a["map_num"])
                if idx is None:
                    raise ValidationError(
                        f"Molecule '{name}' atom {i} references missing map "
                        f"number {a['map_num']}"
                    )
                atom = mol.GetAtomWithIdx(idx)
                atom_type = self._atom_type(atom, a["element"])
                # Molecule charges come from the .fdefn charge table for all
                # force fields (legacy behavior; GAFF molecules need
                # user-supplied charges — see the GAFF charges warning).
                charge = self.typer.charge_dict.get(atom_type, 0.0)
                rows.append({
                    "atom_id": f"{a['element']}{i + 1}",
                    "atom_type": atom_type,
                    "charge": charge,
                    "coords": a["coords"],
                })

            self._write_monomer_lt(name, {"bonds": entry["bonds"]}, rows)
            lt_file = f"{name}.lt"
            written.append(lt_file)
            units.append(UnitSpec(
                id=name,
                kind="molecule",
                lt_file=lt_file,
                count=entry["count"],
                radius=DEFAULT_MOLECULE_RADIUS,
                anchors={},
                monomer_files=[lt_file],
            ))
        return written

    def _add_single_monomer_units(
        self, units: List[UnitSpec], dop1_counts: Dict[str, int]
    ) -> None:
        """DOP=1 polymer chains pack as molecule-like units per variant."""
        for variant_name, count in dop1_counts.items():
            lt_file = f"{variant_name}.lt"
            units.append(UnitSpec(
                id=variant_name,
                kind="molecule",
                lt_file=lt_file,
                count=count,
                radius=DEFAULT_MOLECULE_RADIUS,
                anchors={},
                monomer_files=[lt_file],
            ))
