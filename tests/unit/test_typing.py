"""Tests for UnitTyper (stage 2) and multi-FF typing from one geometry."""

import json
from pathlib import Path

import pytest
from rdkit import Chem

from AutoPoly.exceptions import ValidationError
from AutoPoly.geometry import GeometryBuilder, GeometryConfig
from AutoPoly.molecule import Molecule
from AutoPoly.monomer_generator import SMARTSTyper
from AutoPoly.polymer import Polymer
from AutoPoly.system import System
from AutoPoly.typing import UnitTyper
from AutoPoly.units import UnitLibrary

PE_SEQUENCE = ["CC[*]", "[*]CC[*]", "[*]CC"]


def _build_geometry(tmp_path, models=None, name="pe"):
    system = System(out=str(tmp_path / "out"))
    if models is None:
        models = [Polymer(chain_num=2, sequence=PE_SEQUENCE, tacticity="syndiotactic")]
    config = GeometryConfig(use_mc_chain_growth=False)
    return GeometryBuilder(system, name, config).build(models)


def _lt_atom_types(lt_path):
    """Parse @atom types out of a .lt file's Data Atoms block."""
    types = []
    in_block = False
    for line in Path(lt_path).read_text().splitlines():
        stripped = line.strip()
        if stripped == 'write("Data Atoms") {':
            in_block = True
            continue
        if in_block and stripped == "}":
            break
        if in_block and stripped:
            parts = stripped.split()
            types.append((parts[0].split(":")[1], parts[2], float(parts[3])))
    return types


class TestTyping:
    def test_build_dir_layout_and_manifest(self, tmp_path):
        geom = _build_geometry(tmp_path)
        library = UnitTyper(geom.dir, "oplsaa").type()

        build_dir = Path(geom.dir).parent / "build" / "oplsaa"
        assert (build_dir / "units.json").is_file()
        # 2 chains -> 2 polymer units
        polymers = library.polymer_units
        assert [u.id for u in polymers] == ["poly_1", "poly_2"]
        for unit in polymers:
            assert unit.topology == "linear"
            assert unit.n_monomers == 3
            assert unit.radius > 0
            assert (build_dir / unit.lt_file).is_file()
            assert unit.anchors["head"].startswith("monomer[0]/")
            assert unit.anchors["tail"].startswith("monomer[2]/")
        # Manifest validates on disk
        library.validate(build_dir)

    def test_chain_level_typing_matches_direct_typer(self, tmp_path):
        """Types written into variant .lt files must match SMARTSTyper applied
        to the same whole chain (chain-level typing, joined by map numbers)."""
        geom = _build_geometry(tmp_path)
        UnitTyper(geom.dir, "oplsaa").type()

        data = json.loads((Path(geom.dir) / "geometry.json").read_text())
        chain_mol = Chem.MolFromSmiles(
            data["chain_graphs"]["model_0"]["smiles_mapped"], sanitize=False
        )
        chain_mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(chain_mol)
        typer = SMARTSTyper("oplsaa", verbose=False)
        typer.assign_atom_types(chain_mol)
        ref_types = {
            a.GetAtomMapNum(): a.GetProp("AtomType")
            for a in chain_mol.GetAtoms() if a.GetAtomMapNum() > 0
        }

        build_dir = Path(geom.dir).parent / "build" / "oplsaa"
        variant = data["variants"]["monomer_0_1i"]
        lt_rows = _lt_atom_types(build_dir / "monomer_0_1i.lt")
        assert len(lt_rows) == len(variant["atoms"])
        for row, atom in zip(lt_rows, variant["atoms"]):
            assert row[1] == ref_types[atom["map_num"]]

    def test_gaff_gasteiger_charges(self, tmp_path):
        """GAFF typing writes Gasteiger charges computed on the full chain."""
        geom = _build_geometry(tmp_path)
        UnitTyper(geom.dir, "gaff").type()

        build_dir = Path(geom.dir).parent / "build" / "gaff"
        rows = _lt_atom_types(build_dir / "monomer_0_1i.lt")
        charges = [r[2] for r in rows]
        # Gasteiger charges on a real chain are non-trivial and roughly balanced
        assert any(abs(q) > 1e-6 for q in charges)
        assert sum(charges) == pytest.approx(0.0, abs=0.05)

    def test_oplsaa_charges_from_charge_dict(self, tmp_path):
        geom = _build_geometry(tmp_path)
        UnitTyper(geom.dir, "oplsaa").type()
        build_dir = Path(geom.dir).parent / "build" / "oplsaa"
        rows = _lt_atom_types(build_dir / "monomer_0_2re.lt")
        # OPLS-AA alkane charges from the .fdefn table (methyl C is non-zero)
        assert any(abs(q) > 1e-6 for _, _, q in rows)

    def test_multi_ff_typing_from_one_geometry(self, tmp_path):
        """One geometry can be typed under multiple force fields."""
        geom = _build_geometry(tmp_path)
        lib_oplsaa = UnitTyper(geom.dir, "oplsaa").type()
        lib_gaff2 = UnitTyper(geom.dir, "gaff2").type()

        base = Path(geom.dir).parent / "build"
        assert (base / "oplsaa" / "units.json").is_file()
        assert (base / "gaff2" / "units.json").is_file()
        # Same unit structure, different typing
        assert [u.id for u in lib_oplsaa.units] == [u.id for u in lib_gaff2.units]
        types_o = {r[1] for r in _lt_atom_types(base / "oplsaa" / "monomer_0_1i.lt")}
        types_g = {r[1] for r in _lt_atom_types(base / "gaff2" / "monomer_0_1i.lt")}
        assert types_o != types_g  # different force fields -> different types

    def test_poly_lt_contains_stored_transforms(self, tmp_path):
        geom = _build_geometry(tmp_path)
        UnitTyper(geom.dir, "oplsaa").type()
        build_dir = Path(geom.dir).parent / "build" / "oplsaa"
        poly_lt = (build_dir / "poly_1.lt").read_text()

        assert 'import "oplsaa.lt"' in poly_lt
        assert "poly_1 inherits OPLSAA {" in poly_lt
        assert "monomer[0] = new monomer_0_0le" in poly_lt
        assert "write('Data Bond List')" in poly_lt
        # Deterministic placement: monomer 1 rotated 90 deg about x
        assert ".rot(90.0000,1.0000,0.0000,0.0000)" in poly_lt

    def test_molecule_typing(self, tmp_path):
        system = System(out=str(tmp_path / "out"))
        water = Molecule(Count=4, Smiles="O", Name="water")
        geom = GeometryBuilder(
            system, "sol", GeometryConfig(use_mc_chain_growth=False)
        ).build([water])
        library = UnitTyper(geom.dir, "oplsaa").type()

        water_units = [u for u in library.units if u.id == "water"]
        assert len(water_units) == 1
        assert water_units[0].kind == "molecule"
        assert water_units[0].count == 4
        build_dir = Path(geom.dir).parent / "build" / "oplsaa"
        assert (build_dir / "water.lt").is_file()

    def test_invalid_force_field_raises(self, tmp_path):
        geom = _build_geometry(tmp_path)
        with pytest.raises(ValidationError, match="Invalid force_field"):
            UnitTyper(geom.dir, "charmm")

    def test_missing_geometry_raises(self, tmp_path):
        with pytest.raises(Exception, match="not found"):
            UnitTyper(str(tmp_path / "nope"), "oplsaa")

    def test_manifest_round_trip_from_disk(self, tmp_path):
        geom = _build_geometry(tmp_path)
        UnitTyper(geom.dir, "oplsaa").type()
        build_dir = Path(geom.dir).parent / "build" / "oplsaa"
        loaded = UnitLibrary.load(build_dir)
        assert loaded.force_field == "oplsaa"
        assert len(loaded.polymer_units) == 2
        loaded.validate(build_dir)
