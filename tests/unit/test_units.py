"""Tests for the units.json manifest contract (AutoPoly.units)."""

import json

import pytest

from AutoPoly.exceptions import ValidationError
from AutoPoly.units import (
    FORMAT_VERSION,
    MANIFEST_FILENAME,
    UnitLibrary,
    UnitSpec,
)


def _polymer_unit(**overrides):
    data = dict(
        id="poly_1",
        kind="polymer",
        lt_file="poly_1.lt",
        count=1,
        topology="linear",
        n_monomers=50,
        radius=18.2,
        anchors={"head": "monomer[0]/C1", "tail": "monomer[49]/C2"},
        monomer_files=["monomer_0_0le.lt", "monomer_0_1i.lt"],
    )
    data.update(overrides)
    return UnitSpec(**data)


class TestUnitSpec:
    def test_defaults(self):
        unit = UnitSpec(id="water", kind="molecule", lt_file="water.lt", count=100)
        assert unit.radius == 3.0
        assert unit.anchors == {}
        assert unit.monomer_files == []

    def test_class_name_is_lt_stem(self):
        assert _polymer_unit().class_name == "poly_1"
        assert UnitSpec(id="w", kind="molecule", lt_file="water.lt").class_name == "water"

    def test_invalid_kind_raises(self):
        with pytest.raises(ValidationError, match="Invalid unit kind"):
            UnitSpec(id="x", kind="substrate", lt_file="x.lt")

    def test_invalid_count_raises(self):
        with pytest.raises(ValidationError, match="invalid count"):
            UnitSpec(id="x", kind="molecule", lt_file="x.lt", count=0)

    def test_round_trip(self):
        unit = _polymer_unit()
        assert UnitSpec.from_dict(unit.to_dict()) == unit


class TestUnitLibrary:
    def _library(self):
        return UnitLibrary(
            force_field="oplsaa",
            units=[
                _polymer_unit(),
                UnitSpec(
                    id="water", kind="molecule", lt_file="water.lt",
                    count=100, monomer_files=["water.lt"],
                ),
            ],
            monomer_files=["monomer_0_0le.lt", "monomer_0_1i.lt", "water.lt"],
            build_config={"use_mc_chain_growth": True},
            geometry_source="../geometry/geometry.json",
        )

    def _write_lt_files(self, directory, names):
        for name in names:
            (directory / name).write_text("# fake lt\n")

    def test_save_load_round_trip(self, tmp_path):
        library = self._library()
        path = library.save(tmp_path)
        assert path == tmp_path / MANIFEST_FILENAME

        loaded = UnitLibrary.load(tmp_path)
        assert loaded.force_field == "oplsaa"
        assert loaded.units == library.units
        assert loaded.monomer_files == library.monomer_files
        assert loaded.build_config == library.build_config
        assert loaded.format_version == FORMAT_VERSION

    def test_load_from_file_path(self, tmp_path):
        library = self._library()
        path = library.save(tmp_path)
        loaded = UnitLibrary.load(path)
        assert loaded.units == library.units

    def test_load_missing_manifest_raises(self, tmp_path):
        with pytest.raises(ValidationError, match="not found"):
            UnitLibrary.load(tmp_path)

    def test_load_bad_json_raises(self, tmp_path):
        (tmp_path / MANIFEST_FILENAME).write_text("{not json")
        with pytest.raises(ValidationError, match="Invalid units.json"):
            UnitLibrary.load(tmp_path)

    def test_load_unsupported_version_raises(self, tmp_path):
        data = self._library().to_dict()
        data["format_version"] = FORMAT_VERSION + 1
        (tmp_path / MANIFEST_FILENAME).write_text(json.dumps(data))
        with pytest.raises(ValidationError, match="format_version"):
            UnitLibrary.load(tmp_path)

    def test_validate_ok(self, tmp_path):
        library = self._library()
        self._write_lt_files(
            tmp_path,
            ["poly_1.lt", "monomer_0_0le.lt", "monomer_0_1i.lt", "water.lt"],
        )
        library.validate(tmp_path)  # no raise

    def test_validate_lists_missing_files(self, tmp_path):
        library = self._library()
        self._write_lt_files(tmp_path, ["poly_1.lt", "water.lt"])
        with pytest.raises(ValidationError, match="monomer_0_0le.lt"):
            library.validate(tmp_path)

    def test_convenience_views(self):
        library = self._library()
        assert [u.id for u in library.polymer_units] == ["poly_1"]
        assert [u.id for u in library.molecule_units] == ["water"]
        # 50 monomers + 100 waters
        assert library.total_particle_count() == 150
        assert library.max_chain_length() == 50

    def test_empty_library_views(self):
        library = UnitLibrary(force_field="gaff")
        assert library.total_particle_count() == 1  # clamped to >= 1
        assert library.max_chain_length() == 1
