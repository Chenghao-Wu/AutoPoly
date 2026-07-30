"""Tests for GeometryBuilder (stage 1) and geometry.json."""

import json

import numpy as np
import pytest

from AutoPoly.exceptions import GenerationError, ValidationError
from AutoPoly.geometry import (
    GEOMETRY_FILENAME,
    GeometryBuilder,
    GeometryConfig,
)
from AutoPoly.molecule import Molecule
from AutoPoly.polymer import Polymer
from AutoPoly.system import System

PE_SEQUENCE = ["CC[*]", "[*]CC[*]", "[*]CC"]
ABA_SEQUENCE = [
    "CC[*]",
    "[*]CC[*]",
    "[*]CC(c1ccccc1)[*]",
    "[*]CC(c1ccccc1)[*]",
    "[*]CC(c1ccccc1)[*]",
    "[*]CC[*]",
    "[*]CC",
]


def _det_config(**overrides):
    """Deterministic (non-MC) geometry config."""
    return GeometryConfig(use_mc_chain_growth=False, **overrides)


class TestLinearGeometry:
    def test_geometry_json_schema(self, tmp_path):
        system = System(out=str(tmp_path / "out"))
        poly = Polymer(chain_num=2, sequence=PE_SEQUENCE, tacticity="atactic")
        result = GeometryBuilder(system, "pe", _det_config()).build([poly])

        assert result.json_path.is_file()
        data = json.loads(result.json_path.read_text())

        assert data["format_version"] == 1
        assert "autopoly_version" in data["created_with"]
        assert data["mc_config"]["use_mc_chain_growth"] is False

        # One chain graph per polymer model
        assert set(data["chain_graphs"]) == {"model_0"}
        graph = data["chain_graphs"]["model_0"]
        assert graph["atom_count"] > 0
        assert graph["smiles_mapped"]

        # first/middle/last variants + their _T1 mirrors
        variant_names = set(data["variants"])
        assert "monomer_0_0le" in variant_names
        assert "monomer_0_1i" in variant_names
        assert "monomer_0_2re" in variant_names
        assert "monomer_0_0le_T1" in variant_names
        assert data["variants"]["monomer_0_0le_T1"]["t1_variant_of"] == "monomer_0_0le"
        assert data["variants"]["monomer_0_0le"]["t1_variant_of"] is None

        # Two chains, one per model chain, with per-position placements
        assert len(data["chains"]) == 2
        for chain in data["chains"]:
            assert chain["topology"] == "linear"
            assert chain["n_monomers"] == 3
            assert len(chain["placements"]) == 3
            assert chain["radius"] > 0
            for placement in chain["placements"]:
                assert placement["variant"] in variant_names

    def test_variant_map_numbers_are_chain_atoms(self, tmp_path):
        """Every variant atom map number must exist in the chain graph."""
        from AutoPoly.typing import mol_from_mapped_smiles

        system = System(out=str(tmp_path / "out"))
        poly = Polymer(chain_num=1, sequence=ABA_SEQUENCE, tacticity="atactic")
        result = GeometryBuilder(system, "aba", _det_config()).build([poly])
        data = json.loads(result.json_path.read_text())

        chain_mol = mol_from_mapped_smiles(
            data["chain_graphs"]["model_0"]["smiles_mapped"]
        )
        chain_maps = {
            a.GetAtomMapNum() for a in chain_mol.GetAtoms() if a.GetAtomMapNum() > 0
        }
        assert len(chain_maps) == data["chain_graphs"]["model_0"]["atom_count"]

        for name, variant in data["variants"].items():
            for atom in variant["atoms"]:
                assert atom["map_num"] in chain_maps, (
                    f"{name}: map {atom['map_num']} not in chain graph"
                )
            # Within one variant, map numbers are unique
            maps = [a["map_num"] for a in variant["atoms"]]
            assert len(maps) == len(set(maps))

    def test_copolymer_keeps_distinct_middle_variants(self, tmp_path):
        """Bug B regression: chemically distinct 'middle' monomers must not
        collapse into one variant file."""
        system = System(out=str(tmp_path / "out"))
        poly = Polymer(chain_num=1, sequence=ABA_SEQUENCE, tacticity="atactic")
        result = GeometryBuilder(system, "aba", _det_config()).build([poly])
        data = json.loads(result.json_path.read_text())

        middle_variants = [
            name for name, v in data["variants"].items()
            if v["variant_type"] == "middle" and v["t1_variant_of"] is None
        ]
        # PE middle + PS middle = 2 distinct middle chemistries
        assert len(middle_variants) == 2

    def test_deterministic_seed_reproducibility(self, tmp_path):
        """Deterministic mode: two builds produce identical geometry.json."""
        outputs = []
        for run in (1, 2):
            system = System(out=str(tmp_path / f"run{run}"))
            poly = Polymer(chain_num=1, sequence=PE_SEQUENCE, tacticity="syndiotactic")
            result = GeometryBuilder(system, "pe", _det_config()).build([poly])
            outputs.append(result.json_path.read_text())
        assert outputs[0] == outputs[1]

    def test_mc_seed_reproducibility(self, tmp_path):
        """MC chain growth with a fixed seed produces identical placements."""
        placements = []
        for run in (1, 2):
            system = System(out=str(tmp_path / f"run{run}"))
            poly = Polymer(chain_num=1, sequence=PE_SEQUENCE, tacticity="syndiotactic")
            config = GeometryConfig(use_mc_chain_growth=True, rng_seed=1234)
            result = GeometryBuilder(system, "pe", config).build([poly])
            data = json.loads(result.json_path.read_text())
            placements.append(data["chains"][0]["placements"])
        assert placements[0] == placements[1]

    def test_mc_placements_have_transforms(self, tmp_path):
        system = System(out=str(tmp_path / "out"))
        poly = Polymer(chain_num=1, sequence=PE_SEQUENCE, tacticity="syndiotactic")
        config = GeometryConfig(use_mc_chain_growth=True, rng_seed=7)
        result = GeometryBuilder(system, "pe", config).build([poly])
        data = json.loads(result.json_path.read_text())

        placements = data["chains"][0]["placements"]
        assert len(placements) == 3
        for p in placements:
            # MC-grown monomers always carry a translation
            assert p["translation"] is not None


class TestRingGeometry:
    def test_ring_placements(self, tmp_path):
        system = System(out=str(tmp_path / "out"))
        poly = Polymer(chain_num=1, sequence=PE_SEQUENCE, topology="ring")
        result = GeometryBuilder(system, "ring", _det_config(offset=4.0)).build([poly])
        data = json.loads(result.json_path.read_text())

        chain = data["chains"][0]
        assert chain["topology"] == "ring"

        # Ring: monomers on a circle of radius offset * n / (2*pi)
        n = len(chain["placements"])
        expected_radius = 4.0 * n / (2 * np.pi)
        for p in chain["placements"]:
            x, y, z = p["translation"]
            assert np.hypot(x, y) == pytest.approx(expected_radius, abs=1e-6)
            assert z == 0.0
            assert p["rotation"] is not None

    def test_ring_uses_middle_variants(self, tmp_path):
        system = System(out=str(tmp_path / "out"))
        poly = Polymer(chain_num=1, sequence=PE_SEQUENCE, topology="ring")
        result = GeometryBuilder(system, "ring", _det_config()).build([poly])
        data = json.loads(result.json_path.read_text())

        for p in data["chains"][0]["placements"]:
            name = p["variant"].replace("_T1", "")
            # Ring positions use middle-type variants (internal 'i' naming)
            assert data["variants"][name]["variant_type"] in ("ring", "middle")
            assert "i" in name.rsplit("_", 1)[-1]


class TestMoleculeGeometry:
    def test_molecule_entry(self, tmp_path):
        system = System(out=str(tmp_path / "out"))
        water = Molecule(Count=5, Smiles="O", Name="water")
        result = GeometryBuilder(system, "sol", _det_config()).build([water])
        data = json.loads(result.json_path.read_text())

        assert len(data["molecules"]) == 1
        entry = data["molecules"][0]
        assert entry["name"] == "water"
        assert entry["count"] == 5
        assert len(entry["atoms"]) == 3  # O + 2 H
        assert len(entry["bonds"]) == 2
        # Map numbers are 1..N join keys for typing
        assert sorted(a["map_num"] for a in entry["atoms"]) == [1, 2, 3]

    def test_mixed_polymer_and_molecule(self, tmp_path):
        system = System(out=str(tmp_path / "out"))
        poly = Polymer(chain_num=1, sequence=PE_SEQUENCE)
        water = Molecule(Count=2, Smiles="O", Name="water")
        result = GeometryBuilder(system, "mix", _det_config()).build([poly, water])
        data = json.loads(result.json_path.read_text())
        assert len(data["chains"]) == 1
        assert len(data["molecules"]) == 1


class TestSingleMonomer:
    def test_dop1_chain_has_bare_placement(self, tmp_path):
        system = System(out=str(tmp_path / "out"))
        poly = Polymer(chain_num=3, sequence=["CCCC"])  # butane, no wildcards
        result = GeometryBuilder(system, "single", _det_config()).build([poly])
        data = json.loads(result.json_path.read_text())

        assert len(data["chains"]) == 3
        for chain in data["chains"]:
            assert chain["n_monomers"] == 1
            placement = chain["placements"][0]
            assert placement["rotation"] is None
            assert placement["translation"] is None
            assert "single" in placement["variant"]


class TestGeometryErrors:
    def test_invalid_smiles_raises(self, tmp_path):
        system = System(out=str(tmp_path / "out"))
        poly = Polymer(chain_num=1, sequence=["CC[*]", "[*]CC"])  # too short ok
        poly.sequence = ["CC[*]", "not_a_smiles[*]"]
        with pytest.raises(GenerationError):
            GeometryBuilder(system, "bad", _det_config()).build([poly])

    def test_load_missing_geometry(self, tmp_path):
        with pytest.raises(ValidationError, match="not found"):
            GeometryBuilder.load(str(tmp_path))

    def test_load_unsupported_version(self, tmp_path):
        (tmp_path / GEOMETRY_FILENAME).write_text(
            json.dumps({"format_version": 999})
        )
        with pytest.raises(ValidationError, match="format_version"):
            GeometryBuilder.load(str(tmp_path))
