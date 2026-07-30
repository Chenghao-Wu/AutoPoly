#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integration tests for the three-stage pipeline (Geometry -> Typing -> Packing).

These exercise the stage-level API end to end (including moltemplate runs)
and complement the facade-level regression tests in
test_copolymer_end_to_end.py:

- seeded runs must be fully deterministic (byte-identical system.data)
- ring topologies must produce closed chains with valid bonds
- mixed polymer + small-molecule systems must pack and type correctly
- one geometry must be reusable across multiple force fields
"""

import hashlib
from pathlib import Path

import numpy as np
import pytest

from AutoPoly import System, Polymer, Molecule
from AutoPoly.geometry import GeometryBuilder, GeometryConfig
from AutoPoly.packer import BoxPacker
from AutoPoly.typing import UnitTyper

PE_SEQUENCE = ["CC[*]", "[*]CC[*]", "[*]CC"]


def _parse_counts(data_path):
    """Extract atom/bond counts and box bounds from a LAMMPS data file."""
    counts = {}
    bounds = {}
    for line in Path(data_path).read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[1] in ("atoms", "bonds", "angles", "dihedrals"):
            counts[parts[1]] = int(parts[0])
        if len(parts) >= 4 and parts[-1] in ("xlo", "ylo", "zlo"):
            bounds[parts[-1]] = (float(parts[0]), float(parts[1]))
    return counts, bounds


def _run_pipeline(out_dir, models, ff="oplsaa", geom_seed=99, pack_seed=77,
                  strategy="mc_random"):
    """Run the three stages with fixed seeds; return the project dir."""
    system = System(out=str(out_dir))
    geom = GeometryBuilder(
        system, "proj",
        GeometryConfig(use_mc_chain_growth=True, rng_seed=geom_seed),
    ).build(models)
    units = UnitTyper(geom.dir, ff).type()
    BoxPacker(system, "proj", strategy=strategy, rng_seed=pack_seed).pack(units)
    return Path(out_dir) / "proj"


@pytest.mark.integration
class TestPipelineStages:
    def test_seeded_pe_run_is_deterministic(self, tmp_path):
        """Same seeds -> byte-identical system.data (full-pipeline determinism)."""
        def make():
            return Polymer(chain_num=2, sequence=PE_SEQUENCE,
                           tacticity="syndiotactic")

        proj1 = _run_pipeline(tmp_path / "r1", [make()])
        proj2 = _run_pipeline(tmp_path / "r2", [make()])

        data1 = (proj1 / "system.data").read_bytes()
        data2 = (proj2 / "system.data").read_bytes()
        assert hashlib.md5(data1).hexdigest() == hashlib.md5(data2).hexdigest()

        # Structural invariants: 2 chains x (2C + 6H) x 3 monomers
        counts, bounds = _parse_counts(proj1 / "system.data")
        # first/middle/last = 7/6/7 atoms per monomer -> 20 per chain
        assert counts["atoms"] == 2 * 20
        # intra-monomer bonds: 7-1 + 6-1 + 7-1 = 17, plus 2 inter-monomer
        assert counts["bonds"] == 2 * (17 + 2)
        for lo, hi in bounds.values():
            assert lo == pytest.approx(-hi)

    def test_aba_triblock_stage_level(self, tmp_path):
        """ABA triblock through the stage API: phenyls survive (Bug B guard)."""
        aba = [
            "CC[*]", "[*]CC[*]",
            "[*]CC(c1ccccc1)[*]", "[*]CC(c1ccccc1)[*]", "[*]CC(c1ccccc1)[*]",
            "[*]CC[*]", "[*]CC",
        ]
        poly = Polymer(chain_num=2, sequence=aba, tacticity="atactic")
        proj = _run_pipeline(tmp_path / "aba", [poly])

        counts, _ = _parse_counts(proj / "system.data")
        # 74 atoms per chain (see test_copolymer_end_to_end)
        assert counts["atoms"] == 2 * 74

        # Intermediate artifacts exist with the new layout
        assert (proj / "geometry" / "geometry.json").is_file()
        assert (proj / "build" / "oplsaa" / "units.json").is_file()
        assert (proj / "moltemplate" / "input").is_dir()  # post-processed

    def test_ring_polymer_end_to_end(self, tmp_path):
        """Ring topology: closure bond present, valid moltemplate output."""
        poly = Polymer(chain_num=2, sequence=PE_SEQUENCE, topology="ring")
        proj = _run_pipeline(tmp_path / "ring", [poly])

        counts, _ = _parse_counts(proj / "system.data")
        # Ring uses middle variants for all positions: 6 atoms per monomer
        assert counts["atoms"] == 2 * 3 * 6
        # Per chain: 3 * 5 intra + 3 inter (incl. closure) = 18 bonds
        assert counts["bonds"] == 2 * 18

        # Anchors reference the first/last monomer connection atoms
        import json
        units = json.loads(
            (proj / "build" / "oplsaa" / "units.json").read_text()
        )
        poly_units = [u for u in units["units"] if u["kind"] == "polymer"]
        assert len(poly_units) == 2
        for unit in poly_units:
            assert unit["topology"] == "ring"
            assert unit["anchors"]["head"].startswith("monomer[0]/")
            assert unit["anchors"]["tail"].startswith("monomer[2]/")

    def test_solution_system_gaff(self, tmp_path):
        """Polymer + small molecule under GAFF (Gasteiger chain charges)."""
        poly = Polymer(chain_num=2, sequence=PE_SEQUENCE, tacticity="syndiotactic")
        ethanol = Molecule(Count=3, Smiles="CCO", Name="ethanol")
        proj = _run_pipeline(tmp_path / "sol", [poly, ethanol], ff="gaff")

        counts, _ = _parse_counts(proj / "system.data")
        # 2 chains * 20 atoms + 3 ethanol * 9 atoms
        assert counts["atoms"] == 2 * 20 + 3 * 9

    def test_multi_ff_from_one_geometry(self, tmp_path):
        """One geometry typed under two force fields packs identically in size."""
        system = System(out=str(tmp_path / "multi"))
        poly = Polymer(chain_num=2, sequence=PE_SEQUENCE, tacticity="syndiotactic")
        geom = GeometryBuilder(
            system, "multi",
            GeometryConfig(use_mc_chain_growth=True, rng_seed=99),
        ).build([poly])

        for ff in ("oplsaa", "gaff2"):
            units = UnitTyper(geom.dir, ff).type()
            packer = BoxPacker(system, "multi", strategy="grid",
                               run_moltemplate=False)
            packer.pack(units)
            moltemplate_dir = Path(geom.dir).parent / "moltemplate"
            assert (moltemplate_dir / "system.lt").is_file()
            # Per-FF typed files live in separate build dirs
            build_dir = Path(geom.dir).parent / "build" / ff
            assert (build_dir / "poly_1.lt").is_file()
            assert (build_dir / "units.json").is_file()

    def test_dop1_single_monomer_end_to_end(self, tmp_path):
        """DOP=1 polymer chains pack as molecule-like units."""
        butane = Polymer(chain_num=3, sequence=["CCCC"])
        proj = _run_pipeline(tmp_path / "dop1", [butane])

        counts, _ = _parse_counts(proj / "system.data")
        assert counts["atoms"] == 3 * 14  # C4H10 per chain
