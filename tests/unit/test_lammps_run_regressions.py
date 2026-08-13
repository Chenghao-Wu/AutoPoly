#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Regression tests for issues exposed by running the built
film-on-substrate systems through LAMMPS:

1. random_mc fallback grid must not silently stack instances (it used to
   wrap the index, producing identical positions and atom clashes).
2. system.in.init hybrid styles with zero records in system.data must be
   patched to 'none' (LAMMPS aborts on unused hybrid sub-styles).
3. Substrate runs must emit a ready-to-run in.run (frozen slab, film
   temperature compute).
"""

import pytest

from AutoPoly import System, Polymer, generate
from AutoPoly.core.exceptions import ValidationError
from AutoPoly.packing import SubstrateSpec
from AutoPoly.packing.random_mc import RandomMCStrategy
from AutoPoly.pipeline.geometry import GeometryConfig

PE_SEQUENCE = ["CC[*]", "[*]CC[*]", "[*]CC"]


class _Unit:
    def __init__(self, radius=12.0):
        self.id = "poly_1"
        self.class_name = "poly_1"
        self.radius = radius
        self.count = 1
        self.kind = "polymer"
        self.role = "film"


class TestRandomMCFallback:
    def test_fallback_raises_when_grid_full(self):
        unit = _Unit(radius=12.0)
        with pytest.raises(ValidationError, match="fallback grid is full"):
            # half_box 20: grid_per_dim = int((40-24)/30) = 1 -> capacity 1
            RandomMCStrategy._fallback_item(unit, "polymer_2", 1, 20.0)

    def test_fallback_ok_within_capacity(self):
        unit = _Unit(radius=12.0)
        item = RandomMCStrategy._fallback_item(unit, "polymer_1", 0, 20.0)
        assert item.instance_name == "polymer_1"


@pytest.mark.integration
class TestRunScriptAndStylePatching:
    def test_run_script_and_improper_patch(self, tmp_path):
        system = System(out=str(tmp_path / "out"))
        film = Polymer(chain_num=1, sequence=PE_SEQUENCE)
        spec = SubstrateSpec(builder="alpha_quartz", thickness=12.0,
                             gap=3.0)
        generate(
            system, "pe_q", [film], force_field="gaff",
            substrate=spec,
            box_dims=(30.0, 34.0, 50.0),
            geometry_config=GeometryConfig(use_mc_chain_growth=False),
            rng_seed=3,
        )
        proj = tmp_path / "out" / "pe_q"

        # in.run generated with frozen slab + film temperature compute
        run = (proj / "in.run").read_text()
        assert "group           slab molecule 1" in run
        assert "compute         tfilm film temp" in run
        assert "fix             freeze slab setforce" in run
        assert "thermo_modify   temp tfilm" in run

        # PE has no impropers: hybrid cvff patched to none
        init = (proj / "system.in.init").read_text()
        assert "improper_style  none" in init
        assert "improper_style  hybrid cvff" not in init
        # bonds/angles/dihedrals exist for PE: styles untouched
        assert "bond_style      hybrid harmonic" in init
