#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""End-to-end test: polymer film on a built beta-cristobalite(111) slab."""

import pytest

from AutoPoly import System, Polymer, generate
from AutoPoly.packing import SubstrateSpec
from AutoPoly.pipeline.geometry import GeometryConfig

PE_SEQUENCE = ["CC[*]", "[*]CC[*]", "[*]CC"]

LX = 4 * 10.125769
LY = 3 * 17.538347


def _boundary(content):
    values = {}
    for line in content.splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[2:] in (["xlo", "xhi"],
                                             ["ylo", "yhi"], ["zlo", "zhi"]):
            values[parts[2][0]] = (float(parts[0]), float(parts[1]))
    return values


@pytest.mark.integration
class TestCristobaliteSubstrateEndToEnd:
    def test_pe_film_on_cristobalite(self, tmp_path):
        system = System(out=str(tmp_path / "out"))
        film = Polymer(chain_num=2, sequence=PE_SEQUENCE,
                       tacticity="syndiotactic")
        spec = SubstrateSpec(
            builder="beta_cristobalite",
            thickness=13.0, gap=3.0,
            slab_ff="clayff",
        )
        result = generate(
            system, "pe_on_crist", [film], force_field="gaff",
            substrate=spec,
            box_dims=(40.0, 53.0, 40.0),
            geometry_config=GeometryConfig(use_mc_chain_growth=False),
            rng_seed=3,
        )

        names = [r.instance_name for r in result.records]
        assert sum(n == "substrate" for n in names) == 1
        assert sum(n.startswith("polymer_") for n in names) == 2

        data = (tmp_path / "out" / "pe_on_crist" / "system.data").read_text()
        bounds = _boundary(data)
        assert bounds["x"] == pytest.approx((-LX / 2, LX / 2), abs=1e-3)
        assert bounds["y"] == pytest.approx((-LY / 2, LY / 2), abs=1e-3)

        moltemplate_dir = tmp_path / "out" / "pe_on_crist" / "moltemplate"
        system_lt = (moltemplate_dir / "input" / "system.lt").read_text()
        assert 'import "beta_cristobalite_slab.lt"' in system_lt
        assert "substrate = new CristobaliteSlab" in system_lt

        slab_lt = (moltemplate_dir / "input"
                   / "beta_cristobalite_slab.lt").read_text()
        assert "CristobaliteSlab {" in slab_lt
        assert "@atom:cff_st" in slab_lt  # clayff selected
