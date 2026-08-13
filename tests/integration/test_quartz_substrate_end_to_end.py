#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""End-to-end test: polymer film on a built alpha-quartz(0001) slab.

Runs the full pipeline including moltemplate for a small PE film on a
built-in hydroxylated quartz substrate, and checks the generated LAMMPS
data file (box snapped to quartz surface cells, slab instantiated once).
"""

import pytest

from AutoPoly import System, Polymer, generate
from AutoPoly.packing import SubstrateSpec
from AutoPoly.pipeline.geometry import GeometryConfig

PE_SEQUENCE = ["CC[*]", "[*]CC[*]", "[*]CC"]

#: quartz surface cell (a, a*sqrt(3)) times (nx, ny)
LX = 6 * 4.9019
LY = 4 * 4.9019 * (3.0 ** 0.5)


def _boundary(content):
    values = {}
    for line in content.splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[2:] in (["xlo", "xhi"],
                                             ["ylo", "yhi"], ["zlo", "zhi"]):
            values[parts[2][0]] = (float(parts[0]), float(parts[1]))
    return values


@pytest.mark.integration
class TestQuartzSubstrateEndToEnd:
    def test_pe_film_on_quartz(self, tmp_path):
        system = System(out=str(tmp_path / "out"))
        film = Polymer(chain_num=2, sequence=PE_SEQUENCE,
                       tacticity="syndiotactic")
        spec = SubstrateSpec(
            builder="alpha_quartz",
            thickness=12.0, gap=3.0,
            slab_ff="interface",
        )
        result = generate(
            system, "pe_on_quartz", [film], force_field="gaff",
            substrate=spec,
            box_dims=(30.0, 34.0, 40.0),  # lateral snapped to cells
            geometry_config=GeometryConfig(use_mc_chain_growth=False),
            rng_seed=3,
        )

        names = [r.instance_name for r in result.records]
        assert sum(n == "substrate" for n in names) == 1
        assert sum(n.startswith("polymer_") for n in names) == 2

        data = (tmp_path / "out" / "pe_on_quartz" / "system.data").read_text()
        bounds = _boundary(data)
        assert bounds["x"] == pytest.approx((-LX / 2, LX / 2), abs=1e-3)
        assert bounds["y"] == pytest.approx((-LY / 2, LY / 2), abs=1e-3)
        assert bounds["z"] == (-20.0, 20.0)

        moltemplate_dir = tmp_path / "out" / "pe_on_quartz" / "moltemplate"
        system_lt = (moltemplate_dir / "input" / "system.lt").read_text()
        assert 'import "alpha_quartz_slab.lt"' in system_lt
        assert "substrate = new QuartzSlab" in system_lt

        slab_lt = (moltemplate_dir / "input" / "alpha_quartz_slab.lt").read_text()
        assert "QuartzSlab {" in slab_lt
        assert "@atom:i15_sc4" in slab_lt

        # slab atoms present in the data file with their charges
        assert "1.1" in data  # Si charge column appears

    def test_builder_requires_lateral_dims(self, tmp_path):
        from AutoPoly.core.exceptions import WorkflowError
        system = System(out=str(tmp_path / "out2"))
        film = Polymer(chain_num=1, sequence=PE_SEQUENCE)
        spec = SubstrateSpec(builder="alpha_quartz", thickness=12.0)
        with pytest.raises(WorkflowError, match="lateral box_dims"):
            generate(
                system, "pe_on_quartz_bad", [film], force_field="gaff",
                substrate=spec,
                box_dims=(None, None, 40.0),
                geometry_config=GeometryConfig(use_mc_chain_growth=False),
                rng_seed=3,
            )
