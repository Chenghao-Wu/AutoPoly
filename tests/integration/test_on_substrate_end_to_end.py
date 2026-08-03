"""End-to-end integration test: polymer film on a physical substrate slab.

Runs the full three-stage pipeline including moltemplate for a small
PE film on an ordered water slab, and checks the generated LAMMPS data
file boundaries and instance layout.
"""

import pytest

from AutoPoly import System, Polymer, Molecule, generate
from AutoPoly.packing import SubstrateSpec, CutAbove
from AutoPoly.pipeline.geometry import GeometryConfig

PE_SEQUENCE = ["CC[*]", "[*]CC[*]", "[*]CC"]


def _boundary(content):
    """Parse (xlo, xhi, ylo, yhi, zlo, zhi) from system.data."""
    values = {}
    for line in content.splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[2:] in (["xlo", "xhi"],
                                             ["ylo", "yhi"], ["zlo", "zhi"]):
            values[parts[2][0]] = (float(parts[0]), float(parts[1]))
    return values


@pytest.mark.integration
class TestOnSubstrateEndToEnd:
    def test_pe_film_on_slab_moltemplate(self, tmp_path):
        system = System(out=str(tmp_path / "out"))
        film = Polymer(chain_num=2, sequence=PE_SEQUENCE,
                       tacticity="syndiotactic")
        spec = SubstrateSpec(
            model=Molecule(Count=12, Smiles="CCO", Name="etoh_sub"),
            thickness=8.0, packing="grid", gap=3.0,
        )
        result = generate(
            system, "pe_on_slab", [film], force_field="gaff",
            substrate=spec,
            box_dims=(40.0, 40.0, 50.0),
            geometry_config=GeometryConfig(use_mc_chain_growth=False),
            rng_seed=3,
        )

        # 12 slab + 2 film chains
        names = [r.instance_name for r in result.records]
        assert sum(n.startswith("substrate_") for n in names) == 12
        assert sum(n.startswith("polymer_") for n in names) == 2

        moltemplate_dir = tmp_path / "out" / "pe_on_slab" / "moltemplate"
        # post-processing copies the final data file to the project dir
        data = (tmp_path / "out" / "pe_on_slab" / "system.data").read_text()
        bounds = _boundary(data)
        assert bounds["x"] == (-20.0, 20.0)
        assert bounds["z"] == (-25.0, 25.0)

        # post-processing moves .lt inputs into moltemplate/input/
        system_lt = (moltemplate_dir / "input" / "system.lt").read_text()
        assert "substrate_1 = new etoh_sub" in system_lt
        assert 'import "etoh_sub.lt"' in system_lt

    def test_film_thickness_cut(self, tmp_path):
        """CutAbove at the film bottom removes every film chain."""
        system = System(out=str(tmp_path / "out"))
        film = Polymer(chain_num=2, sequence=PE_SEQUENCE)
        spec = SubstrateSpec(
            model=Molecule(Count=6, Smiles="CCO", Name="etoh_sub"),
            thickness=8.0, gap=3.0,
        )
        # film region starts at z = -25 + 8 + 3 = -14; chains sit above it
        result = generate(
            system, "pe_cut", [film], force_field="gaff",
            substrate=spec,
            box_dims=(40.0, 40.0, 50.0),
            subtract=[CutAbove(z=-14.0)],
            geometry_config=GeometryConfig(use_mc_chain_growth=False),
            rng_seed=3,
        )
        names = [r.instance_name for r in result.records]
        assert all(n.startswith("substrate_") for n in names)
        assert len(names) == 6
