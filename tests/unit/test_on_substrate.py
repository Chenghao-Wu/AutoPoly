"""Tests for the on_substrate packing strategy and workflow integration."""

import re

import pytest

from AutoPoly.core.exceptions import ValidationError, WorkflowError
from AutoPoly.models.molecule import Molecule
from AutoPoly.models.polymer import Polymer
from AutoPoly.packing import (
    BoxSpec,
    PackingContext,
    SubstrateSpec,
    get_strategy,
    registered_strategies,
)
from AutoPoly.packing.on_substrate import OnSubstrateStrategy
from AutoPoly.packing.regions import CutAbove, CutBelow, Cylinder
from AutoPoly.pipeline.units import UnitLibrary, UnitSpec
from AutoPoly.pipeline.workflow import (
    _resolve_strategy,
    _resolve_substrate_models,
)

PE_SEQUENCE = ["CC[*]", "[*]CC[*]", "[*]CC"]


def _move_z(lt_command):
    """Extract the z coordinate from a trailing .move(x,y,z)."""
    coords = lt_command.split(".move(")[1].rstrip(")").split(",")
    return float(coords[2])


def _film_library():
    """One film polymer + one film molecule species (no substrate units)."""
    return UnitLibrary(
        force_field="oplsaa",
        units=[
            UnitSpec(id="poly_1", kind="polymer", lt_file="poly_1.lt",
                     topology="linear", n_monomers=10, radius=5.0),
            UnitSpec(id="water", kind="molecule", lt_file="water.lt",
                     count=3, radius=3.0),
        ],
        monomer_files=["monomer_0_0le.lt", "water.lt"],
        build_config={"offset": 4.0},
    )


def _slab_library():
    """Film units plus a model substrate molecule unit."""
    library = _film_library()
    library.units.append(
        UnitSpec(id="sio2", kind="molecule", lt_file="sio2.lt",
                 count=8, radius=3.0, role="substrate")
    )
    library.monomer_files.append("sio2.lt")
    return library


def _spec(**kwargs):
    defaults = dict(
        model=Molecule(Count=8, Smiles="O=[Si]=O", Name="sio2"),
        thickness=10.0, packing="grid", gap=3.0,
    )
    defaults.update(kwargs)
    return SubstrateSpec(**defaults)


def _ctx(library, spec, **kwargs):
    box = kwargs.pop("box", BoxSpec(box_dims=(40.0, 40.0, 60.0)))
    return PackingContext(
        units=library, box=box, substrate=spec, rng_seed=7, **kwargs
    )


class TestSubstrateSpecValidation:
    def test_exactly_one_source(self):
        with pytest.raises(ValidationError, match="exactly one source"):
            SubstrateSpec(thickness=10.0)
        with pytest.raises(ValidationError, match="exactly one source"):
            SubstrateSpec(model=Molecule(Count=1, Smiles="O", Name="w"),
                          lt_file="slab.lt", class_name="Slab")

    def test_external_requires_class_name(self):
        with pytest.raises(ValidationError, match="class_name"):
            SubstrateSpec(lt_file="slab.lt")

    def test_numeric_validation(self):
        with pytest.raises(ValidationError, match="thickness"):
            _spec(thickness=0.0)
        with pytest.raises(ValidationError, match="gap"):
            _spec(gap=-1.0)
        with pytest.raises(ValidationError, match="packing"):
            _spec(packing="hexagonal")
        with pytest.raises(ValidationError, match="count"):
            _spec(count=-3)
        with pytest.raises(ValidationError, match="density"):
            _spec(density=0.0)


class TestModelSubstrateLayout:
    def test_registered(self):
        assert "on_substrate" in registered_strategies()

    def test_substrate_below_film_with_gap(self):
        spec = _spec()
        result = get_strategy("on_substrate").place(
            _ctx(_slab_library(), spec)
        )
        # 8 substrate + 1 polymer + 3 waters = 12 instances
        assert len(result.records) == 12

        (xlo, xhi), (ylo, yhi), (zlo, zhi) = result.box_bounds
        assert (xlo, xhi) == (-20.0, 20.0)
        assert (zlo, zhi) == (-30.0, 30.0)
        slab_top = zlo + spec.thickness          # -20
        film_zmin = slab_top + spec.gap          # -17

        for record in result.records:
            z = _move_z(record.lt_command)
            if record.instance_name.startswith("substrate_"):
                assert zlo <= z <= slab_top
            else:
                assert film_zmin < z <= zhi

    def test_instance_names(self):
        result = get_strategy("on_substrate").place(_ctx(_slab_library(), _spec()))
        names = [r.instance_name for r in result.records]
        assert names[:8] == [f"substrate_{i}" for i in range(1, 9)]
        assert "polymer_1" in names
        assert "molecule_3" in names

    def test_grid_capacity_error(self):
        library = _slab_library()
        library.units[-1].count = 10000
        with pytest.raises(ValidationError, match="grid capacity"):
            get_strategy("on_substrate").place(_ctx(library, _spec()))

    def test_mc_packed_substrate(self):
        spec = _spec(packing="mc")
        result = get_strategy("on_substrate").place(
            _ctx(_slab_library(), spec)
        )
        assert len(result.records) == 12
        (zlo, zhi) = result.box_bounds[2]
        slab_top = zlo + spec.thickness
        for record in result.records:
            if record.instance_name.startswith("substrate_"):
                assert zlo <= _move_z(record.lt_command) <= slab_top

    def test_mc_substrate_too_dense_errors(self):
        library = _slab_library()
        library.units[-1].count = 500
        library.units[-1].radius = 4.0
        spec = _spec(packing="mc")
        with pytest.raises(ValidationError, match="slab too dense"):
            get_strategy("on_substrate").place(
                _ctx(library, spec, mc_max_attempts=50)
            )

    def test_missing_substrate_units_detected(self):
        with pytest.raises(ValidationError, match="role='substrate'"):
            get_strategy("on_substrate").place(_ctx(_film_library(), _spec()))

    def test_requires_spec(self):
        ctx = PackingContext(units=_film_library())
        with pytest.raises(ValidationError, match="SubstrateSpec"):
            OnSubstrateStrategy().place(ctx)


class TestExternalSubstrate:
    def _external_spec(self, **kwargs):
        defaults = dict(lt_file="/tmp/au_slab.lt", class_name="AuSlab",
                        thickness=8.0, gap=3.0)
        defaults.update(kwargs)
        return SubstrateSpec(**defaults)

    def test_single_slab_instance(self):
        spec = self._external_spec()
        result = get_strategy("on_substrate").place(
            _ctx(_film_library(), spec)
        )
        # 1 slab + 1 polymer + 3 waters
        assert len(result.records) == 5
        slab = result.records[0]
        assert slab.instance_name == "substrate"
        assert slab.lt_command.startswith("substrate = new AuSlab.move(")
        zlo = result.box_bounds[2][0]
        assert _move_z(slab.lt_command) == pytest.approx(
            zlo + spec.thickness / 2.0
        )

    def test_film_respects_gap(self):
        spec = self._external_spec()
        result = get_strategy("on_substrate").place(
            _ctx(_film_library(), spec)
        )
        zlo, zhi = result.box_bounds[2]
        film_zmin = zlo + spec.thickness + spec.gap
        for record in result.records[1:]:
            assert film_zmin < _move_z(record.lt_command) <= zhi


class TestBoxDimResolution:
    def test_auto_lz_from_film_density(self):
        # 10 chains x 100 monomers = 1000 particles -> t_film = 1000/(0.085*40*40)
        library = UnitLibrary(
            force_field="oplsaa",
            units=[
                UnitSpec(id="poly_1", kind="polymer", lt_file="poly_1.lt",
                         topology="linear", n_monomers=100, radius=3.0,
                         count=10),
                UnitSpec(id="sio2", kind="molecule", lt_file="sio2.lt",
                         count=4, radius=3.0, role="substrate"),
            ],
        )
        spec = _spec()
        ctx = PackingContext(
            units=library,
            box=BoxSpec(box_dims=(40.0, 40.0, None)),
            substrate=spec, rng_seed=3,
        )
        result = get_strategy("on_substrate").place(ctx)
        zlo, zhi = result.box_bounds[2]
        t_film = 1000 / (0.085 * 40 * 40)
        expected_lz = spec.thickness + spec.gap + t_film + spec.vacuum
        assert (zhi - zlo) == pytest.approx(expected_lz)

    def test_lz_smaller_than_slab_errors(self):
        spec = _spec()
        ctx = PackingContext(
            units=_slab_library(),
            box=BoxSpec(box_dims=(40.0, 40.0, 5.0)),
            substrate=spec,
        )
        with pytest.raises(ValidationError, match="thickness"):
            get_strategy("on_substrate").place(ctx)

    def test_film_region_thinner_than_chain_errors(self):
        spec = _spec(thickness=25.0, gap=2.0)  # film region = 33-27 = 6 < 2*5
        ctx = PackingContext(
            units=_slab_library(),
            box=BoxSpec(box_dims=(40.0, 40.0, 33.0)),
            substrate=spec,
        )
        with pytest.raises(ValidationError, match="Film region"):
            get_strategy("on_substrate").place(ctx)


class TestSubtract:
    def test_cutabove_removes_all_film_only(self):
        spec = _spec()
        library = _slab_library()
        zlo = -30.0
        film_zmin = zlo + spec.thickness + spec.gap
        ctx = _ctx(library, spec,
                   subtract=[CutAbove(z=film_zmin)])  # film centers all above
        result = get_strategy("on_substrate").place(ctx)
        names = [r.instance_name for r in result.records]
        assert all(n.startswith("substrate_") for n in names)
        assert len(names) == 8

    def test_cutbelow_nothing_in_film_range(self):
        spec = _spec()
        zlo = -30.0
        film_zmin = zlo + spec.thickness + spec.gap
        ctx = _ctx(_slab_library(), spec,
                   subtract=[CutBelow(z=film_zmin)])  # film centers all above
        result = get_strategy("on_substrate").place(ctx)
        assert len(result.records) == 12

    def test_carve_substrate_role(self):
        spec = _spec()
        zlo = -30.0
        ctx = _ctx(_slab_library(), spec,
                   subtract=[CutBelow(z=zlo + spec.thickness - 1.0,
                                      apply_to="substrate")])
        result = get_strategy("on_substrate").place(ctx)
        # grid slab sits at z=-25 (below the cut plane): all 8 removed
        names = [r.instance_name for r in result.records]
        assert not any(n.startswith("substrate_") for n in names)
        assert len(names) == 4

    def test_cylinder_removes_center_instances(self):
        spec = _spec()
        ctx = _ctx(_slab_library(), spec,
                   subtract=[Cylinder(axis="z", center=(0, 0), radius=100.0)])
        result = get_strategy("on_substrate").place(ctx)
        # whole film inside a huge cylinder, substrate untouched
        names = [r.instance_name for r in result.records]
        assert len(names) == 8
        assert all(n.startswith("substrate_") for n in names)

    def test_mc_random_supports_subtract(self):
        ctx = PackingContext(
            units=_film_library(),
            box=BoxSpec(requested_box_size=80.0),
            rng_seed=11,
            subtract=[CutAbove(z=1000.0)],  # removes nothing
        )
        result = get_strategy("mc_random").place(ctx)
        assert len(result.records) == 4  # 1 polymer + 3 waters

        ctx.subtract = [CutBelow(z=1000.0)]  # removes everything
        result = get_strategy("mc_random").place(ctx)
        assert result.records == []

    def test_grid_rejects_subtract_and_substrate(self):
        with pytest.raises(ValidationError, match="subtract"):
            get_strategy("grid").place(
                PackingContext(units=_film_library(),
                               subtract=[CutAbove(z=0.0)])
            )
        with pytest.raises(ValidationError, match="substrate"):
            get_strategy("grid").place(
                PackingContext(units=_film_library(), substrate=_spec())
            )


class TestWorkflowHelpers:
    def test_substrate_auto_selects_strategy(self):
        assert _resolve_strategy("mc_random", _spec(), None) == "on_substrate"

    def test_substrate_conflicts_with_explicit_strategy(self):
        with pytest.raises(ValidationError, match="conflicts"):
            _resolve_strategy("grid", _spec(), None)

    def test_explicit_on_substrate_kept(self):
        assert _resolve_strategy("on_substrate", _spec(), None) == "on_substrate"

    def test_subtract_rejected_for_grid(self):
        with pytest.raises(ValidationError, match="subtract"):
            _resolve_strategy("grid", None, [CutAbove(z=0.0)])

    def test_no_substrate_no_models(self):
        assert _resolve_substrate_models(None, None) == []

    def test_external_substrate_no_models(self):
        spec = SubstrateSpec(lt_file="s.lt", class_name="S", thickness=5.0)
        assert _resolve_substrate_models(spec, None) == []

    def test_auto_count_molecule(self):
        model = Molecule(Count=1, Smiles="O", Name="water")
        spec = _spec(model=model, count="auto", density=0.1)
        models = _resolve_substrate_models(spec, (10.0, 20.0, None))
        assert len(models) == 1
        # 0.1 * 10 * 20 * 10 = 200
        assert models[0].Count == 200
        assert models[0].Smiles == "O"

    def test_auto_count_needs_lateral_dims(self):
        spec = _spec(count="auto")
        with pytest.raises(ValidationError, match="box_dims"):
            _resolve_substrate_models(spec, (None, 20.0, None))

    def test_auto_count_rejects_polymer(self):
        poly = Polymer(chain_num=2, sequence=PE_SEQUENCE)
        spec = _spec(model=poly, count="auto")
        with pytest.raises(ValidationError, match="Molecule"):
            _resolve_substrate_models(spec, (10.0, 10.0, None))

    def test_int_count_override_polymer(self):
        poly = Polymer(chain_num=2, sequence=PE_SEQUENCE)
        spec = _spec(model=poly, count=5)
        models = _resolve_substrate_models(spec, None)
        assert models[0].chain_num == 5
        assert len(models[0].sequenceSet) == 5
        # original untouched
        assert poly.chain_num == 2


class TestGenerateEndToEnd:
    def test_generate_with_model_substrate(self, tmp_path):
        from AutoPoly import System, generate
        from AutoPoly.pipeline.geometry import GeometryConfig

        system = System(out=str(tmp_path / "out"))
        film = Polymer(chain_num=2, sequence=PE_SEQUENCE,
                       tacticity="syndiotactic")
        spec = SubstrateSpec(
            model=Molecule(Count=4, Smiles="O", Name="water_sub"),
            thickness=8.0, packing="grid", gap=3.0,
        )
        result = generate(
            system, "film_on_water", [film],
            substrate=spec,
            box_dims=(40.0, 40.0, 60.0),
            geometry_config=GeometryConfig(use_mc_chain_growth=False),
            rng_seed=3,
            run_moltemplate=False,
        )
        names = [r.instance_name for r in result.records]
        assert names.count("polymer_1") == 1
        assert names.count("polymer_2") == 1
        assert sum(1 for n in names if n.startswith("substrate_")) == 4

        system_lt = (
            tmp_path / "out" / "film_on_water" / "moltemplate" / "system.lt"
        )
        content = system_lt.read_text()
        assert 'import "water_sub.lt"' in content
        assert "substrate_1 = new water_sub" in content

    def test_generate_with_external_substrate(self, tmp_path):
        from AutoPoly import System, generate
        from AutoPoly.pipeline.geometry import GeometryConfig

        # Minimal external slab class
        slab = tmp_path / "slab.lt"
        slab.write_text(
            'Slab inherits OPLSAA {\n'
            '  write("Data Atoms") {\n'
            '    $atom:SI $mol:... @atom:si 0.0  0.0 0.0 0.0\n'
            "  }\n"
            "}\n"
        )

        system = System(out=str(tmp_path / "out"))
        film = Polymer(chain_num=1, sequence=PE_SEQUENCE)
        spec = SubstrateSpec(lt_file=str(slab), class_name="Slab",
                             thickness=8.0, gap=3.0)
        result = generate(
            system, "film_on_slab", [film],
            substrate=spec,
            box_dims=(40.0, 40.0, 60.0),
            geometry_config=GeometryConfig(use_mc_chain_growth=False),
            rng_seed=3,
            run_moltemplate=False,
        )
        names = [r.instance_name for r in result.records]
        assert "substrate" in names
        assert "polymer_1" in names

        moltemplate_dir = tmp_path / "out" / "film_on_slab" / "moltemplate"
        assert (moltemplate_dir / "slab.lt").is_file()
        content = (moltemplate_dir / "system.lt").read_text()
        assert 'import "slab.lt"' in content
        assert "substrate = new Slab.move(" in content

    def test_generate_missing_external_slab_errors(self, tmp_path):
        from AutoPoly import System, generate
        from AutoPoly.pipeline.geometry import GeometryConfig

        system = System(out=str(tmp_path / "out"))
        film = Polymer(chain_num=1, sequence=PE_SEQUENCE)
        spec = SubstrateSpec(lt_file=str(tmp_path / "nope.lt"),
                             class_name="Slab", thickness=8.0)
        with pytest.raises(WorkflowError, match="not found"):
            generate(
                system, "bad_slab", [film],
                substrate=spec,
                box_dims=(40.0, 40.0, 60.0),
                geometry_config=GeometryConfig(use_mc_chain_growth=False),
                run_moltemplate=False,
            )
