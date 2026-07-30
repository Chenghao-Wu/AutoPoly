"""Tests for packing strategies, the registry, and BoxPacker (stage 3)."""

import numpy as np
import pytest

from AutoPoly.core.exceptions import ValidationError
from AutoPoly.pipeline.geometry import GeometryBuilder, GeometryConfig
from AutoPoly.pipeline.packer import BoxPacker
from AutoPoly.packing import (
    PlacementStrategy,
    PackingContext,
    get_strategy,
    register_strategy,
    registered_strategies,
)
from AutoPoly.packing.base import BoxSpec
from AutoPoly.packing.random_mc import compute_auto_box_size
from AutoPoly.models.polymer import Polymer
from AutoPoly.core.system import System
from AutoPoly.pipeline.typing import UnitTyper
from AutoPoly.pipeline.units import UnitLibrary, UnitSpec

PE_SEQUENCE = ["CC[*]", "[*]CC[*]", "[*]CC"]


def _library(**kwargs):
    defaults = dict(
        force_field="oplsaa",
        units=[
            UnitSpec(id="poly_1", kind="polymer", lt_file="poly_1.lt",
                     topology="linear", n_monomers=10, radius=12.0),
            UnitSpec(id="poly_2", kind="polymer", lt_file="poly_2.lt",
                     topology="linear", n_monomers=10, radius=12.0),
            UnitSpec(id="water", kind="molecule", lt_file="water.lt", count=5,
                     radius=3.0),
        ],
        monomer_files=["monomer_0_0le.lt", "water.lt"],
        build_config={"offset": 4.0},
    )
    defaults.update(kwargs)
    return UnitLibrary(**defaults)


def _ctx(library, **kwargs):
    return PackingContext(units=library, box=BoxSpec(), **kwargs)


class TestRegistry:
    def test_builtins_registered(self):
        assert "grid" in registered_strategies()
        assert "mc_random" in registered_strategies()

    def test_get_strategy_returns_instance(self):
        strategy = get_strategy("grid")
        assert isinstance(strategy, PlacementStrategy)
        assert strategy.name == "grid"

    def test_unknown_strategy_raises_with_names(self):
        with pytest.raises(ValidationError, match="grid"):
            get_strategy("does_not_exist")

    def test_register_custom_strategy(self):
        class NullStrategy(PlacementStrategy):
            name = "null_test"

            def place(self, ctx):
                from AutoPoly.packing.base import PlacementResult, symmetric_bounds
                return PlacementResult(records=[], box_bounds=symmetric_bounds(10.0))

        register_strategy("null_test", NullStrategy)
        assert "null_test" in registered_strategies()
        result = get_strategy("null_test").place(_ctx(_library()))
        assert result.records == []

    def test_reserved_names_rejected(self):
        class Bad(PlacementStrategy):
            name = "grafted_surface"

            def place(self, ctx):
                raise NotImplementedError

        with pytest.raises(ValidationError, match="reserved"):
            register_strategy("grafted_surface", Bad)

    def test_non_strategy_class_rejected(self):
        with pytest.raises(ValidationError, match="PlacementStrategy"):
            register_strategy("not_a_strategy", dict)


class TestGridStrategy:
    def test_grid_placements(self):
        library = _library()
        result = get_strategy("grid").place(_ctx(library))

        # 2 polymers + 5 waters = 7 instances
        assert len(result.records) == 7
        names = [r.instance_name for r in result.records]
        assert names == [
            "polymer_1", "polymer_2",
            "molecule_1", "molecule_2", "molecule_3", "molecule_4", "molecule_5",
        ]
        # Polymers instantiate their poly class
        assert result.records[0].lt_command.startswith("polymer_1 = new poly_1")
        # Grid layout stays in the z=0 plane
        for record in result.records:
            assert record.lt_command.rstrip().endswith(")")
            coords = record.lt_command.split(".move(")[1].rstrip(")").split(",")
            assert float(coords[2]) == 0.0

    def test_grid_box_bounds_symmetric(self):
        result = get_strategy("grid").place(_ctx(_library()))
        for lo, hi in result.box_bounds:
            assert lo == -hi
            assert hi > 0

    def test_grid_requested_box_size(self):
        library = _library()
        ctx = PackingContext(
            units=library, box=BoxSpec(requested_box_size=100.0)
        )
        result = get_strategy("grid").place(ctx)
        assert result.box_bounds[0] == (-50.0, 50.0)


class TestRandomMCStrategy:
    def test_mc_placements_and_box(self):
        library = _library()
        result = get_strategy("mc_random").place(_ctx(library, rng_seed=42))

        assert len(result.records) == 7
        for lo, hi in result.box_bounds:
            assert lo == -hi
        assert result.records[0].lt_command.startswith("polymer_1 = new poly_1")

    def test_mc_seed_reproducibility(self):
        library = _library()
        r1 = get_strategy("mc_random").place(_ctx(library, rng_seed=123))
        r2 = get_strategy("mc_random").place(_ctx(library, rng_seed=123))
        assert [r.lt_command for r in r1.records] == [r.lt_command for r in r2.records]

    def test_box_size_formula(self):
        """Box = max(density, chain-length, packing-volume) estimates."""
        library = _library()
        ctx = _ctx(library)
        box = compute_auto_box_size(ctx)

        # density: 2*10 monomers + 5 waters = 25 particles at 0.085/A^3
        from AutoPoly.mc import calculate_box_size
        density_box = calculate_box_size(25, 0.085)
        # chain length: 10^0.6 * 4.0 * 3.0
        chain_box = 10**0.6 * 4.0 * 3.0
        # packing: (2*(4/3)pi*12^3 + 5*(4/3)pi*27)/0.3)^(1/3)
        volume = 2 * (4/3) * np.pi * 12.0**3 + 5 * (4/3) * np.pi * 3.0**3
        packing_box = (volume / 0.3) ** (1/3)

        assert box == pytest.approx(max(density_box, chain_box, packing_box))

    def test_mc_requested_box_size(self):
        library = _library()
        ctx = PackingContext(
            units=library, box=BoxSpec(requested_box_size=200.0), rng_seed=1
        )
        result = get_strategy("mc_random").place(ctx)
        assert result.box_bounds[0] == (-100.0, 100.0)


class TestBoxPacker:
    def _build_units(self, tmp_path, ff="oplsaa"):
        system = System(out=str(tmp_path / "out"))
        poly = Polymer(chain_num=2, sequence=PE_SEQUENCE, tacticity="syndiotactic")
        config = GeometryConfig(use_mc_chain_growth=False)
        geom = GeometryBuilder(system, "proj", config).build([poly])
        library = UnitTyper(geom.dir, ff).type()
        return system, library

    def test_system_lt_layout(self, tmp_path):
        system, library = self._build_units(tmp_path)
        result = BoxPacker(
            system, "proj", strategy="mc_random",
            rng_seed=5, run_moltemplate=False,
        ).pack(library)

        moltemplate_dir = tmp_path / "out" / "proj" / "moltemplate"
        system_lt = moltemplate_dir / "system.lt"
        assert system_lt.is_file()
        content = system_lt.read_text()
        assert 'import "oplsaa.lt"' in content
        assert 'import "poly_1.lt"' in content
        assert 'import "poly_2.lt"' in content
        assert 'write_once("Data Boundary")' in content
        # FF subset + copied monomer files present
        assert (moltemplate_dir / "oplsaa.lt").is_file()
        for f in library.monomer_files:
            assert (moltemplate_dir / f).is_file()
        assert len(result.records) == 2

    def test_pack_from_build_dir_path(self, tmp_path):
        system, library = self._build_units(tmp_path)
        BoxPacker(
            system, "proj", strategy="grid", run_moltemplate=False
        ).pack(library.source_dir)

        moltemplate_dir = tmp_path / "out" / "proj" / "moltemplate"
        assert (moltemplate_dir / "system.lt").is_file()

    def test_unknown_strategy_raises(self, tmp_path):
        system = System(out=str(tmp_path / "out"))
        with pytest.raises(ValidationError, match="Unknown packing strategy"):
            BoxPacker(system, "proj", strategy="nope")

    def test_grid_strategy_pack(self, tmp_path):
        system, library = self._build_units(tmp_path)
        result = BoxPacker(
            system, "proj", strategy="grid", run_moltemplate=False
        ).pack(library)
        names = [r.instance_name for r in result.records]
        assert names == ["polymer_1", "polymer_2"]

    def test_missing_monomer_file_detected(self, tmp_path):
        system, library = self._build_units(tmp_path)
        # Remove one referenced monomer file from the build dir
        import os
        victim = os.path.join(library.source_dir, library.monomer_files[0])
        os.remove(victim)
        with pytest.raises(ValidationError, match="missing"):
            BoxPacker(system, "proj", run_moltemplate=False).pack(library)
