"""Tests for carve regions (whole-instance subtract) and the filter pass."""

import pytest

from AutoPoly.core.exceptions import ValidationError
from AutoPoly.packing.regions import (
    BoxRegion,
    CutAbove,
    CutBelow,
    Cylinder,
    PlacedItem,
    apply_carve_regions,
)


def _item(name, z=0.0, x=0.0, y=0.0, role="film", radius=1.0):
    return PlacedItem(
        unit_id="u1", instance_name=name, role=role,
        position=(x, y, z), radius=radius,
        payload=("command", f"{name} = new u1.move({x},{y},{z})"),
    )


class TestCutAbove:
    def test_center_semantics(self):
        cut = CutAbove(z=10.0)
        assert cut.contains((0, 0, 11.0))
        assert not cut.contains((0, 0, 9.0))

    def test_conservative_inflates_by_radius(self):
        cut = CutAbove(z=10.0, conservative=True)
        # center below plane, but sphere touches it
        assert cut.contains((0, 0, 9.0), radius=2.0)
        assert not cut.contains((0, 0, 7.0), radius=2.0)

    def test_nonconservative_ignores_radius(self):
        cut = CutAbove(z=10.0)
        assert not cut.contains((0, 0, 9.0), radius=5.0)


class TestCutBelow:
    def test_center_semantics(self):
        cut = CutBelow(z=-5.0)
        assert cut.contains((0, 0, -6.0))
        assert not cut.contains((0, 0, 0.0))

    def test_conservative(self):
        cut = CutBelow(z=-5.0, conservative=True)
        assert cut.contains((0, 0, -4.0), radius=2.0)
        assert not cut.contains((0, 0, -2.0), radius=2.0)


class TestCylinder:
    def test_z_axis(self):
        cyl = Cylinder(axis="z", center=(0, 0), radius=5.0)
        assert cyl.contains((3, 3.9, 100.0))    # radial distance < 5
        assert not cyl.contains((3, 4.1, 0.0))  # just outside
        assert cyl.contains((0, 0, -999.0))     # infinite along axis

    def test_x_axis_uses_yz_plane(self):
        cyl = Cylinder(axis="x", center=(0, 0), radius=2.0)
        assert cyl.contains((50.0, 1.0, 1.0))
        assert not cyl.contains((0.0, 3.0, 0.0))

    def test_conservative(self):
        cyl = Cylinder(axis="z", center=(0, 0), radius=5.0, conservative=True)
        assert cyl.contains((6.5, 0, 0), radius=2.0)
        assert not cyl.contains((8.0, 0, 0), radius=2.0)

    def test_validation(self):
        with pytest.raises(ValidationError, match="axis"):
            Cylinder(axis="w")
        with pytest.raises(ValidationError, match="radius"):
            Cylinder(radius=0.0)


class TestBoxRegion:
    def test_bounded_axes(self):
        box = BoxRegion(x=(-10, 0), z=(5, 15))
        assert box.contains((-5, 999.0, 10.0))   # y unbounded
        assert not box.contains((5, 0.0, 10.0))
        assert not box.contains((-5, 0.0, 20.0))

    def test_conservative_expands_ranges(self):
        box = BoxRegion(x=(0, 10), conservative=True)
        assert box.contains((-1.5, 0, 0), radius=2.0)
        assert not box.contains((-3.0, 0, 0), radius=2.0)

    def test_validation(self):
        with pytest.raises(ValidationError, match="lo < hi"):
            BoxRegion(x=(5, -5))
        with pytest.raises(ValidationError, match="at least one bounded"):
            BoxRegion()


class TestApplyCarveRegions:
    def test_no_regions_is_noop(self):
        items = [_item("a"), _item("b")]
        assert apply_carve_regions(items, None) == items
        assert apply_carve_regions(items, []) == items

    def test_role_filtering(self):
        items = [
            _item("film_1", z=20.0, role="film"),
            _item("sub_1", z=20.0, role="substrate"),
        ]
        kept = apply_carve_regions(items, [CutAbove(z=10.0)])  # film default
        assert [i.instance_name for i in kept] == ["sub_1"]

        kept = apply_carve_regions(
            items, [CutAbove(z=10.0, apply_to="substrate")]
        )
        assert [i.instance_name for i in kept] == ["film_1"]

        kept = apply_carve_regions(items, [CutAbove(z=10.0, apply_to="all")])
        assert kept == []

    def test_regions_apply_in_order_first_hit_wins(self):
        items = [_item("a", z=20.0), _item("b", x=3.0)]
        regions = [CutAbove(z=10.0), Cylinder(axis="z", center=(0, 0), radius=5.0)]
        kept = apply_carve_regions(items, regions)
        assert kept == []

    def test_apply_to_validation(self):
        with pytest.raises(ValidationError, match="apply_to"):
            CutAbove(z=0.0, apply_to="everything")
