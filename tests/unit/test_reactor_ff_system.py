"""Unit tests for force-field tables and system.data parsing."""
import pytest
from pathlib import Path

from AutoPoly.reactor.fftables import ForceFieldTables, load_force_field_tables
from AutoPoly.reactor.system_data import SystemData, parse_in_init

FF_DIR = (Path(__file__).parent.parent.parent /
          "AutoPoly" / "extern" / "moltemplate" / "force_fields")


@pytest.fixture(scope="module")
def gaff2():
    return ForceFieldTables(FF_DIR / "gaff2.lt")


def test_gaff2_bond_lookup(gaff2):
    name, coeff = gaff2.bond_for("c3", "oh")
    assert name == "c3-oh"
    assert "harmonic" in coeff


def test_gaff2_bond_order_insensitive(gaff2):
    fwd = gaff2.bond_for("c3", "oh")
    rev = gaff2.bond_for("oh", "c3")
    assert fwd is not None and fwd == rev


def test_gaff2_angle_lookup(gaff2):
    name, coeff = gaff2.angle_for("c3", "c3", "os")
    assert name == "c3-c3-os"
    assert "harmonic" in coeff


def test_gaff2_dihedral_prefers_exact_over_wildcard(gaff2):
    """Exact definitions must beat generic X-c3-c3-X wildcards."""
    name, _ = gaff2.dihedral_for("hc", "c3", "c3", "hc")
    assert name == "hc-c3-c3-hc"


def test_gaff2_dihedral_wildcard_fallback(gaff2):
    """A combination with no exact def falls back to the wildcard def."""
    name, _ = gaff2.dihedral_for("h1", "c3", "c3", "h1")
    assert name.startswith("X-")


def test_gaff2_pair_coeff_clean(gaff2):
    """Pair coeff text must not retain the @atom tokens."""
    coeff = gaff2.pair_for("o", "o")
    assert coeff is not None
    assert "@atom" not in coeff
    assert "lj/charmm/coul/long" in coeff


def test_gaff2_mass(gaff2):
    assert gaff2.mass_for("c3") == pytest.approx(12.01)


def test_load_force_field_tables_registry():
    tables = load_force_field_tables("gaff2")
    assert isinstance(tables, ForceFieldTables)


# ---------------------------------------------------------------------------
# system.data parsing (uses the checked-in example output)
# ---------------------------------------------------------------------------
PEO_DATA = Path("/home/zhenghaowu/autopoly_dev/peo_staged/peo/system.data")
PEO_INIT = Path("/home/zhenghaowu/autopoly_dev/peo_staged/peo/system.in.init")


@pytest.mark.skipif(not PEO_DATA.is_file(), reason="example system.data not present")
def test_system_data_atom_types():
    sd = SystemData(PEO_DATA)
    assert sd.atom_type_names[1] == "c3"
    assert sd.atom_type_id("os") == 6
    assert sd.max_atom_type() == 6


@pytest.mark.skipif(not PEO_DATA.is_file(), reason="example system.data not present")
def test_system_data_bond_by_example():
    sd = SystemData(PEO_DATA)
    # c3-c3 bond exists in the PEO data with a numeric type id
    assert sd.lookup("bond", ("c3", "c3")) == 1
    # order-insensitive
    assert sd.lookup("bond", ("os", "c3")) == sd.lookup("bond", ("c3", "os"))
    # a bond that does not exist returns None
    assert sd.lookup("bond", ("c", "c")) is None


@pytest.mark.skipif(not PEO_INIT.is_file(), reason="example system.in.init not present")
def test_parse_in_init():
    styles = parse_in_init(PEO_INIT)
    assert styles["atom_style"] == "full"
    assert styles["units"] == "real"
    assert "bond_style" in styles
