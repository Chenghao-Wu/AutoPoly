"""End-to-end reactor tests: template construction and script writing."""
import pytest
from pathlib import Path

from AutoPoly.reactor import Reactor
from AutoPoly.reactor.templates import TemplateBuilder

PEO_DIR = Path("/home/zhenghaowu/autopoly_dev/peo_staged/peo")
EG = ("OCCO", "eg")
ADIPIC = ("O=C(O)CCCCC(=O)O", "adipic")

pytestmark = pytest.mark.skipif(
    not (PEO_DIR / "system.data").is_file(),
    reason="example system.data not present",
)


def _build(tmp_path):
    reactor = Reactor(PEO_DIR, monomers=[EG, ADIPIC], force_field="gaff2")
    result = reactor.build(output_dir=tmp_path / "reactor")
    return reactor, result


def test_reactor_builds_template_triplet(tmp_path):
    _, result = _build(tmp_path)
    assert result.reactions
    ts = result.template_sets[0]
    assert ts.pre_file.is_file()
    assert ts.post_file.is_file()
    assert ts.map_file.is_file()
    assert ts.map_file_delete_ids is not None and ts.map_file_delete_ids.is_file()


def test_template_molecule_file_sections(tmp_path):
    _, result = _build(tmp_path)
    content = result.template_sets[0].pre_file.read_text()
    for section in ("atoms", "bonds", "angles", "dihedrals", "Types",
                    "Charges", "Coords", "Bonds"):
        assert section in content


def test_map_file_structure(tmp_path):
    _, result = _build(tmp_path)
    content = result.template_sets[0].map_file_delete_ids.read_text()
    assert "edgeIDs" in content
    assert "equivalences" in content
    assert "deleteIDs" in content
    assert "InitiatorIDs" in content
    assert "EdgeIDs" in content
    assert "Equivalences" in content
    assert "DeleteIDs" in content


def test_equivalences_cover_template_atoms(tmp_path):
    reactor, result = _build(tmp_path)
    meta = reactor.reaction_metadata[0]
    n_template = len(meta.template_atoms)
    content = result.template_sets[0].map_file.read_text()
    assert f"{n_template} equivalences" in content


def test_bond_react_script_contents(tmp_path):
    _, result = _build(tmp_path)
    script = result.script.read_text()
    assert "fix rxns all bond/react" in script
    assert "molecule mol_pre_1" in script
    assert "molecule mol_post_1" in script
    assert "extra/bond/per/atom" in script
    assert "statted_grp_REACT nvt" in script
    assert "rescale_charges yes" in script


def test_new_types_introduced_for_ester(tmp_path):
    """The ester product needs types (o, c, ow, hw) absent from PEO data."""
    _, result = _build(tmp_path)
    new_types = set()
    for a in result.assignments:
        new_types.update(a.new_atom_types)
    assert {"o", "c"} <= new_types


def test_no_reactions_raises(tmp_path):
    from AutoPoly.core.exceptions import GenerationError
    reactor = Reactor(
        PEO_DIR,
        monomers=[("c1ccccc1", "benzene"), ("CCO", "ethanol")],
        force_field="gaff2",
    )
    with pytest.raises(GenerationError):
        reactor.build(output_dir=tmp_path / "reactor")
