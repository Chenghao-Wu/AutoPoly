"""Unit tests for reaction detection and reaction mapping in the reactor."""
import pytest
from rdkit import Chem

from AutoPoly.reactor.detector import detect_reactions
from AutoPoly.reactor.functional_groups import detect_monomer_roles
from AutoPoly.reactor.mapper import prepare_reactions


EG = ("OCCO", "eg")
ADIPIC = ("O=C(O)CCCCC(=O)O", "adipic")
GLYCINE = ("NCC(=O)O", "glycine")
HMDA = ("NCCCCCCN", "hmda")


def _instances(monomers):
    return detect_reactions(detect_monomer_roles(monomers))


def test_diol_diacid_single_reaction():
    instances = _instances([EG, ADIPIC])
    assert len(instances) == 1
    assert "Polyesterification" in instances[0].reaction_name
    assert instances[0].delete_atom is True


def test_no_reaction_for_unreactive_pair():
    instances = _instances([("c1ccccc1", "benzene"), ("CCO", "ethanol")])
    assert instances == []


def test_amino_acid_diamine_diacid_reactions():
    instances = _instances([HMDA, ADIPIC])
    names = {i.reaction_name for i in instances}
    assert any("Polyamidation" in n for n in names)


def test_self_condensation_ab_monomer():
    """An AB monomer (amino acid) self-condenses (A + A via two groups)."""
    instances = _instances([GLYCINE])
    names = {i.reaction_name for i in instances}
    assert any("Polyamidation" in n for n in names)


def test_mapper_eg_adipic_water_byproduct():
    instances = _instances([EG, ADIPIC])
    metas = prepare_reactions(instances)
    assert len(metas) == 1
    m = metas[0]

    # Two initiators on opposite molecules (carbonyl C and hydroxyl O)
    assert len(m.initiators) == 2
    elements = sorted(
        m.reactant_mol.GetAtomWithIdx(i).GetSymbol() for i in m.initiators
    )
    assert elements == ["C", "O"]

    # Byproduct = water (O + 2H)
    assert len(m.byproduct_indices) == 3
    bp_elements = sorted(
        m.reactant_mol.GetAtomWithIdx(i).GetSymbol() for i in m.byproduct_indices
    )
    assert bp_elements == ["H", "H", "O"]

    # Product contains the ester + a free water fragment
    frags = Chem.GetMolFrags(m.product_mol)
    assert len(frags) == 2

    # Full 1-to-1 mapping covering every atom
    assert len(m.reactant_to_product) == m.reactant_mol.GetNumAtoms()
    assert len(m.product_to_reactant) == m.product_mol.GetNumAtoms()

    # Template subset includes initiators and edge atoms are outside first shell
    assert set(m.initiators) <= set(m.template_atoms)
    assert set(m.edge_atoms) <= set(m.template_atoms)


def test_mapper_mapping_is_bidirectional():
    instances = _instances([EG, ADIPIC])
    m = prepare_reactions(instances)[0]
    for r, p in m.reactant_to_product.items():
        assert m.product_to_reactant[p] == r
