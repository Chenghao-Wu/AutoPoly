"""Unit tests for functional-group detection in the reactor."""
import pytest
from rdkit import Chem

from AutoPoly.reactor.functional_groups import (
    detect_functional_groups,
    detect_monomer_roles,
)


def _names(fgs):
    return {fg.fg_name for fg in fgs}


def test_diol_detected():
    mol = Chem.MolFromSmiles("OCCO")
    names = _names(detect_functional_groups(mol))
    assert "diol" in names


def test_single_oh_not_diol():
    """A mono-ol (A1 chain stopper) must NOT register as a diol."""
    mol = Chem.MolFromSmiles("CCCO")
    names = _names(detect_functional_groups(mol))
    assert "diol" not in names


def test_diacid_detected():
    mol = Chem.MolFromSmiles("O=C(O)CCCCC(=O)O")  # adipic acid
    names = _names(detect_functional_groups(mol))
    assert "di_carboxylic_acid" in names


def test_amino_acid_di_different():
    mol = Chem.MolFromSmiles("NCC(=O)O")  # glycine
    names = _names(detect_functional_groups(mol))
    assert "amino_acid" in names


def test_hydroxy_acid_di_different():
    mol = Chem.MolFromSmiles("OCCC(O)=O")  # lactic-like hydroxy acid
    names = _names(detect_functional_groups(mol))
    assert "hydroxy_carboxylic_acid" in names


def test_diisocyanate_detected():
    mol = Chem.MolFromSmiles("O=C=N-C-N=C=O")
    names = _names(detect_functional_groups(mol))
    assert "di_isocyanate" in names


def test_detect_monomer_roles_skips_unreactive():
    roles = detect_monomer_roles([("OCCO", "eg"), ("c1ccccc1", "benzene")])
    names = {r.name for r in roles}
    assert "eg" in names
    assert "benzene" not in names


def test_unparseable_smiles_skipped():
    roles = detect_monomer_roles([("not_a_smiles!!!", "bad"), ("OCCO", "eg")])
    assert {r.name for r in roles} == {"eg"}
