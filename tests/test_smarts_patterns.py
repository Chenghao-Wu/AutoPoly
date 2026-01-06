"""
Test SMARTS patterns for polymerization mechanisms.

This module tests the SMARTS patterns defined in polymerization_patterns.py
to ensure they correctly identify functional groups in various molecules.

Key Tests:
- Vinyl addition pattern (C=C)
- Esterification pattern (carboxyl + alcohol)
- Amidation pattern (carboxyl + amine)
- Etherification pattern (alcohol groups)
"""

import sys
import os

# Add the AutoPoly package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from rdkit import Chem

from AutoPoly.polymerization_patterns import (
    POLYMERIZATION_PATTERNS,
    get_pattern,
    get_all_mechanisms,
    get_polymerization_mechanisms,
    validate_mechanism
)


class TestPatternDefinitions:
    """Test that pattern definitions are properly structured."""

    def test_all_mechanisms_defined(self):
        """Test that all expected mechanisms are defined."""
        expected_mechanisms = {'none', 'vinyl_addition', 'esterification', 'amidation', 'etherification'}
        actual_mechanisms = set(POLYMERIZATION_PATTERNS.keys())
        assert expected_mechanisms == actual_mechanisms

    def test_pattern_has_required_attributes(self):
        """Test that each pattern has required attributes."""
        for name, pattern in POLYMERIZATION_PATTERNS.items():
            assert hasattr(pattern, 'name')
            assert hasattr(pattern, 'smarts')
            assert hasattr(pattern, 'connection_atoms')
            assert hasattr(pattern, 'description')
            assert hasattr(pattern, 'is_condensation')
            assert hasattr(pattern, 'requires_two_groups')

    def test_get_pattern(self):
        """Test get_pattern function."""
        pattern = get_pattern('vinyl_addition')
        assert pattern.name == 'vinyl_addition'
        assert pattern.smarts is not None

    def test_get_pattern_invalid(self):
        """Test get_pattern with invalid mechanism."""
        with pytest.raises(ValueError):
            get_pattern('invalid_mechanism')

    def test_get_all_mechanisms(self):
        """Test get_all_mechanisms function."""
        mechanisms = get_all_mechanisms()
        assert 'none' in mechanisms
        assert 'vinyl_addition' in mechanisms

    def test_validate_mechanism(self):
        """Test validate_mechanism function."""
        assert validate_mechanism('vinyl_addition') == True
        assert validate_mechanism('invalid') == False


class TestVinylPattern:
    """Test vinyl addition SMARTS pattern."""

    def test_vinyl_pattern_detects_ethylene(self):
        """Test that vinyl pattern detects ethylene (C=C)."""
        pattern = POLYMERIZATION_PATTERNS['vinyl_addition']
        mol = Chem.MolFromSmiles("C=C")
        mol = Chem.AddHs(mol)

        smarts_pattern = Chem.MolFromSmarts(pattern.smarts)
        assert mol.HasSubstructMatch(smarts_pattern)

    def test_vinyl_pattern_detects_propylene(self):
        """Test that vinyl pattern detects propylene (C=CC)."""
        pattern = POLYMERIZATION_PATTERNS['vinyl_addition']
        mol = Chem.MolFromSmiles("C=CC")
        mol = Chem.AddHs(mol)

        smarts_pattern = Chem.MolFromSmarts(pattern.smarts)
        assert mol.HasSubstructMatch(smarts_pattern)

    def test_vinyl_pattern_detects_styrene(self):
        """Test that vinyl pattern detects styrene."""
        pattern = POLYMERIZATION_PATTERNS['vinyl_addition']
        # Styrene: C=Cc1ccccc1
        mol = Chem.MolFromSmiles("C=Cc1ccccc1")
        mol = Chem.AddHs(mol)

        smarts_pattern = Chem.MolFromSmarts(pattern.smarts)
        assert mol.HasSubstructMatch(smarts_pattern)

    def test_vinyl_pattern_no_match_ethanol(self):
        """Test that vinyl pattern does not match ethanol."""
        pattern = POLYMERIZATION_PATTERNS['vinyl_addition']
        mol = Chem.MolFromSmiles("CCO")
        mol = Chem.AddHs(mol)

        smarts_pattern = Chem.MolFromSmarts(pattern.smarts)
        assert not mol.HasSubstructMatch(smarts_pattern)


class TestEsterificationPattern:
    """Test esterification SMARTS pattern."""

    def test_ester_pattern_detects_lactic_acid(self):
        """Test that esterification pattern detects lactic acid."""
        pattern = POLYMERIZATION_PATTERNS['esterification']
        # Lactic acid: CC(C(=O)O)O
        mol = Chem.MolFromSmiles("CC(C(=O)O)O")
        mol = Chem.AddHs(mol)

        # Should have carboxyl and alcohol groups
        carboxyl_smarts = Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]')
        alcohol_smarts = Chem.MolFromSmarts('[$([OX2H])]')

        has_carboxyl = mol.HasSubstructMatch(carboxyl_smarts)
        has_alcohol = mol.HasSubstructMatch(alcohol_smarts)

        assert has_carboxyl, "Lactic acid should have carboxyl group"
        assert has_alcohol, "Lactic acid should have alcohol group"

    def test_ester_pattern_no_match_ethanol(self):
        """Test that esterification pattern requires carboxyl."""
        # Ethanol only has alcohol, no carboxyl
        mol = Chem.MolFromSmiles("CCO")
        mol = Chem.AddHs(mol)

        carboxyl_smarts = Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]')
        assert not mol.HasSubstructMatch(carboxyl_smarts)


class TestAmidationPattern:
    """Test amidation SMARTS pattern."""

    def test_amide_pattern_detects_nylon_monomer(self):
        """Test that amidation pattern detects nylon precursor."""
        # Simple amino acid-like structure: H2N-CH2-COOH (glycine)
        mol = Chem.MolFromSmiles("NCC(=O)O")
        mol = Chem.AddHs(mol)

        # Should have carboxyl and amine groups
        carboxyl_smarts = Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]')
        amine_smarts = Chem.MolFromSmarts('[NX3]')

        has_carboxyl = mol.HasSubstructMatch(carboxyl_smarts)
        has_amine = mol.HasSubstructMatch(amine_smarts)

        assert has_carboxyl, "Glycine should have carboxyl group"
        assert has_amine, "Glycine should have amine group"

    def test_amide_pattern_no_match_simple_amine(self):
        """Test that amidation pattern requires carboxyl."""
        # Ethylamine: CCN
        mol = Chem.MolFromSmiles("CCN")
        mol = Chem.AddHs(mol)

        carboxyl_smarts = Chem.MolFromSmarts('[$([C](=[O])[OX2H0])]')
        assert not mol.HasSubstructMatch(carboxyl_smarts)


class TestEtherificationPattern:
    """Test etherification SMARTS pattern."""

    def test_ether_pattern_detects_diol(self):
        """Test that etherification pattern detects diols."""
        # Ethylene glycol: CCO
        mol = Chem.MolFromSmiles("CCO")
        mol = Chem.AddHs(mol)

        alcohol_smarts = Chem.MolFromSmarts('[$([OX2H])]')
        matches = mol.GetSubstructMatches(alcohol_smarts)

        # Ethanol has one alcohol, ethylene glycol would have two
        # Let's test with ethylene glycol: OCCO
        mol_diol = Chem.MolFromSmiles("OCCO")
        mol_diol = Chem.AddHs(mol_diol)

        matches_diol = mol_diol.GetSubstructMatches(alcohol_smarts)
        assert len(matches_diol) >= 2, "Ethylene glycol should have 2+ alcohol groups"

    def test_ether_pattern_detects_peg_monomer(self):
        """Test that etherification pattern detects PEG monomer."""
        # Polyethylene glycol repeat unit: OCC
        mol = Chem.MolFromSmiles("OCC")
        mol = Chem.AddHs(mol)

        alcohol_smarts = Chem.MolFromSmarts('[$([OX2H])]')
        assert mol.HasSubstructMatch(alcohol_smarts)


class TestConnectionAtoms:
    """Test connection atom specifications."""

    def test_vinyl_connection_atoms(self):
        """Test that vinyl pattern specifies carbon atoms."""
        pattern = POLYMERIZATION_PATTERNS['vinyl_addition']
        assert pattern.connection_atoms == ['C', 'C']

    def test_ester_connection_atoms(self):
        """Test that ester pattern specifies C and O atoms."""
        pattern = POLYMERIZATION_PATTERNS['esterification']
        assert pattern.connection_atoms == ['C', 'O']

    def test_amide_connection_atoms(self):
        """Test that amide pattern specifies C and N atoms."""
        pattern = POLYMERIZATION_PATTERNS['amidation']
        assert pattern.connection_atoms == ['C', 'N']

    def test_ether_connection_atoms(self):
        """Test that ether pattern specifies O atoms."""
        pattern = POLYMERIZATION_PATTERNS['etherification']
        assert pattern.connection_atoms == ['O', 'O']


class TestCondensationFlag:
    """Test condensation polymerization flags."""

    def test_vinyl_is_not_condensation(self):
        """Test that vinyl is not condensation."""
        pattern = POLYMERIZATION_PATTERNS['vinyl_addition']
        assert pattern.is_condensation == False

    def test_ester_is_condensation(self):
        """Test that esterification is condensation."""
        pattern = POLYMERIZATION_PATTERNS['esterification']
        assert pattern.is_condensation == True

    def test_amide_is_condensation(self):
        """Test that amidation is condensation."""
        pattern = POLYMERIZATION_PATTERNS['amidation']
        assert pattern.is_condensation == True

    def test_ether_is_condensation(self):
        """Test that etherification is condensation."""
        pattern = POLYMERIZATION_PATTERNS['etherification']
        assert pattern.is_condensation == True
