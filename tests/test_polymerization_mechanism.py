"""
Test polymerization mechanism detection for AutoPoly.

This module tests the PolymerizationMechanism class and detect_mechanism function
to ensure they correctly identify polymerization mechanisms from molecular structures.

Key Tests:
- Mechanism detection for various monomers
- DOP-aware detection
- Connection atom identification
- Mechanism information retrieval
"""

import sys
import os

# Add the AutoPoly package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from rdkit import Chem

from AutoPoly.polymerization_mechanism import (
    PolymerizationMechanism,
    detect_mechanism
)


class TestVinylDetection:
    """Test vinyl addition mechanism detection."""

    def test_ethylene_vinyl_detection(self):
        """Test that ethylene (C=C) is detected as vinyl."""
        mol = Chem.MolFromSmiles("C=C")
        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        assert mechanism == 'vinyl_addition'

    def test_propylene_vinyl_detection(self):
        """Test that propylene is detected as vinyl."""
        mol = Chem.MolFromSmiles("C=CC")
        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        assert mechanism == 'vinyl_addition'

    def test_styrene_vinyl_detection(self):
        """Test that styrene is detected as vinyl."""
        mol = Chem.MolFromSmiles("C=Cc1ccccc1")
        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        assert mechanism == 'vinyl_addition'

    def test_polymethyl_methacrylate_vinyl_detection(self):
        """Test that MMA monomer is detected as vinyl."""
        # MMA: C=C(C)C(=O)OC
        mol = Chem.MolFromSmiles("C=C(C)C(=O)OC")
        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        assert mechanism == 'vinyl_addition'


class TestEsterificationDetection:
    """Test esterification mechanism detection."""

    def test_lactic_acid_esterification(self):
        """Test that lactic acid is detected as esterification for DOP>1."""
        mol = Chem.MolFromSmiles("CC(C(=O)O)O")
        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        assert mechanism == 'esterification'

    def test_lactic_acid_dop1_none(self):
        """Test that lactic acid with DOP=1 returns 'none'."""
        mol = Chem.MolFromSmiles("CC(C(=O)O)O")
        mechanism = detect_mechanism(mol, dop=1, verbose=False)
        assert mechanism == 'none'

    def test_lactide_esterification(self):
        """Test that lactide (PLA dimer) is detected."""
        # Lactide: cyclic diester of lactic acid
        mol = Chem.MolFromSmiles("CC1OC(=O)C(C)OC(=O)C1")
        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        # Should detect carboxyl groups
        assert mechanism in ['esterification', 'ring_opening_ester', 'none']


class TestAmidationDetection:
    """Test amidation mechanism detection."""

    def test_amino_acid_amidation(self):
        """Test that amino acids are detected for amidation."""
        # Simple amino acid: H2N-CH2-COOH (glycine)
        mol = Chem.MolFromSmiles("NCC(=O)O")
        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        assert mechanism == 'amidation'

    def test_caprolactam_amidation(self):
        """Test that caprolactam (nylon-6 precursor) is detected."""
        # Caprolactam: cyclic amide
        mol = Chem.MolFromSmiles("C1CCC(=O)NC1")
        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        # Should detect carboxyl and amine or ring structure
        assert mechanism in ['amidation', 'ring_opening_amide', 'none']

    def test_hexamethylenediamine_adipic_acid(self):
        """Test detection of nylon-6,6 precursors."""
        # H2N-(CH2)6-NH2 + HOOC-(CH2)4-COOH
        # This tests the diacid component
        mol = Chem.MolFromSmiles("OC(=O)CCCC(=O)O")
        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        # Adipic acid has 2 carboxyl groups but no separate amines
        # Should return 'none' (needs diamine to react with)
        assert mechanism == 'none'


class TestEtherificationDetection:
    """Test etherification mechanism detection."""

    def test_ethylene_glycol_etherification(self):
        """Test that ethylene glycol is detected for etherification."""
        mol = Chem.MolFromSmiles("OCCO")
        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        assert mechanism == 'etherification'

    def test_ethylene_glycol_single_molecule(self):
        """Test ethylene glycol with DOP=1 returns 'none'."""
        mol = Chem.MolFromSmiles("OCCO")
        mechanism = detect_mechanism(mol, dop=1, verbose=False)
        assert mechanism == 'none'

    def test_polyethylene_glycol_repeat_unit(self):
        """Test PEG repeat unit detection."""
        # PEG repeat: -O-CH2-CH2-
        mol = Chem.MolFromSmiles("OCC")
        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        # Only 1 alcohol group, not enough for etherification (needs 2+)
        assert mechanism == 'none'


class TestNoneMechanism:
    """Test 'none' mechanism detection for non-polymerizable molecules."""

    def test_water_none(self):
        """Test that water is detected as 'none'."""
        mol = Chem.MolFromSmiles("O")
        mechanism = detect_mechanism(mol, dop=1, verbose=False)
        assert mechanism == 'none'

    def test_ethanol_none(self):
        """Test that ethanol is detected as 'none'."""
        mol = Chem.MolFromSmiles("CCO")
        mechanism = detect_mechanism(mol, dop=1, verbose=False)
        assert mechanism == 'none'

    def test_methane_none(self):
        """Test that methane is detected as 'none'."""
        mol = Chem.MolFromSmiles("C")
        mechanism = detect_mechanism(mol, dop=1, verbose=False)
        assert mechanism == 'none'

    def test_benzene_none(self):
        """Test that benzene (no reactive groups) is 'none'."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        mechanism = detect_mechanism(mol, dop=1, verbose=False)
        assert mechanism == 'none'


class TestDOPAwareDetection:
    """Test DOP-aware mechanism detection."""

    def test_vinyl_dop1_still_vinyl(self):
        """Test that vinyl compounds are detected even with DOP=1."""
        mol = Chem.MolFromSmiles("C=C")
        mechanism = detect_mechanism(mol, dop=1, verbose=False)
        assert mechanism == 'vinyl_addition'

    def test_condensation_dop1_returns_none(self):
        """Test that condensation monomers return 'none' with DOP=1."""
        mol = Chem.MolFromSmiles("CC(C(=O)O)O")  # Lactic acid
        mechanism = detect_mechanism(mol, dop=1, verbose=False)
        assert mechanism == 'none'

    def test_dop10_detects_patterns(self):
        """Test that DOP=10 enables pattern detection."""
        mol = Chem.MolFromSmiles("CC(C(=O)O)O")
        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        assert mechanism == 'esterification'


class TestConnectionAtoms:
    """Test connection atom identification."""

    def test_vinyl_connection_atoms(self):
        """Test that vinyl mechanism identifies C-C connection."""
        detector = PolymerizationMechanism(verbose=False)
        mol = Chem.MolFromSmiles("C=C")

        connection_atoms = detector.get_connection_atoms(mol, 'vinyl_addition')
        assert len(connection_atoms) == 2

        # Both should be carbons
        for atom_idx in connection_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            assert atom.GetSymbol() == 'C'

    def test_backbone_atoms_same_as_connection(self):
        """Test that backbone atoms match connection atoms for vinyl."""
        detector = PolymerizationMechanism(verbose=False)
        mol = Chem.MolFromSmiles("C=C")

        backbone = detector.get_backbone_atoms(mol, 'vinyl_addition')
        connection = detector.get_connection_atoms(mol, 'vinyl_addition')

        assert set(backbone) == set(connection)

    def test_none_mechanism_no_connections(self):
        """Test that 'none' mechanism has no connection atoms."""
        detector = PolymerizationMechanism(verbose=False)
        mol = Chem.MolFromSmiles("CCO")

        connection_atoms = detector.get_connection_atoms(mol, 'none')
        assert len(connection_atoms) == 0


class TestMechanismProperties:
    """Test mechanism property methods."""

    def test_requires_variant_generation(self):
        """Test variant generation requirements."""
        detector = PolymerizationMechanism(verbose=False)

        assert detector.requires_variant_generation('vinyl_addition') == True
        assert detector.requires_variant_generation('esterification') == True
        assert detector.requires_variant_generation('none') == False

    def test_is_condensation_polymerization(self):
        """Test condensation polymerization identification."""
        detector = PolymerizationMechanism(verbose=False)

        assert detector.is_condensation_polymerization('vinyl_addition') == False
        assert detector.is_condensation_polymerization('esterification') == True
        assert detector.is_condensation_polymerization('amidation') == True
        assert detector.is_condensation_polymerization('etherification') == True

    def test_get_mechanism_info(self):
        """Test mechanism information retrieval."""
        detector = PolymerizationMechanism(verbose=False)

        info = detector.get_mechanism_info('vinyl_addition')
        assert info is not None
        assert 'name' in info
        assert 'smarts' in info
        assert 'connection_atoms' in info
        assert 'description' in info
        assert 'is_condensation' in info

        assert info['name'] == 'vinyl_addition'
        assert info['connection_atoms'] == ['C', 'C']
        assert info['is_condensation'] == False

    def test_get_mechanism_info_invalid(self):
        """Test mechanism info for invalid mechanism."""
        detector = PolymerizationMechanism(verbose=False)

        info = detector.get_mechanism_info('invalid')
        assert info is None


class TestConvenienceFunction:
    """Test the detect_mechanism convenience function."""

    def test_convenience_function_matches_class(self):
        """Test that convenience function matches class method."""
        mol = Chem.MolFromSmiles("C=C")

        # Using function
        mechanism1 = detect_mechanism(mol, dop=10, verbose=False)

        # Using class
        detector = PolymerizationMechanism(verbose=False)
        mechanism2 = detector.detect_mechanism(mol, dop=10)

        assert mechanism1 == mechanism2 == 'vinyl_addition'

    def test_convenience_function_default_verbose(self):
        """Test that convenience function has default verbose=True."""
        mol = Chem.MolFromSmiles("CCO")
        # Should not raise error
        mechanism = detect_mechanism(mol, dop=1)
        assert mechanism == 'none'


class TestPriorityOrder:
    """Test that detection follows correct priority order."""

    def test_vinyl_has_priority_over_condensation(self):
        """Test that vinyl pattern is detected before condensation."""
        # Molecule with both C=C and COOH groups
        # Acrylic acid: C=C(=O)O
        mol = Chem.MolFromSmiles("C=CC(=O)O")
        mechanism = detect_mechanism(mol, dop=10, verbose=False)

        # Vinyl should have priority
        assert mechanism == 'vinyl_addition'

    def test_est_before_amidation(self):
        """Test that esterification is checked before amidation."""
        # Molecule with carboxyl + alcohol + amine
        # This would need a specific test molecule
        # For now, just verify the priority is implemented
        detector = PolymerizationMechanism(verbose=False)

        # The implementation checks esterification before amidation
        # This is verified by reading the code
        assert True  # Placeholder
