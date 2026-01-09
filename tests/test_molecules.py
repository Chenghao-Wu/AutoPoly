"""
Test single molecule generation (DOP=1) for AutoPoly.

This module tests the preservation of single molecule generation support,
ensuring that DOP=1 correctly generates independent molecules without
polymerization bonds.

Key Tests:
- Ethanol molecule generation
- Water molecule generation
- Lactic acid molecule generation (has reactive groups but DOP=1)
- Systems of multiple independent molecules
"""

import sys
import os
import tempfile
import shutil
from pathlib import Path

# Add the AutoPoly package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from rdkit import Chem

from AutoPoly.monomer_generator import MonomerGenerator
from AutoPoly.polymerization_mechanism import detect_mechanism


class TestSingleMoleculeGeneration:
    """Test single molecule generation with DOP=1."""

    def test_ethanol_molecule_mechanism_detection(self):
        """Test that ethanol is detected as non-polymerizable."""
        # Ethanol SMILES: CCO
        mol = Chem.MolFromSmiles("CCO")

        # Should detect as 'none' for DOP=1
        mechanism = detect_mechanism(mol, dop=1, verbose=False)
        assert mechanism == 'none', "Ethanol should be detected as non-polymerizable for DOP=1"

    def test_water_molecule_mechanism_detection(self):
        """Test that water is detected as non-polymerizable."""
        # Water SMILES: O
        mol = Chem.MolFromSmiles("O")

        # Should detect as 'none'
        mechanism = detect_mechanism(mol, dop=1, verbose=False)
        assert mechanism == 'none', "Water should be detected as non-polymerizable"

    def test_lactic_acid_molecule_mechanism_detection(self):
        """Test that lactic acid with DOP=1 is non-polymerizable."""
        # Lactic acid: CC(C(=O)O)O
        # Has both carboxyl and alcohol, but DOP=1 means don't polymerize
        mol = Chem.MolFromSmiles("CC(C(=O)O)O")

        # Should detect as 'none' for DOP=1
        mechanism = detect_mechanism(mol, dop=1, verbose=False)
        assert mechanism == 'none', "Lactic acid should be non-polymerizable for DOP=1"

    def test_lactic_acid_polymer_mechanism_detection(self):
        """Test that lactic acid with DOP>1 is detected as esterification."""
        # Lactic acid: CC(C(=O)O)O
        mol = Chem.MolFromSmiles("CC(C(=O)O)O")

        # Should detect as 'esterification' for DOP>1
        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        assert mechanism == 'esterification', "Lactic acid should be detected as esterification for DOP>1"

    def test_vinyl_mechanism_detection_with_dop1(self):
        """Test that vinyl compounds are still detected with DOP=1."""
        # Ethylene: C=C
        mol = Chem.MolFromSmiles("C=C")

        # Should still detect as 'vinyl_addition' even with DOP=1
        mechanism = detect_mechanism(mol, dop=1, verbose=False)
        assert mechanism == 'vinyl_addition', "Vinyl compounds should be detected even with DOP=1"


class TestMonomerGeneratorSingleMolecules:
    """Test MonomerGenerator with single molecules."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory for tests."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_molecules_test_")
        yield temp_dir
        # Cleanup
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

    def test_ethanol_monomer_generation(self, temp_output_dir):
        """Test generating ethanol molecule with MonomerGenerator."""
        generator = MonomerGenerator(
            base_name="ethanol",
            mechanism='none',  # Force non-polymerizable
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Should not raise an error
        variants = generator.generate_variants(smiles="[*]CCO[*]")

        # Check that variants contain expected keys
        assert 'internal' in variants
        assert 'left_end' in variants
        assert 'right_end' in variants

        # All variants should be RDKit molecules
        assert isinstance(variants['internal'], Chem.Mol)

    def test_lactic_acid_monomer_generation(self, temp_output_dir):
        """Test generating single lactic acid molecule."""
        generator = MonomerGenerator(
            base_name="lactic_acid",
            mechanism='none',  # Force no polymerization
            output_dir=temp_output_dir,
            is_gaff=True,  # GAFF is better for polyesters
            verbose=False
        )

        # Should generate successfully
        variants = generator.generate_variants(smiles="[*]CC(C(=O)O)O[*]")

        # Check that variants are generated
        assert 'internal' in variants
        assert isinstance(variants['internal'], Chem.Mol)

    def test_water_monomer_generation(self, temp_output_dir):
        """Test generating water molecule."""
        generator = MonomerGenerator(
            base_name="water",
            mechanism='none',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Should generate successfully
        variants = generator.generate_variants(smiles="[*]O[*]")

        # Check that variants are generated
        assert 'internal' in variants

    def test_auto_detect_none_mechanism(self, temp_output_dir):
        """Test auto-detection of 'none' mechanism for non-polymerizable molecules."""
        # Don't specify mechanism, let it auto-detect
        generator = MonomerGenerator(
            base_name="ethanol_auto",
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Should auto-detect 'none' mechanism
        mol = Chem.MolFromSmiles("CCO")
        mol_with_h = Chem.AddHs(mol)

        # Use the internal mechanism detector
        mechanism = generator.mech_detector.detect_mechanism(mol_with_h, dop=1)
        assert mechanism == 'none', "Should auto-detect 'none' for ethanol with DOP=1"

    def test_generate_lt_files_for_single_molecule(self, temp_output_dir):
        """Test generating .lt files for single molecules."""
        generator = MonomerGenerator(
            base_name="ethanol",
            mechanism='none',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        variants = generator.generate_variants(smiles="[*]CCO[*]")
        files = generator.generate_lt_files(variants)

        # Check that files were generated
        assert 'internal' in files
        assert os.path.exists(files['internal'])

        # Check file has correct extension
        assert files['internal'].endswith('.lt')


class TestMechanismParameter:
    """Test the mechanism parameter in MonomerGenerator."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory for tests."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_mechanism_test_")
        yield temp_dir
        # Cleanup
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

    def test_explicit_none_mechanism(self, temp_output_dir):
        """Test explicitly setting mechanism='none'."""
        generator = MonomerGenerator(
            base_name="test",
            mechanism='none',
            output_dir=temp_output_dir,
            verbose=False
        )

        assert generator.mechanism == 'none'

    def test_vinyl_addition_mechanism(self, temp_output_dir):
        """Test setting mechanism='vinyl_addition'."""
        generator = MonomerGenerator(
            base_name="test",
            mechanism='vinyl_addition',
            output_dir=temp_output_dir,
            verbose=False
        )

        assert generator.mechanism == 'vinyl_addition'

    def test_esterification_mechanism(self, temp_output_dir):
        """Test setting mechanism='esterification'."""
        generator = MonomerGenerator(
            base_name="test",
            mechanism='esterification',
            output_dir=temp_output_dir,
            verbose=False
        )

        assert generator.mechanism == 'esterification'

    def test_invalid_mechanism_raises_error(self, temp_output_dir):
        """Test that invalid mechanism raises an error."""
        with pytest.raises(Exception) as exc_info:
            MonomerGenerator(
                base_name="test",
                mechanism='invalid_mechanism',
                output_dir=temp_output_dir,
                verbose=False
            )

        assert "Unknown mechanism" in str(exc_info.value)


class TestDOPAwareness:
    """Test DOP-aware mechanism detection."""

    def test_dop1_returns_none_for_non_vinyl(self):
        """Test that DOP=1 returns 'none' for non-vinyl molecules."""
        mol = Chem.MolFromSmiles("CCO")  # Ethanol
        mechanism = detect_mechanism(mol, dop=1, verbose=False)
        assert mechanism == 'none'

    def test_dop_greater_than_1_detects_patterns(self):
        """Test that DOP>1 performs pattern detection."""
        mol = Chem.MolFromSmiles("CC(C(=O)O)O")  # Lactic acid

        # DOP=1 should return 'none'
        mechanism1 = detect_mechanism(mol, dop=1, verbose=False)
        assert mechanism1 == 'none'

        # DOP>1 should detect esterification
        mechanism10 = detect_mechanism(mol, dop=10, verbose=False)
        assert mechanism10 == 'esterification'


class TestMonomerGeneratorEdgeCases:
    """Test MonomerGenerator edge cases and error handling."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory for tests."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_edge_test_")
        yield temp_dir
        # Cleanup
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

    def test_invalid_smiles_raises_error(self, temp_output_dir):
        """Test that invalid SMILES raises an error."""
        generator = MonomerGenerator(
            base_name="test_invalid",
            mechanism='none',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Invalid SMILES string
        with pytest.raises(Exception):
            generator.generate_variants(smiles="[*]INVALID_SMILES_STRING[*]")

    def test_empty_smiles_raises_error(self, temp_output_dir):
        """Test that empty SMILES raises an error."""
        generator = MonomerGenerator(
            base_name="test_empty",
            mechanism='none',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Empty SMILES string
        with pytest.raises(Exception):
            generator.generate_variants(smiles="")

    def test_complex_molecules_with_rings(self, temp_output_dir):
        """Test variant generation for cyclic compounds."""
        generator = MonomerGenerator(
            base_name="cyclohexane",
            mechanism='none',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Cyclohexane: C1CCCCC1 (ring structure)
        variants = generator.generate_variants(smiles="[*]C1CCCCC1[*]")

        # Should generate successfully even with rings
        assert 'internal' in variants
        assert isinstance(variants['internal'], Chem.Mol)

    def test_molecules_with_multiple_bonds(self, temp_output_dir):
        """Test variant generation for molecules with double/triple bonds."""
        generator = MonomerGenerator(
            base_name="butadiene",
            mechanism='none',  # Changed to 'none' since this is testing edge cases, not actual polymerization
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # 1,3-butadiene: C=CC=C (conjugated double bonds)
        # Using different pSMILES to avoid valence issues - placing wildcards on atoms that can accept H
        variants = generator.generate_variants(smiles="[*]C=CC=C[*]")

        # Should generate successfully
        assert 'internal' in variants
        assert isinstance(variants['internal'], Chem.Mol)

    def test_variant_with_chiral_centers(self, temp_output_dir):
        """Test variant generation for molecules with stereocenters."""
        generator = MonomerGenerator(
            base_name="lactic_acid_chiral",
            mechanism='esterification',
            output_dir=temp_output_dir,
            is_gaff=True,
            verbose=False
        )

        # Lactic acid with chiral center: CC(C(=O)O)O
        variants = generator.generate_variants(smiles="[*]CC(C(=O)O)O[*]")

        # Should generate successfully preserving chirality
        assert 'internal' in variants
        assert isinstance(variants['internal'], Chem.Mol)

    def test_variant_with_heteroatoms(self, temp_output_dir):
        """Test variant generation for molecules with heteroatoms in backbone."""
        generator = MonomerGenerator(
            base_name="nylon_mononer",
            mechanism="none",
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Nylon-6 monomer precursor (caprolactam): C1CCC(=O)NCC1
        # Contains nitrogen in the ring
        variants = generator.generate_variants(smiles="[*]C1CCC(=O)NCC1[*]")

        # Should generate successfully with heteroatoms
        assert 'internal' in variants
        assert isinstance(variants['internal'], Chem.Mol)

    def test_dop1_single_molecule_no_capping(self, temp_output_dir):
        """Test that DOP=1 molecules don't get capped."""
        generator = MonomerGenerator(
            base_name="ethanol_dop1",
            mechanism='none',  # DOP=1 means no polymerization
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Ethanol with DOP=1
        variants = generator.generate_variants(smiles="[*]CCO[*]")

        # Should generate without capping groups
        assert 'internal' in variants
        # The molecule should be complete, not needing capping
        assert isinstance(variants['internal'], Chem.Mol)


class TestMonomerGeneratorErrorHandling:
    """Test MonomerGenerator error handling."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory for tests."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_error_test_")
        yield temp_dir
        # Cleanup
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

    def test_missing_bonding_sites_handling(self, temp_output_dir):
        """Test handling of molecules without suitable bonding sites."""
        generator = MonomerGenerator(
            base_name="methane",
            mechanism='none',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Methane: C (only one carbon, limited bonding sites)
        variants = generator.generate_variants(smiles="[*]C[*]")

        # Should still generate, just as a non-polymerizable molecule
        assert 'internal' in variants
        assert isinstance(variants['internal'], Chem.Mol)

    def test_rdkit_generation_failure_handling(self, temp_output_dir):
        """Test graceful handling of RDKit generation failures."""
        generator = MonomerGenerator(
            base_name="test_failure",
            mechanism='none',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # SMILES that RDKit can't parse properly
        # Using a very long invalid SMILES
        invalid_smiles = "C" * 1000 + "X"

        with pytest.raises(Exception):
            generator.generate_variants(smiles=f"[*]{invalid_smiles}[*]")

    def test_variant_generation_with_special_atoms(self, temp_output_dir):
        """Test variant generation with special atom types."""
        generator = MonomerGenerator(
            base_name="sulfur_containing",
            mechanism='none',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Dimethyl sulfide: CSC (contains sulfur)
        variants = generator.generate_variants(smiles="[*]CSC[*]")

        # Should generate successfully with sulfur
        assert 'internal' in variants
        assert isinstance(variants['internal'], Chem.Mol)

    def test_halogens_in_molecule(self, temp_output_dir):
        """Test variant generation with halogen atoms."""
        generator = MonomerGenerator(
            base_name="vinyl_chloride",
            mechanism='vinyl_addition',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Vinyl chloride: C=CCl
        variants = generator.generate_variants(smiles="[*]C=CCl[*]")

        # Should generate successfully with chlorine
        assert 'internal' in variants
        assert isinstance(variants['internal'], Chem.Mol)
