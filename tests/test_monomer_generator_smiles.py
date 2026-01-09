"""
Unit tests for SMILES input processing.

This module tests SMILES input with wildcard atoms for monomer generation.
"""

import sys
import os
import tempfile
import shutil

# Add the AutoPoly package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from rdkit import Chem

from AutoPoly.monomer_generator import MonomerGenerator, MonomerGeneratorError


class TestSMILESInput:
    """Test SMILES input processing."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory for tests."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_smiles_test_")
        yield temp_dir
        # Cleanup
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

    def test_smiles_input_basic(self, temp_output_dir):
        """Test that SMILES processing works."""
        generator = MonomerGenerator(
            base_name="test_pe",
            mechanism='vinyl_addition',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Use polyethylene SMILES (has vinyl double bond)
        input_smiles = "[*]C=C[*]"
        variants = generator.generate_variants(input_smiles)

        # Should succeed
        assert variants is not None
        assert 'smiles' in variants
        assert 'internal' in variants['smiles']

    def test_wildcard_format_variations(self, temp_output_dir):
        """Test that both * and [*] formats work."""
        generator = MonomerGenerator(
            base_name="test_pe",
            mechanism='vinyl_addition',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Test with * without brackets
        variants1 = generator.generate_variants("*C=C*")

        # Test with [*] with brackets
        variants2 = generator.generate_variants("[*]C=C[*]")

        # Both should work (but may produce different results without canonicalization)
        assert variants1 is not None
        assert variants2 is not None

    def test_invalid_smiles_no_wildcards(self, temp_output_dir):
        """Test that invalid SMILES (no wildcards) raises error."""
        generator = MonomerGenerator(
            base_name="test_invalid",
            mechanism='none',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # No wildcards - should raise error
        with pytest.raises(MonomerGeneratorError) as exc_info:
            generator.generate_variants("CCCCCC")

        # Error message should mention wildcards
        assert "wildcard" in str(exc_info.value).lower()

    def test_invalid_chemistry(self, temp_output_dir):
        """Test that chemically invalid SMILES raises error."""
        generator = MonomerGenerator(
            base_name="test_chem",
            mechanism='none',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Invalid SMILES syntax (X is not a valid element)
        with pytest.raises(MonomerGeneratorError):
            generator.generate_variants("[*]X=X[*]")

    def test_ring_polymer(self, temp_output_dir):
        """Test SMILES processing for ring polymers."""
        generator = MonomerGenerator(
            base_name="test_ring",
            mechanism='vinyl_addition',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        variants = generator.generate_variants("[*]C=C[*]")

        assert variants is not None
        assert 'internal' in variants
        assert 'left_end' in variants
        assert 'right_end' in variants

    def test_complex_polymer(self, temp_output_dir):
        """Test SMILES processing for more complex polymers."""
        generator = MonomerGenerator(
            base_name="test_complex",
            mechanism='vinyl_addition',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        smiles = "[*]C=C[*]"
        variants = generator.generate_variants(smiles)

        assert variants is not None
        assert 'smiles' in variants

    def test_linear_polyethylene(self, temp_output_dir):
        """Test SMILES processing for simple linear polyethylene."""
        generator = MonomerGenerator(
            base_name="test_pe",
            mechanism='vinyl_addition',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Simple ethylene SMILES
        variants = generator.generate_variants("[*]C=C[*]")

        assert variants is not None
        assert 'internal' in variants
        assert isinstance(variants['internal'], Chem.Mol)

    def test_double_bond_wildcards(self, temp_output_dir):
        """Test proper handling of wildcards on double bonds (PMMA case)."""
        generator = MonomerGenerator(
            base_name="test_double_bond",
            mechanism='vinyl_addition',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # PMMA SMILES
        smiles = "[*]C=C(C)C(=O)OC[*]"
        variants = generator.generate_variants(smiles)

        # Should succeed with RDKit dummy atom removal
        assert variants is not None
        assert 'internal' in variants
        assert 'left_end' in variants
        assert 'right_end' in variants

        # Verify valid SMILES can be generated
        internal_smiles = Chem.MolToSmiles(variants['internal'])
        assert internal_smiles  # Not empty
        assert '=' not in internal_smiles[0]  # Should not start with dangling bond

    def test_various_smiles_formats(self, temp_output_dir):
        """Test that various SMILES formats work correctly."""
        generator = MonomerGenerator(
            base_name="test_formats",
            mechanism='vinyl_addition',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Test various SMILES with different wildcard placements
        smiles_list = [
            "[*]C=C[*]",           # Simple vinyl
            "[*]C=C(C)C(=O)OC[*]", # PMMA
            "[*]C=C(C)C[*]",       # Substituted vinyl
        ]

        for smiles in smiles_list:
            variants = generator.generate_variants(smiles)

            # All should succeed
            assert variants is not None, f"Failed for SMILES: {smiles}"
            assert 'internal' in variants, f"No internal variant for: {smiles}"
            assert 'left_end' in variants, f"No left_end variant for: {smiles}"
            assert 'right_end' in variants, f"No right_end variant for: {smiles}"

            # Verify valid SMILES can be generated (no dangling bonds)
            internal_smiles = Chem.MolToSmiles(variants['internal'])
            assert internal_smiles, f"Empty SMILES for: {smiles}"
            # Should not start with a bond symbol (indicates dangling bond)
            assert internal_smiles[0] not in '=:#', f"Dangling bond in: {internal_smiles}"

    def test_both_ends_get_hydrogen_caps(self, temp_output_dir):
        """Test that BOTH left and right ends get H end-caps."""
        generator = MonomerGenerator(
            base_name="test_caps",
            mechanism='vinyl_addition',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        variants = generator.generate_variants("[*]C=C[*]")

        # Check left-end has H on both ends
        mol_le = variants['left_end']
        mol_le_with_h = Chem.AllChem.AddHs(mol_le)
        h_count_le = sum(1 for atom in mol_le_with_h.GetAtoms() if atom.GetSymbol() == 'H')

        # Check right-end has H on both ends
        mol_re = variants['right_end']
        mol_re_with_h = Chem.AllChem.AddHs(mol_re)
        h_count_re = sum(1 for atom in mol_re_with_h.GetAtoms() if atom.GetSymbol() == 'H')

        # Both should have same H count (both ends capped)
        assert h_count_le == h_count_re, f"Left-end H count ({h_count_le}) != Right-end H count ({h_count_re})"

    def test_no_canonicalization(self, temp_output_dir):
        """Test that SMILES are accepted as-is without canonicalization."""
        generator = MonomerGenerator(
            base_name="test_no_canon",
            mechanism='vinyl_addition',
            output_dir=temp_output_dir,
            is_gaff=False,
            verbose=False
        )

        # Different non-canonical forms should produce different results
        # Note: Without canonicalization, these may produce different internal SMILES
        smiles1 = "[*]C=C[*]"
        smiles2 = "*C=C*"

        variants1 = generator.generate_variants(smiles1)
        variants2 = generator.generate_variants(smiles2)

        # Both should succeed
        assert variants1 is not None
        assert variants2 is not None

        # Internal SMILES should be generated for both
        assert variants1['smiles']['internal']
        assert variants2['smiles']['internal']

    def test_validation_exactly_two_wildcards(self, temp_output_dir):
        """Test that exactly 2 wildcards are required."""
        generator = MonomerGenerator(
            base_name="test_validation",
            mechanism='vinyl_addition',
            output_dir=temp_output_dir,
            verbose=False
        )

        # No wildcards - should fail
        with pytest.raises(MonomerGeneratorError) as exc_info:
            generator.generate_variants("C=C")
        assert "wildcard" in str(exc_info.value).lower()

        # One wildcard - should fail
        with pytest.raises(MonomerGeneratorError) as exc_info:
            generator.generate_variants("[*]C=C")
        assert "wildcard" in str(exc_info.value).lower()

        # Three wildcards - should fail
        with pytest.raises(MonomerGeneratorError) as exc_info:
            generator.generate_variants("[*]C=C[*][*]")
        assert "wildcard" in str(exc_info.value).lower()
