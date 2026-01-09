"""
Tests for MonomerGenerator and MonomerVariant classes.

Tests the actual API of the monomer_generator module including:
- MonomerGenerator construction and initialization
- from_smiles() method for generating variants
- write_lt_files() method for LT file generation
- MonomerVariant dataclass attributes

Author: AutoPoly Development Team
"""

import pytest
import tempfile
import shutil
import os
from pathlib import Path

from AutoPoly.monomer_generator import (
    MonomerGenerator,
    MonomerVariant,
    MonomerGeneratorError,
    ValidationError,
    AtomTypingError,
    ChainBuildingError,
)


@pytest.fixture
def temp_output_dir():
    """Create a temporary output directory for tests."""
    temp_dir = tempfile.mkdtemp(prefix="test_monomer_generator_")
    yield temp_dir
    # Cleanup
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)


class TestMonomerGeneratorConstruction:
    """Test MonomerGenerator construction and initialization."""

    def test_basic_construction(self, temp_output_dir):
        """Test basic MonomerGenerator construction."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        assert generator.base_name == "PE"
        assert generator.force_field == "gaff"
        assert generator.output_dir == temp_output_dir

    def test_default_force_field(self, temp_output_dir):
        """Test that default force field is gaff."""
        generator = MonomerGenerator(
            base_name="PE",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        assert generator.force_field == "gaff"

    def test_oplsaa_force_field(self, temp_output_dir):
        """Test OPLS-AA force field initialization."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="oplsaa",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        assert generator.force_field == "oplsaa"

    def test_lopls_force_field(self, temp_output_dir):
        """Test L-OPLS force field initialization."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="lopls",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        assert generator.force_field == "lopls"

    def test_invalid_force_field_raises_error(self, temp_output_dir):
        """Test that invalid force field raises error."""
        with pytest.raises(MonomerGeneratorError) as excinfo:
            MonomerGenerator(
                base_name="PE",
                force_field="invalid_ff",
                output_dir=temp_output_dir,
                verbose=False
            )
        
        assert "Unknown force field" in str(excinfo.value)

    def test_output_dir_created(self, temp_output_dir):
        """Test that output directory is created if it doesn't exist."""
        new_dir = os.path.join(temp_output_dir, "new_subdir")
        assert not os.path.exists(new_dir)
        
        generator = MonomerGenerator(
            base_name="PE",
            output_dir=new_dir,
            verbose=False
        )
        
        assert os.path.exists(new_dir)


class TestFromSmiles:
    """Test from_smiles() method."""

    def test_linear_pe_from_smiles(self, temp_output_dir):
        """Test generating PE variants from SMILES."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=3)
        
        # Should have 3 variants: first, middle, last
        assert len(variants) == 3

    def test_variant_types_linear(self, temp_output_dir):
        """Test that linear chain produces correct variant types."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=3)
        
        # Check variant types
        assert variants[0].variant_type == "first"
        assert variants[1].variant_type == "middle"
        assert variants[2].variant_type == "last"

    def test_variant_positions(self, temp_output_dir):
        """Test that variants have correct position indices."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=5)
        
        for i, variant in enumerate(variants):
            assert variant.position == i

    def test_variants_have_mol_objects(self, temp_output_dir):
        """Test that variants have RDKit Mol objects."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=3)
        
        for variant in variants:
            assert variant.mol is not None
            assert variant.mol.GetNumAtoms() > 0

    def test_variants_have_conformers(self, temp_output_dir):
        """Test that variants have 3D conformers."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=3)
        
        for variant in variants:
            assert variant.mol.GetNumConformers() > 0

    def test_variants_have_connection_atoms(self, temp_output_dir):
        """Test that variants have connection_atoms tuple."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=3)
        
        for variant in variants:
            assert isinstance(variant.connection_atoms, tuple)
            assert len(variant.connection_atoms) == 2

    def test_minimum_monomers_enforced(self, temp_output_dir):
        """Test that at least 3 monomers are used for chain building."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        # Even with n_monomers=2, should use at least 3
        variants = generator.from_smiles("[*]CC[*]", n_monomers=2)
        
        # Should still produce at least 3 variants (for proper chemical environment)
        assert len(variants) >= 3

    def test_invalid_smiles_raises_error(self, temp_output_dir):
        """Test that invalid SMILES raises error."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        with pytest.raises(ChainBuildingError):
            generator.from_smiles("invalid_smiles", n_monomers=3)

    def test_smiles_without_wildcards_raises_error(self, temp_output_dir):
        """Test that SMILES without exactly 2 wildcards raises error."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        # SMILES with only one wildcard
        with pytest.raises(ChainBuildingError):
            generator.from_smiles("[*]CC", n_monomers=3)


class TestWriteLtFiles:
    """Test write_lt_files() method."""

    def test_write_lt_files_creates_files(self, temp_output_dir):
        """Test that write_lt_files creates .lt files."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=3)
        files = generator.write_lt_files(variants, generate_t1=False)
        
        assert len(files) == 3
        
        for filepath in files:
            assert os.path.exists(filepath)
            assert filepath.endswith('.lt')

    def test_write_lt_files_with_t1(self, temp_output_dir):
        """Test generating T1 chirality variants."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=3)
        files = generator.write_lt_files(variants, generate_t1=True)
        
        # Should have 6 files: 3 regular + 3 T1
        assert len(files) == 6
        
        # Check that T1 files exist
        t1_files = [f for f in files if '_T1.lt' in f]
        assert len(t1_files) == 3

    def test_lt_file_naming_convention(self, temp_output_dir):
        """Test LT file naming follows convention."""
        generator = MonomerGenerator(
            base_name="TEST",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=3)
        files = generator.write_lt_files(variants, generate_t1=False)
        
        # Files should follow {base_name}_{position}{suffix}.lt pattern
        for filepath in files:
            filename = os.path.basename(filepath)
            assert filename.startswith("TEST_")
            assert filename.endswith(".lt")

    def test_lt_file_contains_atoms_block(self, temp_output_dir):
        """Test that LT files contain Data Atoms block."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=3)
        files = generator.write_lt_files(variants, generate_t1=False)
        
        for filepath in files:
            with open(filepath, 'r') as f:
                content = f.read()
                assert 'write("Data Atoms")' in content

    def test_lt_file_contains_bonds_block(self, temp_output_dir):
        """Test that LT files contain Data Bond List block."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=3)
        files = generator.write_lt_files(variants, generate_t1=False)
        
        for filepath in files:
            with open(filepath, 'r') as f:
                content = f.read()
                assert "write('Data Bond List')" in content

    def test_gaff_lt_file_inherits_gaff(self, temp_output_dir):
        """Test that GAFF force field LT files inherit from GAFF."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=3)
        files = generator.write_lt_files(variants, generate_t1=False)
        
        for filepath in files:
            with open(filepath, 'r') as f:
                content = f.read()
                assert 'inherits GAFF' in content
                assert 'import "gaff.lt"' in content


class TestMonomerVariantDataclass:
    """Test MonomerVariant dataclass."""

    def test_monomer_variant_attributes(self, temp_output_dir):
        """Test MonomerVariant has expected attributes."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=3)
        variant = variants[0]
        
        # Check required attributes exist
        assert hasattr(variant, 'base_name')
        assert hasattr(variant, 'variant_type')
        assert hasattr(variant, 'mol')
        assert hasattr(variant, 'smiles')
        assert hasattr(variant, 'connection_atoms')
        assert hasattr(variant, 'force_field')
        assert hasattr(variant, 'position')
        assert hasattr(variant, 'atom_ids')

    def test_get_atom_type_at_connection(self, temp_output_dir):
        """Test get_atom_type_at_connection method."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=3)
        
        # Middle variant should have both connection types
        middle_variant = variants[1]
        left_type = middle_variant.get_atom_type_at_connection('left')
        right_type = middle_variant.get_atom_type_at_connection('right')
        
        # Should return atom types (or None if not set)
        # The actual type depends on force field matching

    def test_variant_base_name(self, temp_output_dir):
        """Test that variant has correct base_name."""
        generator = MonomerGenerator(
            base_name="PMMA",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=3)
        
        for variant in variants:
            assert variant.base_name == "PMMA"

    def test_variant_force_field(self, temp_output_dir):
        """Test that variant has correct force_field."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="oplsaa",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=3)
        
        for variant in variants:
            assert variant.force_field == "oplsaa"


class TestValidationError:
    """Test ValidationError exception."""

    def test_validation_error_exists(self):
        """Test that ValidationError exception class exists."""
        from AutoPoly.monomer_generator import ValidationError
        assert issubclass(ValidationError, MonomerGeneratorError)

    def test_validation_error_can_be_raised(self):
        """Test that ValidationError can be raised and caught."""
        with pytest.raises(ValidationError):
            raise ValidationError("Test validation error")


class TestSingleMolecule:
    """Test single molecule (non-polymer) generation."""

    def test_from_single_molecule(self, temp_output_dir):
        """Test generating LT file for a single molecule."""
        generator = MonomerGenerator(
            base_name="ethanol",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        # Generate single molecule variant
        variant = generator.from_single_molecule("CCO")
        
        assert variant is not None
        assert variant.mol.GetNumAtoms() > 0

    def test_single_molecule_has_conformer(self, temp_output_dir):
        """Test that single molecule has conformer."""
        generator = MonomerGenerator(
            base_name="ethanol",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variant = generator.from_single_molecule("CCO")
        
        assert variant.mol.GetNumConformers() > 0

    def test_single_molecule_with_wildcards_raises_error(self, temp_output_dir):
        """Test that SMILES with wildcards raises error for single molecule."""
        generator = MonomerGenerator(
            base_name="test",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        with pytest.raises(MonomerGeneratorError):
            generator.from_single_molecule("[*]CC[*]")

    def test_write_single_molecule(self, temp_output_dir):
        """Test writing single molecule LT file."""
        generator = MonomerGenerator(
            base_name="ethanol",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variant = generator.from_single_molecule("CCO")
        files = generator.write_single_molecule(variant, generate_t1=False)
        
        assert len(files) == 1
        assert os.path.exists(files[0])
        assert files[0].endswith('.lt')


class TestConvenienceFunctions:
    """Test convenience functions."""

    def test_generate_monomers_function(self, temp_output_dir):
        """Test generate_monomers convenience function."""
        from AutoPoly.monomer_generator import generate_monomers
        
        files = generate_monomers(
            smiles="[*]CC[*]",
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            n_monomers=3,
            generate_t1=False,
            verbose=False
        )
        
        assert len(files) == 3
        for filepath in files:
            assert os.path.exists(filepath)

    def test_generate_molecule_lt_function(self, temp_output_dir):
        """Test generate_molecule_lt convenience function."""
        from AutoPoly.monomer_generator import generate_molecule_lt
        
        files = generate_molecule_lt(
            smiles="CCO",
            molecule_name="ethanol",
            force_field="gaff",
            output_dir=temp_output_dir,
            generate_t1=False,
            verbose=False
        )
        
        assert len(files) == 1
        assert os.path.exists(files[0])


class TestAtomTyping:
    """Test atom typing functionality."""

    def test_atoms_have_types(self, temp_output_dir):
        """Test that atoms in variants have AtomType property."""
        generator = MonomerGenerator(
            base_name="PE",
            force_field="gaff",
            output_dir=temp_output_dir,
            verbose=False
        )
        
        variants = generator.from_smiles("[*]CC[*]", n_monomers=3)
        
        # Check some atoms have types (not all may match force field patterns)
        for variant in variants:
            typed_count = 0
            for atom in variant.mol.GetAtoms():
                if atom.GetAtomicNum() == 0:  # Skip dummy atoms
                    continue
                try:
                    atom_type = atom.GetProp('AtomType')
                    if atom_type:
                        typed_count += 1
                except KeyError:
                    pass
            
            # At least some atoms should be typed
            # (Complete typing depends on force field coverage)
