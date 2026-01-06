"""
Integration tests for Polylactic Acid (PLA) generation in AutoPoly.

These tests verify the complete workflow for generating PLA monomers and polymers
using the esterification mechanism with heteroatom (C, O) backbone support.

Key Tests:
- PLA monomer generation with esterification mechanism
- Mechanism detection for lactic acid
- Atom type assignment for ester linkages
- .lt file generation for PLA variants
"""

import os
import sys
import tempfile
import shutil
from pathlib import Path

# Add the AutoPoly package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

import pytest
from rdkit import Chem

from AutoPoly.monomer_generator import MonomerGenerator
from AutoPoly.polymerization_mechanism import detect_mechanism


@pytest.mark.slow
@pytest.mark.integration
class TestPLAMonomerGeneration:
    """Test PLA monomer generation workflow."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory for tests."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_pla_test_")
        yield temp_dir
        # Cleanup
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

    def test_lactic_acid_mechanism_detection(self):
        """Test that lactic acid is detected as esterification."""
        # Lactic acid: CC(C(=O)O)O
        mol = Chem.MolFromSmiles("CC(C(=O)O)O")

        # Should detect esterification for DOP>1
        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        assert mechanism == 'esterification'

    def test_lactic_acid_dop1_none(self):
        """Test that lactic acid returns 'none' for DOP=1."""
        mol = Chem.MolFromSmiles("CC(C(=O)O)O")
        mechanism = detect_mechanism(mol, dop=1, verbose=False)
        assert mechanism == 'none'

    def test_pla_monomer_generator_initialization(self, temp_output_dir):
        """Test MonomerGenerator initialization for PLA."""
        generator = MonomerGenerator(
            base_name="PLA",
            mechanism='esterification',
            output_dir=temp_output_dir,
            is_gaff=True,  # GAFF is better for polyesters
            verbose=False
        )

        assert generator.base_name == "PLA"
        assert generator.mechanism == 'esterification'
        assert generator.is_gaff == True

    def test_pla_variant_generation(self, temp_output_dir):
        """Test generating PLA monomer variants."""
        generator = MonomerGenerator(
            base_name="PLA",
            mechanism='esterification',
            output_dir=temp_output_dir,
            is_gaff=True,
            verbose=False
        )

        # Generate variants from lactic acid
        # Lactic acid dimer (lactide): CC1OC(=O)C(C)OC(=O)C1
        # Or use monomeric lactic acid: CC(C(=O)O)O
        variants = generator.generate_variants(smiles="CC(C(=O)O)O")

        # Check that all variants were generated
        assert 'internal' in variants
        assert 'left_end' in variants
        assert 'right_end' in variants
        assert 'smiles' in variants
        assert 'bonding_info' in variants

        # Check that variants are RDKit molecules
        assert isinstance(variants['internal'], Chem.Mol)
        assert isinstance(variants['left_end'], Chem.Mol)
        assert isinstance(variants['right_end'], Chem.Mol)

    def test_pla_lt_file_generation(self, temp_output_dir):
        """Test generating .lt files for PLA."""
        generator = MonomerGenerator(
            base_name="PLA",
            mechanism='esterification',
            output_dir=temp_output_dir,
            is_gaff=True,
            verbose=False
        )

        variants = generator.generate_variants(smiles="CC(C(=O)O)O")
        files = generator.generate_lt_files(variants, generate_t1=False)

        # Check that files were generated
        assert 'internal' in files
        assert 'left_end' in files
        assert 'right_end' in files

        # Check that files exist
        assert os.path.exists(files['internal'])
        assert os.path.exists(files['left_end'])
        assert os.path.exists(files['right_end'])

        # Check file extensions
        assert files['internal'].endswith('.lt')
        assert files['left_end'].endswith('.lt')
        assert files['right_end'].endswith('.lt')

    def test_pla_file_contents(self, temp_output_dir):
        """Test that generated .lt files have correct structure."""
        generator = MonomerGenerator(
            base_name="PLA",
            mechanism='esterification',
            output_dir=temp_output_dir,
            is_gaff=True,
            verbose=False
        )

        variants = generator.generate_variants(smiles="CC(C(=O)O)O")
        files = generator.generate_lt_files(variants, generate_t1=False)

        # Read the internal variant file
        with open(files['internal'], 'r') as f:
            content = f.read()

        # Check for GAFF import
        assert 'import "gaff.lt"' in content

        # Check for PLA class definition
        assert 'PLAi inherits GAFF' in content

        # Check for atom definitions
        assert 'write("Data Atoms")' in content

        # Check for bond definitions
        assert 'write(\'Data Bond List\')' in content


@pytest.mark.slow
@pytest.mark.integration
class TestPLAHeteroatomBackbone:
    """Test that PLA uses heteroatom (C, O) backbone correctly."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory for tests."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_pla_heteroatom_test_")
        yield temp_dir
        # Cleanup
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

    def test_pla_has_carbon_oxygen_backbone(self, temp_output_dir):
        """Test that PLA backbone contains C and O atoms."""
        generator = MonomerGenerator(
            base_name="PLA",
            mechanism='esterification',
            output_dir=temp_output_dir,
            is_gaff=True,
            verbose=False
        )

        mol = Chem.MolFromSmiles("CC(C(=O)O)O")
        mol_with_h = Chem.AddHs(mol)

        # Count atoms
        carbons = [atom for atom in mol_with_h.GetAtoms() if atom.GetSymbol() == 'C']
        oxygens = [atom for atom in mol_with_h.GetAtoms() if atom.GetSymbol() == 'O']

        # Should have multiple carbons and oxygens
        assert len(carbons) >= 3
        assert len(oxygens) >= 2

    def test_pla_connection_atoms(self, temp_output_dir):
        """Test that PLA connection atoms are C and O."""
        from AutoPoly.polymerization_patterns import get_pattern

        pattern = get_pattern('esterification')
        assert pattern.connection_atoms == ['C', 'O']


@pytest.mark.slow
@pytest.mark.integration
class TestPLAAtomTyping:
    """Test atom typing for PLA monomers."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory for tests."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_pla_typing_test_")
        yield temp_dir
        # Cleanup
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

    def test_pla_gaff_atom_typing(self, temp_output_dir):
        """Test that PLA monomers get GAFF atom types."""
        generator = MonomerGenerator(
            base_name="PLA",
            mechanism='esterification',
            output_dir=temp_output_dir,
            is_gaff=True,
            verbose=False
        )

        variants = generator.generate_variants(smiles="CC(C(=O)O)O")

        # Check that atoms have AtomType property
        for variant_name in ['internal', 'left_end', 'right_end']:
            mol = variants[variant_name]
            for atom in mol.GetAtoms():
                try:
                    atom_type = atom.GetProp('AtomType')
                    # GAFF types typically start with @atom:
                    assert atom_type.startswith('@atom:')
                except KeyError:
                    pytest.fail(f"Atom {atom.GetIdx()} in {variant_name} has no AtomType")

    def test_pla_carbon_oxygen_types(self, temp_output_dir):
        """Test that C and O atoms have appropriate types."""
        generator = MonomerGenerator(
            base_name="PLA",
            mechanism='esterification',
            output_dir=temp_output_dir,
            is_gaff=True,
            verbose=False
        )

        variants = generator.generate_variants(smiles="CC(C(=O)O)O")
        mol = variants['internal']

        # Check carbon atoms
        carbons = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']
        for atom in carbons:
            atom_type = atom.GetProp('AtomType')
            assert atom_type is not None

        # Check oxygen atoms
        oxygens = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'O']
        for atom in oxygens:
            atom_type = atom.GetProp('AtomType')
            assert atom_type is not None


@pytest.mark.slow
@pytest.mark.integration
class TestPLAComparisonWithVinyl:
    """Test that PLA generation differs from vinyl polymers."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory for tests."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_pla_comparison_test_")
        yield temp_dir
        # Cleanup
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

    def test_pla_mechanism_differs_from_pe(self, temp_output_dir):
        """Test that PLA uses esterification, not vinyl_addition."""
        pla_mol = Chem.MolFromSmiles("CC(C(=O)O)O")
        pe_mol = Chem.MolFromSmiles("C=C")

        pla_mechanism = detect_mechanism(pla_mol, dop=10, verbose=False)
        pe_mechanism = detect_mechanism(pe_mol, dop=10, verbose=False)

        assert pla_mechanism == 'esterification'
        assert pe_mechanism == 'vinyl_addition'
        assert pla_mechanism != pe_mechanism

    def test_pla_connection_atoms_vs_pe(self, temp_output_dir):
        """Test that PLA connection atoms differ from PE."""
        from AutoPoly.polymerization_patterns import get_pattern

        pla_pattern = get_pattern('esterification')
        pe_pattern = get_pattern('vinyl_addition')

        # PLA: C and O
        assert pla_pattern.connection_atoms == ['C', 'O']

        # PE: C and C
        assert pe_pattern.connection_atoms == ['C', 'C']

        # Should be different
        assert pla_pattern.connection_atoms != pe_pattern.connection_atoms


@pytest.mark.slow
@pytest.mark.integration
class TestPLADOP1Support:
    """Test that PLA respects DOP=1 single molecule generation."""

    def test_pla_dop1_generates_single_molecule(self):
        """Test that lactic acid with DOP=1 generates single molecule."""
        mol = Chem.MolFromSmiles("CC(C(=O)O)O")
        mechanism = detect_mechanism(mol, dop=1, verbose=False)

        # Should return 'none' for DOP=1
        assert mechanism == 'none'

    def test_pla_dop10_generates_polymer(self):
        """Test that lactic acid with DOP>1 generates polymer."""
        mol = Chem.MolFromSmiles("CC(C(=O)O)O")
        mechanism = detect_mechanism(mol, dop=10, verbose=False)

        # Should return 'esterification' for DOP>1
        assert mechanism == 'esterification'
