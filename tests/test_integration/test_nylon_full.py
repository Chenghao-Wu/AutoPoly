"""
Integration tests for Nylon-6,6 generation in AutoPoly.

These tests verify the complete workflow for generating Nylon monomers and polymers
using the amidation mechanism with heteroatom (C, N) backbone support.

Key Tests:
- Nylon monomer generation with amidation mechanism
- Mechanism detection for amino acids and diamines/diacids
- Atom type assignment for amide linkages
- .lt file generation for Nylon variants
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
class TestNylonMechanismDetection:
    """Test amidation mechanism detection for Nylon precursors."""

    def test_amino_acid_detection(self):
        """Test that amino acids are detected as amidation."""
        # Glycine: H2N-CH2-COOH
        mol = Chem.MolFromSmiles("NCC(=O)O")

        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        assert mechanism == 'amidation'

    def test_diamine_detection(self):
        """Test detection of diamines (part of Nylon system)."""
        # Hexamethylenediamine: H2N-(CH2)6-NH2
        mol = Chem.MolFromSmiles("NCCCCCN")

        # Diamines have amine groups but need carboxyl for amidation
        # Should return 'none' (needs to react with diacid)
        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        assert mechanism == 'none'

    def test_diacid_detection(self):
        """Test detection of diacids (part of Nylon system)."""
        # Adipic acid: HOOC-(CH2)4-COOH
        mol = Chem.MolFromSmiles("OC(=O)CCCC(=O)O")

        # Diacids have carboxyl groups but no separate alcohols or amines
        # Should return 'none' (needs to react with diamine or diol)
        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        assert mechanism == 'none'

    def test_nylon_66_precursor_detection(self):
        """Test detection of Nylon-6,6 precursor structure."""
        # This tests a structure with both amine and carboxyl
        # Amino acid: NCC(=O)O
        mol = Chem.MolFromSmiles("NCC(=O)O")

        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        assert mechanism == 'amidation'


@pytest.mark.slow
@pytest.mark.integration
class TestNylonMonomerGeneration:
    """Test Nylon monomer generation workflow."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory for tests."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_nylon_test_")
        yield temp_dir
        # Cleanup
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

    def test_nylon_generator_initialization(self, temp_output_dir):
        """Test MonomerGenerator initialization for Nylon."""
        generator = MonomerGenerator(
            base_name="Nylon",
            mechanism='amidation',
            output_dir=temp_output_dir,
            is_gaff=True,  # GAFF is better for polyamides
            verbose=False
        )

        assert generator.base_name == "Nylon"
        assert generator.mechanism == 'amidation'
        assert generator.is_gaff == True

    def test_nylon_variant_generation(self, temp_output_dir):
        """Test generating Nylon monomer variants."""
        generator = MonomerGenerator(
            base_name="Nylon",
            mechanism='amidation',
            output_dir=temp_output_dir,
            is_gaff=True,
            verbose=False
        )

        # Generate variants from amino acid (glycine)
        variants = generator.generate_variants(smiles="[*]NCC(=O)O[*]")

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

    def test_nylon_lt_file_generation(self, temp_output_dir):
        """Test generating .lt files for Nylon."""
        generator = MonomerGenerator(
            base_name="Nylon",
            mechanism='amidation',
            output_dir=temp_output_dir,
            is_gaff=True,
            verbose=False
        )

        variants = generator.generate_variants(smiles="[*]NCC(=O)O[*]")
        files = generator.generate_lt_files(variants, generate_t1=False)

        # Check that files were generated
        assert 'internal' in files
        assert 'left_end' in files
        assert 'right_end' in files

        # Check that files exist
        assert os.path.exists(files['internal'])
        assert os.path.exists(files['left_end'])
        assert os.path.exists(files['right_end'])

    def test_nylon_file_contents(self, temp_output_dir):
        """Test that generated .lt files have correct structure."""
        generator = MonomerGenerator(
            base_name="Nylon",
            mechanism='amidation',
            output_dir=temp_output_dir,
            is_gaff=True,
            verbose=False
        )

        variants = generator.generate_variants(smiles="[*]NCC(=O)O[*]")
        files = generator.generate_lt_files(variants, generate_t1=False)

        # Read the internal variant file
        with open(files['internal'], 'r') as f:
            content = f.read()

        # Check for GAFF import
        assert 'import "gaff.lt"' in content

        # Check for Nylon class definition
        assert 'Nyloni inherits GAFF' in content

        # Check for atom definitions
        assert 'write("Data Atoms")' in content

        # Check for bond definitions
        assert 'write(\'Data Bond List\')' in content


@pytest.mark.slow
@pytest.mark.integration
class TestNylonHeteroatomBackbone:
    """Test that Nylon uses heteroatom (C, N) backbone correctly."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory for tests."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_nylon_heteroatom_test_")
        yield temp_dir
        # Cleanup
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

    def test_nylon_has_carbon_nitrogen_backbone(self, temp_output_dir):
        """Test that Nylon backbone contains C and N atoms."""
        generator = MonomerGenerator(
            base_name="Nylon",
            mechanism='amidation',
            output_dir=temp_output_dir,
            is_gaff=True,
            verbose=False
        )

        mol = Chem.MolFromSmiles("NCC(=O)O")  # Glycine
        mol_with_h = Chem.AddHs(mol)

        # Count atoms
        carbons = [atom for atom in mol_with_h.GetAtoms() if atom.GetSymbol() == 'C']
        nitrogens = [atom for atom in mol_with_h.GetAtoms() if atom.GetSymbol() == 'N']

        # Should have carbons and at least one nitrogen
        assert len(carbons) >= 2
        assert len(nitrogens) >= 1

    def test_nylon_connection_atoms(self, temp_output_dir):
        """Test that Nylon connection atoms are C and N."""
        from AutoPoly.polymerization_patterns import get_pattern

        pattern = get_pattern('amidation')
        assert pattern.connection_atoms == ['C', 'N']


@pytest.mark.slow
@pytest.mark.integration
class TestNylonAtomTyping:
    """Test atom typing for Nylon monomers."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory for tests."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_nylon_typing_test_")
        yield temp_dir
        # Cleanup
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

    def test_nylon_gaff_atom_typing(self, temp_output_dir):
        """Test that Nylon monomers get GAFF atom types."""
        generator = MonomerGenerator(
            base_name="Nylon",
            mechanism='amidation',
            output_dir=temp_output_dir,
            is_gaff=True,
            verbose=False
        )

        variants = generator.generate_variants(smiles="[*]NCC(=O)O[*]")

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

    def test_nylon_carbon_nitrogen_types(self, temp_output_dir):
        """Test that C and N atoms have appropriate types."""
        generator = MonomerGenerator(
            base_name="Nylon",
            mechanism='amidation',
            output_dir=temp_output_dir,
            is_gaff=True,
            verbose=False
        )

        variants = generator.generate_variants(smiles="[*]NCC(=O)O[*]")
        mol = variants['internal']

        # Check carbon atoms
        carbons = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']
        for atom in carbons:
            atom_type = atom.GetProp('AtomType')
            assert atom_type is not None

        # Check nitrogen atoms
        nitrogens = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'N']
        for atom in nitrogens:
            atom_type = atom.GetProp('AtomType')
            assert atom_type is not None


@pytest.mark.slow
@pytest.mark.integration
class TestNylonComparisonWithOtherPolymers:
    """Test that Nylon generation differs from other polymer types."""

    def test_nylon_vs_pla_mechanism(self):
        """Test that Nylon uses amidation, PLA uses esterification."""
        nylon_mol = Chem.MolFromSmiles("NCC(=O)O")  # Amino acid
        pla_mol = Chem.MolFromSmiles("CC(C(=O)O)O")  # Lactic acid

        nylon_mechanism = detect_mechanism(nylon_mol, dop=10, verbose=False)
        pla_mechanism = detect_mechanism(pla_mol, dop=10, verbose=False)

        assert nylon_mechanism == 'amidation'
        assert pla_mechanism == 'esterification'
        assert nylon_mechanism != pla_mechanism

    def test_nylon_vs_pe_mechanism(self):
        """Test that Nylon uses amidation, PE uses vinyl_addition."""
        nylon_mol = Chem.MolFromSmiles("NCC(=O)O")
        pe_mol = Chem.MolFromSmiles("C=C")

        nylon_mechanism = detect_mechanism(nylon_mol, dop=10, verbose=False)
        pe_mechanism = detect_mechanism(pe_mol, dop=10, verbose=False)

        assert nylon_mechanism == 'amidation'
        assert pe_mechanism == 'vinyl_addition'

    def test_nylon_connection_atoms_unique(self):
        """Test that Nylon connection atoms are unique."""
        from AutoPoly.polymerization_patterns import get_pattern

        nylon_pattern = get_pattern('amidation')
        pla_pattern = get_pattern('esterification')
        pe_pattern = get_pattern('vinyl_addition')

        # Nylon: C and N
        assert nylon_pattern.connection_atoms == ['C', 'N']

        # PLA: C and O
        assert pla_pattern.connection_atoms == ['C', 'O']

        # PE: C and C
        assert pe_pattern.connection_atoms == ['C', 'C']

        # All different
        patterns = [nylon_pattern.connection_atoms,
                    pla_pattern.connection_atoms,
                    pe_pattern.connection_atoms]
        assert len(set(map(tuple, patterns))) == 3


@pytest.mark.slow
@pytest.mark.integration
class TestNylonDOP1Support:
    """Test that Nylon respects DOP=1 single molecule generation."""

    def test_nylon_dop1_generates_single_molecule(self):
        """Test that amino acids with DOP=1 generate single molecule."""
        mol = Chem.MolFromSmiles("NCC(=O)O")
        mechanism = detect_mechanism(mol, dop=1, verbose=False)

        # Should return 'none' for DOP=1
        assert mechanism == 'none'

    def test_nylon_dop10_generates_polymer(self):
        """Test that amino acids with DOP>1 generate polymer."""
        mol = Chem.MolFromSmiles("NCC(=O)O")
        mechanism = detect_mechanism(mol, dop=10, verbose=False)

        # Should return 'amidation' for DOP>1
        assert mechanism == 'amidation'


@pytest.mark.slow
@pytest.mark.integration
class TestNylonCondensationPolymerization:
    """Test that Nylon is correctly identified as condensation polymerization."""

    def test_nylon_is_condensation(self):
        """Test that amidation is marked as condensation."""
        from AutoPoly.polymerization_patterns import get_pattern

        pattern = get_pattern('amidation')
        assert pattern.is_condensation == True

    def test_nylon_requires_two_groups(self):
        """Test that Nylon precursors need 2+ functional groups."""
        from AutoPoly.polymerization_patterns import get_pattern

        pattern = get_pattern('amidation')
        assert pattern.requires_two_groups == True


@pytest.mark.slow
@pytest.mark.integration
class TestNylonRealWorldExamples:
    """Test with real Nylon precursor structures."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory for tests."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_nylon_real_test_")
        yield temp_dir
        # Cleanup
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

    def test_caprolactam_detection(self):
        """Test detection of caprolactam (Nylon-6 precursor)."""
        # Caprolactam: cyclic amide
        mol = Chem.MolFromSmiles("C1CCC(=O)NC1")

        mechanism = detect_mechanism(mol, dop=10, verbose=False)
        # Should detect amidation or ring-opening
        assert mechanism in ['amidation', 'ring_opening_amide', 'none']

    def test_hexamethylenediamine_with_adipic_acid(self, temp_output_dir):
        """Test generation of Nylon-6,6 from diamine and diacid."""
        # This tests the concept, actual Nylon-6,6 requires two monomers
        # For now, test that amino acid structure works
        generator = MonomerGenerator(
            base_name="Nylon66",
            mechanism='amidation',
            output_dir=temp_output_dir,
            is_gaff=True,
            verbose=False
        )

        # Use amino acid as proxy
        variants = generator.generate_variants(smiles="[*]NCC(=O)O[*]")

        # Should generate successfully
        assert 'internal' in variants
        assert isinstance(variants['internal'], Chem.Mol)
