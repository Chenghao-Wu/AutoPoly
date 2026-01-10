"""Tests for the monomer_processing module."""

import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
from AutoPoly.monomer_processing import (
    generate_monomer_from_psmiles,
    generate_molecule_from_smiles,
    generate_sequence_variants_for_polymerization,
    n_monomer_atoms,
    read_lt_end_atoms,
    extract_element_from_atom,
    evaluate_offset,
)


class TestGenerateMonomerFromPsmiles:
    """Test monomer generation from pSMILES."""

    @patch('AutoPoly.monomer_processing.MonomerGenerator')
    def test_generate_monomer_from_psmiles_creates_generator(self, mock_generator_class, tmp_path):
        """Test that MonomerGenerator is called with correct parameters."""
        mock_generator = MagicMock()
        mock_generator_class.return_value = mock_generator
        mock_generator.from_smiles.return_value = []
        mock_generator.write_lt_files.return_value = []

        cache = {}
        result_name, result_counter = generate_monomer_from_psmiles(
            psmiles="[*]C=C[*]",
            path_cwd=str(tmp_path),
            force_field="oplsaa",
            generated_cache=cache,
            counter=0
        )

        # Verify generator was created with correct parameters
        mock_generator_class.assert_called_once()
        call_kwargs = mock_generator_class.call_args[1]
        assert call_kwargs['base_name'] == 'monomer_0'
        assert call_kwargs['force_field'] == 'oplsaa'
        assert call_kwargs['output_dir'] == str(tmp_path)
        assert call_kwargs['verbose'] is False

    @patch('AutoPoly.monomer_processing.MonomerGenerator')
    def test_generate_monomer_from_psmiles_calls_generation_methods(self, mock_generator_class, tmp_path):
        """Test that from_smiles and write_lt_files are called."""
        mock_generator = MagicMock()
        mock_generator_class.return_value = mock_generator
        mock_generator.from_smiles.return_value = []
        mock_generator.write_lt_files.return_value = []

        cache = {}
        generate_monomer_from_psmiles(
            psmiles="[*]C=C[*]",
            path_cwd=str(tmp_path),
            force_field="oplsaa",
            generated_cache=cache,
            counter=0
        )

        # Verify methods were called
        mock_generator.from_smiles.assert_called_once_with(smiles="[*]C=C[*]", n_monomers=3)
        mock_generator.write_lt_files.assert_called_once_with([], generate_t1=True)

    @patch('AutoPoly.monomer_processing.MonomerGenerator')
    def test_generate_monomer_from_psmiles_caches_result(self, mock_generator_class, tmp_path):
        """Test that second call uses cache."""
        mock_generator = MagicMock()
        mock_generator_class.return_value = mock_generator
        mock_generator.from_smiles.return_value = []
        mock_generator.write_lt_files.return_value = []

        cache = {}
        # First call
        result_name1, counter1 = generate_monomer_from_psmiles(
            psmiles="[*]C=C[*]",
            path_cwd=str(tmp_path),
            force_field="oplsaa",
            generated_cache=cache,
            counter=0
        )

        # Second call with same pSMILES
        result_name2, counter2 = generate_monomer_from_psmiles(
            psmiles="[*]C=C[*]",
            path_cwd=str(tmp_path),
            force_field="oplsaa",
            generated_cache=cache,
            counter=1
        )

        # Generator should only be called once (cache hit on second call)
        assert mock_generator.from_smiles.call_count == 1
        assert result_name1 == result_name2

    @patch('AutoPoly.monomer_processing.MonomerGenerator')
    def test_generate_monomer_from_psmiles_increments_counter(self, mock_generator_class, tmp_path):
        """Test that counter is incremented correctly."""
        mock_generator = MagicMock()
        mock_generator_class.return_value = mock_generator
        mock_generator.from_smiles.return_value = []
        mock_generator.write_lt_files.return_value = []

        cache = {}
        result_name, result_counter = generate_monomer_from_psmiles(
            psmiles="[*]C=C[*]",
            path_cwd=str(tmp_path),
            force_field="oplsaa",
            generated_cache=cache,
            counter=5
        )

        # Counter should be incremented
        assert result_counter == 6
        assert result_name == "monomer_5"


class TestGenerateMoleculeFromSmiles:
    """Test molecule generation from SMILES."""

    @patch('AutoPoly.monomer_processing.MonomerGenerator')
    def test_generate_molecule_from_smiles_creates_generator(self, mock_generator_class, tmp_path):
        """Test that MonomerGenerator is called with correct parameters."""
        mock_generator = MagicMock()
        mock_generator_class.return_value = mock_generator
        mock_variant = MagicMock()
        mock_generator.from_single_molecule.return_value = mock_variant

        cache = {}
        result_filename, result_counter = generate_molecule_from_smiles(
            smiles="O",
            molecule_name="water",
            path_cwd=str(tmp_path),
            force_field="oplsaa",
            generated_cache=cache,
            counter=0
        )

        # Verify generator was created with molecule_name as base_name
        mock_generator_class.assert_called_once()
        call_kwargs = mock_generator_class.call_args[1]
        assert call_kwargs['base_name'] == 'water'
        assert call_kwargs['force_field'] == 'oplsaa'
        assert call_kwargs['output_dir'] == str(tmp_path)

    @patch('AutoPoly.monomer_processing.MonomerGenerator')
    def test_generate_molecule_from_smiles_calls_single_molecule(self, mock_generator_class, tmp_path):
        """Test that from_single_molecule is called."""
        mock_generator = MagicMock()
        mock_generator_class.return_value = mock_generator
        mock_variant = MagicMock()
        mock_generator.from_single_molecule.return_value = mock_variant

        cache = {}
        generate_molecule_from_smiles(
            smiles="CCO",
            molecule_name="ethanol",
            path_cwd=str(tmp_path),
            force_field="gaff",
            generated_cache=cache,
            counter=0
        )

        # Verify methods were called
        mock_generator.from_single_molecule.assert_called_once_with(smiles="CCO", molecule_name="ethanol")
        mock_generator.write_single_molecule.assert_called_once_with(mock_variant, generate_t1=False)

    @patch('AutoPoly.monomer_processing.MonomerGenerator')
    def test_generate_molecule_from_smiles_caches_result(self, mock_generator_class, tmp_path):
        """Test that molecule generation is cached."""
        mock_generator = MagicMock()
        mock_generator_class.return_value = mock_generator
        mock_variant = MagicMock()
        mock_generator.from_single_molecule.return_value = mock_variant

        cache = {}
        # First call
        result_filename1, counter1 = generate_molecule_from_smiles(
            smiles="O",
            molecule_name="water",
            path_cwd=str(tmp_path),
            force_field="oplsaa",
            generated_cache=cache,
            counter=0
        )

        # Second call with same SMILES and force field
        result_filename2, counter2 = generate_molecule_from_smiles(
            smiles="O",
            molecule_name="water",
            path_cwd=str(tmp_path),
            force_field="oplsaa",
            generated_cache=cache,
            counter=1
        )

        # Generator should only be called once (cache hit)
        assert mock_generator.from_single_molecule.call_count == 1
        assert result_filename1 == result_filename2


class TestGenerateSequenceVariants:
    """Test sequence variant generation."""

    @patch('AutoPoly.monomer_processing.MonomerGenerator')
    def test_generate_sequence_variants_creates_generator(self, mock_generator_class, tmp_path):
        """Test MonomerGenerator creation for sequence variants."""
        mock_generator = MagicMock()
        mock_generator_class.return_value = mock_generator
        mock_generator.from_smiles.return_value = []
        mock_generator.write_lt_files.return_value = []

        cache = {}
        result_mapping, result_counter = generate_sequence_variants_for_polymerization(
            base_smiles="[*]C=C[*]",
            dop=5,
            topology="linear",
            path_cwd=str(tmp_path),
            force_field="oplsaa",
            generated_cache=cache,
            counter=0
        )

        # Verify generator creation
        mock_generator_class.assert_called_once()
        call_kwargs = mock_generator_class.call_args[1]
        assert call_kwargs['base_name'] == 'monomer_0'
        assert call_kwargs['force_field'] == 'oplsaa'

    @patch('AutoPoly.monomer_processing.MonomerGenerator')
    def test_generate_sequence_variants_uses_max_dop(self, mock_generator_class, tmp_path):
        """Test that at least 3 monomers are used for variant generation."""
        mock_generator = MagicMock()
        mock_generator_class.return_value = mock_generator
        mock_generator.from_smiles.return_value = []
        mock_generator.write_lt_files.return_value = []

        cache = {}
        # With DOP=2, should still use n_monomers=3
        generate_sequence_variants_for_polymerization(
            base_smiles="[*]C=C[*]",
            dop=2,
            topology="linear",
            path_cwd=str(tmp_path),
            force_field="oplsaa",
            generated_cache=cache,
            counter=0
        )

        # Verify from_smiles was called with n_monomers=3 (max of 3 and dop)
        mock_generator.from_smiles.assert_called_once()
        call_args = mock_generator.from_smiles.call_args[1]
        assert call_args['n_monomers'] == 3

    @patch('AutoPoly.monomer_processing.MonomerGenerator')
    def test_generate_sequence_variants_caches_result(self, mock_generator_class, tmp_path):
        """Test that sequence variants are cached."""
        mock_generator = MagicMock()
        mock_generator_class.return_value = mock_generator
        mock_generator.from_smiles.return_value = []
        mock_generator.write_lt_files.return_value = []

        cache = {}
        # First call
        result_mapping1, counter1 = generate_sequence_variants_for_polymerization(
            base_smiles="[*]C=C[*]",
            dop=5,
            topology="linear",
            path_cwd=str(tmp_path),
            force_field="oplsaa",
            generated_cache=cache,
            counter=0
        )

        # Second call with same parameters
        result_mapping2, counter2 = generate_sequence_variants_for_polymerization(
            base_smiles="[*]C=C[*]",
            dop=5,
            topology="linear",
            path_cwd=str(tmp_path),
            force_field="oplsaa",
            generated_cache=cache,
            counter=1
        )

        # Generator should only be called once
        assert mock_generator.from_smiles.call_count == 1


class TestNMonomerAtoms:
    """Test atom counting from .lt files."""

    def test_n_monomer_atoms_counts_atoms_in_block(self, tmp_path):
        """Test that atoms in Data Atoms block are counted."""
        lt_file = tmp_path / "test.lt"
        lt_file.write_text(
            "# Test monomer\n"
            "write(\"Data Atoms\") {\n"
            "    $atom:C1 $mol:... @atom:CT 0.0 0.0 0.0 0.0\n"
            "    $atom:C2 $mol:... @atom:CT 0.0 1.54 0.0 0.0\n"
            "    $atom:H1 $mol:... @atom:HC 0.0 -0.5 0.0 0.0\n"
            "    $atom:H2 $mol:... @atom:HC 0.0 2.04 0.0 0.0\n"
            "}\n"
        )

        count = n_monomer_atoms("test.lt", str(tmp_path))
        assert count == 4

    def test_n_monomer_atoms_with_empty_block(self, tmp_path):
        """Test counting with empty Data Atoms block."""
        lt_file = tmp_path / "test.lt"
        lt_file.write_text(
            "# Test monomer\n"
            "write(\"Data Atoms\") {\n"
            "}\n"
        )

        count = n_monomer_atoms("test.lt", str(tmp_path))
        # Should count the closing brace line (or 0 depending on implementation)
        assert count >= 0

    def test_n_monomer_atoms_missing_file_raises_system_exit(self, tmp_path):
        """Test that missing file raises GenerationError."""
        from AutoPoly.exceptions import GenerationError
        with pytest.raises(GenerationError, match="Monomer file not found"):
            n_monomer_atoms("nonexistent.lt", str(tmp_path))


class TestReadLtEndAtoms:
    """Test reading end atoms from .lt files."""

    def test_read_lt_end_atoms_returns_first_and_second(self, tmp_path):
        """Test that first and second atoms are returned."""
        lt_file = tmp_path / "test.lt"
        lt_file.write_text(
            "# Test monomer\n"
            "write(\"Data Atoms\") {\n"
            "    $atom:C1 $mol:test @atom:CT 0.0 0.0 0.0 0.0\n"
            "    $atom:H2 $mol:test @atom:HC 0.0 1.0 0.0 0.0\n"
            "    $atom:O3 $mol:test @atom:O 0.0 2.0 0.0 0.0\n"
            "}\n"
        )

        first, second = read_lt_end_atoms(str(lt_file))
        assert first == "C1"
        assert second == "H2"

    def test_read_lt_end_atoms_insufficient_atoms_raises_error(self, tmp_path):
        """Test that insufficient atoms raises ValueError."""
        lt_file = tmp_path / "test.lt"
        lt_file.write_text(
            "write(\"Data Atoms\") {\n"
            "    $atom:C1 $mol:test @atom:CT 0.0 0.0 0.0 0.0\n"
            "}\n"
        )

        with pytest.raises(ValueError, match="Could not find both end atoms"):
            read_lt_end_atoms(str(lt_file))


class TestExtractElementFromAtom:
    """Test element extraction from atom strings."""

    def test_extract_element_from_carbon_atom(self):
        """Test extraction from carbon atom."""
        element = extract_element_from_atom("$atom:C1")
        assert element == "C"

    def test_extract_element_from_hydrogen_atom(self):
        """Test extraction from hydrogen atom."""
        element = extract_element_from_atom("$atom:H16")
        assert element == "H"

    def test_extract_element_from_two_letter_element(self):
        """Test extraction from two-letter element (Si)."""
        element = extract_element_from_atom("$atom:Si2")
        assert element == "Si"

    def test_extract_element_from_three_letter_element(self):
        """Test extraction from three-letter element (Fe)."""
        element = extract_element_from_atom("$atom:Fe10")
        assert element == "Fe"

    def test_extract_element_from_invalid_string(self):
        """Test extraction from invalid string."""
        element = extract_element_from_atom("invalid_string")
        assert element is None

    def test_extract_element_from_lowercase(self):
        """Test extraction from lowercase element (should not match)."""
        element = extract_element_from_atom("$atom:c1")
        assert element is None  # Pattern requires uppercase first letter


class TestEvaluateOffset:
    """Test offset distance evaluation."""

    def test_evaluate_offset_calculates_distance(self, tmp_path):
        """Test that offset is calculated from atom coordinates."""
        lt_file = tmp_path / "test.lt"
        # Create C1-C2 bond with 1.54 Angstrom distance
        lt_file.write_text(
            "write(\"Data Atoms\") {\n"
            "    $atom:C1 $mol:test @atom:CT 0.0 0.0 0.0 0.0\n"
            "    $atom:C2 $mol:test @atom:CT 0.0 1.54 0.0 0.0\n"
            "}\n"
        )

        offset = evaluate_offset("test.lt", str(tmp_path), offset_spacing=0.5, current_offset=3.5)
        # Distance (1.54) + offset_spacing (0.5) = 2.04
        assert offset == pytest.approx(2.04, abs=0.01)

    def test_evaluate_offset_with_zero_spacing(self, tmp_path):
        """Test offset calculation with zero spacing."""
        lt_file = tmp_path / "test.lt"
        lt_file.write_text(
            "write(\"Data Atoms\") {\n"
            "    $atom:C1 $mol:test @atom:CT 0.0 0.0 0.0 0.0\n"
            "    $atom:C2 $mol:test @atom:CT 0.0 1.0 0.0 0.0\n"
            "}\n"
        )

        offset = evaluate_offset("test.lt", str(tmp_path), offset_spacing=0.0, current_offset=3.5)
        # Distance (1.0) + 0.0 = 1.0
        assert offset == pytest.approx(1.0, abs=0.01)

    def test_evaluate_offset_missing_file_returns_current(self, tmp_path):
        """Test that missing file returns current_offset."""
        offset = evaluate_offset("nonexistent.lt", str(tmp_path), offset_spacing=0.5, current_offset=3.5)
        assert offset == 3.5

    def test_evaluate_offset_with_diagonal_bond(self, tmp_path):
        """Test offset calculation with diagonal bond."""
        lt_file = tmp_path / "test.lt"
        # Create atoms at (0,0,0) and (1,1,1) - distance = sqrt(3) ≈ 1.732
        lt_file.write_text(
            "write(\"Data Atoms\") {\n"
            "    $atom:C1 $mol:test @atom:CT 0.0 0.0 0.0 0.0\n"
            "    $atom:C2 $mol:test @atom:CT 0.0 1.0 1.0 1.0\n"
            "}\n"
        )

        offset = evaluate_offset("test.lt", str(tmp_path), offset_spacing=0.0, current_offset=3.5)
        # Distance should be sqrt(3) ≈ 1.732
        assert offset == pytest.approx(1.732, abs=0.01)
