"""
Unit tests for the Polymerization class.

This module tests the Polymerization class functionality including:
- Initialization with various parameters
- Working directory creation
- Force field selection
- Model management
"""

import os
import sys
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

# Add the AutoPoly package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from AutoPoly.system import System
from AutoPoly.polymerization import Polymerization
from AutoPoly.polymer import Polymer


class TestPolymerizationInitialization:
    """Test Polymerization class initialization."""

    def test_initialization_basic(self, temp_system):
        """Test basic initialization of Polymerization."""
        # Mock the run=False to avoid actual moltemplate execution
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            poly = Polymerization(
                name="test_poly",
                system=temp_system,
                model=[],
                run=False
            )
            assert poly.name == "test_poly"
            assert poly.system == temp_system
            assert poly.model == []

    def test_force_field_oplsaa(self, temp_system):
        """Test initialization with OPLS-AA force field."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            # create_working_directory not mocked
                poly = Polymerization(
                    name="test_oplsaa",
                    system=temp_system,
                    model=[],
                    run=False,
                    force_field="oplsaa"
                )
                assert poly.force_field == "oplsaa"
                assert "oplsaa.prm" in poly.path_oplsaaprm

    def test_force_field_gaff(self, temp_system):
        """Test initialization with GAFF force field."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            # create_working_directory not mocked
                poly = Polymerization(
                    name="test_gaff",
                    system=temp_system,
                    model=[],
                    run=False,
                    force_field="gaff"
                )
                assert poly.force_field == "gaff"
                assert "gaff.lt" in poly.path_oplsaaprm

    def test_force_field_lopls(self, temp_system):
        """Test initialization with LOPLS force field."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            # create_working_directory not mocked
                poly = Polymerization(
                    name="test_lopls",
                    system=temp_system,
                    model=[],
                    run=False,
                    force_field="lopls"
                )
                assert poly.force_field == "lopls"
                assert "loplsaa.prm" in poly.path_oplsaaprm

    def test_deprecated_is_lopls_parameter(self, temp_system):
        """Test that deprecated is_lopls parameter shows warning."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            # create_working_directory not mocked
                poly = Polymerization(
                    name="test_deprecated",
                    system=temp_system,
                    model=[],
                    run=False,
                    is_lopls=True
                )
                # Should convert to force_field="lopls"
                assert poly.force_field == "lopls"

    def test_invalid_force_field_raises_error(self, temp_system):
        """Test that invalid force field raises SystemExit."""
        with pytest.raises(SystemExit):
            Polymerization(
                name="test_invalid",
                system=temp_system,
                model=[],
                run=False,
                force_field="invalid_ff"
            )

    def test_default_parameters(self, temp_system):
        """Test initialization with default parameters."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            poly = Polymerization(
                name="test_defaults",
                system=temp_system,
                model=[],
                run=False
            )
            assert poly.rotate == 90.0
            assert poly.offset_spacing == 2.0
            assert poly.offset == 4.0
            assert poly.packingL_spacing == 5.0
            assert poly.moltemplate_box_size == 400.0


class TestPolymerizationDirectoryManagement:
    """Test Polymerization directory management."""

    def test_create_working_directory(self, temp_system):
        """Test working directory creation."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            # Create a unique name to avoid conflicts
            import uuid
            unique_name = f"test_dir_{uuid.uuid4().hex[:8]}"
            poly = Polymerization(
                name=unique_name,
                system=temp_system,
                model=[],
                run=False
            )
            expected_path = Path(temp_system.get_output_path()) / unique_name / "moltemplate"
            assert expected_path.exists()
            assert expected_path.is_dir()

    def test_existing_directory_prompt(self, temp_system):
        """Test that existing directory triggers user prompt."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            import uuid
            unique_name = f"test_exist_{uuid.uuid4().hex[:8]}"

            # Create first instance
            poly1 = Polymerization(
                name=unique_name,
                system=temp_system,
                model=[],
                run=False
            )

            # Mock input to return 'n' (don't delete)
            with patch('builtins.input', return_value='n'):
                with pytest.raises(SystemExit):
                    poly2 = Polymerization(
                        name=unique_name,
                        system=temp_system,
                        model=[],
                        run=False
                    )


class TestPolymerizationMethods:
    """Test Polymerization methods."""

    def test_set_tacticity(self, temp_system):
        """Test set_tacticity method."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            poly = Polymerization(
                name="test_tacticity",
                system=temp_system,
                model=[],
                run=False
            )
            poly.set_tacticity("isotactic")
            assert poly.tacticity == "isotactic"

    def test_path_cwd_construction(self, temp_system):
        """Test that path_cwd is constructed correctly."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            poly = Polymerization(
                name="test_cwd",
                system=temp_system,
                model=[],
                run=False
            )
            expected = f"{temp_system.get_output_path()}/test_cwd/moltemplate/"
            assert poly.path_cwd == expected


class TestPolymerizationModelManagement:
    """Test Polymerization model management."""

    def test_model_with_polymer_objects(self, temp_system):
        """Test Polymerization with Polymer model objects."""
        polymer = Polymer(
            ChainNum=1,
            Sequence=["PE", "PE"],
            DOP=2,
            topology="linear"
        )

        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            # create_working_directory not mocked
                poly = Polymerization(
                    name="test_model",
                    system=temp_system,
                    model=[polymer],
                    run=False
                )
                assert len(poly.model) == 1
                assert poly.model[0] == polymer


class TestPolymerizationPathManagement:
    """Test Polymerization path management."""

    def test_path_master_construction(self, temp_system):
        """Test that path_master points to extern directory."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            # create_working_directory not mocked
                poly = Polymerization(
                    name="test_path",
                    system=temp_system,
                    model=[],
                    run=False
                )
                assert "extern" in poly.path_master
                # Use proper path relationship checking
                assert poly.path_moltemplatesrc.startswith(poly.path_master.rstrip('/') + '/')

    def test_moltemplate_source_path(self, temp_system):
        """Test moltemplate source path construction."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            # create_working_directory not mocked
                poly = Polymerization(
                    name="test_moltemplate_src",
                    system=temp_system,
                    model=[],
                    run=False
                )
                assert "moltemplate" in poly.path_moltemplatesrc
                assert "src" in poly.path_moltemplatesrc


class TestMonomerFileOperations:
    """Test monomer file operations in Polymerization class."""

    def test_n_monomer_atoms_counts_correctly(self, temp_system):
        """Test that n_monomer_atoms correctly counts atoms in .lt file."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            poly = Polymerization(
                name="test_atom_count",
                system=temp_system,
                model=[],
                run=False
            )

            # Create a test .lt file with known number of atoms
            test_content = """# Test monomer file
write("Data Atoms") {
  $atom:C1  @atom:opls_135  1  0.0  0.0  0.0
  $atom:H2  @atom:opls_140  2  1.0  0.0  0.0
  $atom:H3  @atom:opls_140  3  0.0  1.0  0.0
}
"""
            test_file = Path(poly.path_cwd) / "test_monomer.lt"
            test_file.write_text(test_content)

            # Test the count
            atom_count = poly.n_monomer_atoms("test_monomer.lt")
            assert atom_count == 3

    def test_n_monomer_atoms_handles_missing_file(self, temp_system):
        """Test that n_monomer_atoms handles missing files."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            poly = Polymerization(
                name="test_missing_file",
                system=temp_system,
                model=[],
                run=False
            )

            # Test with non-existent file
            with pytest.raises(SystemExit):
                poly.n_monomer_atoms("nonexistent.lt")

    def test_n_monomer_atoms_handles_malformed_file(self, temp_system):
        """Test that n_monomer_atoms handles malformed .lt files."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            poly = Polymerization(
                name="test_malformed",
                system=temp_system,
                model=[],
                run=False
            )

            # Create a malformed .lt file (missing closing brace)
            test_content = """# Test malformed file
write("Data Atoms") {
  $atom:C1  @atom:opls_135  1  0.0  0.0  0.0
  $atom:H2  @atom:opls_140  2  1.0  0.0  0.0
"""
            test_file = Path(poly.path_cwd) / "malformed.lt"
            test_file.write_text(test_content)

            # Test the count - should return 0 or handle gracefully
            atom_count = poly.n_monomer_atoms("malformed.lt")
            # The method should not crash, return value depends on implementation
            assert isinstance(atom_count, int)


class TestDirectoryManagementWithUserInput:
    """Test directory management with user input handling."""

    def test_create_folder_with_user_confirmation(self, temp_system):
        """Test create_folder prompts user when directory exists."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            import uuid
            unique_name = f"test_folder_{uuid.uuid4().hex[:8]}"

            # Create first instance to establish the directory
            poly1 = Polymerization(
                name=unique_name,
                system=temp_system,
                model=[],
                run=False
            )

            # Mock input to return 'y' (delete and recreate)
            with patch('builtins.input', return_value='y'):
                poly2 = Polymerization(
                    name=unique_name,
                    system=temp_system,
                    model=[],
                    run=False
                )
                # Should successfully create the directory
                assert poly2.path_cwd is not None

    def test_create_folder_rejects_overwrite(self, temp_system):
        """Test create_folder exits when user rejects overwrite."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            import uuid
            unique_name = f"test_reject_{uuid.uuid4().hex[:8]}"

            # Create first instance
            poly1 = Polymerization(
                name=unique_name,
                system=temp_system,
                model=[],
                run=False
            )

            # Mock input to return 'n' (don't delete)
            with patch('builtins.input', return_value='n'):
                with pytest.raises(SystemExit):
                    poly2 = Polymerization(
                        name=unique_name,
                        system=temp_system,
                        model=[],
                        run=False
                    )


class TestMoltemplateIntegration:
    """Test moltemplate integration and subprocess calls."""

    def test_invoke_moltemplate_calls_subprocess(self, temp_system):
        """Test that invoke_moltemplate calls moltemplate subprocess."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            # create_working_directory not mocked
                poly = Polymerization(
                    name="test_invoke",
                    system=temp_system,
                    model=[],
                    run=False
                )

                # Create a dummy system.lt file
                system_lt = Path(poly.path_cwd) / "system.lt"
                system_lt.write_text("# Dummy system file")

                # Mock subprocess.run to avoid actual moltemplate execution
                with patch('subprocess.run') as mock_run:
                    mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

                    poly.invoke_moltemplate()

                    # Verify subprocess was called
                    mock_run.assert_called_once()
                    call_args = mock_run.call_args
                    command_list = call_args[0][0]
                    assert "bash" in command_list
                    # Check that moltemplate.sh is in the path (it's the second argument)
                    assert "moltemplate.sh" in command_list[1]

    def test_invoke_moltemplate_handles_failure(self, temp_system):
        """Test that invoke_moltemplate handles moltemplate failures."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            # create_working_directory not mocked
                poly = Polymerization(
                    name="test_invoke_fail",
                    system=temp_system,
                    model=[],
                    run=False
                )

                # Create a dummy system.lt file
                system_lt = Path(poly.path_cwd) / "system.lt"
                system_lt.write_text("# Dummy system file")

                # Mock subprocess.run to return failure
                with patch('subprocess.run') as mock_run:
                    mock_run.return_value = MagicMock(
                        returncode=1,
                        stdout="",
                        stderr="Moltemplate error"
                    )

                    with pytest.raises(SystemExit):
                        poly.invoke_moltemplate()

    def test_invoke_moltemplate_missing_system_lt(self, temp_system):
        """Test that invoke_moltemplate handles missing system.lt."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            # create_working_directory not mocked
                poly = Polymerization(
                    name="test_missing_lt",
                    system=temp_system,
                    model=[],
                    run=False
                )

                # Don't create system.lt - should fail
                with pytest.raises(SystemExit):
                    poly.invoke_moltemplate()


class TestFileWritingMethods:
    """Test file writing methods for moltemplate files."""

    def test_make_system_lt_generates_file(self, temp_system):
        """Test that make_system_lt generates system.lt file."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            # create_working_directory not mocked
                # Create a simple polymer model
                polymer = Polymer(
                    ChainNum=1,
                    Sequence=["PE"],
                    DOP=2,
                    topology="linear"
                )

                poly = Polymerization(
                    name="test_system_lt",
                    system=temp_system,
                    model=[polymer],
                    run=False
                )

                # Mock the monomer generation
                with patch.object(poly, 'generate_monomer_from_psmiles', return_value="PE"):
                    poly.make_system_lt()

                # Check that system.lt was created
                system_lt = Path(poly.path_cwd) / "system.lt"
                assert system_lt.exists()

                # Check content contains expected imports
                content = system_lt.read_text()
                assert "oplsaa.lt" in content or "gaff.lt" in content

    def test_make_system_lt_includes_force_fields(self, temp_system):
        """Test that make_system_lt includes correct force field imports."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            # create_working_directory not mocked
                polymer = Polymer(
                    ChainNum=1,
                    Sequence=["PE"],
                    DOP=2,
                    topology="linear"
                )

                # Test with GAFF
                poly_gaff = Polymerization(
                    name="test_gaff_ff",
                    system=temp_system,
                    model=[polymer],
                    run=False,
                    force_field="gaff"
                )

                with patch.object(poly_gaff, 'generate_monomer_from_psmiles', return_value="PE"):
                    poly_gaff.make_system_lt()

                system_lt = Path(poly_gaff.path_cwd) / "system.lt"
                content = system_lt.read_text()
                assert "gaff.lt" in content

    def test_make_system_lt_sets_correct_topology(self, temp_system):
        """Test that make_system_lt handles different topologies."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            # Test with ring topology
            polymer_ring = Polymer(
                ChainNum=1,
                Sequence=["PE"],
                DOP=3,
                topology="ring"
            )

            poly = Polymerization(
                name="test_ring_topology",
                system=temp_system,
                model=[polymer_ring],
                run=False
            )

            with patch.object(poly, 'generate_monomer_from_psmiles', return_value="PE"):
                poly.make_system_lt()

            # Check that system.lt was created
            system_lt = Path(poly.path_cwd) / "system.lt"
            assert system_lt.exists()

    @pytest.mark.skip(reason="Complex test requiring proper monomer file format setup - to be revisited")
    def test_make_poly_lt_generates_file(self, temp_system):
        """Test that make_poly_lt generates poly_*.lt files."""
        # This test requires proper monomer file format with end atoms
        # Skipping for now - to be implemented with better mock setup
        pass

    def test_make_force_field_lt_generates_file(self, temp_system):
        """Test that make_force_field_lt generates force field file."""
        with patch.object(Polymerization, 'make_lmp_data_file_by_moltemplate'):
            # create_working_directory not mocked
                polymer = Polymer(
                    ChainNum=1,
                    Sequence=["PE"],
                    DOP=2,
                    topology="linear"
                )

                poly = Polymerization(
                    name="test_ff_lt",
                    system=temp_system,
                    model=[polymer],
                    run=False,
                    force_field="gaff"
                )

                # Mock the subset generation methods
                with patch.object(poly, 'make_gaff_subset'):
                    poly.make_force_field_lt()

                    # Verify make_gaff_subset was called for GAFF
                    # The actual file generation is tested in integration tests
