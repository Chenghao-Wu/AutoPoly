"""Tests for file_management module."""

import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
from AutoPoly.core.file_management import create_working_directory
from AutoPoly.core.file_management import get_rid_of_lj_cut_coul_long
from AutoPoly.core.file_management import mv_files
from AutoPoly.core.exceptions import WorkflowError


class TestCreateWorkingDirectory:
    """Test working directory creation."""

    def test_create_working_directory_creates_structure(self, tmp_path):
        """Test that create_working_directory creates correct directory structure."""
        # Create a mock system object
        mock_system = MagicMock()
        mock_system.get_folder_path.return_value = str(tmp_path / "test_project")

        result_path = create_working_directory(mock_system, "simulation")

        # Should create project_dir/moltemplate/
        assert result_path.exists()
        assert result_path.is_dir()
        assert result_path.name == "moltemplate"

        # Verify parent directory structure
        polymer_path = result_path.parent
        assert polymer_path.name == "simulation"

    def test_create_working_directory_with_existing_dir_prompts_user(self, monkeypatch, tmp_path):
        """Test that existing directory prompts user for confirmation."""
        # Create existing directory
        existing_dir = tmp_path / "test_project"
        existing_dir.mkdir(parents=True)

        mock_system = MagicMock()
        mock_system.get_folder_path.return_value = str(tmp_path)

        # Mock user input to return 'y' (yes)
        monkeypatch.setattr('builtins.input', lambda x: 'y')

        result_path = create_working_directory(mock_system, "test_project")

        # Should create moltemplate directory
        assert result_path.exists()
        assert result_path.is_dir()

    def test_create_working_directory_user_declines_exits(self, monkeypatch, tmp_path):
        """Test that declining to overwrite exits with WorkflowError."""
        # Create existing directory
        existing_dir = tmp_path / "test_project"
        existing_dir.mkdir(parents=True)

        mock_system = MagicMock()
        mock_system.get_folder_path.return_value = str(tmp_path)

        # Mock user input to return 'n' (no)
        monkeypatch.setattr('builtins.input', lambda x: 'n')

        with pytest.raises(WorkflowError, match="Directory exists"):
            create_working_directory(mock_system, "test_project")


class TestGetRidOfLjCutCoulLong:
    """Test removal of LJ/cut/coul/long from settings file."""

    def test_remove_lj_cut_coul_long_from_settings(self, tmp_path):
        """Test that LJ/cut/coul/long is removed from pair_coeff lines."""
        settings_file = tmp_path / "system.in.settings"
        settings_file.write_text(
            "pair_style lj/cut/coul/long 10.0\n"
            "pair_coeff 1 1 lj/cut/coul/long 0.0 1.0 1.0\n"
            "pair_coeff 2 2 lj/cut/coul/long 0.0 1.0 1.0\n"
        )

        get_rid_of_lj_cut_coul_long(str(tmp_path))

        content = settings_file.read_text()
        # pair_style line should remain unchanged
        assert "pair_style lj/cut/coul/long 10.0" in content
        # pair_coeff lines should have lj/cut/coul/long removed
        assert "pair_coeff 1 1 0.0 1.0 1.0" in content
        assert "pair_coeff 2 2 0.0 1.0 1.0" in content
        # No lj/cut/coul/long should remain in pair_coeff lines
        lines = content.split('\n')
        for line in lines:
            if line.strip().startswith('pair_coeff'):
                assert 'lj/cut/coul/long' not in line

    def test_remove_lj_cut_coul_long_preserves_other_settings(self, tmp_path):
        """Test that other pair_style settings are preserved."""
        settings_file = tmp_path / "system.in.settings"
        settings_file.write_text(
            "pair_style lj/cut/coul/long 10.0\n"
            "pair_coeff * * lj/cut/coul/long 0.0 1.0 1.0\n"
            "bond_style harmonic\n"
            "angle_style harmonic\n"
        )

        get_rid_of_lj_cut_coul_long(str(tmp_path))

        content = settings_file.read_text()
        # Other settings should be preserved
        assert "bond_style harmonic" in content
        assert "angle_style harmonic" in content

    def test_settings_file_not_found_raises_system_exit(self, tmp_path):
        """Test that missing settings file raises WorkflowError."""
        # Don't create settings file
        with pytest.raises(WorkflowError, match="system.in.settings does not exist"):
            get_rid_of_lj_cut_coul_long(str(tmp_path))

    def test_remove_multiple_lj_cut_coul_long_occurrences(self, tmp_path):
        """Test that multiple occurrences in pair_coeff lines are all removed."""
        settings_file = tmp_path / "system.in.settings"
        settings_file.write_text(
            "pair_style lj/cut/coul/long 10.0\n"
            "pair_coeff 1 1 lj/cut/coul/long 0.0 1.0 1.0\n"
            "pair_coeff 1 2 lj/cut/coul/long 0.0 1.0 1.0\n"
            "pair_coeff 2 2 lj/cut/coul/long 0.0 1.0 1.0\n"
        )

        get_rid_of_lj_cut_coul_long(str(tmp_path))

        content = settings_file.read_text()
        # pair_style should remain
        assert "pair_style lj/cut/coul/long 10.0" in content
        # All pair_coeff occurrences should be removed
        lines = content.split('\n')
        for line in lines:
            if line.strip().startswith('pair_coeff'):
                assert 'lj/cut/coul/long' not in line


class TestMoveFiles:
    """Test file movement and organization."""

    def test_move_files_creates_output_and_input_dirs(self, tmp_path):
        """Test that mv_files creates output/ and input/ directories."""
        # Create working directory structure
        work_dir = tmp_path / "simulation" / "moltemplate"
        work_dir.mkdir(parents=True)

        # Create system.lt (minimal requirement)
        (work_dir / "system.lt").write_text("# Test system file\n")

        mv_files(str(work_dir))

        # Should create output/ and input/ subdirectories in the working dir
        output_dir = work_dir / "output"
        input_dir = work_dir / "input"
        assert output_dir.exists()
        assert input_dir.exists()
        assert output_dir.is_dir()
        assert input_dir.is_dir()

    def test_move_files_copies_system_files_to_parent(self, tmp_path):
        """Test that system.data and system.in are copied to parent."""
        work_dir = tmp_path / "simulation" / "moltemplate"
        work_dir.mkdir(parents=True)

        # Create system.lt and output files
        (work_dir / "system.lt").write_text("# Test\n")
        (work_dir / "system.data").write_text("Test data\n")
        (work_dir / "system.in").write_text("Test input\n")
        (work_dir / "system.in.settings").write_text("Test settings\n")

        mv_files(str(work_dir))

        # Verify files copied to parent (also check they still exist in moltemplate due to copy2)
        parent_dir = tmp_path / "simulation"
        assert (parent_dir / "system.data").exists()
        assert (parent_dir / "system.in").exists()
        assert (parent_dir / "system.in.settings").exists()

    def test_move_files_moves_in_files_to_output(self, tmp_path):
        """Test that system.in* files are moved to output/."""
        work_dir = tmp_path / "simulation" / "moltemplate"
        work_dir.mkdir(parents=True)

        (work_dir / "system.lt").write_text("# Test\n")
        (work_dir / "system.in").write_text("Test input\n")
        (work_dir / "system.init").write_text("Test init\n")
        (work_dir / "output_ttree").mkdir()

        mv_files(str(work_dir))

        # Verify .in* files in output/
        final_output = work_dir / "output"
        assert (final_output / "system.in").exists()
        assert (final_output / "system.init").exists()
        assert (final_output / "output_ttree").exists()

    def test_move_files_moves_lt_files_to_input(self, tmp_path):
        """Test that .lt files are moved to input/."""
        work_dir = tmp_path / "simulation" / "moltemplate"
        work_dir.mkdir(parents=True)

        (work_dir / "system.lt").write_text("# Test\n")
        (work_dir / "PE.lt").write_text("PE monomer\n")
        (work_dir / "PS.lt").write_text("PS monomer\n")
        (work_dir / "oplsaa.prm").write_text("Force field\n")

        mv_files(str(work_dir))

        # Verify .lt and .prm files in input/
        input_dir = work_dir / "input"
        assert (input_dir / "PE.lt").exists()
        assert (input_dir / "PS.lt").exists()
        assert (input_dir / "oplsaa.prm").exists()

        # Verify they're no longer in moltemplate/
        assert not (work_dir / "PE.lt").exists()
        assert not (work_dir / "PS.lt").exists()

    def test_move_files_handles_missing_files_gracefully(self, tmp_path):
        """Test that mv_files doesn't fail if some files are missing."""
        work_dir = tmp_path / "simulation" / "moltemplate"
        work_dir.mkdir(parents=True)

        # Only create system.lt, no output files
        (work_dir / "system.lt").write_text("# Test\n")

        # Should not raise error
        mv_files(str(work_dir))

        # Verify directories are still created
        output_dir = work_dir / "output"
        input_dir = work_dir / "input"
        assert output_dir.exists()
        assert input_dir.exists()

    def test_move_files_moves_data_files(self, tmp_path):
        """Test that system.data file is moved to output/."""
        work_dir = tmp_path / "simulation" / "moltemplate"
        work_dir.mkdir(parents=True)

        (work_dir / "system.lt").write_text("# Test\n")
        # Create a data file that matches system*data pattern (must end with "data")
        (work_dir / "system.data").write_text("Test data\n")

        mv_files(str(work_dir))

        # Verify system.data is copied to parent and also moved to output
        output_dir = work_dir / "output"
        parent_dir = tmp_path / "simulation"
        assert (output_dir / "system.data").exists()
        assert (parent_dir / "system.data").exists()  # Copied first, then moved
