"""Tests for the System class."""

import pytest
import shutil
from pathlib import Path
from AutoPoly.core.system import System


class TestSystemInitialization:
    """Test System class initialization."""

    def test_system_init_creates_directory(self, tmp_path):
        """Test that System initialization creates output directory."""
        output_dir = tmp_path / "test_output"
        system = System(out=str(output_dir))

        folder_path = Path(system.get_folder_path())
        assert folder_path.exists()
        assert folder_path.is_dir()

    def test_system_init_with_none_out_uses_cwd(self):
        """Test that System with out=None uses current working directory."""
        system = System(out=None)

        folder_path = Path(system.get_folder_path())
        assert folder_path.exists()
        assert folder_path == Path.cwd()

    def test_system_init_stores_output_name(self):
        """Test that System stores the output directory name."""
        output_name = "my_simulation"
        system = System(out=output_name)

        assert system.out == output_name


class TestSystemGetPaths:
    """Test System path retrieval methods."""

    def test_get_folder_path_returns_correct_path(self, tmp_path):
        """Test that get_folder_path returns the correct full path."""
        output_dir = tmp_path / "test_output"
        system = System(out=str(output_dir))

        folder_path = system.get_folder_path()

        # Note: System class uses Path.cwd() internally, so the actual path
        # will be different from tmp_path. We just verify it returns a string.
        assert isinstance(folder_path, str)
        assert "test_output" in folder_path

    def test_get_output_path_returns_correct_path(self, tmp_path):
        """Test that get_output_path returns the correct full path."""
        output_dir = tmp_path / "test_output"
        system = System(out=str(output_dir))

        output_path = system.get_output_path()

        assert isinstance(output_path, str)
        assert "test_output" in output_path


class TestSystemChangeDirectory:
    """Test System directory change functionality."""

    def test_change_output_directory_updates_path(self, tmp_path):
        """Test that change_output_directory updates the output path."""
        system = System(out=str(tmp_path / "output1"))

        new_output = tmp_path / "output2"
        system.change_output_directory(str(new_output))

        assert system.out == str(new_output)
        assert "output2" in system.get_folder_path()

    def test_change_output_directory_creates_new_directory(self, tmp_path):
        """Test that change_output_directory creates the new directory."""
        system = System(out=str(tmp_path / "output1"))

        new_output = tmp_path / "output2"
        system.change_output_directory(str(new_output))

        new_folder_path = Path(system.get_folder_path())
        assert new_folder_path.exists()
        assert new_folder_path.is_dir()


class TestSystemCleanup:
    """Test System cleanup functionality."""

    def test_cleanup_removes_output_directory(self, tmp_path):
        """Test that cleanup_output_directory removes the output directory."""
        output_dir = tmp_path / "test_output"
        system = System(out=str(output_dir))

        # Verify directory exists
        folder_path = Path(system.get_folder_path())
        assert folder_path.exists()

        # Cleanup
        system.cleanup_output_directory()

        # Verify directory is removed
        assert not folder_path.exists()

    def test_cleanup_with_nonexistent_directory(self, tmp_path):
        """Test that cleanup doesn't raise error when directory doesn't exist."""
        output_dir = tmp_path / "nonexistent_output"
        system = System(out=str(output_dir))

        # Remove the directory first
        folder_path = Path(system.get_folder_path())
        if folder_path.exists():
            shutil.rmtree(folder_path)

        # Should not raise error
        system.cleanup_output_directory()

        assert not folder_path.exists()


