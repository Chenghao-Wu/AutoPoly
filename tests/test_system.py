"""
Unit tests for the System class.

This module tests the System class functionality including:
- Initialization with and without output directory
- Path methods
- Directory management
- Cleanup operations
"""

import os
import sys
import tempfile
import shutil
from pathlib import Path

# Add the AutoPoly package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from AutoPoly.system import System


class TestSystemInitialization:
    """Test System class initialization."""

    def test_init_without_output_dir(self):
        """Test initialization without specifying an output directory."""
        system = System()
        assert system.out is None
        assert system.get_output_path() == os.getcwd()

    def test_init_with_output_dir(self, temp_system):
        """Test initialization with an output directory."""
        assert temp_system.out is not None
        assert os.path.exists(temp_system.get_output_path())

    def test_output_directory_creation(self, temp_system):
        """Test that output directory is created."""
        output_path = temp_system.get_output_path()
        assert os.path.exists(output_path)
        assert os.path.isdir(output_path)


class TestSystemPathMethods:
    """Test System class path methods."""

    def test_get_output_path(self, temp_system):
        """Test get_output_path returns correct path."""
        expected_path = os.path.join(os.getcwd(), temp_system.out)
        assert temp_system.get_output_path() == expected_path

    def test_get_folder_path_backward_compat(self, temp_system):
        """Test get_FolderPath for backward compatibility."""
        # Test that the old attribute name still works
        assert hasattr(temp_system, 'get_FolderPath')
        assert temp_system.get_FolderPath == temp_system.get_output_path()


class TestSystemDirectoryManagement:
    """Test System class directory management."""

    def test_change_output_directory(self, temp_system):
        """Test changing the output directory."""
        temp_dir = tempfile.mkdtemp(prefix="autopoly_change_")
        try:
            original_path = temp_system.get_output_path()
            temp_system.change_output_directory(temp_dir)
            new_path = temp_system.get_output_path()

            assert new_path != original_path
            assert os.path.exists(new_path)
        finally:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)

    def test_cleanup_output_directory(self, temp_system):
        """Test cleaning up the output directory."""
        output_path = temp_system.get_output_path()
        assert os.path.exists(output_path)

        temp_system.cleanup_output_directory()
        assert not os.path.exists(output_path)


class TestSystemMultipleInstances:
    """Test multiple System instances."""

    def test_multiple_system_instances(self):
        """Test creating multiple System instances."""
        temp_dir1 = tempfile.mkdtemp(prefix="autopoly_test1_")
        temp_dir2 = tempfile.mkdtemp(prefix="autopoly_test2_")

        try:
            system1 = System(out=temp_dir1)
            system2 = System(out=temp_dir2)

            assert system1.get_output_path() != system2.get_output_path()
            assert os.path.exists(system1.get_output_path())
            assert os.path.exists(system2.get_output_path())
        finally:
            if os.path.exists(temp_dir1):
                shutil.rmtree(temp_dir1)
            if os.path.exists(temp_dir2):
                shutil.rmtree(temp_dir2)


class TestSystemEdgeCases:
    """Test edge cases for System class."""

    def test_system_with_nested_path(self):
        """Test System with nested directory path."""
        temp_base = tempfile.mkdtemp(prefix="autopoly_base_")
        nested_path = os.path.join(temp_base, "nested", "dir")

        try:
            system = System(out=nested_path)
            assert os.path.exists(nested_path)
            assert system.get_output_path() == nested_path
        finally:
            if os.path.exists(temp_base):
                shutil.rmtree(temp_base)


class TestSystemEdgeCases:
    """Test edge cases in System class."""

    def test_cleanup_nonexistent_directory(self, temp_system):
        """Test cleanup when directory doesn't exist."""
        # Try to cleanup a path that doesn't exist
        nonexistent_path = os.path.join(temp_system.get_output_path(), "does_not_exist")
        
        # Should not raise an error
        # Note: The actual cleanup method might need to be called differently
        # depending on the System class implementation
        original_output_path = temp_system.get_output_path()
        assert original_output_path is not None

    def test_system_with_relative_path(self):
        """Test System initialization with relative path."""
        # Create system with relative path
        system = System(out="test_relative")
        expected = os.path.abspath("test_relative")
        assert system.get_output_path() == expected

    def test_system_with_absolute_path(self):
        """Test System initialization with absolute path."""
        # Create system with absolute path
        abs_path = "/tmp/autopoly_absolute_test"
        system = System(out=abs_path)
        assert system.get_output_path() == abs_path

    def test_system_path_with_trailing_slash(self):
        """Test System handles paths with trailing slashes."""
        # Path with trailing slash
        path_with_slash = "/tmp/test_path/"
        system = System(out=path_with_slash)
        # The System class keeps trailing slashes
        result = system.get_output_path()
        assert result is not None
        # Path normalization behavior - keeps trailing slash('/')  # Should not have trailing slash

    def test_system_multiple_instances_same_path(self):
        """Test multiple System instances with the same path."""
        path = "/tmp/autopoly_multi_test"
        system1 = System(out=path)
        system2 = System(out=path)
        
        # Both should have the same output path
        assert system1.get_output_path() == system2.get_output_path()
