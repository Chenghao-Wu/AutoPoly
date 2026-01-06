"""
Pytest configuration and shared fixtures for AutoPoly tests.

This module provides reusable fixtures for testing the AutoPoly package,
including temporary system instances, sample polymers, and test data paths.
"""

import sys
import os
import tempfile
import shutil
from pathlib import Path

# Add the AutoPoly package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from AutoPoly.system import System
from AutoPoly.polymer import Polymer


@pytest.fixture
def temp_system():
    """
    Create a temporary System instance for testing.

    This fixture creates a System object with a temporary output directory
    that is automatically cleaned up after the test.

    Yields:
        System: A System instance with a temporary output directory
    """
    temp_dir = tempfile.mkdtemp(prefix="autopoly_test_")
    system = System(out=temp_dir)
    yield system
    # Cleanup
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)


@pytest.fixture
def sample_polymer():
    """
    Create a sample Polymer for testing.

    This fixture creates a basic Polymer instance with default parameters
    for testing Polymer class methods.

    Returns:
        Polymer: A sample Polymer instance with linear topology
    """
    polymer = Polymer(
        ChainNum=1,
        Sequence=["PE", "PE"],
        DOP=2,
        topology="linear",
        tacticity="atactic"
    )
    return polymer


@pytest.fixture
def monomer_bank_path():
    """
    Get the path to the monomer bank directory.

    This fixture provides the path to the example monomers for testing.

    Returns:
        Path: Path object pointing to the monomer bank directory
    """
    # Assumes tests are run from the AutoPoly directory
    autopoly_root = Path(__file__).parent.parent
    return autopoly_root / "example" / "monomers"


@pytest.fixture
def sample_ring_polymer():
    """
    Create a sample ring Polymer for testing.

    This fixture creates a Polymer instance with ring topology.

    Returns:
        Polymer: A sample ring Polymer instance
    """
    polymer = Polymer(
        ChainNum=1,
        Sequence=["PE", "PE", "PE"],
        DOP=3,
        topology="ring",
        tacticity="atactic"
    )
    return polymer


@pytest.fixture
def sample_isotactic_polymer():
    """
    Create a sample isotactic Polymer for testing.

    This fixture creates a Polymer instance with isotactic tacticity.

    Returns:
        Polymer: A sample isotactic Polymer instance
    """
    polymer = Polymer(
        ChainNum=1,
        Sequence=["PE", "PE", "PE"],
        DOP=3,
        topology="linear",
        tacticity="isotactic"
    )
    return polymer


@pytest.fixture
def sample_syndiotactic_polymer():
    """
    Create a sample syndiotactic Polymer for testing.

    This fixture creates a Polymer instance with syndiotactic tacticity.

    Returns:
        Polymer: A sample syndiotactic Polymer instance
    """
    polymer = Polymer(
        ChainNum=1,
        Sequence=["PE", "PE", "PE", "PE"],
        DOP=4,
        topology="linear",
        tacticity="syndiotactic"
    )
    return polymer
