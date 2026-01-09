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


@pytest.fixture
def sample_monomer_lt_file():
    """
    Create a sample monomer .lt file for testing.

    This fixture creates a temporary .lt file with properly formatted
    moltemplate atom data for testing monomer processing functions.

    Returns:
        Path: Path to the created .lt file
    """
    fd, temp_path = tempfile.mkstemp(suffix='.lt')
    temp_file = Path(temp_path)

    lt_content = """# Test monomer file
write("Data Atoms") {
  $atom:C1  @atom:opls_135  1  0.0  0.0  0.0
  $atom:H2  @atom:opls_140  2  1.0  0.0  0.0
  $atom:H3  @atom:opls_140  3  0.0  1.0  0.0
  $atom:C4  @atom:opls_135  4  1.54  0.0  0.0
}
"""
    temp_file.write_text(lt_content)

    yield temp_file

    # Cleanup
    os.close(fd)
    temp_file.unlink()


@pytest.fixture
def sample_settings_file():
    """
    Create a sample system.in.settings file for testing.

    This fixture creates a temporary LAMMPS settings file with
    pair_coeff entries for testing the get_rid_of_lj_cut_coul_long function.

    Returns:
        Path: Path to the created settings file
    """
    fd, temp_path = tempfile.mkstemp(suffix='.in.settings')
    temp_file = Path(temp_path)

    settings_content = """# LAMMPS settings file
pair_style hybrid/overlay lj/cut/coul/long 10.0 10.0 coul/long 10.0
pair_coeff * * lj/cut/coul/long 0.0 0.0
pair_coeff 1 2 lj/cut/coul/long 0.1 2.5
pair_coeff 2 3 lj/cut/coul/long 0.2 3.0
"""
    temp_file.write_text(settings_content)

    yield temp_file

    # Cleanup
    os.close(fd)
    temp_file.unlink()


@pytest.fixture
def mock_rdkit_mol():
    """
    Create a mock RDKit molecule object for testing.

    This fixture creates a simple RDKit molecule (ethylene) for
    testing mechanism detection and pattern matching functions
    that require RDKit Mol objects.

    Returns:
        Chem.Mol: RDKit molecule object
    """
    from rdkit import Chem
    return Chem.MolFromSmiles("C=C")


@pytest.fixture
def mock_rdkit_mol_with_h():
    """
    Create a mock RDKit molecule with explicit hydrogens.

    Returns:
        Chem.Mol: RDKit molecule object with hydrogens
    """
    from rdkit import Chem
    mol = Chem.MolFromSmiles("CCO")
    return Chem.AddHs(mol)


@pytest.fixture
def sample_psmiles_list():
    """
    Provide a list of sample pSMILES strings for testing.

    Returns:
        list: List of pSMILES strings representing various monomers
    """
    return [
        "[*]C=C[*]",              # Polyethylene
        "[*]C=C(C)C(=O)OC[*]",    # PMMA
        "[*]C=C(C)c1ccccc1[*]",   # Polystyrene
        "[*]C(C)(C)C(=O)O[*]",    # PLA (esterification)
        "[*]C(=O)N[*]",           # Nylon (amidation)
    ]


@pytest.fixture
def temp_monomer_dir():
    """
    Create a temporary directory for monomer file operations.

    This fixture creates a temporary directory that can be used
    for testing monomer generation and file operations without
    polluting the actual monomer bank.

    Returns:
        Path: Path to temporary directory
    """
    temp_dir = Path(tempfile.mkdtemp(prefix="autopoly_monomer_test_"))

    yield temp_dir

    # Cleanup
    shutil.rmtree(temp_dir)
