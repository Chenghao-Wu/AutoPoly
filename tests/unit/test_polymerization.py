"""Tests for the Polymerization class."""

import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
from AutoPoly.polymerization import Polymerization
from AutoPoly.exceptions import ValidationError


class TestPolymerizationInitialization:
    """Test Polymerization class initialization."""

    @patch('AutoPoly.polymerization.ForceFieldManager')
    @patch('AutoPoly.polymerization.WorkflowManager')
    @patch('AutoPoly.polymerization.create_working_directory')
    def test_polymerization_init_with_gaff(self, mock_create_dir, mock_workflow_class, mock_ff_class, tmp_path):
        """Test initialization with gaff force field."""
        mock_system = MagicMock()
        mock_system.get_folder_path.return_value = str(tmp_path)
        mock_create_dir.return_value = str(tmp_path / "test" / "moltemplate")
        mock_workflow = MagicMock()
        mock_workflow_class.return_value = mock_workflow

        poly = Polymerization(name="test", system=mock_system, model=None, run=False, force_field="gaff")

        assert poly.force_field == "gaff"
        assert "gaff.lt" in poly.path_oplsaaprm
        assert poly.gaff_analyzer is not None

    @patch('AutoPoly.polymerization.ForceFieldManager')
    @patch('AutoPoly.polymerization.WorkflowManager')
    @patch('AutoPoly.polymerization.create_working_directory')
    def test_polymerization_init_with_invalid_force_field_raises_system_exit(
        self, mock_create_dir, mock_workflow_class, mock_ff_class, tmp_path
    ):
        """Test that invalid force field raises ValidationError."""
        mock_system = MagicMock()
        mock_system.get_folder_path.return_value = str(tmp_path)

        with pytest.raises(ValidationError, match="Invalid force_field"):
            Polymerization(name="test", system=mock_system, model=None, run=False, force_field="invalid")

    @patch('AutoPoly.polymerization.ForceFieldManager')
    @patch('AutoPoly.polymerization.WorkflowManager')
    @patch('AutoPoly.polymerization.create_working_directory')
    def test_polymerization_init_creates_working_directory(
        self, mock_create_dir, mock_workflow_class, mock_ff_class, tmp_path
    ):
        """Test that working directory is created."""
        mock_system = MagicMock()
        mock_system.get_folder_path.return_value = str(tmp_path)
        mock_create_dir.return_value = str(tmp_path / "test" / "moltemplate")
        mock_workflow = MagicMock()
        mock_workflow_class.return_value = mock_workflow

        Polymerization(name="test", system=mock_system, model=None, run=False, force_field="oplsaa")

        # Verify create_working_directory was called
        mock_create_dir.assert_called_once()


class TestPolymerizationDelegation:
    """Test Polymerization delegation methods."""

    @patch('AutoPoly.polymerization.ForceFieldManager')
    @patch('AutoPoly.polymerization.WorkflowManager')
    @patch('AutoPoly.polymerization.create_working_directory')
    @patch('AutoPoly.polymerization.monomer_processing.generate_monomer_from_psmiles')
    def test_generate_monomer_from_psmiles_delegates_to_monomer_processing(
        self, mock_generate, mock_create_dir, mock_workflow_class, mock_ff_class, tmp_path
    ):
        """Test that generate_monomer_from_psmiles delegates correctly."""
        mock_system = MagicMock()
        mock_system.get_folder_path.return_value = str(tmp_path)
        mock_create_dir.return_value = str(tmp_path / "test" / "moltemplate")
        mock_workflow = MagicMock()
        mock_workflow_class.return_value = mock_workflow
        mock_generate.return_value = ("monomer_0", 1)

        poly = Polymerization(name="test", system=mock_system, model=None, run=False, force_field="oplsaa")

        result_name, result_counter = poly.generate_monomer_from_psmiles("[*]C=C[*]")

        # Verify delegation with correct parameters
        mock_generate.assert_called_once()
        call_args = mock_generate.call_args[0]
        assert call_args[0] == "[*]C=C[*]"  # psmiles
        assert call_args[1] == poly.path_cwd  # path_cwd
        assert call_args[2] == "oplsaa"  # force_field
        # Verify counter was updated
        assert poly._smiles_to_name_counter == 1
        assert result_name == "monomer_0"

    @patch('AutoPoly.polymerization.ForceFieldManager')
    @patch('AutoPoly.polymerization.WorkflowManager')
    @patch('AutoPoly.polymerization.create_working_directory')
    @patch('AutoPoly.polymerization.monomer_processing.generate_molecule_from_smiles')
    def test_generate_molecule_from_smiles_delegates_to_monomer_processing(
        self, mock_generate, mock_create_dir, mock_workflow_class, mock_ff_class, tmp_path
    ):
        """Test that generate_molecule_from_smiles delegates correctly."""
        mock_system = MagicMock()
        mock_system.get_folder_path.return_value = str(tmp_path)
        mock_create_dir.return_value = str(tmp_path / "test" / "moltemplate")
        mock_workflow = MagicMock()
        mock_workflow_class.return_value = mock_workflow
        mock_generate.return_value = ("water.lt", 1)

        poly = Polymerization(name="test", system=mock_system, model=None, run=False, force_field="oplsaa")

        result_filename, result_counter = poly.generate_molecule_from_smiles("O", "water")

        # Verify delegation with correct parameters
        mock_generate.assert_called_once()
        call_args = mock_generate.call_args[0]
        assert call_args[0] == "O"  # smiles
        assert call_args[1] == "water"  # molecule_name
        assert call_args[2] == poly.path_cwd  # path_cwd
        assert call_args[3] == "oplsaa"  # force_field
        # Verify counter was updated
        assert poly._smiles_to_name_counter == 1
        assert result_filename == "water.lt"

    @patch('AutoPoly.polymerization.ForceFieldManager')
    @patch('AutoPoly.polymerization.WorkflowManager')
    @patch('AutoPoly.polymerization.create_working_directory')
    @patch('AutoPoly.polymerization.monomer_processing.n_monomer_atoms')
    def test_n_monomer_atoms_delegates_to_monomer_processing(
        self, mock_n_atoms, mock_create_dir, mock_workflow_class, mock_ff_class, tmp_path
    ):
        """Test that n_monomer_atoms delegates correctly."""
        mock_system = MagicMock()
        mock_system.get_folder_path.return_value = str(tmp_path)
        mock_create_dir.return_value = str(tmp_path / "test" / "moltemplate")
        mock_workflow = MagicMock()
        mock_workflow_class.return_value = mock_workflow
        mock_n_atoms.return_value = 5

        poly = Polymerization(name="test", system=mock_system, model=None, run=False, force_field="oplsaa")

        result = poly.n_monomer_atoms("test.lt")

        # Verify delegation
        mock_n_atoms.assert_called_once_with("test.lt", poly.path_cwd)
        assert result == 5

    @patch('AutoPoly.polymerization.ForceFieldManager')
    @patch('AutoPoly.polymerization.WorkflowManager')
    @patch('AutoPoly.polymerization.create_working_directory')
    @patch('AutoPoly.polymerization.monomer_processing.extract_element_from_atom')
    def test_extract_element_from_atom_delegates_to_monomer_processing(
        self, mock_extract, mock_create_dir, mock_workflow_class, mock_ff_class, tmp_path
    ):
        """Test that extract_element_from_atom delegates correctly."""
        mock_system = MagicMock()
        mock_system.get_folder_path.return_value = str(tmp_path)
        mock_create_dir.return_value = str(tmp_path / "test" / "moltemplate")
        mock_workflow = MagicMock()
        mock_workflow_class.return_value = mock_workflow
        mock_extract.return_value = "C"

        poly = Polymerization(name="test", system=mock_system, model=None, run=False, force_field="oplsaa")

        result = poly.extract_element_from_atom("$atom:C1")

        # Verify delegation
        mock_extract.assert_called_once_with("$atom:C1")
        assert result == "C"

    @patch('AutoPoly.polymerization.ForceFieldManager')
    @patch('AutoPoly.polymerization.WorkflowManager')
    @patch('AutoPoly.polymerization.create_working_directory')
    @patch('AutoPoly.polymerization.monomer_processing.read_lt_end_atoms')
    def test_read_lt_end_atoms_delegates_to_monomer_processing(
        self, mock_read, mock_create_dir, mock_workflow_class, mock_ff_class, tmp_path
    ):
        """Test that read_lt_end_atoms delegates correctly."""
        mock_system = MagicMock()
        mock_system.get_folder_path.return_value = str(tmp_path)
        mock_create_dir.return_value = str(tmp_path / "test" / "moltemplate")
        mock_workflow = MagicMock()
        mock_workflow_class.return_value = mock_workflow
        mock_read.return_value = ("C1", "C2")

        poly = Polymerization(name="test", system=mock_system, model=None, run=False, force_field="oplsaa")

        result = poly.read_lt_end_atoms("test.lt")

        # Verify delegation
        mock_read.assert_called_once_with("test.lt")
        assert result == ("C1", "C2")


class TestPolymerizationFileManagement:
    """Test Polymerization file management methods."""

    @patch('AutoPoly.polymerization.ForceFieldManager')
    @patch('AutoPoly.polymerization.WorkflowManager')
    @patch('AutoPoly.polymerization.create_working_directory')
    @patch('AutoPoly.polymerization.get_rid_of_lj_cut_coul_long')
    def test_get_rid_of_lj_cut_coul_long_delegates(
        self, mock_get_rid, mock_create_dir, mock_workflow_class, mock_ff_class, tmp_path
    ):
        """Test that get_rid_of_lj_cut_coul_long delegates correctly."""
        mock_system = MagicMock()
        mock_system.get_folder_path.return_value = str(tmp_path)
        mock_create_dir.return_value = str(tmp_path / "test" / "moltemplate")
        mock_workflow = MagicMock()
        mock_workflow_class.return_value = mock_workflow

        poly = Polymerization(name="test", system=mock_system, model=None, run=False, force_field="oplsaa")
        poly.get_rid_of_lj_cut_coul_long()

        # Verify delegation
        mock_get_rid.assert_called_once_with(poly.path_cwd)

    @patch('AutoPoly.polymerization.ForceFieldManager')
    @patch('AutoPoly.polymerization.WorkflowManager')
    @patch('AutoPoly.polymerization.create_working_directory')
    @patch('AutoPoly.polymerization.mv_files')
    def test_mv_files_delegates(self, mock_mv, mock_create_dir, mock_workflow_class, mock_ff_class, tmp_path):
        """Test that mv_files delegates correctly."""
        mock_system = MagicMock()
        mock_system.get_folder_path.return_value = str(tmp_path)
        mock_create_dir.return_value = str(tmp_path / "test" / "moltemplate")
        mock_workflow = MagicMock()
        mock_workflow_class.return_value = mock_workflow

        poly = Polymerization(name="test", system=mock_system, model=None, run=False, force_field="oplsaa")
        poly.mv_files()

        # Verify delegation
        mock_mv.assert_called_once_with(poly.path_cwd)

    @patch('AutoPoly.polymerization.ForceFieldManager')
    @patch('AutoPoly.polymerization.WorkflowManager')
    @patch('AutoPoly.polymerization.create_working_directory')
    @patch('AutoPoly.polymerization.monomer_processing.evaluate_offset')
    def test_evaluate_offset_delegates_and_updates(
        self, mock_eval, mock_create_dir, mock_workflow_class, mock_ff_class, tmp_path
    ):
        """Test that evaluate_offset delegates and updates self.offset."""
        mock_system = MagicMock()
        mock_system.get_folder_path.return_value = str(tmp_path)
        mock_create_dir.return_value = str(tmp_path / "test" / "moltemplate")
        mock_workflow = MagicMock()
        mock_workflow_class.return_value = mock_workflow
        mock_eval.return_value = 3.5

        poly = Polymerization(name="test", system=mock_system, model=None, run=False, force_field="oplsaa")
        initial_offset = poly.offset
        poly.evaluate_offset("test.lt")

        # Verify delegation with correct parameters
        mock_eval.assert_called_once()
        call_args = mock_eval.call_args[0]
        assert call_args[0] == "test.lt"
        assert call_args[1] == poly.path_cwd
        assert call_args[2] == poly.offset_spacing
        assert call_args[3] == initial_offset
        # Verify offset was updated
        assert poly.offset == 3.5


class TestPolymerizationWorkflow:
    """Test Polymerization workflow methods."""

    @patch('AutoPoly.polymerization.ForceFieldManager')
    @patch('AutoPoly.polymerization.WorkflowManager')
    @patch('AutoPoly.polymerization.create_working_directory')
    def test_make_lmp_data_file_delegates_to_workflow(
        self, mock_create_dir, mock_workflow_class, mock_ff_class, tmp_path
    ):
        """Test that make_lmp_data_file delegates to WorkflowManager."""
        mock_system = MagicMock()
        mock_system.get_folder_path.return_value = str(tmp_path)
        mock_create_dir.return_value = str(tmp_path / "test" / "moltemplate")
        mock_workflow = MagicMock()
        mock_workflow_class.return_value = mock_workflow

        poly = Polymerization(name="test", system=mock_system, model=None, run=False, force_field="oplsaa")
        poly.make_lmp_data_file_by_moltemplate()

        # Verify delegation
        mock_workflow.make_lmp_data_file_by_moltemplate.assert_called_once()


class TestPolymerizationSetters:
    """Test Polymerization setter methods."""

    @patch('AutoPoly.polymerization.ForceFieldManager')
    @patch('AutoPoly.polymerization.WorkflowManager')
    @patch('AutoPoly.polymerization.create_working_directory')
    def test_set_tacticity(self, mock_create_dir, mock_workflow_class, mock_ff_class, tmp_path):
        """Test setting tacticity."""
        mock_system = MagicMock()
        mock_system.get_folder_path.return_value = str(tmp_path)
        mock_create_dir.return_value = str(tmp_path / "test" / "moltemplate")
        mock_workflow = MagicMock()
        mock_workflow_class.return_value = mock_workflow

        poly = Polymerization(name="test", system=mock_system, model=None, run=False, force_field="oplsaa")
        poly.set_tacticity("isotactic")

        assert poly.tacticity == "isotactic"
