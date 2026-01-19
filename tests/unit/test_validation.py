"""
Tests for SMILES validation.

This module tests the validation functions for SMILES strings and other
user inputs to prevent injection attacks and ensure data integrity.
"""

import pytest
from AutoPoly.validation import validate_smiles, validate_smiles_list
from AutoPoly.exceptions import ValidationError


class TestSMILESValidation:
    """Test SMILES validation."""

    def test_valid_smiles_with_wildcard(self):
        """Test valid SMILES with connection points."""
        assert validate_smiles("[*]CC[*]") is True

    def test_valid_smiles_multiple_wildcards(self):
        """Test valid SMILES with multiple connection points."""
        assert validate_smiles("[*]CC([*])C(=O)OC") is True

    def test_valid_smiles_no_wildcard(self):
        """Test valid SMILES without wildcards (with allow_wildcards=False)."""
        assert validate_smiles("CCO", allow_wildcards=False) is True

    def test_invalid_smiles(self):
        """Test that invalid SMILES raises ValidationError."""
        with pytest.raises(ValidationError, match="Invalid SMILES"):
            validate_smiles("INVALID_SMILES")

    def test_missing_wildcard_with_default(self):
        """Test that SMILES without connection points raises ValidationError (default)."""
        with pytest.raises(ValidationError, match="must contain \\[\\*\\]"):
            validate_smiles("CC")

    def test_wildcard_with_allow_false(self):
        """Test that wildcards are rejected when allow_wildcards=False."""
        with pytest.raises(ValidationError, match="should not contain"):
            validate_smiles("[*]CC[*]", allow_wildcards=False)

    def test_empty_smiles_raises_error(self):
        """Test that empty SMILES raises ValidationError."""
        with pytest.raises(ValidationError, match="must be a non-empty string"):
            validate_smiles("")

    def test_non_string_smiles_raises_error(self):
        """Test that non-string SMILES raises ValidationError."""
        with pytest.raises(ValidationError, match="must be a non-empty string"):
            validate_smiles(123)

    def test_validate_smiles_list_valid(self):
        """Test validating a list of valid SMILES strings."""
        assert validate_smiles_list(["[*]CC[*]", "[*]C=C[*]", "[*]CC[*]"]) is True

    def test_validate_smiles_list_empty(self):
        """Test that empty list raises ValidationError."""
        with pytest.raises(ValidationError, match="must be a non-empty list"):
            validate_smiles_list([])

    def test_validate_smiles_list_invalid_element(self):
        """Test that invalid SMILES in list raises ValidationError."""
        with pytest.raises(ValidationError, match="position 1 is invalid"):
            validate_smiles_list(["[*]CC[*]", "invalid", "[*]CC[*]"])

    def test_complex_valid_smiles(self):
        """Test validation of complex SMILES."""
        assert validate_smiles("[*]CC([*])C(=O)OC") is True
