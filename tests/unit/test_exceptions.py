"""
Tests for custom exception hierarchy.

This module tests the custom exceptions used throughout AutoPoly for proper
error handling and testing.
"""

import pytest
from AutoPoly import Polymer
from AutoPoly.core.exceptions import ValidationError, AutoPolyError


class TestExceptions:
    """Test custom exception hierarchy."""

    def test_zero_chain_num_raises_validation_error(self):
        """Test that zero chain_num raises ValidationError."""
        with pytest.raises(ValidationError, match="chain_num must be greater than 0"):
            Polymer(chain_num=0, sequence=["CC[*]", "[*]CC"])

    def test_invalid_topology_raises_value_error(self):
        """Test that invalid topology raises ValueError."""
        with pytest.raises(ValueError, match="topology must be"):
            Polymer(
                chain_num=1,
                sequence=["[*]CC[*]"],
                topology="invalid"
            )

    def test_invalid_tacticity_raises_value_error(self):
        """Test that invalid tacticity raises ValueError."""
        with pytest.raises(ValueError, match="tacticity must be"):
            Polymer(
                chain_num=1,
                sequence=["[*]CC[*]"],
                tacticity="invalid"
            )

    def test_empty_sequence_raises_validation_error(self):
        """Test that empty sequence raises ValidationError."""
        with pytest.raises(ValidationError, match="sequence cannot be empty"):
            Polymer(chain_num=1, sequence=[])

    def test_sequence_too_long_raises_validation_error(self):
        """Test that excessive sequence length raises ValidationError."""
        with pytest.raises(ValidationError, match="exceeds maximum"):
            Polymer(
                chain_num=1,
                sequence=["[*]CC[*]"] * 15000
            )

    def test_validation_error_is_auto_poly_error(self):
        """Test that ValidationError inherits from AutoPolyError."""
        try:
            Polymer(chain_num=0, sequence=["[*]CC[*]"])
            assert False, "Should have raised ValidationError"
        except AutoPolyError:
            # This is expected
            assert True

    def test_validation_error_catch_as_base_exception(self):
        """Test catching ValidationError as AutoPolyError."""
        with pytest.raises(AutoPolyError):
            Polymer(chain_num=0, sequence=["[*]CC[*]"])
