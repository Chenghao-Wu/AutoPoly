"""Pytest configuration and fixtures for AutoPoly tests."""

import pytest


# pytest-randomly will automatically seed all tests for reproducibility
# No manual seeding needed - the plugin handles it
# If a test fails, pytest-randomly shows the seed: FAILED - Random seed: 1234567890
# Reproduce with: pytest --randomly-seed=1234567890

# Add custom fixtures here as needed (when you see code duplication 3+ times)
# For now, use built-in pytest fixtures: tmp_path, mocker, monkeypatch
