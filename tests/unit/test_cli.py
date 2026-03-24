# -*- coding: utf-8 -*-
"""Tests for AutoPoly CLI."""
import json
import subprocess
import sys

import pytest


class TestCLI:
    def _run(self, *args):
        """Run autopoly CLI command via python -m."""
        result = subprocess.run(
            [sys.executable, "-m", "AutoPoly.cli"] + list(args),
            capture_output=True, text=True, timeout=30,
        )
        return result

    def test_info_returns_json(self):
        result = self._run("info")
        assert result.returncode == 0
        data = json.loads(result.stdout)
        assert "force_fields" in data
        assert "examples" in data

    def test_no_command_exits_1(self):
        result = self._run()
        assert result.returncode == 1

    def test_describe_smiles(self):
        result = self._run("describe", "CC")
        assert result.returncode == 0
        data = json.loads(result.stdout)
        assert "atoms" in data
        assert data["has_wildcards"] is False

    def test_describe_wildcard(self):
        result = self._run("describe", "[*]CC[*]")
        assert result.returncode == 0
        data = json.loads(result.stdout)
        assert data["wildcard_count"] == 2

    def test_validate_stdin(self, tmp_path):
        config = json.dumps({
            "type": "bead_spring",
            "name": "test",
            "n_chains": 5,
            "bead_types": [{"name": "A"}],
            "sequence": [["A", 10]],
        })
        result = subprocess.run(
            [sys.executable, "-m", "AutoPoly.cli", "validate", "--stdin"],
            capture_output=True, text=True, timeout=30,
            input=config,
        )
        assert result.returncode == 0
        data = json.loads(result.stdout)
        assert data["success"] is True

    def test_validate_invalid_config_exits_1(self):
        config = json.dumps({"type": "atomistic"})
        result = subprocess.run(
            [sys.executable, "-m", "AutoPoly.cli", "validate", "--stdin"],
            capture_output=True, text=True, timeout=30,
            input=config,
        )
        assert result.returncode == 1
        data = json.loads(result.stdout)
        assert data["success"] is False

    def test_validate_file(self, tmp_path):
        config_file = tmp_path / "config.json"
        config_file.write_text(json.dumps({
            "type": "bead_spring",
            "name": "test",
            "n_chains": 5,
            "bead_types": [{"name": "A"}],
            "sequence": [["A", 10]],
        }))
        result = self._run("validate", str(config_file))
        assert result.returncode == 0
        data = json.loads(result.stdout)
        assert data["success"] is True

    def test_generate_bead_spring(self, tmp_path):
        config_file = tmp_path / "config.json"
        config_file.write_text(json.dumps({
            "type": "bead_spring",
            "name": "cli_test",
            "output_dir": str(tmp_path),
            "n_chains": 2,
            "bead_types": [{"name": "A"}],
            "sequence": [["A", 5]],
            "generation_method": "geometric",
        }))
        result = self._run("generate", str(config_file))
        assert result.returncode == 0
        data = json.loads(result.stdout)
        assert data["success"] is True
        assert data["files_created"]
