# -*- coding: utf-8 -*-
"""Tests for AutoPoly agent API."""
import json
import os

import pytest

from AutoPoly import agent
from AutoPoly.results import GenerationResult, ValidationResult


class TestInfo:
    def test_returns_dict(self):
        result = agent.info()
        assert isinstance(result, dict)

    def test_has_required_keys(self):
        result = agent.info()
        assert "force_fields" in result
        assert "topologies" in result
        assert "tacticities" in result
        assert "bead_spring_options" in result
        assert "limits" in result
        assert "examples" in result

    def test_force_fields_have_descriptions(self):
        result = agent.info()
        ffs = result["force_fields"]
        assert "oplsaa" in ffs
        assert "gaff2" in ffs
        assert isinstance(ffs["oplsaa"], str)

    def test_limits_are_integers(self):
        result = agent.info()
        limits = result["limits"]
        assert isinstance(limits["max_dop"], int)
        assert limits["max_dop"] == 10000

    def test_examples_are_valid_configs(self):
        result = agent.info()
        for name, config in result["examples"].items():
            assert "type" in config
            assert config["type"] in ("atomistic", "bead_spring")

    def test_json_serializable(self):
        result = agent.info()
        json_str = json.dumps(result)
        assert json_str  # no exception


class TestValidate:
    def test_missing_type(self):
        vr = agent.validate({})
        assert isinstance(vr, ValidationResult)
        assert vr.success is False
        assert any("type" in e for e in vr.errors)

    def test_unknown_type(self):
        vr = agent.validate({"type": "quantum"})
        assert vr.success is False
        assert any("quantum" in e for e in vr.errors)

    def test_atomistic_minimal_valid(self):
        vr = agent.validate({
            "type": "atomistic",
            "name": "test",
            "polymers": [{
                "chain_num": 2,
                "sequence": ["CC[*]", "[*]CC[*]", "[*]CC"],
            }],
        })
        assert vr.success is True
        assert vr.errors == []

    def test_atomistic_missing_name(self):
        vr = agent.validate({
            "type": "atomistic",
            "polymers": [{"chain_num": 1, "sequence": ["CC[*]", "[*]CC"]}],
        })
        assert vr.success is False
        assert any("name" in e for e in vr.errors)

    def test_atomistic_bad_force_field(self):
        vr = agent.validate({
            "type": "atomistic",
            "name": "test",
            "force_field": "amber",
            "polymers": [{"chain_num": 1, "sequence": ["CC[*]", "[*]CC"]}],
        })
        assert vr.success is False
        assert any("amber" in e for e in vr.errors)
        # should have a suggestion
        assert len(vr.suggestions) > 0

    def test_atomistic_no_polymers_or_molecules(self):
        vr = agent.validate({
            "type": "atomistic",
            "name": "test",
        })
        assert vr.success is False
        assert any("polymers" in e or "molecules" in e for e in vr.errors)

    def test_atomistic_bad_wildcard_first(self):
        vr = agent.validate({
            "type": "atomistic",
            "name": "test",
            "polymers": [{
                "chain_num": 1,
                "sequence": ["[*]CC[*]", "[*]CC[*]", "[*]CC"],
            }],
        })
        assert vr.success is False
        assert any("First monomer" in e for e in vr.errors)

    def test_atomistic_bad_wildcard_last(self):
        vr = agent.validate({
            "type": "atomistic",
            "name": "test",
            "polymers": [{
                "chain_num": 1,
                "sequence": ["CC[*]", "[*]CC[*]", "[*]CC[*]"],
            }],
        })
        assert vr.success is False
        assert any("Last monomer" in e for e in vr.errors)

    def test_atomistic_molecule_with_wildcard(self):
        vr = agent.validate({
            "type": "atomistic",
            "name": "test",
            "molecules": [{"count": 10, "smiles": "[*]CC"}],
        })
        assert vr.success is False
        assert any("wildcard" in e.lower() for e in vr.errors)

    def test_bead_spring_minimal_valid(self):
        vr = agent.validate({
            "type": "bead_spring",
            "name": "test",
            "n_chains": 10,
            "bead_types": [{"name": "A"}],
            "sequence": [["A", 50]],
        })
        assert vr.success is True

    def test_bead_spring_missing_n_chains(self):
        vr = agent.validate({
            "type": "bead_spring",
            "name": "test",
            "bead_types": [{"name": "A"}],
            "sequence": [["A", 50]],
        })
        assert vr.success is False
        assert any("n_chains" in e for e in vr.errors)

    def test_bead_spring_bad_topology(self):
        vr = agent.validate({
            "type": "bead_spring",
            "name": "test",
            "n_chains": 10,
            "bead_types": [{"name": "A"}],
            "sequence": [["A", 50]],
            "topology": "star",
        })
        assert vr.success is False
        assert any("star" in e for e in vr.errors)

    def test_bead_spring_bad_bond_style(self):
        vr = agent.validate({
            "type": "bead_spring",
            "name": "test",
            "n_chains": 10,
            "bead_types": [{"name": "A"}],
            "sequence": [["A", 50]],
            "bond_style": "morse",
        })
        assert vr.success is False
        assert any("morse" in e for e in vr.errors)

    def test_never_raises(self):
        # Even with garbage input, validate should not raise
        vr = agent.validate({"type": "atomistic", "polymers": "not a list"})
        # Should either succeed or fail gracefully
        assert isinstance(vr, ValidationResult)


class TestDescribeSmiles:
    def test_simple_smiles(self):
        result = agent.describe_smiles("CC")
        assert "atoms" in result
        assert "elements" in result
        assert result["has_wildcards"] is False

    def test_wildcard_smiles(self):
        result = agent.describe_smiles("[*]CC[*]")
        assert result["has_wildcards"] is True
        assert result["wildcard_count"] == 2

    def test_invalid_smiles(self):
        result = agent.describe_smiles("not_a_smiles_at_all_xyz")
        assert "error" in result


class TestSuggestForceField:
    def test_organic(self):
        suggestions = agent.suggest_force_field(["CC", "CCO"])
        assert "oplsaa" in suggestions
        assert "gaff2" in suggestions

    def test_always_includes_dreiding(self):
        suggestions = agent.suggest_force_field(["CC"])
        assert "dreiding" in suggestions


class TestGenerate:
    def test_bead_spring_generates(self, tmp_path):
        config = {
            "type": "bead_spring",
            "name": "test_bs",
            "output_dir": str(tmp_path),
            "n_chains": 2,
            "bead_types": [{"name": "A"}],
            "sequence": [["A", 10]],
            "generation_method": "geometric",
        }
        result = agent.generate(config)
        assert isinstance(result, GenerationResult)
        assert result.success is True
        assert result.output_dir is not None
        assert result.files_created  # at least one file

    def test_bead_spring_with_options(self, tmp_path):
        config = {
            "type": "bead_spring",
            "name": "test_opts",
            "output_dir": str(tmp_path),
            "n_chains": 3,
            "bead_types": [{"name": "A"}, {"name": "B", "epsilon": 1.5}],
            "sequence": [["A", 5], ["B", 5]],
            "topology": "linear",
            "bond_style": "fene",
            "pair_style": "wca",
            "generation_method": "saw",
        }
        result = agent.generate(config)
        assert result.success is True
        assert result.metadata["bond_style"] == "fene"

    def test_invalid_config_returns_error(self):
        result = agent.generate({"type": "bead_spring"})
        assert isinstance(result, GenerationResult)
        assert result.success is False
        assert len(result.errors) > 0

    def test_generate_json_serializable(self, tmp_path):
        config = {
            "type": "bead_spring",
            "name": "test_json",
            "output_dir": str(tmp_path),
            "n_chains": 2,
            "bead_types": [{"name": "A"}],
            "sequence": [["A", 5]],
            "generation_method": "geometric",
        }
        result = agent.generate(config)
        json_str = str(result)
        parsed = json.loads(json_str)
        assert parsed["success"] is True
