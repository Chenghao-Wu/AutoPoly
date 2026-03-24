# -*- coding: utf-8 -*-
"""LangChain integration tests for AutoPoly tools.py."""
import json

import pytest

try:
    from langchain_core.tools import BaseTool
    HAS_LANGCHAIN = True
except ImportError:
    HAS_LANGCHAIN = False

pytestmark = pytest.mark.skipif(not HAS_LANGCHAIN, reason="langchain-core not installed")


@pytest.fixture(scope="module")
def tools():
    from AutoPoly.tools import get_autopoly_tools
    return get_autopoly_tools()


@pytest.fixture(scope="module")
def tool_map(tools):
    return {t.name: t for t in tools}


# ---------------------------------------------------------------------------
# TestToolsLoading
# ---------------------------------------------------------------------------

class TestToolsLoading:
    def test_get_autopoly_tools_returns_4_tools(self, tools):
        assert len(tools) == 4

    def test_tool_names(self, tool_map):
        expected = {"autopoly_info", "autopoly_generate_atomistic",
                    "autopoly_generate_bead_spring", "autopoly_describe_smiles"}
        assert set(tool_map.keys()) == expected

    def test_tools_have_descriptions(self, tools):
        for t in tools:
            assert t.description, f"Tool {t.name} has no description"

    def test_tools_are_langchain_base_tools(self, tools):
        for t in tools:
            assert isinstance(t, BaseTool)


# ---------------------------------------------------------------------------
# TestToolSchemas
# ---------------------------------------------------------------------------

class TestToolSchemas:
    def test_atomistic_schema_has_required_fields(self, tool_map):
        schema = tool_map["autopoly_generate_atomistic"].args_schema
        field_names = set(schema.model_fields.keys())
        for f in ("name", "force_field", "polymers", "molecules", "output_dir", "placement_method"):
            assert f in field_names, f"Missing field '{f}' in atomistic schema"

    def test_bead_spring_schema_has_required_fields(self, tool_map):
        schema = tool_map["autopoly_generate_bead_spring"].args_schema
        field_names = set(schema.model_fields.keys())
        for f in ("name", "n_chains", "bead_types", "sequence", "topology",
                   "bond_style", "pair_style", "generation_method", "density"):
            assert f in field_names, f"Missing field '{f}' in bead-spring schema"

    def test_schema_field_descriptions_present(self, tool_map):
        for tool_name in ("autopoly_generate_atomistic", "autopoly_generate_bead_spring"):
            schema = tool_map[tool_name].args_schema
            for name, field_info in schema.model_fields.items():
                assert field_info.description, (
                    f"Tool '{tool_name}' field '{name}' has no description"
                )


# ---------------------------------------------------------------------------
# TestToolInvocation
# ---------------------------------------------------------------------------

class TestToolInvocation:
    def test_info_tool_invoke(self, tool_map):
        result = tool_map["autopoly_info"].invoke({})
        data = json.loads(result)
        assert "force_fields" in data
        assert "examples" in data
        assert "limits" in data

    def test_describe_smiles_invoke(self, tool_map):
        result = tool_map["autopoly_describe_smiles"].invoke({"smiles": "CC"})
        data = json.loads(result)
        assert "atoms" in data
        assert data["has_wildcards"] is False

    def test_describe_smiles_wildcard_invoke(self, tool_map):
        result = tool_map["autopoly_describe_smiles"].invoke({"smiles": "[*]CC[*]"})
        data = json.loads(result)
        assert data["wildcard_count"] == 2
        assert data["has_wildcards"] is True

    def test_bead_spring_generate_invoke(self, tool_map, tmp_path):
        result = tool_map["autopoly_generate_bead_spring"].invoke({
            "name": "test_bs",
            "output_dir": str(tmp_path),
            "n_chains": 2,
            "bead_types": [{"name": "A", "mass": 1.0, "epsilon": 1.0, "sigma": 1.0}],
            "sequence": [["A", 10]],
            "topology": "linear",
            "bond_style": "harmonic",
            "pair_style": "lj",
            "generation_method": "saw",
        })
        data = json.loads(result)
        assert data["success"] is True
        assert len(data["files_created"]) > 0

    def test_bead_spring_generate_invalid_invoke(self, tool_map, tmp_path):
        result = tool_map["autopoly_generate_bead_spring"].invoke({
            "name": "bad",
            "output_dir": str(tmp_path),
            "n_chains": 2,
            "bead_types": [{"name": "A", "mass": 1.0, "epsilon": 1.0, "sigma": 1.0}],
            "sequence": [["A", 10]],
            "topology": "invalid_topology",
            "bond_style": "harmonic",
            "pair_style": "lj",
            "generation_method": "saw",
        })
        data = json.loads(result)
        assert data["success"] is False
        assert len(data["errors"]) > 0

    def test_atomistic_generate_no_polymers_invoke(self, tool_map, tmp_path):
        result = tool_map["autopoly_generate_atomistic"].invoke({
            "name": "empty",
            "output_dir": str(tmp_path),
            "force_field": "oplsaa",
        })
        data = json.loads(result)
        assert data["success"] is False
        assert any("polymers" in e or "molecules" in e for e in data["errors"])
