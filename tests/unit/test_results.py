# -*- coding: utf-8 -*-
"""Tests for AutoPoly result types."""
import json

from AutoPoly.core.results import GenerationResult, Result, ValidationResult


class TestResult:
    def test_success_result(self):
        r = Result(success=True)
        assert r.success is True
        assert r.errors == []
        assert r.warnings == []

    def test_failure_result(self):
        r = Result(success=False, errors=["something broke"])
        assert r.success is False
        assert r.errors == ["something broke"]

    def test_to_dict(self):
        r = Result(success=True, warnings=["heads up"])
        d = r.to_dict()
        assert d == {"success": True, "errors": [], "warnings": ["heads up"]}

    def test_json_serialization(self):
        r = Result(success=True)
        parsed = json.loads(str(r))
        assert parsed["success"] is True


class TestValidationResult:
    def test_with_suggestions(self):
        vr = ValidationResult(
            success=False,
            errors=["bad field"],
            suggestions=["try this instead"],
        )
        d = vr.to_dict()
        assert d["suggestions"] == ["try this instead"]
        assert "suggestions" in json.loads(str(vr))


class TestGenerationResult:
    def test_success_with_files(self):
        gr = GenerationResult(
            success=True,
            output_dir="/tmp/out",
            data_file="system.data",
            files_created=["system.data", "system.in"],
            metadata={"force_field": "oplsaa"},
        )
        d = gr.to_dict()
        assert d["output_dir"] == "/tmp/out"
        assert d["data_file"] == "system.data"
        assert len(d["files_created"]) == 2
        assert d["metadata"]["force_field"] == "oplsaa"

    def test_failure_no_files(self):
        gr = GenerationResult(success=False, errors=["bad config"])
        d = gr.to_dict()
        assert d["output_dir"] is None
        assert d["files_created"] == []
