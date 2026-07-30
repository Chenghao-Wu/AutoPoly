# -*- coding: utf-8 -*-
"""
Result types for the AutoPoly agent API.

All result types serialize to JSON via __str__ for easy agent consumption.
"""
import dataclasses
import json
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class Result:
    """Base result type."""
    success: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)

    def __str__(self) -> str:
        return json.dumps(self.to_dict(), indent=2)


@dataclass
class ValidationResult(Result):
    """Result of config validation with actionable fix suggestions."""
    suggestions: List[str] = field(default_factory=list)


@dataclass
class GenerationResult(Result):
    """Result of polymer system generation."""
    output_dir: Optional[str] = None
    data_file: Optional[str] = None
    files_created: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
