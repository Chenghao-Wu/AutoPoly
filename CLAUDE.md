# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

AutoPoly generates polymer structures and prepares them for molecular dynamics simulations with LAMMPS. It transforms SMILES strings into complete LAMMPS input files through a pipeline: complement SMILES → monomer `.lt` templates → polymer chains → moltemplate → LAMMPS data files.

## Build & Test Commands

```bash
pip install -e .                          # Install in dev mode
pip install -e ".[dev]"                   # Install with test dependencies
pytest tests/                             # Run all tests
pytest tests/unit/                        # Unit tests only
pytest tests/integration/                 # Integration tests only
pytest tests/ -m "not slow"              # Skip slow tests
pytest tests/unit/test_polymer.py         # Single test file
pytest tests/unit/test_polymer.py::TestPolymer::test_method  # Single test
pytest --cov=AutoPoly                    # With coverage report
```

Dependencies: `numpy`, `rdkit` (>=2022.09.1). External: `moltemplate` (bundled in `extern/`), LAMMPS (user-installed).

## Architecture

### Pipeline (3-step user API, multi-stage internal)

```
User: System(out=dir) → Polymer/Molecule → Polymerization(force_field=...)
                                              ↓
Internal: MonomerGenerator (SMILES→.lt) → ForceFieldManager (atom typing, params)
          → WorkflowManager (chain placement, system.lt) → moltemplate.sh → LAMMPS files
```

### Key Modules

- **`polymerization.py`** — Public entry point. Orchestrates the entire workflow by delegating to WorkflowManager.
- **`workflow.py`** (~1000 lines) — Core pipeline logic. Processes models, generates `.lt` files, invokes moltemplate, organizes output. Has two code paths: deterministic grid placement (`make_poly_lt`/`make_system_lt`) and MC placement (`make_poly_lt_mc`/`make_system_lt_mc`).
- **`monomer_generator.py`** (~2000 lines) — Converts complement SMILES to moltemplate `.lt` files. Generates 6 variants per monomer (first/middle/last × 2 chirality). Uses RDKit for geometry optimization and atom typing via SMARTS patterns.
- **`force_field.py`** (~2700 lines) — `ForceFieldManager` handles all 6 force fields. Atom typing, parameter lookup, dihedral modification. GAFF/GAFF2 use `gaff_analysis.py` to filter parameters to used atom types.
- **`polymer.py`** / **`molecule.py`** — Data classes. Both implement `get_sequenceSet()`, `get_mer_set()` so WorkflowManager handles them uniformly. Polymer uses complement SMILES (wildcards `[*]`); Molecule uses regular SMILES.
- **`mc/`** — Monte Carlo module: `chain_growth.py` (SAW algorithm), `collision.py` (cell-linked list spatial hashing), `placement.py` (random molecular placement).
- **`bead_spring.py`** — Coarse-grained bead-spring models with direct LAMMPS data file generation (bypasses moltemplate).
- **`extern/`** — Bundled moltemplate with 52 force field `.lt` files in `extern/moltemplate/force_fields/`. `rdlt.py` handles RDKit-based `.lt` operations.

### Complement SMILES Convention

AutoPoly's unique monomer specification: wildcard `[*]` position determines monomer role.
- First: `"CC[*]"` (1 wildcard, right) — chain start
- Middle: `"[*]CC[*]"` (2 wildcards) — chain interior
- Last: `"[*]CC"` (1 wildcard, left) — chain end
- Ring: all middle variants

### Force Fields

Six supported: `"oplsaa"`, `"lopls"`, `"gaff"`, `"gaff2"`, `"dreiding"`, `"compass"`. Files live in `AutoPoly/extern/moltemplate/force_fields/`. Atom typing uses manually-built SMARTS patterns (exploring BESMARTS for data-driven approach).

### Resource Limits (conf.py)

MAX_DOP = 10,000; MAX_SEQUENCE_LENGTH = 10,000; MAX_UNIQUE_MONOMERS = 100.

## Testing Conventions

- Tests use `pytest` with `tmp_path` fixture for file operations
- Markers: `@pytest.mark.slow`, `@pytest.mark.integration`
- CI runs Python 3.10/3.11/3.12 via GitHub Actions (`.github/workflows/test.yml`)
- Contract tests in `tests/integration/test_contracts.py` verify RDKit compatibility
- `conftest.py` uses `pytest-randomly` for reproducible test ordering
