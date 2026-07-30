# Installation

## Requirements

- **Python** 3.10 or later
- **RDKit** ≥ 2022.09.1 and **NumPy** (installed automatically as dependencies)
- **LAMMPS** — required to *run* the simulations AutoPoly prepares (not needed for file generation)

[Moltemplate](https://moltemplate.org/) is bundled with AutoPoly under `AutoPoly/extern/` — no separate install is needed for the standard workflow.

## Install AutoPoly

Clone the repository and install in editable mode:

```bash
git clone https://github.com/WuGroup-XJTLU/AutoPoly.git
cd AutoPoly
pip install -e .
```

### Optional extras

| Extra | Install | Provides |
|---|---|---|
| `dev` | `pip install -e ".[dev]"` | Test suite: pytest, pytest-cov, pytest-mock |
| `agent` | `pip install -e ".[agent]"` | LangChain tool wrappers for the [agent API](../guides/agent.md) |
| `docs` | `pip install -e ".[docs]"` | Build this documentation site locally with `mkdocs serve` |

## Verify the installation

```python
import AutoPoly
print(AutoPoly.__version__)   # 1.0.0

from AutoPoly import System, Polymer, Polymerization
print("All components available")
```

You can also check the command-line interface:

```bash
autopoly info
```

This prints a JSON summary of force fields, topologies, limits, and example configs.

## Building the documentation locally

```bash
pip install -e ".[docs]"
mkdocs serve     # live preview at http://127.0.0.1:8000
mkdocs build     # static site in ./site/
```

## Troubleshooting installation

Having problems? See the [Installation Issues](../guides/troubleshooting.md#installation-issues) section of the Troubleshooting guide.
