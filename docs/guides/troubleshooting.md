# Troubleshooting

Common issues and solutions for AutoPoly.

## API Errors

### TypeError: unexpected keyword argument 'ChainNum' / 'DOP' / 'Sequence'

**Cause:** v0.x calling convention used with v1.0+ AutoPoly.

**Solution:** v1.0 uses snake_case, and DOP is always `len(sequence)`:

```python
# OLD (v0.x)
poly = Polymer(ChainNum=10, Sequence=["PE"], DOP=50)

# NEW (v1.0+)
poly = Polymer(
    chain_num=10,
    sequence=["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"],   # DOP = 50
)
print(poly.dop)   # 50
```

## Sequence Errors

### ValidationError: sequence cannot be empty

Provide at least one monomer: `sequence=["[*]CC[*]"] * 10`.

### ValidationError: sequence length exceeds maximum

Sequences are capped at 10 000 entries (a safety limit against memory exhaustion). Split the work or reduce DOP:

```python
sequence = ["[*]CC[*]"] * 5000   # within the limit
```

### ValidationError: number of unique monomers exceeds maximum

A sequence may contain at most 100 **unique** monomer types. Reuse monomer strings instead of generating near-duplicate variants:

```python
pe = "[*]CC[*]"
ps = "[*]CC([*])c1ccccc1"
sequence = [pe] * 50 + [ps] * 50   # only 2 unique types
```

## SMILES Errors

### ValidationError: invalid SMILES in sequence

Check, in order:

1. **Wildcard count** — first = 1, middle = 2, last = 1. See [Complement SMILES](complement-smiles.md).
2. **Brackets** — wildcards must be `[*]`, never a bare `*`.
3. **SMILES syntax** — unbalanced parentheses, invalid valences.

Quick validation helper:

```python
from rdkit import Chem

def valid(smiles: str) -> bool:
    return Chem.MolFromSmiles(smiles.replace("[*]", "C")) is not None

assert valid("[*]CC[*]")
```

### Ring polymer failures

Ring polymers must use **only middle variants** (2 wildcards everywhere) and `topology="ring"`:

```python
poly = Polymer(chain_num=10, sequence=["[*]CC[*]"] * 50, topology="ring")
```

## Force Field Errors

### ValidationError: invalid force_field

Valid values: `"oplsaa"`, `"lopls"`, `"gaff"`, `"gaff2"`, `"dreiding"`, `"compass"`. Watch for typos like `"opls"`.

### GAFF simulations give poor energies

GAFF/GAFF2 runs use automatically assigned **Gasteiger charges**, which are approximate. For production, compute AM1-BCC charges with Antechamber (or RESP with your QM package) and edit them into `system.in.charges`. See the [Force Fields guide](force-fields.md).

## Moltemplate Errors

### Monomer .lt file generation failed

Causes and fixes:

- **Invalid SMILES** — validate as shown above
- **Unsupported functional groups** — simplify the structure, or try a more general force field (`"dreiding"`)
- **Exotic chemistries** — check that your elements are covered by the chosen force field's atom types

### Moltemplate execution failed

Inspect the `.lt` files in the output's `moltemplate/` directory — the terminal output names the file and line of the syntax error.

## Installation Issues

### ImportError: no module named 'AutoPoly'

AutoPoly isn't installed (or the wrong environment is active). Install in editable mode from the repo root:

```bash
cd /path/to/AutoPoly
pip install -e .
```

Verify:

```python
import AutoPoly
print(AutoPoly.__version__)
```

### Dependency conflicts

Use a fresh virtual environment:

```bash
python -m venv autopoly_env
source autopoly_env/bin/activate    # Windows: autopoly_env\Scripts\activate
pip install -e /path/to/AutoPoly
```

### RDKit-related import warnings

`System` and `BeadSpringPolymer` work without RDKit, but `Polymer`, `Molecule`, `generate`, and `MonomerGenerator` require it. Install RDKit (`pip install rdkit` or `conda install -c conda-forge rdkit`) for full functionality.

## Performance Issues

### Generation takes too long

- **Start small** — 2–5 chains, DOP 10–20 — then scale up
- **Use bead-spring** for very large systems: `BeadSpringPolymer` bypasses moltemplate entirely (see the [Bead-Spring guide](bead-spring.md))
- **Tune MC placement** — lower `monomer_density` or raise `mc_max_attempts` if placement retries dominate (see [MC Placement](mc-placement.md))

### Memory errors

Reduce `chain_num` or sequence length, respect the 10 000-entry sequence limit, and make sure you run 64-bit Python:

```bash
python -c "import sys; print(sys.maxsize > 2**32)"   # should print True
```

## File Permission Errors

`PermissionError` writing to the output directory — check permissions, pre-create the directory, or write elsewhere:

```python
import os
home = os.path.expanduser("~")
system = System(out=os.path.join(home, "autopoly_output"))
```

## Quick Diagnostics

```python
import AutoPoly, sys, shutil
print(f"AutoPoly version: {AutoPoly.__version__}")
print(f"Python version: {sys.version}")
print(f"Moltemplate found: {shutil.which('moltemplate.sh') is not None}")

try:
    import rdkit
    print("RDKit: OK")
except ImportError:
    print("RDKit: NOT FOUND")
```

## Getting Help

If your issue isn't covered:

1. Check the [API reference](../reference/index.md)
2. Review the [tutorials](../tutorials/index.md) for working code
3. Open an issue on [GitHub](https://github.com/WuGroup-XJTLU/AutoPoly/issues) with the full traceback, minimal reproducing code, AutoPoly version, Python version, and OS
