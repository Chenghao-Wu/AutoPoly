# Small Molecules & Mixtures

The `Molecule` class builds boxes of small molecules — solvents, additives, pure liquids — using **regular SMILES** (no `[*]` wildcards; those belong to polymers). Two scripts cover the ground:

- `example_molecules.py` — water, water + ethanol, and PE chains in water, selectable with `--example N`
- `example_benzene_system.py` — the simplest possible case: 100 benzene molecules

You will learn:

- The `Molecule(Count=..., Smiles=..., Name=...)` calling convention (note the capitalized parameters)
- How to build **mixtures** by listing several molecules in `model`
- Why GAFF/GAFF2 is recommended for small organic molecules, especially aromatics

## Molecules and mixtures

```python
--8<-- "examples/example_molecules.py"
```

Run all three cases, or pick one:

```bash
cd examples
python example_molecules.py              # run all examples
python example_molecules.py --example 2  # water + ethanol only
```

## A pure liquid: benzene

```python
--8<-- "examples/example_benzene_system.py"
```

```bash
python example_benzene_system.py
```

## Next

- Combine molecules with a polymer chain — the [Polymer Solution tutorial](peo-solution.md)
- Understand why `Molecule` rejects wildcards — [Polymers vs Molecules](../guides/polymers-vs-molecules.md)
