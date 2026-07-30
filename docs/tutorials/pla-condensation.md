# Condensation Polymer (PLA)

Polylactic acid (PLA) is a biodegradable polyester formed by **step-growth condensation** — chemically different from the vinyl (chain-growth) polymers in the other tutorials.

```
Vinyl addition (PE, PP, PS, PMMA):   all-carbon backbone, no byproducts
Condensation (PLA, Nylon, PET):      heteroatom backbone, small-molecule byproducts

PLA:  n HO-CH(CH3)-COOH  ->  [-O-CH(CH3)-C(=O)-]n  +  n H2O
```

You will learn:

- The complement SMILES for an **ester-backbone** repeat unit
- Why GAFF is the recommended force field for oxygen-rich polymers
- That the AutoPoly workflow is identical for condensation and vinyl polymers — only the monomer definition changes

PLA repeat unit:

| Role | Complement SMILES |
|---|---|
| First | `OC(C)C(=O)[*]` |
| Middle | `[*]OC(C)C(=O)[*]` |
| Last | `[*]OC(C)C(=O)O` |

## The script

```python
--8<-- "examples/example_pla_condensation.py"
```

## Run it

```bash
cd examples
python example_pla_condensation.py
```

!!! tip "Validate against experiment"
    PLA: density ~1.24–1.26 g/cm³, Tg ~330 K. See the [Force Fields guide](../guides/force-fields.md) for benchmarking advice.
