# Polymer Solution (PEO + water)

This tutorial mixes the `Polymer` and `Molecule` classes in one system: poly(ethylene oxide) chains solvated in **200 explicit water molecules**. PEO/water is a classic biocompatible polymer solution (hydrogels, drug delivery, battery electrolytes).

You will learn:

- How to pass polymers **and** molecules together in `generate(system, name, [...])`
- The two SMILES conventions side by side: wildcards for the polymer, plain SMILES for the solvent
- Why GAFF is a convenient single force field for mixed organic systems

## The script

```python
--8<-- "examples/example_peo_solution.py"
```

## Run it

```bash
cd examples
python example_peo_solution.py
```

!!! note "Charges and water models"
    Gasteiger charges are assigned automatically. For production runs, consider AM1-BCC/RESP charges and a water-specific model if quantitative aqueous properties matter — see [Force Fields](../guides/force-fields.md#treat-gasteiger-charges-as-a-starting-point).

## Variations to try

- Change the water count to tune concentration
- Swap water for ethanol (`Smiles="CCO"`) or a benzene/THF mixture
- Scale up: more chains, longer PEO — but equilibrate carefully (low initial density, staged NPT)
