# First Polymer: PMMA

**The recommended starting point for new AutoPoly users.** This tutorial walks through the complete workflow for generating linear poly(methyl methacrylate) — PMMA — using complement SMILES notation.

You will learn:

- How the **first / middle / last** wildcard variants define a chain
- How `MonomerGenerator` turns them into moltemplate `.lt` templates (including mirror-tacticity `_T1` variants)
- The full **System → Polymer → Polymerization** sequence
- Where the LAMMPS output files land

The PMMA repeat unit needs all three complement SMILES variants:

| Role | Complement SMILES |
|---|---|
| First | `CC(C)(C(=O)OC)[*]` |
| Middle | `[*]CC([*])(C)C(=O)OC` |
| Last | `[*]CC(C)(C(=O)OC)` |

## The script

```python
--8<-- "examples/example_pmma_linear.py"
```

## Run it

```bash
cd examples
python example_pmma_linear.py
```

Output is written to the directory named in `System(out=...)` — look for `system.data` and the `system.in.*` files, as described in [Output Files](../guides/output-files.md).

## Next

- Vary `chain_num`, DOP, and `tacticity` and see how the output changes
- Move on to the [Block Copolymer tutorial](block-copolymer.md) to mix monomers in one sequence
- Read the [Complement SMILES guide](../guides/complement-smiles.md) for the full format rules
