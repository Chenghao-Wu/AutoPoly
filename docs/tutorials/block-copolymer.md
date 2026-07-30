# Block Copolymer (PE-PS-PE)

Block copolymers are where complement SMILES shines: the `sequence` list gives you **explicit control over the monomer at every position**, so any block arrangement is just a matter of listing monomers in order.

This example builds a **PE-PS-PE ABA triblock**: polyethylene end blocks flanking a polystyrene mid-block.

You will learn:

- How to write a mixed-monomer sequence with clean block boundaries
- Why internal positions always need **two** wildcards, even across a block boundary
- How block lengths follow directly from list repetition counts

## The script

```python
--8<-- "examples/example_block_copolymer.py"
```

## Run it

```bash
cd examples
python example_block_copolymer.py
```

## Variations to try

- Change the block length ratio (e.g. symmetric 10-20-10 vs asymmetric 5-30-5)
- Make a **diblock** (PE-PS) by ending the sequence with a PS `last` variant
- Make an **alternating** copolymer by interleaving single monomers — see [Common Patterns](../guides/complement-smiles.md#common-patterns)
