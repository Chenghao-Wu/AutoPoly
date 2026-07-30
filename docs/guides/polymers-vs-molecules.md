# Polymers vs Molecules

AutoPoly has two model classes. Picking the right one is the first decision in any script.

| Feature | `Polymer` | `Molecule` |
|---------|-----------|------------|
| **Use for** | Polymer chains of any length | Solvents, additives, any small molecule |
| **SMILES** | Complement SMILES with `[*]` wildcards | Regular SMILES — **no** wildcards |
| **Defined by** | `chain_num`, `sequence`, `topology`, `tacticity` | `Count`, `Smiles`, `Name` |
| **Chain bonds** | Bonds created between sequence neighbors | Independent molecules, no inter-molecule bonds |
| **Example** | `Polymer(chain_num=10, sequence=["CC[*]"]+["[*]CC[*]"]*48+["[*]CC"])` | `Molecule(Count=100, Smiles="O", Name="water")` |

!!! warning "Wildcard rule"
    `[*]` wildcards belong to `Polymer` sequences only. A `Molecule` SMILES containing `[*]` is rejected — this is the most common mix-up between the two classes.

## When to use Polymer

Anything with a backbone built from repeated units — homopolymers, block copolymers, alternating copolymers, ring polymers:

```python
from AutoPoly import Polymer

polymer = Polymer(
    chain_num=10,                    # chains in the box
    sequence=["CC[*]"] + ["[*]CC[*]"] * 48 + ["[*]CC"],
    topology="linear",               # or "ring"
    tacticity="atactic",             # or "isotactic", "syndiotactic"
)
```

Key points:

- **DOP is implicit** — `len(sequence)`, never a separate parameter
- **Sequence order is chain order** — mix monomers for copolymers
- **Ring topology** uses only middle variants — see [Complement SMILES](complement-smiles.md#ring-vs-linear-polymers)

## When to use Molecule

Anything that stays a single, unbonded molecule in the box:

```python
from AutoPoly import Molecule

water   = Molecule(Count=100, Smiles="O",   Name="water")
ethanol = Molecule(Count=20,  Smiles="CCO", Name="ethanol")
```

Note the capitalized parameter names (`Count`, `Smiles`, `Name`) — `Molecule` keeps this older calling convention, unlike the snake_case `Polymer`.

## Mixed systems

Pass both classes together in the `models` list to build polymer solutions, melts with additives, or multi-solvent mixtures:

```python
from AutoPoly import System, Polymer, Molecule, generate

system = System(out="solution")

polymer = Polymer(
    chain_num=5,
    sequence=["COC[*]"] + ["[*]COC[*]"] * 48 + ["[*]COC"],   # PEO
)
water = Molecule(Count=200, Smiles="O", Name="water")

generate(
    system,
    "peo_solution",
    [polymer, water],     # polymers and molecules mix freely
    force_field="gaff",
)
```

For force field choice in mixed systems, see the [Force Fields guide](force-fields.md) — GAFF/GAFF2 usually handle both polymers and small organic solvents well. A full worked example is the [Polymer Solution tutorial](../tutorials/peo-solution.md).

## Related

- [Molecule API](../reference/molecule.md) and [Polymer API](../reference/polymer.md) references
- [Small Molecules & Mixtures tutorial](../tutorials/molecules.md)
