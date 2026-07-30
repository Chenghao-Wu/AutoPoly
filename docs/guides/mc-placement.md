# MC Placement & Chain Growth

AutoPoly builds the initial configuration of your simulation box with Monte Carlo (MC) methods: chains are **grown** monomer-by-monomer as self-avoiding walks (SAW), and whole molecules are **placed** with collision detection. The result is a coiled, overlap-free starting structure that equilibrates far faster than a grid or an extended chain.

## Placement methods

Choose with `placement_method`:

| Value | Behavior |
|---|---|
| `"mc_random"` (default) | Monte Carlo placement with SAW chain growth — realistic coiled conformations |
| `"grid"` | Deterministic grid placement — extended chains on a lattice; simple, but far from equilibrium |

```python
Polymerization(
    name="polymer_mc",
    system=system,
    model=[polymer],
    force_field="oplsaa",
    placement_method="mc_random",     # Monte Carlo placement (default)
    use_mc_chain_growth=True,         # SAW chain growth (default)
    mc_max_attempts=10000,            # max placement attempts
    mc_monomer_density=0.085,         # target density (monomers / Å³)
    mc_bond_angle_min=50.0,           # min deflection angle (degrees)
    mc_bond_angle_max=90.0,           # max deflection angle (degrees)
)
```

## How SAW chain growth works

1. Each chain starts from a random seed monomer.
2. Every next monomer is attached at a random orientation within the allowed **deflection angle** window (`mc_bond_angle_min` … `mc_bond_angle_max`). Deflection = 180° − bond angle, so a tetrahedral carbon corresponds to ≈ 70.5°.
3. A cell-linked-list **collision detector** rejects placements that overlap existing atoms — both other chains and, with `mc_intrachain_exclude_neighbors` (default 2), the growing chain itself beyond its bonded neighbors.
4. If no valid position is found within `mc_max_attempts`, growth backtracks and retries.

Set `use_mc_chain_growth=False` to place whole chains rigidly instead of growing them.

## Box sizing

Boxes are sized from SAW scaling — `N^0.6 × bond_length` for a chain of `N` monomers — combined with the target `mc_monomer_density`, rather than the fully extended chain length. This produces compact, realistic boxes at low initial density, leaving room for overlap-free placement before NPT compression.

## Tuning guide

| Symptom | Knob |
|---|---|
| Placement fails / "max attempts exceeded" | Lower `mc_monomer_density` (looser box), raise `mc_max_attempts` |
| Chains too extended or too knotted | Adjust `mc_bond_angle_min/max` toward your chemistry's real bond angles |
| False-positive collisions along a chain | Raise `mc_intrachain_exclude_neighbors` (2 is recommended) |
| Want reproducible grid layout | `placement_method="grid"` |

## The mc module

The placement engine is a standalone subpackage — see the [mc API reference](../reference/mc.md) for `CollisionDetector`, `ChainGrowthMC`, and `MolecularPlacementMC`. The same collision machinery drives [bead-spring](bead-spring.md) generation.

## See also

- [MC Placement tutorial](../tutorials/mc-placement.md) — grid vs MC random vs chain growth, side by side
- [Polymerization API](../reference/polymerization.md) — all `mc_*` parameters
