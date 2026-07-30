# MC Placement Methods

Three ways to fill the simulation box, compared on the same system (PEO, 5 chains of 20 monomers):

1. **Grid placement** — deterministic lattice arrangement; simple, but chains start extended and far from equilibrium
2. **MC random placement** — random positions and orientations with collision detection
3. **MC chain growth** — self-avoiding walk: chains grown monomer-by-monomer into coiled, realistic conformations

You will learn:

- The `placement_method` and `use_mc_chain_growth` switches
- The `mc_*` tuning parameters (attempts, density, bond angles)
- Why SAW-grown coils equilibrate faster than grid-placed extended chains

## The script

```python
--8<-- "examples/example_peo_mc_placement.py"
```

## Run it

```bash
cd examples
python example_peo_mc_placement.py
```

## What to compare

- **Box size** — MC runs size the box from SAW scaling (`N^0.6`), giving more compact boxes than grid placement
- **Initial conformation** — visualize the three `system.data` files: extended lattice vs random coils
- **Time to equilibrate** — grid systems typically need much longer NVT/NPT equilibration

For the theory and tuning guide, see [MC Placement & Chain Growth](../guides/mc-placement.md).
