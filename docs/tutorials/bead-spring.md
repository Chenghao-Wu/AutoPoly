# Bead-Spring Models

Coarse-grained bead-spring polymers, written directly as LAMMPS data files — **no SMILES, no moltemplate**. One script demonstrates the full feature set:

- **Homopolymers and block copolymers** (multiple bead types)
- **Linear and ring topologies**
- **Harmonic and FENE bonds**, plus angle potentials
- **Monte Carlo pre-equilibration** and **SAW generation**
- **Density-based box sizing**

You will learn:

- `BeadType` / `AngleType` / `MCConfig` / `SAWConfig` in action
- The `sequence` formats (`[("A", 50), ("B", 50)]`, `"AABB"`, explicit lists)
- When coarse-grained beats atomistic (large systems, long chains, polymer physics)

## The script

```python
--8<-- "examples/example_bead_spring.py"
```

## Run it

```bash
cd examples
python example_bead_spring.py
```

Each sub-example writes its own self-contained LAMMPS data file. For the concepts behind the knobs, see the [Bead-Spring guide](../guides/bead-spring.md) and the [API reference](../reference/bead-spring.md).
