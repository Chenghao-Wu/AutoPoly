# Bead-Spring (Coarse-Grained)

`BeadSpringPolymer` builds coarse-grained bead-spring polymer models — Kremer–Grest-style chains of Lennard-Jones beads — and writes LAMMPS data files **directly**, with no SMILES, no atom typing, and no moltemplate. Use it for large systems, long chains, and polymer physics studies where chemical detail is not the point.

```python
from AutoPoly import System, BeadSpringPolymer, BeadType

system = System(out="cg_melt")

bsp = BeadSpringPolymer(
    name="homopolymer",
    system=system,
    n_chains=50,
    bead_types=[BeadType(name="A", mass=1.0, epsilon=1.0, sigma=1.0)],
    sequence=[("A", 100)],        # 100 beads of type A per chain
    topology="linear",
    bond_style="harmonic",
    pair_style="lj",
    generation_method="saw",
)
bsp.generate_data_file()
```

## Building blocks

| Class | Purpose |
|---|---|
| `BeadType` | One bead species: `name`, `mass`, LJ `epsilon` and `sigma` |
| `AngleType` | Bending stiffness for a bead triplet (`k`, `theta0`) |
| `SAWConfig` | Tuning for self-avoiding-walk generation (trials, backtracking, ring closure) |
| `MCConfig` | Tuning for Monte Carlo pre-equilibration (moves, temperature, steps) |

## Sequences and topologies

The `sequence` describes one chain and accepts several forms:

```python
sequence=[("A", 100)]              # homopolymer: 100 × A
sequence=[("A", 50), ("B", 50)]    # diblock: A₅₀-b-B₅₀
sequence="AABB"                    # shorthand by bead name
sequence=["A", "A", "B", "B"]      # explicit per-bead list
```

Topologies:

- `"linear"` — chain with two free ends
- `"ring"` — cyclic chain; the generator spends extra trials closing the ring (`SAWConfig.ring_closure_trials`, `ring_closure_tolerance`)

## Potentials

- **Bonds** — `"harmonic"` (E = K(r−r₀)²) or `"fene"` (finitely extensible nonlinear elastic, the standard for bead-spring melts)
- **Pairs** — `"lj"` (full Lennard-Jones, cutoff 2.5σ) or `"wca"` (repulsive Weeks–Chandler–Andersen, cutoff ≈ 1.122σ)
- **Angles** — optional per-triplet stiffness via `AngleType`

## Generation methods

| `generation_method` | Speed | Result |
|---|---|---|
| `"geometric"` | Fastest | Simple geometric placement; may contain overlaps |
| `"saw"` (default) | Fast | Self-avoiding walk — overlap-free coils |
| `"mc"` | Slow | Monte Carlo pre-equilibration (`MCConfig`) — relaxed starting configurations |

The SAW generator shares its collision detector with the atomistic [MC placement engine](mc-placement.md); box size comes from the target bead density (default 0.85 beads/σ³).

## Output

`generate_data_file()` writes a complete LAMMPS data file (masses, pair/bond/angle coefficients, atoms, bonds, angles) under `<System out>/<name>/` — no `system.in.*` split files, since everything is self-contained in reduced units.

## See also

- [Bead-Spring tutorial](../tutorials/bead-spring.md) — homopolymer, diblock with FENE, ring with angles, MC equilibration
- [Bead-Spring API](../reference/bead-spring.md) — full class and dataclass reference
