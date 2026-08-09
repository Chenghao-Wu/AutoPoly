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

## Architectures (graph core)

For anything beyond linear/ring, pass a `BeadArchitecture` built with the
factories in `AutoPoly.models.architectures` (a chain is an explicit graph:
nodes = beads, edges = bonds — everything else is derived):

```python
from AutoPoly.models import architectures as arch

bsp = BeadSpringPolymer(
    name="comb",
    system=system,
    n_chains=20,
    bead_types=[BeadType("A"), BeadType("B")],
    architecture=arch.comb(backbone=[("A", 100)], side=("B", 5), every=10),
    use_angles=True,
)
```

| Factory | Architecture |
|---|---|
| `arch.linear(seq)` | linear chain (same as `topology="linear"`) |
| `arch.ring(seq)` | cyclic chain |
| `arch.star(center, arms)` | star / miktoarm: center bead + f arms (arms may differ) |
| `arch.comb(backbone, side, every, offset=0)` | regular comb: side chains every `every` backbone beads |
| `arch.graft(backbone, {idx: side_seq, ...})` | graft polymer with explicit graft points |
| `arch.tadpole(ring_seq, tail, attach=0)` | ring + linear tail (lariat) |
| `arch.dendrimer(core, branch, branch_factor, generations)` | regularly branched dendrimer |
| `arch.custom(bead_types, bonds)` | arbitrary graph escape hatch |

Every sequence slot (backbone, arm, side chain, tail) accepts the same
formats as `sequence` — strings, explicit lists, or block tuples. Explicit
monomer-level structure is available through `MonomerTemplate` (multi-bead
monomers with named connection points `head`/`tail`/`side`) and
`BeadArchitecture.from_monomers(...)`.

**Copolymer sequence generators** (usable in any sequence slot):

```python
from AutoPoly import block_sequence, alternating_sequence, random_sequence, gradient_sequence

arch.linear(block_sequence([("A", 50), ("B", 50)]))     # diblock
arch.linear(alternating_sequence(["A", "B"], 100))      # ABAB...
arch.linear(random_sequence(["A", "B"], 100, weights=[0.7, 0.3], seed=42))
arch.linear(gradient_sequence("A", "B", 100, seed=1))   # A→B composition gradient
```

**Branch angles** are configurable: angle triplets centered on branch
points (graft points, star centers — beads with degree > 2) are included
by default and get their own canonical triplets (e.g. side–backbone–side
vs. backbone–backbone–backbone), so they can carry separate `AngleType`
stiffness. Pass `include_branch_angles=False` to exclude them.

**Branched MC equilibration** works out of the box: tree-pivot moves
(rotate the subtree beyond any bridge bond) and segment crankshaft moves
(on linear segments) replace the path-based pivot/crankshaft; reptation
is linear-only and is skipped for branched chains.

## Mixtures (multi-species melts)

`BeadSpringSystem` packs several species — any mix of architectures —
into one box and one data file:

```python
from AutoPoly import BeadSpringSystem, BeadType
from AutoPoly.models import architectures as arch

bss = BeadSpringSystem(
    name="ring_linear_blend",
    system=system,
    bead_types=[BeadType("A")],
    bond_style="fene",
    pair_style="wca",
    density=0.85,
)
bss.add_species(arch.ring([("A", 50)]), n_chains=20)
bss.add_species(arch.linear([("A", 50)]), n_chains=20)
bss.generate_data_file()
```

All species share one bead-type table; chains of all species are placed
collision-free with the graph-based SAW generator.

## Moltemplate backend

`generate_moltemplate()` emits moltemplate `.lt` files mirroring the
atomistic pipeline (`bead_spring.lt` force field, `bead_<Type>.lt`
monomer objects, `chains.lt` with one object per chain at the generated
coordinates, `system.lt`) and runs the bundled moltemplate to produce
`system.data` + `system.in.init/settings`:

```python
bsp.generate_moltemplate(run_moltemplate=True)   # -> <out>/<name>/moltemplate/
bss.generate_moltemplate(run_moltemplate=True)   # mixtures too
```

Use it when you want the standard moltemplate output layout or need to
mix coarse-grained chains with other moltemplate objects.

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
