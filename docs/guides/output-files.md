# Output Files

Every `generate` run writes into `<System out>/<name>/`:

```
my_polymer/polyethylene/
├── moltemplate/           # Intermediate files (.lt inputs, moltemplate output)
├── system.data            # LAMMPS data file (topology & coordinates)
├── system.in.init         # Units, atom/bond/angle styles
├── system.in.settings     # Force field parameters
└── system.in.charges      # Atomic charges
```

## The files

| File | Contents |
|---|---|
| `system.data` | The LAMMPS data file: box dimensions, masses, atoms, bonds, angles, dihedrals — full topology and coordinates |
| `system.in.init` | Simulation header: `units`, `atom_style`, `bond_style`, `angle_style`, `dihedral_style`, `pair_style`, boundary conditions |
| `system.in.settings` | Force field parameters: `pair_coeff`, `bond_coeff`, `angle_coeff`, `dihedral_coeff` entries for every type in the system |
| `system.in.charges` | Atomic partial charges (Gasteiger charges for GAFF/GAFF2 — replace with AM1-BCC/RESP for production) |
| `moltemplate/` | Everything moltemplate produced: the monomer `.lt` templates, the assembled `system.lt`, and moltemplate's raw output |

## Running with LAMMPS

Include the pieces from your own input script:

```bash
lmp -in your_run.in
```

```
# your_run.in
include system.in.init      # styles and units
read_data system.data       # topology and coordinates
include system.in.settings  # force field parameters

# ... your minimization / equilibration / production commands
```

!!! tip "Equilibrate in stages"
    Initial boxes are built at low density for overlap-free placement. Minimize → NVT heating → NPT compression → NPT production, with a 0.5–1.0 fs timestep for atomistic systems.

## Bead-spring output

`BeadSpringPolymer` writes a single self-contained LAMMPS data file (including all coefficients in reduced units) instead of the `system.in.*` split — see the [Bead-Spring guide](bead-spring.md).
