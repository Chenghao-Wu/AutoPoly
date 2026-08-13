# Reactive MD with the Reactor

The **Reactor** prepares LAMMPS [`fix bond/react`](https://docs.lammps.org/fix_bond_react.html)
inputs for an AutoPoly system, so monomers packed into a box can **react during
MD** to form polymers (step-growth polymerization, crosslinking, curing, ...).
It is the AutoPoly counterpart of the
[AutoREACTER](https://github.com/NanoCIPHER-Lab/AutoREACTER) workflow.

## What it does

Given a generated system and its monomers, the reactor:

1. **Detects functional groups** on each monomer (diol, diacid, diamine,
   diisocyanate, hydroxy-acid, amino-acid, ...).
2. **Enumerates reactions** from a built-in reaction-SMARTS library
   (polyesterification, polyamidation, polyurethane formation,
   polyanhydride / polythioester formation, ...).
3. **Builds the reaction templates** LAMMPS needs: a pre-reaction and
   post-reaction *molecule template* plus a *map file* (superimpose file)
   for each unique reaction — typed with the **same force field** as the
   system and carrying embedded 3D coordinates.
4. **Writes supplementary force-field parameters** for any types the reaction
   creates that are not yet in `system.data` (e.g. a new cross-monomer ester
   bond, a rehybridized carbonyl, a water byproduct).
5. **Writes `in.bond_react`** — a ready-to-run reaction-stage input script.

## Quick start

```python
from AutoPoly import System, Molecule, generate, Reactor

system = System(out="reactor_out")
eg     = Molecule(Count=20, Smiles="OCCO", Name="eg")
adipic = Molecule(Count=20, Smiles="O=C(O)CCCCC(=O)O", Name="adipic")

# Build the monomer melt (stage 1-3 pipeline)
generate(system, "melt", [eg, adipic], force_field="gaff2")

# Detect and build the reaction templates
reactor = Reactor(system.get_folder_path() + "/melt", monomers=[eg, adipic])
for inst in reactor.detect_reactions():
    print(inst.description)

result = reactor.build(temperature=300.0, nevery=100, rmin=0.0, rmax=3.5)
print(result.reactions)
```

Then run the reaction stage:

```bash
cd reactor_out/melt
lmp -in in.bond_react
```

## The Reactor API

```python
Reactor(
    project_dir,            # AutoPoly project dir (contains system.data + build/<ff>/)
    monomers=None,          # Molecule/Polymer models or (smiles, name) pairs;
                            #   if omitted, read from geometry/geometry.json
    force_field=None,       # inferred from the single build/<ff>/ dir if omitted
    functional_groups=None, # custom functional-group library (dict)
    reaction_library=None,  # custom reaction library (dict)
)
```

- **`reactor.detect_reactions()`** → list of `ReactionInstance` (each has a
  `.description` string). Stored on `reactor.reaction_instances`.
- **`reactor.build(**kwargs)`** → `ReactorResult` with `reactions`,
  `template_files`, `template_sets`, `assignments`, `warnings`, and the
  `script` path. Build keywords: `temperature`, `nevery`, `rmin`, `rmax`,
  `prob`, `run_steps`, `stabilize_steps`, `seed`, `output_dir`.

## Output files

Written under `<project>/reactor/` (and the script at `<project>/in.bond_react`):

```
melt/
├── in.bond_react                          # reaction-stage LAMMPS input
├── reactor/
│   ├── template_pre_1.molecule            # pre-reaction template
│   ├── template_post_1.molecule           # post-reaction template
│   ├── RXN_1.map                          # superimpose map file
│   ├── RXN_1_with_delete_ids.map          # map incl. DeleteIDs (byproduct)
│   └── system.in.settings.reactor         # coeffs for reaction-created types
```

`in.bond_react` includes `system.in.init/.settings/.charges` plus the
supplementary `system.in.settings.reactor`, reads `system.data` with the
extra per-atom topology slots `fix bond/react` needs, declares the templates,
and issues one `react ...` clause per reaction.

!!! note "Reactions with a byproduct (water, HCl, ...)"
    Condensation reactions split off a small molecule. Use the
    `RXN_*_with_delete_ids.map` variant (referenced in the script comments) so
    LAMMPS deletes those atoms, and prefer NPT to absorb the density change.

## Supported reactions

The built-in library covers the common step-growth condensations/additions:

| Reaction | Reactants | Byproduct |
|----------|-----------|-----------|
| Polyesterification | diol + diacid / diacid-halide; hydroxy-acid (AB) | H₂O / HCl |
| Transesterification | diol + diester | alcohol |
| Polyamidation | diamine + diacid / diacid-halide; amino-acid (AB) | H₂O / HCl |
| Polyanhydride | COOH/acid-halide (AB) | HCl |
| Polythioesterification | dithiol + diacid / diacid-halide | H₂O / HCl |
| Polyurethane (polyaddition) | diol + diisocyanate | — |

Custom reactions can be added by passing a `reaction_library` dict (same
schema as `AutoPoly.reactor.reaction_library.REACTIONS`) and, if new
functional groups are needed, a matching `functional_groups` dict.

## How it works (internals)

- Reaction SMARTS atom-map numbers drive template construction: map numbers
  **1 and 2** mark the two *initiator* atoms between which the new bond forms;
  other mapped atoms form the *first shell*.
- A breadth-first walk (4 bonds) around the first shell defines the *template*
  atoms; the outermost shell becomes the map file's *EdgeIDs*.
- Atom types come from the same SMARTS typer used by the pipeline, applied to
  the combined reactant/product molecules (full chemical environment).
- Numeric types **reuse the ids already in `system.data`** whenever a matching
  interaction exists (by example); only genuinely new interactions get fresh
  ids, with coefficients taken from the force-field tables.
