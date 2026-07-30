# Force Fields

AutoPoly supports six force fields. Pass one as the `force_field` string to `Polymerization`:

| Force Field | Value | Best For | Limitations |
|------------|-------|----------|-------------|
| **OPLS-AA** | `"oplsaa"` | General organic polymers, biomolecules | May lack parameters for exotic chemistries |
| **LOPLS** | `"lopls"` | Liquid-phase, better densities | Less extensive than OPLS-AA |
| **GAFF** | `"gaff"` | Small molecules, drug-like compounds | Charges need care (see below) |
| **GAFF2** | `"gaff2"` | Updated GAFF, improved parameters | Newer, less extensively tested |
| **DREIDING** | `"dreiding"` | Generic systems, metals, inorganics | Generic parameters, less accurate for organics |
| **COMPASS** | `"compass"` | Commercial polymers, condensed phases | Class II — needs the LAMMPS CLASS2 package |

## The force fields in detail

### OPLS-AA — `"oplsaa"`

The safe default for organic systems. Extensively parameterized and validated, with wide coverage of organic functional groups and a large literature base. Best for commodity polymers (PE, PP, PS), organic liquids, and biomolecules.

*Citation: Jorgensen et al., J. Am. Chem. Soc. 1996, 118, 11225–11236.*

### LOPLS — `"lopls"`

OPLS re-optimized for the liquid phase: better densities and enthalpies of vaporization for condensed systems. A good choice for polymer melts and solutions where density accuracy matters.

*Citation: Siu et al., J. Chem. Theory Comput. 2012, 8, 1459–1470.*

### GAFF — `"gaff"`

The general AMBER force field — comprehensive coverage of organic and drug-like molecules, and the standard choice for small-molecule solvents and mixed systems. AutoPoly assigns Gasteiger charges automatically; for production runs replace them with AM1-BCC or RESP charges in `system.in.charges`.

*Citation: Wang et al., J. Comput. Chem. 2004, 25, 1157–1174.*

### GAFF2 — `"gaff2"`

GAFF with improved bond, angle, and torsional parameters. Prefer it over GAFF for new projects when no legacy compatibility is needed.

### DREIDING — `"dreiding"`

A generic force field with near-universal element coverage, including metals and inorganics. Ideal for exploratory structure generation and systems lacking specific parameters — but too approximate for quantitative production properties of organics.

*Citation: Mayo et al., J. Phys. Chem. 1990, 94, 8897–8909.*

### COMPASS — `"compass"`

A Class II force field (with cross-terms) validated extensively against commercial polymers — excellent densities and mechanical properties for polyolefins, polyesters, and polyamides. Requires LAMMPS built with the CLASS2 package.

*Citation: Sun, J. Phys. Chem. B 1998, 102, 7338–7364.*

## Choosing: a decision tree

```
System type?
│
├─ Polymer only
│  ├─ Commodity (PE, PP, PS) → OPLS-AA or LOPLS
│  ├─ Commercial/industrial  → COMPASS
│  └─ Novel/exploratory      → OPLS-AA, or DREIDING for a first look
│
├─ Small molecules only
│  ├─ Organic solvents  → GAFF or GAFF2
│  └─ Simple molecules  → OPLS-AA or GAFF
│
├─ Mixed (polymer + solvent)
│  ├─ Organic solvent → GAFF / GAFF2 (handles both)
│  └─ Aqueous         → GAFF (with TIP3P-style water) or OPLS-AA
│
└─ Special systems
   ├─ Metals, MOFs, inorganics → DREIDING
   └─ High-accuracy polymer properties → COMPASS
```

## Common pitfalls

### Don't mix force fields across components

One `Polymerization` call uses one force field for everything in `model`. That is the correct behavior — never combine, say, an OPLS-AA polymer with a GAFF solvent in one box without validated combination rules:

```python
# GOOD: one force field covers all components
Polymerization(
    name="compatible_system",
    system=system,
    model=[polymer, solvent],
    force_field="gaff",
)
```

### Treat Gasteiger charges as a starting point

AutoPoly assigns Gasteiger (partial equalization) charges automatically for GAFF/GAFF2. They are fine for setup and testing, but for production runs compute AM1-BCC charges with Antechamber (or RESP with your QM package) and edit them into `system.in.charges`.

### Don't use DREIDING for final numbers

DREIDING's generic parameters are for exploration. Generate and equilibrate the structure with DREIDING if nothing else covers your chemistry, then rebuild with a specific force field for production.

## Validate before production

1. **Compare with literature** for similar systems.
2. **Benchmark basic properties** (density, energy) against experiment — if off by more than ~10%, reconsider the force field.
3. **Start small** (10–20 chains, DOP ≈ 50) before scaling up.

## See also

- [Force Field Comparison tutorial](../tutorials/force-field-comparison.md) — the same polymer built with all six
- [Gasteiger charges example](../tutorials/index.md) — automatic charge assignment with GAFF
- [Polymerization API](../reference/polymerization.md) — the `force_field` parameter
