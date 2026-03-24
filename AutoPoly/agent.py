# -*- coding: utf-8 -*-
"""
Config-driven agent API for AutoPoly.

Three core functions: info(), validate(), generate()
Two utilities: describe_smiles(), suggest_force_field()

All functions accept plain dicts and return JSON-serializable results.
"""
import os
from typing import List, Optional

from .conf import (
    EXAMPLE_CONFIGS,
    FORCE_FIELD_DESCRIPTIONS,
    MAX_DOP,
    MAX_SEQUENCE_LENGTH,
    MAX_UNIQUE_MONOMERS,
    TACTICITY_DESCRIPTIONS,
    TOPOLOGY_DESCRIPTIONS,
)
from .results import GenerationResult, ValidationResult


def info() -> dict:
    """Return all available options, limits, and example configs.

    An agent should call this first to discover what AutoPoly can do.
    """
    return {
        "force_fields": FORCE_FIELD_DESCRIPTIONS,
        "topologies": TOPOLOGY_DESCRIPTIONS,
        "tacticities": TACTICITY_DESCRIPTIONS,
        "bead_spring_options": {
            "bond_styles": {
                "harmonic": "Harmonic bond: E = K*(r-r0)^2",
                "fene": "FENE bond: finitely extensible nonlinear elastic",
            },
            "pair_styles": {
                "lj": "Full Lennard-Jones (attractive + repulsive), cutoff 2.5 sigma",
                "wca": "Weeks-Chandler-Andersen (repulsive only), cutoff ~1.122 sigma",
            },
            "generation_methods": {
                "geometric": "Simple geometric placement. Fast but may have overlaps.",
                "saw": "Self-Avoiding Random Walk. Fast and overlap-free.",
                "mc": "Monte Carlo equilibration. Slow but produces equilibrated configs.",
            },
            "topologies": {
                "linear": "Linear chain with two free ends.",
                "ring": "Ring (cyclic) chain, no free ends.",
            },
        },
        "limits": {
            "max_dop": MAX_DOP,
            "max_sequence_length": MAX_SEQUENCE_LENGTH,
            "max_unique_monomers": MAX_UNIQUE_MONOMERS,
        },
        "examples": EXAMPLE_CONFIGS,
    }


def validate(config: dict) -> ValidationResult:
    """Validate a config dict without generating any files.

    Returns ValidationResult with success, errors, warnings, and suggestions.
    Never writes to disk, never raises exceptions.
    """
    try:
        config_type = config.get("type")
        if config_type == "atomistic":
            return _validate_atomistic(config)
        elif config_type == "bead_spring":
            return _validate_bead_spring(config)
        elif config_type is None:
            return ValidationResult(
                success=False,
                errors=["Missing required field 'type'. Must be 'atomistic' or 'bead_spring'."],
                suggestions=["Add \"type\": \"atomistic\" or \"type\": \"bead_spring\" to your config."],
            )
        else:
            return ValidationResult(
                success=False,
                errors=[f"Unknown type '{config_type}'. Must be 'atomistic' or 'bead_spring'."],
                suggestions=[
                    "Use \"type\": \"atomistic\" for all-atom polymer systems with force fields.",
                    "Use \"type\": \"bead_spring\" for coarse-grained bead-spring models.",
                ],
            )
    except Exception as e:
        return ValidationResult(
            success=False,
            errors=[f"Validation error: {e}"],
        )


def generate(config: dict) -> GenerationResult:
    """Validate config, then generate a polymer system.

    Returns GenerationResult with file paths and metadata on success,
    or errors on failure. Never raises exceptions.
    """
    vr = validate(config)
    if not vr.success:
        return GenerationResult(
            success=False,
            errors=vr.errors,
            warnings=vr.warnings,
        )

    try:
        config_type = config["type"]
        if config_type == "atomistic":
            return _generate_atomistic(config)
        else:
            return _generate_bead_spring(config)
    except Exception as e:
        return GenerationResult(
            success=False,
            errors=[f"Generation failed: {e}"],
            warnings=vr.warnings,
        )


def describe_smiles(smiles: str) -> dict:
    """Get structural info about a SMILES string.

    Returns dict with atom count, elements, molecular weight, wildcard info.
    Requires RDKit.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError:
        return {"error": "RDKit is required for describe_smiles. Install with: pip install rdkit"}

    wildcard_count = smiles.count("[*]")
    clean = smiles.replace("[*]", "[H]")

    mol = Chem.MolFromSmiles(clean)
    if mol is None:
        return {"error": f"Invalid SMILES: '{smiles}'", "smiles": smiles}

    mol = Chem.AddHs(mol)
    elements = sorted(set(atom.GetSymbol() for atom in mol.GetAtoms()))

    return {
        "smiles": smiles,
        "atoms": mol.GetNumAtoms(),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "elements": elements,
        "molecular_weight": round(Descriptors.ExactMolWt(mol), 4),
        "has_wildcards": wildcard_count > 0,
        "wildcard_count": wildcard_count,
    }


def suggest_force_field(smiles_list: List[str]) -> List[str]:
    """Suggest compatible force fields based on elements in the SMILES.

    Returns list of force field names, most general first.
    """
    elements = set()
    for smiles in smiles_list:
        info = describe_smiles(smiles)
        if "error" in info:
            return list(FORCE_FIELD_DESCRIPTIONS.keys())
        elements.update(info["elements"])

    elements.discard("H")

    organic_elements = {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I"}
    is_organic = elements.issubset(organic_elements)

    suggestions = []
    if is_organic:
        suggestions.extend(["oplsaa", "gaff2", "gaff", "lopls"])
    suggestions.append("dreiding")
    suggestions.append("compass")

    return suggestions


# ---------------------------------------------------------------------------
# Atomistic validation & generation
# ---------------------------------------------------------------------------

def _validate_atomistic(config: dict) -> ValidationResult:
    errors = []
    warnings = []
    suggestions = []

    if "name" not in config:
        errors.append("Missing required field 'name'.")
        suggestions.append("Add a system name, e.g. \"name\": \"my_polymer\".")

    ff = config.get("force_field", "oplsaa")
    valid_ffs = list(FORCE_FIELD_DESCRIPTIONS.keys())
    if ff not in valid_ffs:
        errors.append(f"Unknown force field '{ff}'.")
        _add_did_you_mean(ff, valid_ffs, suggestions, "force field")

    polymers = config.get("polymers", [])
    molecules = config.get("molecules", [])
    if not polymers and not molecules:
        errors.append("At least one of 'polymers' or 'molecules' is required.")
        suggestions.append(
            "Add a polymer: \"polymers\": [{\"chain_num\": 2, "
            "\"sequence\": [\"CC[*]\", \"[*]CC[*]\", \"[*]CC\"]}]"
        )

    for i, p in enumerate(polymers):
        prefix = f"polymers[{i}]"
        if "chain_num" not in p:
            errors.append(f"{prefix}: Missing 'chain_num'.")
        elif not isinstance(p["chain_num"], int) or p["chain_num"] < 1:
            errors.append(f"{prefix}: 'chain_num' must be a positive integer.")

        seq = p.get("sequence", [])
        if not seq:
            errors.append(f"{prefix}: Missing or empty 'sequence'.")
        else:
            if len(seq) > MAX_SEQUENCE_LENGTH:
                errors.append(
                    f"{prefix}: Sequence length {len(seq)} exceeds limit {MAX_SEQUENCE_LENGTH}."
                )
            _validate_smiles_sequence(seq, prefix, errors, warnings, suggestions)

        topo = p.get("topology", "linear")
        if topo not in TOPOLOGY_DESCRIPTIONS:
            errors.append(f"{prefix}: Unknown topology '{topo}'.")
            _add_did_you_mean(topo, list(TOPOLOGY_DESCRIPTIONS.keys()), suggestions, "topology")

        tact = p.get("tacticity", "atactic")
        if tact not in TACTICITY_DESCRIPTIONS:
            errors.append(f"{prefix}: Unknown tacticity '{tact}'.")
            _add_did_you_mean(tact, list(TACTICITY_DESCRIPTIONS.keys()), suggestions, "tacticity")

    for i, m in enumerate(molecules):
        prefix = f"molecules[{i}]"
        if "smiles" not in m:
            errors.append(f"{prefix}: Missing 'smiles'.")
        elif "[*]" in m.get("smiles", ""):
            errors.append(f"{prefix}: Molecule SMILES must not contain wildcards [*]. Use 'polymers' for polymers.")
        if "count" not in m:
            errors.append(f"{prefix}: Missing 'count'.")
        elif not isinstance(m["count"], int) or m["count"] < 1:
            errors.append(f"{prefix}: 'count' must be a positive integer.")

    placement = config.get("placement_method", "mc_random")
    if placement not in ("grid", "mc_random"):
        errors.append(f"Unknown placement_method '{placement}'. Must be 'grid' or 'mc_random'.")

    return ValidationResult(
        success=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        suggestions=suggestions,
    )


def _validate_smiles_sequence(seq, prefix, errors, warnings, suggestions):
    """Validate complement SMILES wildcard patterns."""
    if len(seq) == 1:
        wc = seq[0].count("[*]")
        if wc == 2:
            pass  # ring-capable single monomer
        elif wc > 2:
            errors.append(f"{prefix}: Single monomer has {wc} wildcards, expected 0 or 2.")
        return

    for j, smi in enumerate(seq):
        wc = smi.count("[*]")
        if j == 0:
            if wc != 1:
                errors.append(f"{prefix}.sequence[{j}]: First monomer should have exactly 1 wildcard (e.g. 'CC[*]'), got {wc}.")
                if wc == 0:
                    suggestions.append(f"Add a wildcard to the right: '{smi}[*]'")
                elif wc == 2:
                    suggestions.append(f"First monomer needs 1 wildcard. Use a middle monomer '{smi}' in positions 1..N-2.")
        elif j == len(seq) - 1:
            if wc != 1:
                errors.append(f"{prefix}.sequence[{j}]: Last monomer should have exactly 1 wildcard (e.g. '[*]CC'), got {wc}.")
                if wc == 0:
                    suggestions.append(f"Add a wildcard to the left: '[*]{smi}'")
                elif wc == 2:
                    suggestions.append(f"Last monomer needs 1 wildcard. Use a middle monomer '{smi}' in positions 1..N-2.")
        else:
            if wc != 2:
                errors.append(f"{prefix}.sequence[{j}]: Middle monomer should have exactly 2 wildcards (e.g. '[*]CC[*]'), got {wc}.")

    # Try to validate with RDKit if available
    try:
        from rdkit import Chem
        for j, smi in enumerate(seq):
            clean = smi.replace("[*]", "[H]")
            if Chem.MolFromSmiles(clean) is None:
                errors.append(f"{prefix}.sequence[{j}]: Invalid SMILES '{smi}'.")
    except ImportError:
        warnings.append("RDKit not available; SMILES syntax not validated.")


def _validate_bead_spring(config: dict) -> ValidationResult:
    from .bead_spring import BeadSpringPolymer

    errors = []
    warnings = []
    suggestions = []

    if "name" not in config:
        errors.append("Missing required field 'name'.")

    if "n_chains" not in config:
        errors.append("Missing required field 'n_chains'.")
    elif not isinstance(config["n_chains"], int) or config["n_chains"] < 1:
        errors.append("'n_chains' must be a positive integer.")

    bead_types = config.get("bead_types")
    if not bead_types:
        errors.append("Missing or empty 'bead_types'. Need at least one, e.g. [{\"name\": \"A\"}].")
    else:
        for i, bt in enumerate(bead_types):
            if "name" not in bt:
                errors.append(f"bead_types[{i}]: Missing 'name'.")

    seq = config.get("sequence")
    if not seq:
        errors.append("Missing or empty 'sequence'. E.g. [[\"A\", 50]] or \"AABB\" or [\"A\",\"A\",\"B\"].")

    topo = config.get("topology", "linear")
    if topo not in BeadSpringPolymer.VALID_TOPOLOGIES:
        errors.append(f"Unknown topology '{topo}'.")
        _add_did_you_mean(topo, BeadSpringPolymer.VALID_TOPOLOGIES, suggestions, "topology")

    bs = config.get("bond_style", "harmonic")
    if bs not in BeadSpringPolymer.VALID_BOND_STYLES:
        errors.append(f"Unknown bond_style '{bs}'.")
        _add_did_you_mean(bs, BeadSpringPolymer.VALID_BOND_STYLES, suggestions, "bond_style")

    ps = config.get("pair_style", "lj")
    if ps not in BeadSpringPolymer.VALID_PAIR_STYLES:
        errors.append(f"Unknown pair_style '{ps}'.")
        _add_did_you_mean(ps, BeadSpringPolymer.VALID_PAIR_STYLES, suggestions, "pair_style")

    gm = config.get("generation_method", "saw")
    if gm not in BeadSpringPolymer.VALID_GENERATION_METHODS:
        errors.append(f"Unknown generation_method '{gm}'.")
        _add_did_you_mean(gm, BeadSpringPolymer.VALID_GENERATION_METHODS, suggestions, "generation_method")

    return ValidationResult(
        success=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        suggestions=suggestions,
    )


# ---------------------------------------------------------------------------
# Generation
# ---------------------------------------------------------------------------

def _generate_atomistic(config: dict) -> GenerationResult:
    from .molecule import Molecule
    from .polymer import Polymer
    from .polymerization import Polymerization
    from .system import System

    name = config["name"]
    output_dir = config.get("output_dir", ".")
    ff = config.get("force_field", "oplsaa")
    placement = config.get("placement_method", "mc_random")
    use_mc = config.get("use_mc_chain_growth", True)

    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    sys = System(out=output_dir)

    models = []
    for p in config.get("polymers", []):
        poly = Polymer(
            chain_num=p["chain_num"],
            sequence=p["sequence"],
            topology=p.get("topology", "linear"),
            tacticity=p.get("tacticity", "atactic"),
        )
        models.append(poly)

    for m in config.get("molecules", []):
        mol = Molecule(
            Count=m["count"],
            Smiles=m["smiles"],
            Name=m.get("name"),
        )
        models.append(mol)

    Polymerization(
        name=name,
        system=sys,
        model=models,
        run=True,
        force_field=ff,
        placement_method=placement,
        use_mc_chain_growth=use_mc,
    )

    result_dir = os.path.join(output_dir, name)
    files = _scan_output_files(result_dir)
    data_file = _find_data_file(files)

    return GenerationResult(
        success=True,
        output_dir=result_dir,
        data_file=data_file,
        files_created=files,
        metadata={
            "force_field": ff,
            "placement_method": placement,
            "n_polymers": len(config.get("polymers", [])),
            "n_molecules": len(config.get("molecules", [])),
        },
    )


def _generate_bead_spring(config: dict) -> GenerationResult:
    from .bead_spring import BeadSpringPolymer, BeadType
    from .system import System

    name = config["name"]
    output_dir = config.get("output_dir", ".")

    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    sys = System(out=output_dir)

    bead_types = [
        BeadType(
            name=bt["name"],
            mass=bt.get("mass", 1.0),
            epsilon=bt.get("epsilon", 1.0),
            sigma=bt.get("sigma", 1.0),
        )
        for bt in config["bead_types"]
    ]

    bsp = BeadSpringPolymer(
        name=name,
        system=sys,
        n_chains=config["n_chains"],
        bead_types=bead_types,
        sequence=_normalize_sequence(config["sequence"]),
        topology=config.get("topology", "linear"),
        bond_style=config.get("bond_style", "harmonic"),
        pair_style=config.get("pair_style", "lj"),
        generation_method=config.get("generation_method", "saw"),
        density=config.get("density"),
        bond_length=config.get("bond_length", 1.0),
        k_bond=config.get("k_bond", 30.0),
    )

    bsp.generate_data_file()

    result_dir = os.path.join(output_dir, name)
    files = _scan_output_files(result_dir)
    data_file = _find_data_file(files)

    sys_info = bsp.get_system_info()

    return GenerationResult(
        success=True,
        output_dir=result_dir,
        data_file=data_file,
        files_created=files,
        metadata={
            "topology": config.get("topology", "linear"),
            "bond_style": config.get("bond_style", "harmonic"),
            "pair_style": config.get("pair_style", "lj"),
            "generation_method": config.get("generation_method", "saw"),
            "n_chains": config["n_chains"],
            "n_beads_per_chain": sys_info.get("beads_per_chain"),
            "total_beads": sys_info.get("total_beads"),
            "box_size": sys_info.get("box_size"),
        },
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _normalize_sequence(seq):
    """Convert JSON lists to tuples for BeadSpringPolymer sequence format."""
    if isinstance(seq, str):
        return seq
    return [tuple(item) if isinstance(item, list) and len(item) == 2
            and isinstance(item[1], (int, float)) else item for item in seq]


def _scan_output_files(directory: str) -> List[str]:
    """Scan directory for created files, return relative paths."""
    if not os.path.isdir(directory):
        return []
    files = []
    for root, _, filenames in os.walk(directory):
        for fn in filenames:
            full = os.path.join(root, fn)
            files.append(os.path.relpath(full, directory))
    return sorted(files)


def _find_data_file(files: List[str]) -> Optional[str]:
    """Find the LAMMPS data file among output files."""
    for f in files:
        if f.endswith(".data") or f.endswith(".lammps"):
            return f
    return None


def _add_did_you_mean(value: str, valid: List[str], suggestions: List[str], label: str):
    """Add a 'did you mean' suggestion using simple substring matching."""
    matches = [v for v in valid if value.lower() in v.lower() or v.lower() in value.lower()]
    if matches:
        suggestions.append(f"Unknown {label} '{value}'. Did you mean '{matches[0]}'? Available: {', '.join(valid)}")
    else:
        suggestions.append(f"Unknown {label} '{value}'. Available: {', '.join(valid)}")
