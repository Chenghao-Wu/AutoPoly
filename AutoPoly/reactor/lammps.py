#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LAMMPS Reaction Input Script Writer for the AutoPoly Reactor

Writes a ready-to-run ``in.bond_react`` script that:

1. Includes the AutoPoly ``system.in.init`` / ``system.in.settings`` /
   ``system.in.charges`` (so styles and force-field parameters come from
   the generated system),
2. reads ``system.data`` with the extra per-atom topology slots
   ``fix bond/react`` needs to create bonds,
3. declares the reaction templates (``molecule``) and the ``fix bond/react``
   command with one ``react`` clause per reaction,
4. thermostats the non-reacting group (``<stab_group>_REACT``), and
5. runs the reaction stage with reaction-count thermo output.

Supplementary coefficients for parameter types the reaction introduces
(new cross-monomer bonds, changed atom types, ...) are written to
``system.in.settings.reactor`` and included after the main settings.

Created on 2026-08-04
@author: zwu
"""
from pathlib import Path
from typing import Dict, List, Optional

from ..core.system import logger
from .templates import TemplateSet, TypeAssignment
def write_settings_reactor(
    path: Path,
    assignments: List[TypeAssignment],
    masses_path: Optional[Path] = None,
) -> bool:
    """
    Write supplementary coeff lines for reactor-introduced types.

    New atom types also need ``mass`` commands, which must run BEFORE
    ``read_data``; those are written to ``masses_path`` (default:
    ``system.in.masses.reactor`` next to ``path``) when any new masses exist.

    Returns True when at least one coeff line was written.
    """
    coeff_lines: List[str] = []
    mass_lines: List[str] = []
    seen = set()
    for assignment in assignments:
        for line in assignment.new_coeff_lines:
            if line not in seen:
                seen.add(line)
                coeff_lines.append(line)
        for line in assignment.new_mass_lines:
            if line not in seen:
                seen.add(line)
                mass_lines.append(line)

    if mass_lines:
        if masses_path is None:
            masses_path = Path(path).with_name("system.in.masses.reactor")
        # Convert "id  mass  # name" entries into LAMMPS "mass" commands.
        mass_cmds = []
        for line in mass_lines:
            parts = line.split("#")[0].split()
            if len(parts) >= 2:
                comment = line.split("#", 1)[1].strip() if "#" in line else ""
                suffix = f"  # {comment}" if comment else ""
                mass_cmds.append(f"mass {parts[0]} {parts[1]}{suffix}")
        header = [
            "# Masses for atom types introduced by AutoPoly reactor",
            "# (must be included BEFORE read_data)",
            "",
        ]
        Path(masses_path).write_text("\n".join(header + mass_cmds) + "\n")

    if not coeff_lines:
        return bool(mass_lines)

    lines = [
        "# Supplementary parameters introduced by AutoPoly reactor",
        "# (types created by reactions; not present in the initial system.data)",
        "",
    ]
    lines.extend(coeff_lines)
    Path(path).write_text("\n".join(lines) + "\n")
    return True


def write_bond_react_script(
    path: Path,
    template_sets: List[TemplateSet],
    *,
    temperature: float = 300.0,
    nevery: int = 100,
    rmin: float = 0.0,
    rmax: float = 3.5,
    prob: float = 1.0,
    stabilize_steps: int = 60,
    run_steps: int = 1_000_000,
    thermo_every: int = 100,
    dump_every: int = 1000,
    data_file: str = "system.data",
    init_file: str = "system.in.init",
    settings_file: str = "system.in.settings",
    charges_file: Optional[str] = "system.in.charges",
    reactor_settings_file: Optional[str] = "system.in.settings.reactor",
    reactor_masses_file: Optional[str] = "system.in.masses.reactor",
    box_bounds: Optional[Dict[str, tuple]] = None,
    max_new_types: Optional[Dict[str, int]] = None,
    seed: int = 12345,
    extra_per_atom: int = 50,
) -> Path:
    """
    Write the ``in.bond_react`` LAMMPS input script.

    Args:
        path: Output script path.
        template_sets: Built template triplets (one per reaction).
        temperature: Reaction-stage temperature (K).
        nevery: Reaction attempt frequency (timesteps).
        rmin / rmax: Initiator distance cutoffs (Angstrom).
        prob: Reaction probability once constraints are met.
        stabilize_steps: nve/limit stabilization length after each reaction.
        run_steps: Reaction-stage run length.
        thermo_every / dump_every: Output frequencies.
        data_file / init_file / settings_file / charges_file: file names
            (relative to the run directory) produced by AutoPoly stage 3.
        reactor_settings_file: supplementary settings (None to skip include).
        seed: velocity RNG seed.
        extra_per_atom: extra topology slots per atom for fix bond/react.

    Returns:
        The script path.
    """
    path = Path(path)
    # When the reaction introduces atom types absent from system.data,
    # mass/pair_coeff/*_coeff need the simulation box first, so use the
    # create_box + read_data add/merge pattern; otherwise a plain read_data.
    has_new_types = reactor_masses_file is not None
    lines: List[str] = [
        "# AutoPoly reactor: fix bond/react reaction stage",
        "# Run from the directory containing system.data / system.in.*",
        "",
        "#------------Init (styles)------------",
        f"include {init_file}",
    ]

    if has_new_types:
        if box_bounds is None or max_new_types is None:
            raise ValueError(
                "box_bounds and max_new_types are required when the reaction "
                "introduces new atom types (reactor_masses_file is set)"
            )
        x = box_bounds["x"]; y = box_bounds["y"]; z = box_bounds["z"]
        lines += [
            "",
            "#----- Create box with slots for reaction-created types -----",
            f"region simbox block {x[0]} {x[1]} {y[0]} {y[1]} {z[0]} {z[1]} units box",
            f"create_box {max_new_types.get('atom', 0)} simbox &",
            f"    bond/types {max_new_types.get('bond', 0)} &",
            f"    angle/types {max_new_types.get('angle', 0)} &",
            f"    dihedral/types {max_new_types.get('dihedral', 0)} &",
            f"    improper/types {max_new_types.get('improper', 0)} &",
            f"    extra/bond/per/atom {extra_per_atom} &",
            f"    extra/angle/per/atom {extra_per_atom} &",
            f"    extra/dihedral/per/atom {extra_per_atom} &",
            f"    extra/improper/per/atom {extra_per_atom} &",
            f"    extra/special/per/atom {extra_per_atom}",
            "",
            "# Masses for the new atom types (box now exists)",
            f"include {reactor_masses_file}",
            "",
            "#----- Append the real coordinates -----",
            f"read_data {data_file} add append offset 0 0 0 0 0",
        ]
    else:
        lines += [
            "",
            "#------------Read box (defines simulation box)------------",
            f"read_data {data_file} &",
            f"    extra/bond/per/atom {extra_per_atom} &",
            f"    extra/angle/per/atom {extra_per_atom} &",
            f"    extra/dihedral/per/atom {extra_per_atom} &",
            f"    extra/improper/per/atom {extra_per_atom} &",
            f"    extra/special/per/atom {extra_per_atom}",
        ]

    lines += [
        "",
        "#------------Force-field settings------------",
        f"include {settings_file}",
    ]
    if charges_file:
        lines.append(f"include {charges_file}")
    if reactor_settings_file:
        lines.append(f"include {reactor_settings_file}")
    lines += [
        "",
        "#------------Minimization and velocity------------",
        "# kspace_modify gewald: needed for (near-)uncharged systems so PPPM runs",
        "kspace_modify gewald 0.2",
        "minimize 1.0e-4 1.0e-6 1000 10000",
        f"velocity all create {temperature} {seed} dist gaussian",
        "timestep 1.0",
        f"thermo {thermo_every}",
        "reset_timestep 0",
        "",
        "#------------Reaction templates------------",
    ]

    react_clauses: List[str] = []
    for i, ts in enumerate(template_sets, start=1):
        pre_id = f"mol_pre_{i}"
        post_id = f"mol_post_{i}"
        lines.append(f"molecule {pre_id} {ts.pre_file.name}")
        lines.append(f"molecule {post_id} {ts.post_file.name}")
        if ts.map_file_delete_ids is not None:
            lines.append(
                f"# NOTE: reaction {i} deletes byproduct atoms; use "
                f"{ts.map_file_delete_ids.name} (with NPT) instead of "
                f"{ts.map_file.name} if atom deletion is desired"
            )
        react_clauses.append(
            f"react rxn_{i} all {nevery} {rmin} {rmax} {pre_id} {post_id} "
            f"{ts.map_file.name} prob {prob} {seed} rescale_charges yes"
        )
        lines.append("")

    all_reacts = " &\n                ".join(react_clauses)
    lines += [
        "#------------fix bond/react------------",
        "fix rxns all bond/react stabilization yes statted_grp 0.03 &",
        f"                {all_reacts}",
        "",
        "# Thermostat the NON-reacting group (statted_grp_REACT).",
        "# Reacting atoms are integrated by the fix's internal nve/limit.",
        f"fix nvt_react statted_grp_REACT nvt temp {temperature} {temperature} 100.0",
        "",
        "thermo_style custom step time temp f_rxns[*] press density vol pe ke etotal",
        f"dump traj all xyz {dump_every} reacted.xyz",
        "dump_modify traj types numeric",
        "",
        f"run {run_steps}",
        "write_data reacted.data nofix",
    ]

    path.write_text("\n".join(lines) + "\n")
    logger.info(f"Wrote reaction input script: {path}")
    return path
