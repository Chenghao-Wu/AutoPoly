#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bead-Spring Polymers with Side Groups (Comb / Graft Architectures)
===================================================================

This example demonstrates branched bead-spring models built with the
architecture graph core (AutoPoly.models.architectures):

- Example 1: Comb polymer with single-bead side groups, regularly spaced
             (like short-chain-branched polyethylene: backbone beads + one
             branch bead per graft point), with branch-point angles.
- Example 2: Graft copolymer with oligomeric side chains at explicitly
             chosen graft points (e.g. backbone-g-oligomer).
- Example 3: Comb polymer assembled from an explicitly defined monomer
             (MonomerTemplate: backbone bead + side-group bead with
             head/tail/side connection points).
- Example 4: The direct writer backend (backend="direct") — a lightweight
             alternative that writes polymer.data + in.polymer without
             running moltemplate (useful for very large melts).

All examples use the standard moltemplate backend via generate(): .lt
files are emitted and the bundled moltemplate produces
moltemplate/system.data + system.in.init/settings (the standard AutoPoly
output layout). The direct writer is available through
generate(backend="direct") and is demonstrated in Example 4.

Requires: pip install -e .  (from the AutoPoly repo root)
"""

from AutoPoly import (
    System,
    BeadSpringPolymer,
    BeadType,
    AngleType,
    MonomerTemplate,
)
from AutoPoly.models import architectures as arch
from AutoPoly.models.architectures import BeadArchitecture


def example_1_comb_with_side_groups():
    """
    Example 1: Comb polymer with single-bead side groups
    ----------------------------------------------------
    A backbone of 40 A beads with one B side-group bead grafted every 4
    backbone beads (10 graft points) — the coarse-grained analog of a
    short-chain-branched polymer.

    Branch-point angles (centered on graft beads, degree 3) are included
    and configured separately: the backbone-backbone-side triplet
    (A-A-B) is stiffer than the backbone triplets (A-A-A).
    """
    print("\n" + "=" * 60)
    print("Example 1: Comb polymer with single-bead side groups")
    print("=" * 60)

    system = System(out="bead_spring_comb")

    bead_A = BeadType(name="A", mass=1.0, epsilon=1.0, sigma=1.0)  # backbone
    bead_B = BeadType(name="B", mass=1.0, epsilon=1.0, sigma=0.8)  # side group

    comb = arch.comb(
        backbone=[("A", 40)],   # 40 backbone beads
        side="B",               # one side-group bead per graft point
        every=4,                # graft every 4 backbone beads
    )

    polymer = BeadSpringPolymer(
        name="comb",
        system=system,
        n_chains=10,
        bead_types=[bead_A, bead_B],
        architecture=comb,
        bond_style="fene",          # standard for bead-spring melts
        pair_style="wca",           # purely repulsive (Kremer-Grest style)
        use_angles=True,
        default_k_angle=5.0,        # backbone bending (A-A-A)
        angle_types=[
            # Stiffer angle involving the side group at branch points
            AngleType(("A", "A", "B"), k=20.0, theta0=120.0),
        ],
        # Moderate density keeps SAW insertion reliable for branched chains;
        # compress to melt density (~0.85) with NPT during equilibration
        density=0.4,
    )
    polymer.generate()  # moltemplate backend (default)

    info = polymer.get_system_info()
    print(f"  Architecture: {info['architecture']} (branched: {info['is_branched']})")
    print(f"  Chains: {info['n_chains']}")
    print(f"  Beads per chain: {info['n_beads_per_chain']} "
          f"(40 backbone + 10 side groups)")
    print(f"  Total atoms: {info['total_atoms']}")
    print(f"  Total bonds: {info['total_bonds']}")
    print(f"  Total angles: {info['total_angles']} (incl. branch-point triplets)")
    print(f"  Output: {info['output_path']}/moltemplate")

    return polymer


def example_2_graft_copolymer():
    """
    Example 2: Graft copolymer with oligomeric side chains
    ------------------------------------------------------
    A-g-B graft copolymer: an A backbone with B5 side chains at
    explicitly chosen graft points (arch.graft gives full control over
    grafting positions; side chains can be any sequence, e.g. a random
    copolymer).
    """
    print("\n" + "=" * 60)
    print("Example 2: Graft copolymer with oligomeric side chains")
    print("=" * 60)

    system = System(out="bead_spring_graft")

    bead_A = BeadType(name="A", mass=1.0, epsilon=1.0, sigma=1.0)
    bead_B = BeadType(name="B", mass=1.0, epsilon=1.0, sigma=1.0)

    graft_arch = arch.graft(
        backbone=[("A", 30)],
        grafts={
            5:  ("B", 5),       # B5 side chain at backbone bead 5
            15: ("B", 5),       # ... and at bead 15
            25: ("B", 5),       # ... and at bead 25
        },
    )

    polymer = BeadSpringPolymer(
        name="graft",
        system=system,
        n_chains=10,
        bead_types=[bead_A, bead_B],
        architecture=graft_arch,
        bond_style="fene",
        pair_style="wca",
        density=0.4,
    )
    polymer.generate()  # moltemplate backend (default)

    info = polymer.get_system_info()
    print(f"  Architecture: {info['architecture']}")
    print(f"  Chains: {info['n_chains']}")
    print(f"  Beads per chain: {info['n_beads_per_chain']} "
          f"(30 backbone + 3 x 5 side-chain beads)")
    print(f"  Total atoms: {info['total_atoms']}")
    print(f"  Total bonds: {info['total_bonds']}")
    print(f"  Output: {info['output_path']}/moltemplate")

    return polymer


def example_3_explicit_monomer():
    """
    Example 3: Comb built from an explicitly defined monomer
    --------------------------------------------------------
    The same comb topology, but defined at the monomer level: a
    MonomerTemplate with one backbone bead and one side-group bead,
    connected head-to-tail into a chain. This mirrors how monomers are
    defined in the atomistic pipeline and makes the monomer structure
    (backbone vs side group) explicit and reusable.
    """
    print("\n" + "=" * 60)
    print("Example 3: Comb from an explicit MonomerTemplate")
    print("=" * 60)

    system = System(out="bead_spring_monomer_comb")

    bead_A = BeadType(name="A", mass=1.0, epsilon=1.0, sigma=1.0)  # backbone
    bead_B = BeadType(name="B", mass=1.0, epsilon=1.0, sigma=0.8)  # side group

    # A graft monomer: backbone bead 0 (polymerization points head/tail)
    # with a side-group bead 1 attached to it.
    graft_monomer = MonomerTemplate(
        name="G",
        beads=["A", "B"],               # bead 0 = backbone, bead 1 = side group
        internal_bonds=[(0, 1)],        # backbone--side-group bond
        connections={"head": 0, "tail": 0, "side": 1},
    )

    # Polymerize 20 graft monomers head-to-tail
    n_monomers = 20
    comb = BeadArchitecture.from_monomers(
        templates={"G": graft_monomer},
        instances=["G"] * n_monomers,
        inter_bonds=[
            (i, "tail", i + 1, "head") for i in range(n_monomers - 1)
        ],
        name="comb_from_monomer",
    )

    polymer = BeadSpringPolymer(
        name="monomer_comb",
        system=system,
        n_chains=5,
        bead_types=[bead_A, bead_B],
        architecture=comb,
        bond_style="fene",
        pair_style="wca",
        use_angles=True,
        default_k_angle=5.0,
        angle_types=[
            AngleType(("A", "A", "B"), k=20.0, theta0=120.0),
        ],
        density=0.4,
    )
    polymer.generate()  # moltemplate backend (default)

    info = polymer.get_system_info()
    print(f"  Architecture: {info['architecture']} (branched: {info['is_branched']})")
    print(f"  Monomer: G = 1 backbone bead (A) + 1 side-group bead (B)")
    print(f"  Chains: {info['n_chains']}")
    print(f"  Beads per chain: {info['n_beads_per_chain']} "
          f"(20 monomers x 2 beads)")
    print(f"  Total atoms: {info['total_atoms']}")
    print(f"  Total bonds: {info['total_bonds']}")
    print(f"  Total angles: {info['total_angles']}")
    print(f"  Output: {info['output_path']}/moltemplate")

    return polymer


def example_4_direct_backend():
    """
    Example 4: Direct writer backend (lightweight alternative)
    ----------------------------------------------------------
    generate(backend="direct") writes polymer.data + in.polymer directly,
    without the moltemplate build step. Use it for very large melts where
    the moltemplate backend (default) is slow; the physics (coordinates,
    bonds, angles) is identical.
    """
    print("\n" + "=" * 60)
    print("Example 4: Direct writer backend")
    print("=" * 60)

    system = System(out="bead_spring_comb_direct")

    bead_A = BeadType(name="A", mass=1.0, epsilon=1.0, sigma=1.0)
    bead_B = BeadType(name="B", mass=1.0, epsilon=1.0, sigma=0.8)

    comb = arch.comb(backbone=[("A", 20)], side="B", every=4)

    polymer = BeadSpringPolymer(
        name="comb_direct",
        system=system,
        n_chains=5,
        bead_types=[bead_A, bead_B],
        architecture=comb,
        bond_style="fene",
        pair_style="wca",
        use_angles=True,
        angle_types=[AngleType(("A", "A", "B"), k=20.0, theta0=120.0)],
        density=0.4,
    )
    polymer.generate(backend="direct")  # lightweight; no moltemplate run

    info = polymer.get_system_info()
    print(f"  Total atoms: {info['total_atoms']}")
    print(f"  Total bonds: {info['total_bonds']}")
    print(f"  Output: {info['output_path']} (polymer.data + in.polymer)")

    return polymer


def main():
    """Run all side-group bead-spring examples."""
    print("\n" + "#" * 60)
    print("# Bead-Spring Polymers with Side Groups")
    print("#" * 60)

    example_1_comb_with_side_groups()
    example_2_graft_copolymer()
    example_3_explicit_monomer()
    example_4_direct_backend()

    print("\n" + "=" * 60)
    print("All examples completed successfully!")
    print("=" * 60)
    print("\nGenerated output directories:")
    print("  - bead_spring_comb/          (moltemplate backend)")
    print("  - bead_spring_graft/         (moltemplate backend)")
    print("  - bead_spring_monomer_comb/  (moltemplate backend)")
    print("  - bead_spring_comb_direct/   (direct writer)")
    print("\nTo run a LAMMPS simulation (moltemplate backend, examples 1-3):")
    print("  cd <output_dir>/<name>/moltemplate")
    print("  lmp -in in.polymer")
    print("\nTo run a LAMMPS simulation (direct writer, example 4):")
    print("  cd bead_spring_comb_direct/comb_direct")
    print("  lmp -in in.polymer")


if __name__ == "__main__":
    main()
