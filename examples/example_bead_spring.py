#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bead-Spring Coarse-Grained Polymer Examples
============================================

This example demonstrates the BeadSpringPolymer class for generating
coarse-grained bead-spring polymer models for LAMMPS simulations.

Features demonstrated:
- Homopolymers and block copolymers
- Linear and ring topologies
- Harmonic and FENE bonds
- Angle potentials
- Monte Carlo pre-equilibration
- Density-based box sizing

Output files (for each example):
- polymer.data: LAMMPS data file
- in.polymer: LAMMPS input script

Requires: pip install -e .  (from the AutoPoly repo root)
"""

from AutoPoly import System, BeadSpringPolymer, BeadType, AngleType, MCConfig, SAWConfig


def example_1_homopolymer():
    """
    Example 1: Simple Homopolymer with Harmonic Bonds
    -------------------------------------------------
    Creates a basic homopolymer system with:
    - Single bead type (A)
    - Linear topology
    - Harmonic bonds
    - Density-based box sizing
    """
    print("\n" + "="*60)
    print("Example 1: Simple Homopolymer")
    print("="*60)

    # Create system for output files
    system = System(out="bead_spring_homopolymer")

    # Define a single bead type with LJ parameters
    bead_A = BeadType(
        name="A",
        mass=1.0,       # Reduced mass
        epsilon=1.0,    # LJ energy parameter
        sigma=1.0       # LJ length parameter
    )

    # Create homopolymer: 10 chains, 50 beads each
    polymer = BeadSpringPolymer(
        name="homopolymer",
        system=system,
        n_chains=10,
        bead_types=[bead_A],
        sequence=[("A", 50)],       # 50 A beads per chain
        topology="linear",
        bond_style="harmonic",
        bond_length=1.0,
        k_bond=30.0,                # Bond spring constant
        density=0.85,               # Target density (beads/sigma^3)
    )

    # Generate LAMMPS files
    polymer.generate_data_file()

    # Print system information
    info = polymer.get_system_info()
    print(f"  Chains: {info['n_chains']}")
    print(f"  Beads per chain: {info['n_beads_per_chain']}")
    print(f"  Total atoms: {info['total_atoms']}")
    print(f"  Total bonds: {info['total_bonds']}")
    print(f"  Bond style: {info['bond_style']}")
    print(f"  Output: {info['output_path']}")

    return polymer


def example_2_diblock_fene():
    """
    Example 2: Diblock Copolymer with FENE Bonds
    --------------------------------------------
    Creates a diblock copolymer system with:
    - Two bead types (A, B) with different LJ parameters
    - FENE bond style (finitely extensible nonlinear elastic)
    - Different block lengths (25 A beads + 25 B beads)
    """
    print("\n" + "="*60)
    print("Example 2: Diblock Copolymer with FENE Bonds")
    print("="*60)

    system = System(out="bead_spring_diblock")

    # Define two bead types with different interactions
    bead_A = BeadType(name="A", mass=1.0, epsilon=1.0, sigma=1.0)
    bead_B = BeadType(name="B", mass=1.0, epsilon=1.2, sigma=1.1)  # Slightly different

    # Create diblock copolymer: A25-B25
    polymer = BeadSpringPolymer(
        name="diblock",
        system=system,
        n_chains=8,
        bead_types=[bead_A, bead_B],
        sequence=[("A", 25), ("B", 25)],  # Block pattern: 25 A's then 25 B's
        topology="linear",
        bond_style="fene",
        bond_length=0.97,           # Equilibrium bond length
        k_bond=30.0,                # FENE spring constant
        fene_r0=1.5,                # FENE maximum extension
        density=0.85,
    )

    polymer.generate_data_file()

    info = polymer.get_system_info()
    print(f"  Chains: {info['n_chains']}")
    print(f"  Beads per chain: {info['n_beads_per_chain']}")
    print(f"  Bead types: {info['bead_types']}")
    print(f"  Sequence pattern: A25-B25")
    print(f"  Bond style: {info['bond_style']}")
    print(f"  Output: {info['output_path']}")

    return polymer


def example_3_ring_with_angles():
    """
    Example 3: Ring Polymer with Angle Potentials
    ----------------------------------------------
    Creates a ring polymer system with:
    - Ring topology (closed chains)
    - Angle potentials for chain stiffness
    - Custom angle parameters for different triplets
    """
    print("\n" + "="*60)
    print("Example 3: Ring Polymer with Angle Potentials")
    print("="*60)

    system = System(out="bead_spring_ring")

    # Define two bead types for a ring copolymer
    bead_A = BeadType(name="A", mass=1.0, epsilon=1.0, sigma=1.0)
    bead_B = BeadType(name="B", mass=1.0, epsilon=1.0, sigma=1.0)

    # Define custom angle parameters for specific triplets
    angle_types = [
        AngleType(triplet=("A", "A", "A"), k=20.0, theta0=180.0),  # Stiff A-A-A
        AngleType(triplet=("A", "A", "B"), k=15.0, theta0=180.0),  # Medium at junctions
        AngleType(triplet=("A", "B", "B"), k=10.0, theta0=180.0),  # Flexible B region
    ]

    # Create ring polymer with alternating pattern
    # Using string sequence for explicit pattern
    polymer = BeadSpringPolymer(
        name="ring",
        system=system,
        n_chains=5,
        bead_types=[bead_A, bead_B],
        sequence="AAAABBBBAAAABBBB",    # Explicit string sequence (16 beads)
        topology="ring",                 # Closed ring
        bond_style="harmonic",
        bond_length=1.0,
        k_bond=30.0,
        use_angles=True,                # Enable angle potentials
        default_k_angle=10.0,           # Default angle stiffness
        default_theta0=180.0,           # Default equilibrium angle
        angle_types=angle_types,        # Custom angle parameters
        density=0.5,                    # Lower density for rings
    )

    polymer.generate_data_file()

    info = polymer.get_system_info()
    print(f"  Chains: {info['n_chains']}")
    print(f"  Beads per chain: {info['n_beads_per_chain']}")
    print(f"  Topology: {info['topology']}")
    print(f"  Total angles: {info['total_angles']}")
    print(f"  Angle stiffness: custom per triplet")
    print(f"  Output: {info['output_path']}")

    return polymer


def example_4_equilibrated_melt():
    """
    Example 4: Multi-Chain Melt with MC Equilibration
    -------------------------------------------------
    Creates a polymer melt with Monte Carlo pre-equilibration:
    - Multiple chains at melt density
    - MC equilibration to remove overlaps
    - Custom MC configuration
    - Reports acceptance rate statistics

    MC equilibration uses several move types:
    - Single bead displacement
    - Crankshaft rotation
    - Pivot rotation
    - Reptation (slithering snake)
    - Chain translation
    - Chain rotation
    """
    print("\n" + "="*60)
    print("Example 4: Multi-Chain Melt with MC Equilibration")
    print("="*60)

    system = System(out="bead_spring_melt")

    bead_A = BeadType(name="A", mass=1.0, epsilon=1.0, sigma=1.0)

    # Configure MC equilibration
    # Note: For production runs, use n_steps=10000 or more.
    # Here we use fewer steps for demonstration.
    mc_config = MCConfig(
        density=0.85,               # Melt density
        n_steps=1000,               # Number of MC steps (increase for production)
        temperature=1.0,            # Reduced temperature for Metropolis
        max_displacement=0.5,       # Max single bead displacement
        max_angle=0.3,              # Max rotation angle (radians)
        lj_cutoff=2.5,              # LJ cutoff in sigma units
        bond_k=100.0,               # Bond spring constant for energy
    )

    print("  MC Configuration:")
    print(f"    Steps: {mc_config.n_steps}")
    print(f"    Temperature: {mc_config.temperature}")
    print(f"    Max displacement: {mc_config.max_displacement}")
    print(f"    Max rotation angle: {mc_config.max_angle} rad")

    # Create multi-chain melt with equilibration
    polymer = BeadSpringPolymer(
        name="melt",
        system=system,
        n_chains=50,
        bead_types=[bead_A],
        sequence=[("A", 100)],       # 100-bead chains
        topology="linear",
        bond_style="harmonic",
        bond_length=1.0,
        k_bond=30.0,
        density=0.85,
        equilibrate=True,           # Enable MC equilibration
        mc_config=mc_config,        # Custom MC parameters
    )

    print("\n  Running MC equilibration...")
    polymer.generate_data_file()

    info = polymer.get_system_info()
    print(f"\n  Final system:")
    print(f"    Chains: {info['n_chains']}")
    print(f"    Beads per chain: {info['n_beads_per_chain']}")
    print(f"    Total atoms: {info['total_atoms']}")
    print(f"    Output: {info['output_path']}")

    return polymer


def example_5_saw_generation():
    """
    Example 5: Fast SAW-Based Configuration Generation
    ---------------------------------------------------
    Creates a polymer system using Self-Avoiding Random Walk:
    - Generates overlap-free configurations directly
    - Much faster than MC equilibration (100-1000x for large systems)
    - Suitable for moderate densities (< 0.6 beads/sigma^3)

    SAW algorithm:
    1. Grows chains bead-by-bead
    2. Each new bead placed at bond length distance
    3. Positions filtered by angle and collision constraints
    4. Backtracking when stuck
    """
    print("\n" + "="*60)
    print("Example 5: Fast SAW-Based Configuration Generation")
    print("="*60)

    system = System(out="bead_spring_saw")

    bead_A = BeadType(name="A", mass=1.0, epsilon=1.0, sigma=1.0)

    # Configure SAW algorithm
    saw_config = SAWConfig(
        collision_sigma=1.0,         # Bead diameter for collision detection
        collision_tolerance=0.1,     # Slight overlap tolerance
        n_trials=50,                 # Trial positions per bead
        max_backtrack_depth=10,      # Max beads to remove when stuck
        max_total_backtracks=1000,   # Total backtrack budget per chain
        bond_angle_min=60.0,         # Min bond angle (degrees)
        bond_angle_max=180.0,        # Max bond angle
    )

    print("  SAW Configuration:")
    print(f"    Collision sigma: {saw_config.collision_sigma}")
    print(f"    Trials per bead: {saw_config.n_trials}")
    print(f"    Angle range: {saw_config.bond_angle_min}° - {saw_config.bond_angle_max}°")
    print(f"    Max backtracks: {saw_config.max_total_backtracks}")

    # Create polymer with SAW generation
    polymer = BeadSpringPolymer(
        name="saw_polymer",
        system=system,
        n_chains=20,
        bead_types=[bead_A],
        sequence=[("A", 50)],        # 50-bead chains
        topology="linear",
        bond_style="harmonic",
        bond_length=1.0,
        k_bond=30.0,
        density=0.3,                 # Moderate density for reliable SAW
        generation_method="saw",     # Use SAW instead of MC or geometric
        saw_config=saw_config,
    )

    print("\n  Generating configurations with SAW...")
    polymer.generate_data_file()

    info = polymer.get_system_info()
    print(f"\n  Generated system:")
    print(f"    Chains: {info['n_chains']}")
    print(f"    Beads per chain: {info['n_beads_per_chain']}")
    print(f"    Total atoms: {info['total_atoms']}")
    print(f"    Generation method: SAW (Self-Avoiding Random Walk)")
    print(f"    Output: {info['output_path']}")

    return polymer


def main():
    """Run all bead-spring polymer examples."""
    print("\n" + "#"*60)
    print("# Bead-Spring Coarse-Grained Polymer Examples")
    print("#"*60)

    # Run each example.
    # Example 4 (MC equilibration) is opt-in: it is much slower than the
    # others. Uncomment the call below to run it.
    example_1_homopolymer()
    example_2_diblock_fene()
    example_3_ring_with_angles()
    # example_4_equilibrated_melt()
    example_5_saw_generation()

    print("\n" + "="*60)
    print("All examples completed successfully!")
    print("="*60)
    print("\nGenerated output directories:")
    print("  - bead_spring_homopolymer/")
    print("  - bead_spring_diblock/")
    print("  - bead_spring_ring/")
    print("  - bead_spring_saw/")
    print("  (example 4 would also create bead_spring_melt/)")
    print("\nTo run a LAMMPS simulation:")
    print("  cd <output_dir>/<name>")
    print("  lmp -in in.polymer")


if __name__ == "__main__":
    main()
