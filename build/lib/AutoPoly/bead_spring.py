"""
Bead-Spring Polymer Model Module

This module provides the BeadSpringPolymer class for generating simplified
bead-spring polymer models for coarse-grained molecular dynamics simulations.

The BeadSpringPolymer class handles:
- Generation of bead-spring polymer structures
- Support for linear and ring topologies
- LAMMPS data file generation
- Input script creation for molecular dynamics simulations

This module is useful for:
- Quick polymer structure generation
- Coarse-grained simulations
- Educational purposes
- Testing polymer simulation workflows
"""
import os
import sys
from pathlib import Path
import numpy as np
from typing import Optional
from .logger import setup_logger

logger = setup_logger()

class BeadSpringPolymer:
    """
    Bead-spring polymer model generator for LAMMPS simulations.
    
    This class generates simplified bead-spring polymer models where each
    monomer is represented as a single bead connected by harmonic springs.
    It supports both linear and ring polymer topologies and generates
    complete LAMMPS input files for molecular dynamics simulations.
    
    Attributes:
        name (str): Name for the output files
        system (object): System object containing path information
        path (str): Full path to output directory
        n_chains (int): Number of polymer chains
        n_beads (int): Number of beads per chain
        topology (str): Polymer topology ('linear' or 'ring')
        bond_length (float): Equilibrium bond length
        mass (float): Mass of each bead
        epsilon (float): LJ energy parameter
        sigma (float): LJ distance parameter
    """
    
    def __init__(self, name: str = None, system: object = None, 
                 n_chains: int = 1, n_beads: int = 10, 
                 topology: str = "linear", bond_length: float = 1.0,
                 mass: float = 1.0, epsilon: float = 1.0, sigma: float = 1.0) -> None:
        """
        Initialize bead-spring polymer generator.
        
        Args:
            name (str, optional): Name for the output files. Defaults to None.
            system (object, optional): System object containing path information.
                                     Defaults to None.
            n_chains (int, optional): Number of polymer chains. Defaults to 1.
            n_beads (int, optional): Number of beads per chain. Defaults to 10.
            topology (str, optional): "linear" or "ring". Defaults to "linear".
            bond_length (float, optional): Equilibrium bond length. Defaults to 1.0.
            mass (float, optional): Mass of each bead. Defaults to 1.0.
            epsilon (float, optional): LJ energy parameter. Defaults to 1.0.
            sigma (float, optional): LJ distance parameter. Defaults to 1.0.
        
        Raises:
            ValueError: If topology is not 'linear' or 'ring'
        """
        if topology not in ["linear", "ring"]:
            raise ValueError("Topology must be either 'linear' or 'ring'")
            
        self.name = name
        self.system = system
        self.path = f"{self.system.get_FolderPath}/{self.name}" if system else f"./{name}"
        self.n_chains = n_chains
        self.n_beads = n_beads
        self.topology = topology
        self.bond_length = bond_length
        self.mass = mass
        self.epsilon = epsilon
        self.sigma = sigma
        
        # Create output directory
        Path(self.path).mkdir(parents=True, exist_ok=True)
        logger.info(f"Initialized bead-spring polymer: {n_chains} chains, {n_beads} beads each, {topology} topology")
        
    def generate_data_file(self) -> None:
        """
        Generate LAMMPS data file for bead-spring polymer.
        
        This method creates a complete LAMMPS data file containing:
        - System information (number of atoms, bonds, atom types, bond types)
        - Box dimensions
        - Atom masses
        - Atom coordinates
        - Bond connectivity
        
        The method handles different topologies:
        - Linear: Chains arranged along x-axis with spacing
        - Ring: Chains arranged in 3D grid with random orientations
        """
        n_bonds = self.n_beads - 1 if self.topology == "linear" else self.n_beads
        total_beads = self.n_chains * self.n_beads
        total_bonds = self.n_chains * n_bonds
        
        # Calculate spacing for ring polymers
        radius = self.bond_length * self.n_beads / (2 * np.pi)
        if self.topology == "ring":
            # Increase box size to accommodate rings with spacing
            spacing = radius * 3  # Use 3x radius for good separation
            # Calculate grid arrangement in 3D
            n_per_dim = int(np.ceil(np.cbrt(self.n_chains)))  # Cubic root for 3D arrangement
            box_size = max(n_per_dim * spacing * 2, 50.0)
        else:
            box_size = max(self.n_beads * self.bond_length * 2, 50.0)
        
        with open(f"{self.path}/polymer.data", 'w') as f:
            # Header
            f.write("LAMMPS Bead-Spring Polymer Data File\n\n")
            f.write(f"{total_beads} atoms\n")
            f.write(f"{total_bonds} bonds\n\n")
            f.write("1 atom types\n")
            f.write("1 bond types\n\n")
            
            # Box dimensions
            f.write(f"{-box_size/2:.1f} {box_size/2:.1f} xlo xhi\n")
            f.write(f"{-box_size/2:.1f} {box_size/2:.1f} ylo yhi\n")
            f.write(f"{-box_size/2:.1f} {box_size/2:.1f} zlo zhi\n\n")
            
            # Masses
            f.write("Masses\n\n")
            f.write(f"1 {self.mass:.3f}\n\n")
            
            # Atoms section: atom-ID molecule-ID atom-type x y z
            f.write("Atoms\n\n")
            atom_id = 1
            for chain in range(self.n_chains):
                if self.topology == "ring":
                    # Calculate position in 3D grid
                    ix = chain % n_per_dim
                    iy = (chain // n_per_dim) % n_per_dim
                    iz = chain // (n_per_dim * n_per_dim)
                    
                    # Center of the current ring
                    center_x = (ix - n_per_dim/2 + 0.5) * spacing
                    center_y = (iy - n_per_dim/2 + 0.5) * spacing
                    center_z = (iz - n_per_dim/2 + 0.5) * spacing
                    
                    # Randomly choose rotation angles for the ring
                    theta = np.random.uniform(0, np.pi)  # Rotation around x-axis
                    phi = np.random.uniform(0, 2*np.pi)  # Rotation around z-axis
                    
                    # Create rotation matrices
                    Rx = np.array([[1, 0, 0],
                                 [0, np.cos(theta), -np.sin(theta)],
                                 [0, np.sin(theta), np.cos(theta)]])
                    Rz = np.array([[np.cos(phi), -np.sin(phi), 0],
                                 [np.sin(phi), np.cos(phi), 0],
                                 [0, 0, 1]])
                    R = Rz @ Rx  # Combined rotation matrix
                    
                    # Place beads in a rotated circle
                    for bead in range(self.n_beads):
                        angle = 2 * np.pi * bead / self.n_beads
                        # Initial position in xy-plane
                        pos = np.array([
                            radius * np.cos(angle),
                            radius * np.sin(angle),
                            0.0
                        ])
                        # Apply rotation and translation
                        pos = R @ pos + np.array([center_x, center_y, center_z])
                        # Format: atom-ID molecule-ID atom-type x y z
                        f.write(f"{atom_id} {chain+1} 1 {pos[0]:.3f} {pos[1]:.3f} {pos[2]:.3f}\n")
                        atom_id += 1
                else:
                    # Linear chain along x-axis
                    for bead in range(self.n_beads):
                        x = bead * self.bond_length
                        y = chain * self.bond_length * 2  # Space chains apart
                        z = 0.0
                        # Format: atom-ID molecule-ID atom-type x y z
                        f.write(f"{atom_id} {chain+1} 1 {x:.3f} {y:.3f} {z:.3f}\n")
                        atom_id += 1
            
            # Bonds
            f.write("\nBonds\n\n")
            bond_id = 1
            for chain in range(self.n_chains):
                start_id = chain * self.n_beads + 1
                for bead in range(self.n_beads - 1):
                    f.write(f"{bond_id} 1 {start_id + bead} {start_id + bead + 1}\n")
                    bond_id += 1
                # Add closing bond for ring topology
                if self.topology == "ring":
                    f.write(f"{bond_id} 1 {start_id + self.n_beads - 1} {start_id}\n")
                    bond_id += 1
        
        # Generate LAMMPS input script
        self._generate_input_script()
        logger.info(f"Generated bead-spring polymer files in {self.path}")
    
    def _generate_input_script(self) -> None:
        """
        Generate LAMMPS input script for the bead-spring polymer.
        
        This method creates a complete LAMMPS input script with:
        - Units and atom style settings
        - Force field parameters (LJ + harmonic bonds)
        - Energy minimization
        - Molecular dynamics simulation settings
        - Output and analysis commands
        """
        with open(f"{self.path}/in.polymer", 'w') as f:
            f.write("# LAMMPS input script for bead-spring polymer\n\n")
            f.write("units           lj\n")
            f.write("atom_style      molecular\n")
            f.write("boundary        p p p\n\n")
            
            f.write("read_data      polymer.data\n\n")
            
            # Pair and bond potentials
            f.write(f"pair_style     lj/cut {2.5*self.sigma}\n")
            f.write(f"pair_coeff     1 1 {self.epsilon} {self.sigma}\n")
            f.write("bond_style     harmonic\n")
            f.write(f"bond_coeff     1 30.0 {self.bond_length}\n\n")
            
            # Basic simulation settings
            f.write("neighbor        0.3 bin\n")
            f.write("neigh_modify   every 1 delay 0 check yes\n\n")
            
            # Minimization
            f.write("minimize       1.0e-4 1.0e-6 1000 10000\n")
            f.write("reset_timestep 0\n\n")
            
            # MD settings
            f.write("timestep       0.005\n")
            f.write("fix           1 all nvt temp 1.0 1.0 0.5\n")
            f.write("dump          1 all custom 1000 dump.lammpstrj id type x y z\n")
            f.write("thermo        1000\n")
            f.write("run           100000\n")
    
    def get_system_info(self) -> dict:
        """
        Get comprehensive information about the bead-spring polymer system.
        
        Returns:
            dict: Dictionary containing system properties
        """
        n_bonds = self.n_beads - 1 if self.topology == "linear" else self.n_beads
        total_beads = self.n_chains * self.n_beads
        total_bonds = self.n_chains * n_bonds
        
        return {
            'name': self.name,
            'n_chains': self.n_chains,
            'n_beads_per_chain': self.n_beads,
            'topology': self.topology,
            'total_beads': total_beads,
            'total_bonds': total_bonds,
            'bond_length': self.bond_length,
            'mass': self.mass,
            'epsilon': self.epsilon,
            'sigma': self.sigma,
            'output_path': self.path
        }
    
    def modify_parameters(self, **kwargs) -> None:
        """
        Modify polymer parameters after initialization.
        
        Args:
            **kwargs: Keyword arguments for parameters to modify
                     (n_chains, n_beads, topology, bond_length, mass, epsilon, sigma)
        """
        valid_params = ['n_chains', 'n_beads', 'topology', 'bond_length', 
                       'mass', 'epsilon', 'sigma']
        
        for param, value in kwargs.items():
            if param in valid_params:
                setattr(self, param, value)
                logger.info(f"Modified {param} to {value}")
            else:
                logger.warning(f"Unknown parameter: {param}")
        
        if 'topology' in kwargs:
            if kwargs['topology'] not in ["linear", "ring"]:
                raise ValueError("Topology must be either 'linear' or 'ring'") 