import AutoPoly
from AutoPoly.bead_spring import BeadSpringPolymer

# Define the system
system = AutoPoly.System(out="bead_spring_test")

# Create a linear bead-spring polymer
linear_polymer = BeadSpringPolymer(
    name="linear_polymer",
    system=system,
    n_chains=10,
    n_beads=50,
    topology="linear",
    bond_length=1.0,
    mass=1.0,
    epsilon=1.0,
    sigma=1.0
)
linear_polymer.generate_data_file()

# Create a ring bead-spring polymer
ring_polymer = BeadSpringPolymer(
    name="ring_polymer",
    system=system,
    n_chains=60,
    n_beads=50,
    topology="ring",
    bond_length=1.0,
    mass=1.0,
    epsilon=1.0,
    sigma=1.0
)
ring_polymer.generate_data_file() 