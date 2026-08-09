"""
Integration tests: run the bundled moltemplate on bead-spring .lt files
and validate the generated LAMMPS system.data.
"""

import re

import pytest
import tempfile

from AutoPoly.models.bead_spring import BeadSpringPolymer, BeadType
from AutoPoly.models.bead_spring_system import BeadSpringSystem
from AutoPoly.models import architectures as arch


class MockSystem:
    def __init__(self, path: str):
        self._path = path

    def get_folder_path(self) -> str:
        return self._path


@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def mock_system(temp_dir):
    return MockSystem(temp_dir)


A = BeadType("A", mass=1.0, epsilon=1.0, sigma=1.0)
B = BeadType("B", mass=2.0, epsilon=0.5, sigma=1.0)


def read_system_data_counts(path):
    with open(path) as f:
        txt = f.read()
    counts = {}
    for key in ("atoms", "bonds", "angles", "atom types",
                "bond types", "angle types"):
        m = re.search(rf"(\d+)\s+{key}\b", txt)
        counts[key] = int(m.group(1)) if m else 0
    return counts


@pytest.mark.integration
class TestMoltemplateExecution:
    def test_comb_polymer_end_to_end(self, mock_system, temp_dir):
        comb = arch.comb(backbone=[("A", 8)], side="B", every=3)
        bsp = BeadSpringPolymer(
            name="comb", system=mock_system, n_chains=2,
            bead_types=[A, B], architecture=comb,
            use_angles=True, density=0.1,
        )
        mtd = bsp.generate_moltemplate(run_moltemplate=True)

        assert (mtd / "system.data").exists()
        assert (mtd / "system.in.init").exists()
        assert (mtd / "system.in.settings").exists()

        counts = read_system_data_counts(mtd / "system.data")
        # comb: 8 backbone + 3 grafts = 11 beads, 10 bonds per chain
        assert counts["atoms"] == 2 * 11
        assert counts["bonds"] == 2 * 10
        assert counts["atom types"] == 2
        assert counts["bond types"] == 1
        assert counts["angles"] == 2 * 11

    def test_mixture_end_to_end(self, mock_system, temp_dir):
        bss = BeadSpringSystem(
            name="mix", system=mock_system, bead_types=[A, B],
            use_angles=True, density=0.1,
            bond_style="fene", pair_style="wca",
        )
        bss.add_species(arch.ring([("A", 8)]), 2)
        bss.add_species(arch.star(center="A", arms=[("B", 4)] * 3), 1)
        mtd = bss.generate_moltemplate(run_moltemplate=True)

        counts = read_system_data_counts(mtd / "system.data")
        # rings: 8 beads/8 bonds; star: 13 beads/12 bonds
        assert counts["atoms"] == 2 * 8 + 13
        assert counts["bonds"] == 2 * 8 + 12
        assert counts["angles"] == 2 * 8 + 12
        assert counts["atom types"] == 2

        with open(mtd / "system.in.init") as f:
            init = f.read()
        assert "units lj" in init
        assert "bond_style fene" in init
        assert "special_bonds fene" in init
        assert "pair_style lj/cut 1.12246" in init

        with open(mtd / "system.in.settings") as f:
            settings = f.read()
        assert "bond_coeff 1 30.0 1.5" in settings
        assert "angle_coeff" in settings

    def test_moltemplate_data_positions_match_saw(self, mock_system, temp_dir):
        """Coordinates in system.data should match the generated SAW coords."""
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A], sequence=[("A", 6)], density=0.1,
        )
        mtd = bsp.generate_moltemplate(run_moltemplate=True)
        with open(mtd / "system.data") as f:
            txt = f.read()
        atoms_section = txt.split("Atoms")[1].split("Bonds")[0]
        coords = []
        for line in atoms_section.splitlines():
            tokens = line.split()
            if len(tokens) >= 7 and tokens[0].isdigit():
                coords.append([float(tokens[4]), float(tokens[5]), float(tokens[6])])
        assert len(coords) == 6
        for k in range(6):
            saw_pos = bsp._positions[k]
            for dim in range(3):
                assert coords[k][dim] == pytest.approx(saw_pos[dim], abs=1e-3)
