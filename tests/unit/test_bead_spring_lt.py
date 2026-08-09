"""
Unit tests for the moltemplate backend (.lt emission) of bead-spring models.
These tests only check the generated .lt files; tests that execute
moltemplate live in tests/integration/test_bead_spring_moltemplate.py.
"""

import pytest
import tempfile

from AutoPoly.models.bead_spring import BeadSpringPolymer, BeadType, AngleType
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
B = BeadType("B", mass=2.0, epsilon=0.5, sigma=1.2)


class TestForceFieldLt:
    def test_masses_and_styles(self, mock_system, temp_dir):
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A, B], sequence=[("A", 5)],
            bond_style="fene", pair_style="wca",
            use_angles=True, density=0.1,
        )
        mtd = bsp.generate_moltemplate(run_moltemplate=False)
        content = (mtd / "bead_spring.lt").read_text()
        assert "@atom:A 1.000" in content
        assert "@atom:B 2.000" in content
        assert "units lj" in content
        assert "atom_style full" in content
        assert "bond_style fene" in content
        assert "special_bonds fene" in content
        # WCA cutoff = 2^(1/6) * max_sigma = 1.12246 * 1.2
        assert f"pair_style lj/cut {2 ** (1/6) * 1.2:.5f}" in content
        assert "pair_modify shift yes" in content
        # FENE coeff: K R0 eps sigma (from first bead type)
        assert "bond_coeff @bond:bs1 30.0 1.5000 1.0000 1.0000" in content
        assert "angle_style harmonic" in content

    def test_pair_coeffs_use_mixed_types(self, mock_system, temp_dir):
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A, B], sequence=[("A", 5)],
            density=0.1,
        )
        mtd = bsp.generate_moltemplate(run_moltemplate=False)
        content = (mtd / "bead_spring.lt").read_text()
        # Lorentz-Berthelot: eps_AB = sqrt(1*0.5), sigma_AB = (1+1.2)/2
        assert "pair_coeff @atom:A @atom:B 0.7071 1.1000" in content
        assert "pair_coeff @atom:B @atom:B 0.5000 1.2000" in content

    def test_angle_coeffs_per_triplet(self, mock_system, temp_dir):
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A, B],
            architecture=arch.comb(backbone=[("A", 6)], side="B", every=2),
            use_angles=True,
            angle_types=[AngleType(("A", "A", "B"), k=55.0, theta0=109.5)],
            density=0.1,
        )
        mtd = bsp.generate_moltemplate(run_moltemplate=False)
        content = (mtd / "bead_spring.lt").read_text()
        assert "angle_coeff @angle:at" in content
        assert "55.0000 109.5  # A-A-B" in content


class TestBeadAndChainsLt:
    def test_bead_lt_objects(self, mock_system, temp_dir):
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A, B], sequence=[("A", 3)], density=0.1,
        )
        mtd = bsp.generate_moltemplate(run_moltemplate=False)
        bead_a = (mtd / "bead_A.lt").read_text()
        assert "bead_A inherits BeadSpringFF" in bead_a
        assert "@atom:A" in bead_a
        assert (mtd / "bead_B.lt").exists()

    def test_chains_lt_graph_bonds(self, mock_system, temp_dir):
        graft_arch = arch.graft(["A", "A", "A"], {1: "B"})
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A, B], architecture=graft_arch,
            generation_method="geometric",
        )
        mtd = bsp.generate_moltemplate(run_moltemplate=False)
        content = (mtd / "chains.lt").read_text()
        assert "chain_1 inherits BeadSpringFF" in content
        assert "$atom:bead[1]/b $atom:bead[3]/b" in content  # graft bond
        assert "@bond:bs1" in content
        # 4 beads instantiated
        assert content.count("= new bead_") == 4

    def test_chains_lt_typed_angles(self, mock_system, temp_dir):
        star = arch.star(center="A", arms=[("B", 3)] * 3)
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A, B], architecture=star,
            use_angles=True, generation_method="geometric",
        )
        mtd = bsp.generate_moltemplate(run_moltemplate=False)
        content = (mtd / "chains.lt").read_text()
        assert 'write("Data Angles")' in content
        assert "@angle:at" in content
        # 3-arm star: C(3,2)=3 branch triplets at the center
        # + 2 arm triplets per arm = 3 + 6 = 9
        assert content.count("$angle:a") == 9

    def test_no_angles_omits_angle_sections(self, mock_system, temp_dir):
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A], sequence=[("A", 4)], density=0.1,
        )
        mtd = bsp.generate_moltemplate(run_moltemplate=False)
        chains = (mtd / "chains.lt").read_text()
        ff = (mtd / "bead_spring.lt").read_text()
        assert "Data Angles" not in chains
        assert "angle_style" not in ff
        assert "angle_coeff" not in ff


class TestSystemLt:
    def test_system_lt_instances_and_boundary(self, mock_system, temp_dir):
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=3,
            bead_types=[A], sequence=[("A", 4)],
            box_size=12.0,
        )
        mtd = bsp.generate_moltemplate(run_moltemplate=False)
        content = (mtd / "system.lt").read_text()
        for k in range(3):
            assert f"chains[{k}] = new chain_{k + 1}" in content
        assert "-6.0000  6.0000  xlo xhi" in content

    def test_mixture_system_lt(self, mock_system, temp_dir):
        bss = BeadSpringSystem(
            name="m", system=mock_system, bead_types=[A, B], density=0.1,
        )
        bss.add_species(arch.ring([("A", 6)]), 2)
        bss.add_species(arch.star(center="A", arms=[("B", 3)] * 2), 1)
        mtd = bss.generate_moltemplate(run_moltemplate=False)
        content = (mtd / "system.lt").read_text()
        assert content.count("= new chain_") == 3
        chains = (mtd / "chains.lt").read_text()
        assert "(ring)" in chains
        assert "(star)" in chains

    def test_input_script_written(self, mock_system, temp_dir):
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A], sequence=[("A", 4)], density=0.1,
        )
        mtd = bsp.generate_moltemplate(run_moltemplate=False)
        content = (mtd / "in.polymer").read_text()
        assert "read_data       system.data" in content
        assert "include         system.in.init" in content
        assert "include         system.in.settings" in content
