"""
Unit tests for graph-architecture bead-spring polymers and mixtures:
BeadSpringPolymer with architecture=..., branched MC moves, and
BeadSpringSystem mixtures.
"""

import numpy as np
import pytest
import tempfile
from pathlib import Path

from AutoPoly.models.bead_spring import (
    BeadSpringPolymer, BeadType, AngleType, MCConfig, SAWConfig,
    mc_tree_pivot_move, mc_segment_crankshaft_move, build_chain_graph,
)
from AutoPoly.models.bead_spring_system import BeadSpringSystem
from AutoPoly.models import architectures as arch


class MockSystem:
    """Mock System object for testing."""

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
B = BeadType("B", mass=1.0, epsilon=1.0, sigma=1.0)
C = BeadType("C", mass=2.0, epsilon=0.5, sigma=1.2)


def read_data_file(path):
    """Parse a LAMMPS data file into header counts and sections."""
    with open(path) as f:
        lines = f.read().splitlines()
    header = {}
    for line in lines:
        for key in ("atoms", "bonds", "angles", "atom types",
                    "bond types", "angle types"):
            if line.strip().endswith(key):
                header[key] = int(line.split()[0])
    return header, lines


# =============================================================================
# Construction with architecture=
# =============================================================================

class TestArchitectureConstruction:
    def test_star_construction(self, mock_system):
        star = arch.star(center="A", arms=[("B", 5)] * 3)
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=2,
            bead_types=[A, B], architecture=star,
        )
        assert bsp.n_beads == 16
        assert bsp.topology == "star"
        assert bsp._sequence == star.bead_types

    def test_architecture_and_sequence_mutually_exclusive(self, mock_system):
        with pytest.raises(ValueError, match="not both"):
            BeadSpringPolymer(
                name="t", system=mock_system, n_chains=1,
                bead_types=[A], sequence=[("A", 5)],
                architecture=arch.linear([("A", 5)]),
            )

    def test_either_sequence_or_architecture_required(self, mock_system):
        with pytest.raises(ValueError, match="required"):
            BeadSpringPolymer(
                name="t", system=mock_system, n_chains=1, bead_types=[A],
            )

    def test_unknown_bead_type_in_architecture_raises(self, mock_system):
        star = arch.star(center="A", arms=[("Z", 5)] * 2)
        with pytest.raises(ValueError, match="Unknown bead type"):
            BeadSpringPolymer(
                name="t", system=mock_system, n_chains=1,
                bead_types=[A], architecture=star,
            )

    def test_legacy_ring_equivalent_to_arch_ring(self, mock_system):
        legacy = BeadSpringPolymer(
            name="t1", system=mock_system, n_chains=1,
            bead_types=[A], sequence=[("A", 10)], topology="ring",
        )
        via_arch = BeadSpringPolymer(
            name="t2", system=mock_system, n_chains=1,
            bead_types=[A], architecture=arch.ring([("A", 10)]),
        )
        assert legacy._bonds_local == via_arch._bonds_local
        assert legacy._sequence == via_arch._sequence
        assert legacy.topology == via_arch.topology == "ring"


# =============================================================================
# Branch angles
# =============================================================================

class TestBranchAngles:
    def _comb_polymer(self, mock_system, include_branch):
        comb = arch.comb(backbone=[("A", 10)], side="B", every=2)
        return BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A, B], architecture=comb,
            use_angles=True, include_branch_angles=include_branch,
        )

    def test_branch_angles_included_by_default(self, mock_system):
        bsp = self._comb_polymer(mock_system, True)
        # 15 beads, 14 bonds. Grafts at backbone 0,2,4,6,8; graft bead 0 has
        # degree 2 (backbone end), so 4 branch points (degree 3) contribute
        # C(3,2)=3 triplets each = 12 branch triplets.
        # Non-branch centers: degree-2 beads 0,1,3,5,7 -> 5 triplets.
        n_all = len(bsp._angle_triplets_local)
        assert n_all == 5 + 12

    def test_branch_angles_excluded(self, mock_system):
        bsp = self._comb_polymer(mock_system, False)
        assert len(bsp._angle_triplets_local) == 5

    def test_branch_angle_types_configurable(self, mock_system):
        comb = arch.comb(backbone=[("A", 6)], side="B", every=2)
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A, B], architecture=comb,
            use_angles=True,
            angle_types=[AngleType(("A", "A", "B"), k=99.0, theta0=120.0)],
        )
        k, theta0 = bsp._get_angle_params(("A", "A", "B"))
        assert k == pytest.approx(99.0)
        assert theta0 == pytest.approx(120.0)
        # A-B-A (two side beads cannot meet on one backbone bead here,
        # but A-A-A exists)
        assert ("A", "A", "A") in bsp._angle_type_map
        assert ("A", "A", "B") in bsp._angle_type_map


# =============================================================================
# Data file generation for architectures
# =============================================================================

class TestArchitectureDataFiles:
    def test_star_data_file_counts(self, mock_system, temp_dir):
        star = arch.star(center="A", arms=[("B", 8)] * 4)
        bsp = BeadSpringPolymer(
            name="star", system=mock_system, n_chains=3,
            bead_types=[A, B], architecture=star,
            use_angles=True, density=0.3,
        )
        bsp.generate_data_file()
        header, _ = read_data_file(f"{temp_dir}/star/polymer.data")
        assert header["atoms"] == 3 * 33
        assert header["bonds"] == 3 * 32
        # Angles: center C(4,2)=6 + 7 per arm = 6 + 28 = 34 per chain
        assert header["angles"] == 3 * 34
        assert header["atom types"] == 2

    def test_comb_data_file_bonds_use_graph_edges(self, mock_system, temp_dir):
        comb = arch.graft(["A", "A", "A"], {1: "B"})
        bsp = BeadSpringPolymer(
            name="g", system=mock_system, n_chains=1,
            bead_types=[A, B], architecture=comb,
            generation_method="geometric",
        )
        bsp.generate_data_file()
        _, lines = read_data_file(f"{temp_dir}/g/polymer.data")
        i = lines.index("Bonds")
        bond_lines = [l for l in lines[i + 2:] if l.strip()]
        bond_pairs = {tuple(sorted((int(l.split()[2]), int(l.split()[3]))))
                      for l in bond_lines}
        assert bond_pairs == {(1, 2), (2, 3), (2, 4)}

    def test_tadpole_ring_closure_in_data_file(self, mock_system, temp_dir):
        td = arch.tadpole([("A", 6)], [("B", 3)])
        bsp = BeadSpringPolymer(
            name="td", system=mock_system, n_chains=1,
            bead_types=[A, B], architecture=td,
            generation_method="geometric",
        )
        bsp.generate_data_file()
        _, lines = read_data_file(f"{temp_dir}/td/polymer.data")
        i = lines.index("Bonds")
        bond_lines = [l for l in lines[i + 2:] if l.strip()]
        bond_pairs = {tuple(sorted((int(l.split()[2]), int(l.split()[3]))))
                      for l in bond_lines}
        assert (1, 6) in bond_pairs   # ring closure
        assert (1, 7) in bond_pairs   # tail attachment


# =============================================================================
# SAW generation for branched architectures
# =============================================================================

class TestGraphSAW:
    def test_comb_saw_bond_lengths_and_no_overlaps(self, mock_system):
        comb = arch.comb(backbone=[("A", 12)], side=("B", 2), every=3)
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=3,
            bead_types=[A, B], architecture=comb,
            generation_method="saw", density=0.3,
        )
        assert bsp.saw_generate()
        pos = np.array(bsp._positions)
        # All bonds have length ~= bond_length
        for i, j in bsp._bonds:
            dist = np.linalg.norm(pos[i] - pos[j])
            assert dist == pytest.approx(1.0, abs=0.05)
        # No overlaps (all pairs beyond collision sigma*0.9)
        seen = set()
        for i in range(len(pos)):
            for j in range(i + 1, len(pos)):
                key = (min(i, j), max(i, j))
                if key in seen:
                    continue
                seen.add(key)
                dist = np.linalg.norm(pos[i] - pos[j])
                if dist < 0.9:
                    # bonded neighbors may be at ~1.0; anything below 0.9
                    # is a genuine overlap
                    pytest.fail(f"overlap between beads {i} and {j}: {dist}")

    def test_tadpole_saw_ring_closure(self, mock_system):
        td = arch.tadpole([("A", 8)], [("B", 3)])
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A, B], architecture=td,
            generation_method="saw", box_size=10.0,
        )
        assert bsp.saw_generate()
        pos = np.array(bsp._positions)
        # Ring closure bond (7, 0) must be satisfied
        dist = np.linalg.norm(pos[7] - pos[0])
        assert dist == pytest.approx(1.0, abs=0.25)


# =============================================================================
# Branched MC moves
# =============================================================================

class TestBranchedMCMoves:
    def test_tree_pivot_preserves_bond_lengths(self):
        star = arch.star(center="A", arms=[("B", 5)] * 3)
        rng = np.random.default_rng(0)
        positions = [rng.normal(size=3) for _ in range(star.n_beads)]
        # Enforce bond lengths first
        positions[0] = np.zeros(3)
        order, parents, _ = star.growth_order()
        for bead in order[1:]:
            p = parents[bead]
            d = rng.normal(size=3)
            d /= np.linalg.norm(d)
            positions[bead] = positions[p] + d

        new_positions, valid = mc_tree_pivot_move(
            positions, pivot_idx=0, subtree_indices=[1, 2, 3, 4, 5],
            max_angle=1.0,
        )
        assert valid
        for i, j in star.bonds:
            old_d = np.linalg.norm(positions[i] - positions[j])
            new_d = np.linalg.norm(new_positions[i] - new_positions[j])
            assert new_d == pytest.approx(old_d, abs=1e-8)

    def test_segment_crankshaft_preserves_bonds(self):
        positions = [np.array([float(i), 0.1 * (i % 2), 0.0]) for i in range(8)]
        segment = [1, 2, 3, 4, 5, 6]
        new_positions, valid = mc_segment_crankshaft_move(
            positions, segment, max_angle=0.5,
        )
        assert valid
        # Endpoints of the crankshaft window unchanged structure:
        for a, b in zip(segment[:-1], segment[1:]):
            old_d = np.linalg.norm(positions[a] - positions[b])
            new_d = np.linalg.norm(new_positions[a] - new_positions[b])
            assert new_d == pytest.approx(old_d, abs=1e-8)

    def test_segment_crankshaft_too_short_invalid(self):
        positions = [np.zeros(3) for _ in range(3)]
        _, valid = mc_segment_crankshaft_move(positions, [0, 1, 2])
        assert not valid

    def test_build_chain_graph_linear_not_branched(self):
        g = build_chain_graph(arch.linear([("A", 10)]), start=0)
        assert not g["branched"]

    def test_build_chain_graph_star(self):
        star = arch.star(center="A", arms=[("B", 4)] * 3)
        g = build_chain_graph(star, start=5)
        assert g["branched"]
        # All 12 edges are bridges (tree)
        assert len(g["pivot_edges"]) == 12
        assert (5, 6) in g["edges"]

    def test_build_chain_graph_ring_has_no_bridges(self):
        g = build_chain_graph(arch.ring([("A", 8)]), start=0)
        assert g["branched"]  # not a simple path
        assert g["pivot_edges"] == []
        # But it has a crankshaft-able segment
        assert len(g["segments"]) == 1

    def test_branched_mc_equilibration_runs(self, mock_system):
        star = arch.star(center="A", arms=[("B", 8)] * 3)
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=2,
            bead_types=[A, B], architecture=star,
            generation_method="saw", density=0.3,
        )
        assert bsp.saw_generate()
        bsp.equilibrate(MCConfig(n_steps=100))
        # Most bonds stay near equilibrium after MC (displacement moves
        # allow thermal stretching, and 100 steps in a dense box does not
        # fully equilibrate; tree-pivot and segment crankshaft preserve
        # bond lengths exactly). Minimum-image distance: MC wraps
        # coordinates into the box.
        pos = np.array(bsp._positions)
        box = bsp._calculate_box_size()
        n_ok = 0
        for i, j in bsp._bonds:
            dr = pos[i] - pos[j]
            dr = dr - box * np.round(dr / box)
            if abs(np.linalg.norm(dr) - 1.0) < 0.3:
                n_ok += 1
        assert n_ok >= 0.7 * len(bsp._bonds)


class TestBackendSelection:
    def test_default_backend_is_moltemplate(self, mock_system, temp_dir):
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A], sequence=[("A", 5)], density=0.1,
        )
        assert bsp._backend == "moltemplate"
        bsp.generate(run_moltemplate=False)
        assert (Path(temp_dir) / "t" / "moltemplate" / "system.lt").exists()
        # Direct-writer outputs are NOT produced by default
        assert not (Path(temp_dir) / "t" / "polymer.data").exists()

    def test_generate_with_direct_backend(self, mock_system, temp_dir):
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A], sequence=[("A", 5)], density=0.1,
        )
        bsp.generate(backend="direct")
        assert (Path(temp_dir) / "t" / "polymer.data").exists()
        assert not (Path(temp_dir) / "t" / "moltemplate").exists()

    def test_constructor_direct_backend(self, mock_system, temp_dir):
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A], sequence=[("A", 5)], density=0.1,
            backend="direct",
        )
        bsp.generate()
        assert (Path(temp_dir) / "t" / "polymer.data").exists()

    def test_invalid_backend_raises(self, mock_system):
        with pytest.raises(ValueError, match="Backend"):
            BeadSpringPolymer(
                name="t", system=mock_system, n_chains=1,
                bead_types=[A], sequence=[("A", 5)], backend="banana",
            )

    def test_generate_invalid_backend_raises(self, mock_system):
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A], sequence=[("A", 5)], density=0.1,
        )
        with pytest.raises(ValueError, match="Backend"):
            bsp.generate(backend="banana")

    def test_system_default_backend_is_moltemplate(self, mock_system, temp_dir):
        bss = BeadSpringSystem(
            name="m", system=mock_system, bead_types=[A], density=0.1,
        )
        assert bss._backend == "moltemplate"
        bss.add_species(arch.linear([("A", 4)]), 1)
        bss.generate(run_moltemplate=False)
        assert (Path(temp_dir) / "m" / "moltemplate" / "system.lt").exists()
        assert not (Path(temp_dir) / "m" / "polymer.data").exists()

    def test_system_direct_backend(self, mock_system, temp_dir):
        bss = BeadSpringSystem(
            name="m", system=mock_system, bead_types=[A], density=0.1,
            backend="direct",
        )
        bss.add_species(arch.linear([("A", 4)]), 1)
        bss.generate()
        assert (Path(temp_dir) / "m" / "polymer.data").exists()


class TestSAWRetries:
    def test_system_retries_on_failure(self, mock_system, monkeypatch):
        """Whole-system SAW retries until success (system_retries)."""
        import AutoPoly.models.bead_spring as bs_module

        calls = {"n": 0}
        real_fn = bs_module.saw_generate_graphs

        def flaky(*args, **kwargs):
            calls["n"] += 1
            if calls["n"] < 3:
                return None, None, {
                    "success": False,
                    "chains_completed": 0,
                    "total_backtracks": 0,
                    "failure_reason": "chain_growth",
                }
            return real_fn(*args, **kwargs)

        monkeypatch.setattr(bs_module, "saw_generate_graphs", flaky)

        star = arch.star(center="A", arms=[("B", 5)] * 3)
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A, B], architecture=star,
            box_size=10.0,
            saw_config=SAWConfig(system_retries=5),
        )
        assert bsp.saw_generate()
        assert calls["n"] == 3

    def test_system_retries_exhausted(self, mock_system, monkeypatch):
        """Returns False after all retries are exhausted."""
        import AutoPoly.models.bead_spring as bs_module

        calls = {"n": 0}

        def always_fail(*args, **kwargs):
            calls["n"] += 1
            return None, None, {
                "success": False,
                "chains_completed": 0,
                "total_backtracks": 0,
                "failure_reason": "chain_growth",
            }

        monkeypatch.setattr(bs_module, "saw_generate_graphs", always_fail)

        star = arch.star(center="A", arms=[("B", 5)] * 3)
        bsp = BeadSpringPolymer(
            name="t", system=mock_system, n_chains=1,
            bead_types=[A, B], architecture=star,
            box_size=10.0,
            saw_config=SAWConfig(system_retries=2),
        )
        assert not bsp.saw_generate()
        assert calls["n"] == 2


# =============================================================================
# Mixtures (BeadSpringSystem)
# =============================================================================

class TestBeadSpringSystem:
    def test_add_species_and_counts(self, mock_system):
        bss = BeadSpringSystem(name="m", system=mock_system, bead_types=[A, B])
        bss.add_species(arch.ring([("A", 10)]), 2)
        bss.add_species(arch.comb(backbone=[("A", 8)], side="B", every=2), 3)
        assert bss.n_species == 2
        assert bss.n_chains == 5
        assert bss.total_beads == 2 * 10 + 3 * 12

    def test_add_species_validation(self, mock_system):
        bss = BeadSpringSystem(name="m", system=mock_system, bead_types=[A])
        with pytest.raises(ValueError, match="n_chains"):
            bss.add_species(arch.linear([("A", 5)]), 0)
        with pytest.raises(ValueError, match="Unknown bead type"):
            bss.add_species(arch.linear([("Z", 5)]), 1)

    def test_generate_without_species_raises(self, mock_system):
        bss = BeadSpringSystem(name="m", system=mock_system, bead_types=[A])
        with pytest.raises(ValueError, match="add_species"):
            bss.generate_data_file()

    def test_mixture_data_file_counts(self, mock_system, temp_dir):
        bss = BeadSpringSystem(
            name="blend", system=mock_system, bead_types=[A, B],
            use_angles=True, density=0.2,
        )
        bss.add_species(arch.ring([("A", 8)]), 2)
        bss.add_species(arch.linear([("A", 6), ("B", 2)]), 3)
        bss.generate_data_file()
        header, _ = read_data_file(f"{temp_dir}/blend/polymer.data")
        assert header["atoms"] == 2 * 8 + 3 * 8
        assert header["bonds"] == 2 * 8 + 3 * 7
        # Angles: ring 8/chain, linear 6/chain
        assert header["angles"] == 2 * 8 + 3 * 6
        assert header["atom types"] == 2

    def test_mixture_molecule_ids(self, mock_system, temp_dir):
        bss = BeadSpringSystem(
            name="blend", system=mock_system, bead_types=[A],
            density=0.2,
        )
        bss.add_species(arch.linear([("A", 4)]), 3)
        bss.generate_data_file()
        _, lines = read_data_file(f"{temp_dir}/blend/polymer.data")
        i = lines.index("Atoms  # molecular")
        atom_lines = [l for l in lines[i + 2:] if l.strip() and not l.startswith("Bonds")]
        mol_ids = [int(l.split()[1]) for l in atom_lines]
        assert sorted(set(mol_ids)) == [1, 2, 3]

    def test_mixture_input_script_written(self, mock_system, temp_dir):
        bss = BeadSpringSystem(
            name="blend", system=mock_system, bead_types=[A, B],
            bond_style="fene", pair_style="wca", density=0.2,
        )
        bss.add_species(arch.linear([("A", 5)]), 2)
        bss.generate_data_file()
        with open(f"{temp_dir}/blend/in.polymer") as f:
            content = f.read()
        assert "bond_style      fene" in content
        assert "special_bonds   fene" in content
        assert "pair_coeff      1 2" in content  # cross pair from mixed types

    def test_mixture_shared_type_table(self, mock_system, temp_dir):
        bss = BeadSpringSystem(
            name="blend", system=mock_system, bead_types=[A, B, C],
            density=0.2,
        )
        bss.add_species(arch.linear([("C", 4)]), 1)
        bss.add_species(arch.star(center="A", arms=[("B", 3)] * 2), 1)
        bss.generate_data_file()
        _, lines = read_data_file(f"{temp_dir}/blend/polymer.data")
        i = lines.index("Masses")
        mass_lines = []
        for l in lines[i + 2:]:
            if not l.strip():
                break
            mass_lines.append(l)
        assert len(mass_lines) == 3
        # C beads (type 3) used in first species
        j = lines.index("Atoms  # molecular")
        atom_lines = [l for l in lines[j + 2:] if l.strip() and not l.startswith("Bonds")]
        first_species_types = {int(l.split()[2]) for l in atom_lines[:4]}
        assert first_species_types == {3}

    def test_get_system_info(self, mock_system):
        bss = BeadSpringSystem(name="m", system=mock_system, bead_types=[A, B])
        bss.add_species(arch.ring([("A", 10)]), 2, name="rings")
        bss.add_species(arch.linear(arch.random_sequence(["A", "B"], 12, seed=1)), 1)
        info = bss.get_system_info()
        assert info["n_species"] == 2
        assert info["n_chains"] == 3
        assert info["species"][0]["name"] == "rings"
        assert info["species"][0]["n_bonds_per_chain"] == 10
        assert info["total_bonds"] == 2 * 10 + 11
