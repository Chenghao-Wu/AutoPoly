"""
Unit tests for the polymer architecture graph core
(AutoPoly.models.architectures).
"""

import numpy as np
import pytest

from AutoPoly.models.architectures import (
    BeadArchitecture,
    MonomerTemplate,
    normalize_sequence,
    linear,
    ring,
    star,
    comb,
    graft,
    tadpole,
    dendrimer,
    custom,
    block_sequence,
    alternating_sequence,
    random_sequence,
    gradient_sequence,
)


class TestNormalizeSequence:
    def test_string(self):
        assert normalize_sequence("AABB") == ["A", "A", "B", "B"]

    def test_explicit_list(self):
        assert normalize_sequence(["A", "B"]) == ["A", "B"]

    def test_block_format(self):
        assert normalize_sequence([("A", 3), ("B", 2)]) == ["A"] * 3 + ["B"] * 2

    def test_bare_block_tuple(self):
        assert normalize_sequence(("B", 3)) == ["B", "B", "B"]

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            normalize_sequence([])

    def test_zero_count_block_raises(self):
        with pytest.raises(ValueError):
            normalize_sequence([("A", 0)])


class TestMonomerTemplate:
    def test_valid(self):
        tpl = MonomerTemplate(
            name="G", beads=["A", "B"],
            internal_bonds=[(0, 1)],
            connections={"head": 0, "tail": 0, "side": 1},
        )
        assert tpl.name == "G"

    def test_no_beads_raises(self):
        with pytest.raises(ValueError):
            MonomerTemplate(name="X", beads=[])

    def test_internal_bond_out_of_range_raises(self):
        with pytest.raises(ValueError):
            MonomerTemplate(name="X", beads=["A"], internal_bonds=[(0, 5)])

    def test_self_bond_raises(self):
        with pytest.raises(ValueError):
            MonomerTemplate(name="X", beads=["A", "B"], internal_bonds=[(1, 1)])

    def test_connection_out_of_range_raises(self):
        with pytest.raises(ValueError):
            MonomerTemplate(name="X", beads=["A"], connections={"head": 3})

    def test_single_bead(self):
        tpl = MonomerTemplate.single_bead("A")
        assert tpl.beads == ["A"]
        assert tpl.connections == {"head": 0, "tail": 0}


class TestBeadArchitectureValidation:
    def test_duplicate_bond_raises(self):
        arch = BeadArchitecture(bead_types=["A"] * 3,
                                bonds=[(0, 1), (1, 2), (2, 1)])
        with pytest.raises(ValueError, match="Duplicate"):
            arch.validate()

    def test_self_bond_raises(self):
        arch = BeadArchitecture(bead_types=["A"] * 2, bonds=[(0, 0)])
        with pytest.raises(ValueError, match="Self bond"):
            arch.validate()

    def test_out_of_range_raises(self):
        arch = BeadArchitecture(bead_types=["A"] * 2, bonds=[(0, 7)])
        with pytest.raises(ValueError, match="out of range"):
            arch.validate()

    def test_disconnected_raises(self):
        arch = BeadArchitecture(bead_types=["A"] * 4,
                                bonds=[(0, 1), (2, 3)])
        with pytest.raises(ValueError, match="not connected"):
            arch.validate()

    def test_unknown_bead_type_raises(self):
        arch = BeadArchitecture(bead_types=["A", "Z"], bonds=[(0, 1)])
        with pytest.raises(ValueError, match="Unknown bead type"):
            arch.validate(known_bead_types=["A"])

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            BeadArchitecture(bead_types=[], bonds=[]).validate()


class TestGraphProperties:
    def test_linear_degrees(self):
        arch = linear([("A", 5)])
        assert arch.degrees() == [1, 2, 2, 2, 1]
        assert not arch.is_branched
        assert not arch.is_cyclic
        assert arch.max_degree == 2

    def test_ring_degrees(self):
        arch = ring([("A", 5)])
        assert arch.degrees() == [2] * 5
        assert not arch.is_branched
        assert arch.is_cyclic

    def test_star_degrees(self):
        arch = star(center="A", arms=[("B", 2)] * 3)
        assert arch.degrees()[0] == 3
        assert arch.is_branched
        assert not arch.is_cyclic


class TestAngleTriplets:
    def test_linear_matches_path_enumeration(self):
        arch = linear([("A", 5)])
        assert arch.angle_triplets() == [(0, 1, 2), (1, 2, 3), (2, 3, 4)]

    def test_ring_wrap_arounds(self):
        arch = ring([("A", 4)])
        triplets = set(arch.angle_triplets())
        # Every bead is a center with two neighbors
        expected = {
            (1, 0, 3), (0, 1, 2), (1, 2, 3), (0, 3, 2),
        }
        assert triplets == expected

    def test_branch_triplet_count(self):
        # Star with 3 arms: center degree 3 -> C(3,2)=3 branch triplets,
        # each 2-bead arm adds 1 triplet (0, start, start+1)
        arch = star(center="A", arms=[("B", 2)] * 3)
        all_t = arch.angle_triplets()
        no_branch = arch.angle_triplets(include_branch=False)
        assert len(all_t) == 3 + 3
        assert len(no_branch) == 3

    def test_branch_classification(self):
        arch = star(center="A", arms=[("B", 2)] * 3)
        classified = arch.angle_triplets_classified()
        branch_triplets = [t for t, is_b in classified if is_b]
        assert all(t[1] == 0 for t in branch_triplets)  # centered on center bead
        assert len(branch_triplets) == 3

    def test_comb_graft_point_triplets(self):
        # Backbone A-A-A with side bead B on middle bead (degree 3)
        arch = graft(["A", "A", "A"], {1: "B"})
        all_t = set(arch.angle_triplets())
        # Branch triplets at bead 1: (0,1,2), (0,1,3), (2,1,3)
        assert (0, 1, 2) in all_t
        assert (0, 1, 3) in all_t
        assert (2, 1, 3) in all_t
        no_branch = set(arch.angle_triplets(include_branch=False))
        assert no_branch == set()


class TestGrowthOrder:
    def test_parent_before_child(self):
        arch = comb(backbone=[("A", 8)], side="B", every=3)
        order, parents, closing = arch.growth_order()
        position = {b: i for i, b in enumerate(order)}
        for bead, parent in parents.items():
            if parent is not None:
                assert position[parent] < position[bead]

    def test_ring_has_one_closing_edge(self):
        arch = ring([("A", 6)])
        _, _, closing = arch.growth_order()
        assert len(closing) == 1

    def test_tree_has_no_closing_edges(self):
        arch = star(center="A", arms=[("B", 3)] * 4)
        _, _, closing = arch.growth_order()
        assert closing == []

    def test_tadpole_has_one_closing_edge(self):
        arch = tadpole([("A", 6)], [("B", 3)])
        _, _, closing = arch.growth_order()
        assert len(closing) == 1
        # The closing edge is one of the ring bonds (both endpoints < 6);
        # which one depends on the BFS spanning tree
        u, v = closing[0]
        assert u < 6 and v < 6


class TestSubtreesAndSegments:
    def test_subtree_beyond_edge_star_arm(self):
        arch = star(center="A", arms=[("B", 3)] * 2)
        # Arm 1 = beads 1,2,3
        subtree = arch.subtree_beyond_edge(0, 1)
        assert subtree == [1, 2, 3]

    def test_subtree_beyond_edge_comb_side(self):
        arch = graft(["A", "A", "A"], {1: "B"})
        assert arch.subtree_beyond_edge(1, 3) == [3]

    def test_linear_segments_whole_chain(self):
        arch = linear([("A", 6)])
        segs = arch.linear_segments()
        assert len(segs) == 1
        assert segs[0] == [0, 1, 2, 3, 4, 5]

    def test_ring_segment(self):
        arch = ring([("A", 5)])
        segs = arch.linear_segments()
        assert len(segs) == 1
        # Segment walks the full ring back to start
        assert segs[0][0] == 0 and segs[0][-1] == 0
        assert len(segs[0]) == 6

    def test_comb_segments_cover_backbone_edges(self):
        arch = comb(backbone=[("A", 10)], side="B", every=2)
        segs = arch.linear_segments()
        seg_edges = set()
        for seg in segs:
            for a, b in zip(seg[:-1], seg[1:]):
                seg_edges.add((min(a, b), max(a, b)))
        # Backbone edges up to the last graft point belong to some
        # crankshaft-able segment; the final tail edge (8, 9) forms a
        # length-2 segment which is below the minimum length
        for i in range(8):
            assert (i, i + 1) in seg_edges
        assert (8, 9) not in seg_edges
        # Single-bead side chains at true branch points (degree 3) form
        # length-2 segments which are below the minimum length; the side
        # chain at bead 0 (degree 2) legitimately joins a segment
        branch_side_edges = {(2, 11), (4, 12), (6, 13), (8, 14)}
        assert not (branch_side_edges & seg_edges)


class TestFactories:
    def test_linear(self):
        arch = linear([("A", 4)])
        assert arch.n_beads == 4
        assert arch.n_bonds == 3
        assert arch.name == "linear"

    def test_ring_too_small_raises(self):
        with pytest.raises(ValueError, match="at least 3"):
            ring([("A", 2)])

    def test_star_counts(self):
        arch = star(center="A", arms=[("B", 5), ("C", 3), ("B", 5)])
        assert arch.n_beads == 1 + 5 + 3 + 5
        assert arch.n_bonds == 5 + 3 + 5
        assert arch.bead_types[0] == "A"

    def test_star_too_few_arms_raises(self):
        with pytest.raises(ValueError, match="at least 2 arms"):
            star(center="A", arms=[("B", 5)])

    def test_comb_counts(self):
        arch = comb(backbone=[("A", 10)], side=("B", 2), every=2)
        # 5 graft points (0,2,4,6,8), 2 beads each
        assert arch.n_beads == 10 + 5 * 2
        assert arch.n_bonds == 9 + 5 * 2

    def test_comb_offset(self):
        arch = comb(backbone=[("A", 10)], side="B", every=3, offset=1)
        # graft points at 1, 4, 7
        assert arch.n_beads == 10 + 3

    def test_comb_no_graft_points_raises(self):
        with pytest.raises(ValueError, match="No graft points"):
            comb(backbone=[("A", 5)], side="B", every=3, offset=10)

    def test_graft_explicit(self):
        arch = graft([("A", 6)], {0: "B", 3: ("C", 2), 5: "B"})
        assert arch.n_beads == 6 + 1 + 2 + 1
        assert arch.bead_types[6] == "B"
        assert arch.bead_types[7:9] == ["C", "C"]

    def test_graft_out_of_range_raises(self):
        with pytest.raises(ValueError, match="out of range"):
            graft([("A", 5)], {7: "B"})

    def test_tadpole_counts(self):
        arch = tadpole([("A", 8)], [("B", 4)])
        assert arch.n_beads == 12
        assert arch.n_bonds == 12  # 8 ring + 1 attach + 3 tail
        assert arch.is_cyclic

    def test_tadpole_attach_out_of_range_raises(self):
        with pytest.raises(ValueError, match="out of range"):
            tadpole([("A", 5)], [("B", 2)], attach=9)

    def test_dendrimer_counts(self):
        arch = dendrimer(core="A", branch="B", branch_factor=3, generations=2)
        # 1 + 3 + 9 = 13 beads, tree -> 12 bonds
        assert arch.n_beads == 13
        assert arch.n_bonds == 12
        assert arch.bead_types[0] == "A"

    def test_dendrimer_invalid_params_raise(self):
        with pytest.raises(ValueError):
            dendrimer(core="A", branch="B", branch_factor=1)
        with pytest.raises(ValueError):
            dendrimer(core="A", branch="B", generations=0)

    def test_custom(self):
        arch = custom(["A", "B", "C"], [(0, 1), (1, 2)])
        assert arch.n_beads == 3
        assert arch.name == "custom"

    def test_custom_invalid_raises(self):
        with pytest.raises(ValueError):
            custom(["A", "B"], [(0, 5)])


class TestFromMonomers:
    def _graft_monomer(self):
        return MonomerTemplate(
            name="G", beads=["A", "B"],
            internal_bonds=[(0, 1)],
            connections={"head": 0, "tail": 0, "side": 1},
        )

    def test_expansion(self):
        tpl = self._graft_monomer()
        arch = BeadArchitecture.from_monomers(
            {"G": tpl}, ["G", "G", "G"],
            [(0, "tail", 1, "head"), (1, "tail", 2, "head")],
        )
        assert arch.bead_types == ["A", "B"] * 3
        assert arch.n_bonds == 3 + 2  # internal + inter-monomer
        assert arch.is_branched

    def test_unknown_template_raises(self):
        with pytest.raises(ValueError, match="Unknown monomer template"):
            BeadArchitecture.from_monomers({}, ["G"], [])

    def test_missing_connection_raises(self):
        tpl = self._graft_monomer()
        with pytest.raises(ValueError, match="no connection"):
            BeadArchitecture.from_monomers(
                {"G": tpl}, ["G", "G"], [(0, "tail", 1, "nope")]
            )


class TestSequenceGenerators:
    def test_block(self):
        assert block_sequence([("A", 2), ("B", 1)]) == ["A", "A", "B"]

    def test_alternating(self):
        seq = alternating_sequence(["A", "B"], 6)
        assert seq == ["A", "B"] * 3

    def test_alternating_three_types(self):
        seq = alternating_sequence(["A", "B", "C"], 7)
        assert seq[:4] == ["A", "B", "C", "A"]

    def test_random_seeded_reproducible(self):
        s1 = random_sequence(["A", "B"], 20, seed=42)
        s2 = random_sequence(["A", "B"], 20, seed=42)
        assert s1 == s2
        assert len(s1) == 20
        assert set(s1) <= {"A", "B"}

    def test_random_weights(self):
        seq = random_sequence(["A", "B"], 200, weights=[0.99, 0.01], seed=1)
        assert seq.count("A") > seq.count("B")

    def test_random_bad_weights_raise(self):
        with pytest.raises(ValueError):
            random_sequence(["A", "B"], 10, weights=[1.0])
        with pytest.raises(ValueError):
            random_sequence(["A", "B"], 10, weights=[0.0, 0.0])

    def test_gradient_endpoints(self):
        seq = gradient_sequence("A", "B", 50, seed=3)
        assert len(seq) == 50
        assert seq[0] == "A"   # p(B)=0 at first bead
        assert seq[-1] == "B"  # p(B)=1 at last bead

    def test_gradient_too_short_raises(self):
        with pytest.raises(ValueError):
            gradient_sequence("A", "B", 1)
