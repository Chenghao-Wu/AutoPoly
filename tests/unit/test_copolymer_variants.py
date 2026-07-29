#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Regression tests for block-copolymer variant mapping (Bug B) and junction
bond length in MC chain growth (Bug A).

Bug B: generate_sequence_variants_for_polymerization used to deduplicate
monomer variants by variant_type alone, so two chemically distinct 'middle'
monomers in one sequence collapsed to the first one (wrong chemistry).

Bug A: ChainGrowthMC used to place each incoming monomer's left connection
exactly on the previous right connection, producing zero-length junction
bonds.
"""
import numpy as np
import pytest

from AutoPoly import monomer_processing
from AutoPoly.mc import CollisionDetector, ChainGrowthMC

# ABA triblock: two chemically distinct middle monomers (PE and PS)
ABA_SEQUENCE = [
    "CC[*]",                # first
    "[*]CC[*]",             # middle (ethylene)
    "[*]CC(c1ccccc1)[*]",   # middle (styrene)
    "[*]CC(c1ccccc1)[*]",   # middle (styrene)
    "[*]CC(c1ccccc1)[*]",   # middle (styrene)
    "[*]CC[*]",             # middle (ethylene)
    "[*]CC",                # last
]

PE_SEQUENCE = ["CC[*]"] + ["[*]CC[*]"] * 4 + ["[*]CC"]


@pytest.fixture(scope="module")
def aba_mapping(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("variants")
    mapping, _ = monomer_processing.generate_sequence_variants_for_polymerization(
        smiles_list=ABA_SEQUENCE,
        topology="linear",
        path_cwd=str(tmp),
        force_field="oplsaa",
        generated_cache={},
        counter=0,
    )
    return mapping, tmp


class TestCopolymerVariantMapping:
    def test_by_position_present_and_full_length(self, aba_mapping):
        mapping, _ = aba_mapping
        assert "_by_position" in mapping
        assert len(mapping["_by_position"]) == len(ABA_SEQUENCE)

    def test_distinct_middles_get_distinct_files(self, aba_mapping):
        """The PE middle (pos 1) and PS middles (pos 2-4) must map to
        different .lt files — previously all middles collapsed to one file."""
        mapping, _ = aba_mapping
        by_position = mapping["_by_position"]
        pe_middle_file = by_position[1][0]
        ps_middle_file = by_position[2][0]
        assert pe_middle_file != ps_middle_file
        # identical chemistries share a file
        assert by_position[2][0] == by_position[3][0] == by_position[4][0]
        assert by_position[1][0] == by_position[5][0]

    def test_template_files_exist(self, aba_mapping):
        mapping, tmp = aba_mapping
        for fname, fname_t1 in mapping["_by_position"]:
            assert (tmp / fname).exists(), fname

    def test_styrene_template_has_phenyl(self, aba_mapping):
        """The PS middle template must contain 16 atoms (C8H8); if it has 6,
        the ethylene template was used instead (the bug)."""
        mapping, tmp = aba_mapping
        ps_file = tmp / mapping["_by_position"][2][0]
        n_atoms = (ps_file.read_text()
                   .split('write("Data Atoms")')[1]
                   .split('}')[0].count('$atom:'))
        assert n_atoms == 16

    def test_legacy_keys_still_present(self, aba_mapping):
        mapping, _ = aba_mapping
        for key in ("first", "middle", "last"):
            assert key in mapping
            assert f"{key}_T1" in mapping


@pytest.fixture(scope="module")
def pe_mapping(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("variants_pe")
    mapping, _ = monomer_processing.generate_sequence_variants_for_polymerization(
        smiles_list=PE_SEQUENCE,
        topology="linear",
        path_cwd=str(tmp),
        force_field="oplsaa",
        generated_cache={},
        counter=0,
    )
    return mapping, tmp


class TestChainGrowthJunction:
    def test_junction_bonds_have_physical_length(self, pe_mapping):
        """Grown chains must not contain zero-length junction bonds."""
        mapping, tmp = pe_mapping
        lt_files = [str(tmp / fname) for fname, _ in mapping["_by_position"]]

        bounds = ((-50.0, 50.0), (-50.0, 50.0), (-50.0, 50.0))
        detector = CollisionDetector(bounds, cell_size=5.0)
        mc = ChainGrowthMC(detector, max_attempts=1000,
                           bond_angle_min=50.0, bond_angle_max=90.0)
        placements = mc.grow_chain(lt_files, chain_id=0)

        for prev, cur in zip(placements, placements[1:]):
            d = np.linalg.norm(cur.world_left_conn - prev.world_right_conn)
            assert d == pytest.approx(mc.junction_bond_length, abs=1e-6), (
                f"junction at monomer {cur.monomer_index} has length {d:.4f} A"
            )

    def test_default_junction_length(self):
        detector = CollisionDetector(((-1, 1), (-1, 1), (-1, 1)))
        mc = ChainGrowthMC(detector)
        assert mc.junction_bond_length == pytest.approx(1.54)
