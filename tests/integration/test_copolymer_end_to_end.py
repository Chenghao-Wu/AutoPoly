#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
End-to-end regression tests for copolymer chemistry and junction geometry.

These run the COMPLETE generation pipeline and validate the generated
system.data, because the two v1.0 workflow bugs were invisible at the API
level (Polymer metadata looked correct while the output was wrong):

  * Bug B — positional variant keying: two chemically distinct 'middle'
    monomers in one sequence silently collapsed to the first one.  These
    tests assert the output contains the right atoms (phenyls / branches).
  * Bug A — ChainGrowthMC produced zero-length junction bonds.  These tests
    assert every bond in the output has a physical length.
"""

import numpy as np
import pytest
from collections import defaultdict

from AutoPoly import System, Polymer, generate

# PE2 - PS3 - PE2 triblock (corrected styrene complement SMILES)
ABA_SEQUENCE = [
    "CC[*]",                # ethylene first
    "[*]CC[*]",             # ethylene middle
    "[*]CC(c1ccccc1)[*]",   # styrene middle
    "[*]CC(c1ccccc1)[*]",   # styrene middle
    "[*]CC(c1ccccc1)[*]",   # styrene middle
    "[*]CC[*]",             # ethylene middle
    "[*]CC",                # ethylene last
]

# PE20 -b- PP20 diblock
PE_PP_SEQUENCE = (
    ["CC[*]"] + ["[*]CC[*]"] * 19 + ["[*]CC(C)[*]"] * 19 + ["[*]CC(C)"]
)


def _parse_data_file(path):
    """Parse a LAMMPS data file -> (masses, atoms, bonds).

    atoms: {id: (mol, type, charge, xyz)}; bonds: [(i, j)].
    """
    lines = open(path).read().splitlines()

    def section(name):
        out, on = [], False
        for line in lines:
            if line.strip().startswith(name):
                on = True
                continue
            if on:
                if not line.strip():
                    if out:
                        break
                    continue
                if line.strip() in ("Atoms", "Bonds", "Angles", "Masses"):
                    break
                out.append(line.split("#")[0].split())
        return out

    masses = {int(r[0]): float(r[1]) for r in section("Masses")}
    atoms = {
        int(r[0]): (int(r[1]), int(r[2]), float(r[3]),
                    np.array([float(r[4]), float(r[5]), float(r[6])]))
        for r in section("Atoms")
    }
    bonds = [(int(r[2]), int(r[3])) for r in section("Bonds")]
    return masses, atoms, bonds


def _chains(atoms, bonds):
    adj = defaultdict(set)
    for i, j in bonds:
        adj[i].add(j)
        adj[j].add(i)
    seen, comps = set(), []
    for i in atoms:
        if i in seen:
            continue
        stack, comp = [i], []
        while stack:
            u = stack.pop()
            if u in seen:
                continue
            seen.add(u)
            comp.append(u)
            stack.extend(adj[u] - seen)
        comps.append(comp)
    return comps, adj


def _assert_no_zero_length_bonds(atoms, bonds, min_len=0.9):
    for i, j in bonds:
        d = np.linalg.norm(atoms[i][3] - atoms[j][3])
        assert d > min_len, f"bond {i}-{j} has length {d:.4f} A"


@pytest.mark.integration
@pytest.mark.slow
class TestCopolymerEndToEnd:
    def test_aba_triblock_chemistry(self, tmp_path):
        """ABA (PE-PS-PE) output must contain the phenyl rings — Bug B used
        to silently replace the styrene middles with ethylene."""
        system = System(out=str(tmp_path / "aba"))
        poly = Polymer(chain_num=2, sequence=ABA_SEQUENCE, tacticity="atactic")
        generate(system, "aba", [poly], force_field="oplsaa",
                 box_dims=(40.0, 40.0, 40.0))

        masses, atoms, bonds = _parse_data_file(tmp_path / "aba" / "aba" / "system.data")
        comps, adj = _chains(atoms, bonds)

        def is_carbon(i):
            return masses[atoms[i][1]] > 10

        # 74 atoms per chain: PE first/last (7+7) + 2 PE mid (12) + 3 PS mid (48)
        assert sorted(len(c) for c in comps) == [74, 74]
        for comp in comps:
            nC = sum(1 for i in comp if is_carbon(i))
            nH = sum(1 for i in comp if not is_carbon(i))
            assert (nC, nH) == (32, 42)
            # 3 phenyl ipso carbons + 3 backbone CH branch points per chain;
            # zero means the phenyls were lost (the bug)
            ipso = sum(1 for i in comp if is_carbon(i)
                       and sum(1 for j in adj[i] if is_carbon(j)) == 3)
            assert ipso == 6

        _assert_no_zero_length_bonds(atoms, bonds)

        q = sum(a[2] for a in atoms.values())
        assert q == pytest.approx(0.0, abs=1e-4)

    def test_pe_pp_diblock_chemistry(self, tmp_path):
        """PE-b-PP output must contain the 19 methyl branches per chain —
        Bug B used to silently build linear polyethylene instead."""
        system = System(out=str(tmp_path / "pe_pp"))
        poly = Polymer(chain_num=2, sequence=PE_PP_SEQUENCE)
        generate(system, "pe_pp", [poly], force_field="oplsaa")

        masses, atoms, bonds = _parse_data_file(tmp_path / "pe_pp" / "pe_pp" / "system.data")
        comps, adj = _chains(atoms, bonds)

        def is_carbon(i):
            return masses[atoms[i][1]] > 10

        assert sorted(len(c) for c in comps) == [302, 302]
        for comp in comps:
            nC = sum(1 for i in comp if is_carbon(i))
            nH = sum(1 for i in comp if not is_carbon(i))
            assert (nC, nH) == (100, 202)
            branch = sum(1 for i in comp if is_carbon(i)
                         and sum(1 for j in adj[i] if is_carbon(j)) == 3)
            assert branch == 19

        _assert_no_zero_length_bonds(atoms, bonds)

        q = sum(a[2] for a in atoms.values())
        assert q == pytest.approx(0.0, abs=1e-4)
