#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for the beta-cristobalite(111) substrate slab builder."""

import numpy as np
import pytest

from AutoPoly.core.exceptions import ValidationError
from AutoPoly.packing import SubstrateSpec
from AutoPoly.surfaces import CristobaliteBuilder, SLAB_BUILDERS


def _build(**kw):
    args = dict(lx_target=40.0, ly_target=53.0, thickness=12.0)
    args.update(kw)
    return CristobaliteBuilder(**args).build()


class TestCellSnapping:
    def test_lateral_dims_are_cell_multiples(self):
        slab = _build()
        assert slab.lx == pytest.approx(4 * 10.125769)
        assert slab.ly == pytest.approx(3 * 17.538347)

    def test_minimum_one_cell(self):
        slab = CristobaliteBuilder(lx_target=5.0, ly_target=9.0,
                                   thickness=13.0).build()
        assert slab.lx == pytest.approx(10.125769)
        assert slab.ly == pytest.approx(17.538347)


class TestStructure:
    def test_every_si_four_coordinated_and_neutral(self):
        slab = _build()
        n_si = sum(a.role == "Si" for a in slab.atoms)
        n_ob = sum(a.role == "OB" for a in slab.atoms)
        n_oh = sum(a.role == "OH" for a in slab.atoms)
        n_ho = sum(a.role == "HO" for a in slab.atoms)
        assert 2 * n_ob + n_oh == 4 * n_si
        assert n_ho == n_oh
        assert abs(sum(a.charge for a in slab.atoms)) < 1e-6

    def test_sio2_stoichiometry(self):
        """O:Si ratio between 2:1 (bulk) and the surface-enriched value;
        one H per silanol O."""
        slab = _build()
        n_si = sum(a.role == "Si" for a in slab.atoms)
        n_o = sum(a.role in ("OB", "OH") for a in slab.atoms)
        n_ho = sum(a.role == "HO" for a in slab.atoms)
        n_oh = sum(a.role == "OH" for a in slab.atoms)
        assert n_ho == n_oh
        assert n_o > 2.0 * n_si  # surface O enrichment
        # each silanol O replaces half a bridging O: O = 2*Si + silanols/2
        assert n_o == 2 * n_si + (slab.n_silanol_top
                                  + slab.n_silanol_bottom) // 2

    def test_slab_fits_envelope(self):
        slab = _build(thickness=12.0)
        assert slab.z_extent <= 12.0

    def test_intrinsic_q3_density(self):
        """beta-cristobalite(111): one dangling bond per surface Si,
        ~4.5 silanols/nm^2 (matches the Zhuravlev density)."""
        slab = _build(oh_density=4.6)
        area_nm2 = slab.lx * slab.ly / 100.0
        top = slab.n_silanol_top / area_nm2
        assert top == pytest.approx(4.5, abs=0.15)
        bottom = slab.n_silanol_bottom / area_nm2
        assert bottom == pytest.approx(4.5, abs=0.15)

    def test_reproducible_with_seed(self):
        a = _build(seed=1)
        b = _build(seed=1)
        pa = np.array([x.position for x in a.atoms])
        pb = np.array([x.position for x in b.atoms])
        assert len(a.atoms) == len(b.atoms)
        np.testing.assert_allclose(pa, pb)

    def test_bond_lengths_and_angles(self):
        """All Si-O bonds 1.61 A (lateral minimum image)."""
        slab = CristobaliteBuilder(lx_target=21.0, ly_target=36.0,
                                   thickness=13.0).build()
        si = [a for a in slab.atoms if a.role == "Si"]
        ox = [a for a in slab.atoms if a.role in ("OB", "OH")]
        for s in si[:20]:
            d = np.array([o.position for o in ox]) - s.position
            d[:, 0] -= slab.lx * np.round(d[:, 0] / slab.lx)
            d[:, 1] -= slab.ly * np.round(d[:, 1] / slab.ly)
            dist = np.sort(np.linalg.norm(d, axis=1))
            assert dist[:4] == pytest.approx([1.61] * 4, abs=1e-2)


class TestLtOutput:
    def test_write_lt(self, tmp_path):
        slab = _build()
        path = slab.write_lt(tmp_path)
        text = path.read_text()
        assert f"{slab.class_name} {{" in text
        assert 'write("Data Atoms")' in text
        n_lines = sum(1 for l in text.splitlines()
                      if l.strip().startswith("$atom:"))
        assert n_lines == len(slab.atoms)
        # bond topology with zero force constants (special_bonds
        # exclusion; slab is frozen so no forces are needed)
        assert 'write("Data Bonds")' in text
        bond_lines = [l for l in text.splitlines()
                      if l.strip().startswith("$bond:")]
        n_si = sum(a.role == "Si" for a in slab.atoms)
        assert len(bond_lines) == 4 * n_si + slab.n_silanol_top \
            + slab.n_silanol_bottom
        for l in text.splitlines():
            if "bond_coeff" in l:
                assert " 0.0 " in l
        assert slab.bonds is not None and len(slab.bonds) == len(bond_lines)
        assert "@atom:i15_sc4" in text


class TestRegistry:
    def test_builder_registry(self):
        assert "beta_cristobalite" in SLAB_BUILDERS
        assert "alpha_quartz" in SLAB_BUILDERS

    def test_spec_accepts_cristobalite(self):
        spec = SubstrateSpec(builder="beta_cristobalite", thickness=13.0)
        assert spec.is_builder

    def test_unknown_builder_still_rejected(self):
        with pytest.raises(ValidationError, match="Unknown substrate builder"):
            SubstrateSpec(builder="tridymite")

    def test_thickness_floor(self):
        with pytest.raises(ValidationError, match="thickness"):
            _build(thickness=5.0)
