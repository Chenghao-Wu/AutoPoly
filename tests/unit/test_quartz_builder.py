#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for the alpha-quartz(0001) substrate slab builder."""

import numpy as np
import pytest

from AutoPoly.core.exceptions import ValidationError
from AutoPoly.packing import SubstrateSpec
from AutoPoly.surfaces import QuartzBuilder
from AutoPoly.surfaces.ff_tables import (
    CLAYFF,
    INTERFACE_FF,
    SLAB_ROLES,
)


def _build(**kw):
    args = dict(lx_target=30.0, ly_target=34.0, thickness=12.0)
    args.update(kw)
    return QuartzBuilder(**args).build()


class TestCellSnapping:
    def test_lateral_dims_are_cell_multiples(self):
        slab = _build()
        assert slab.lx == pytest.approx(6 * 4.9019)
        assert slab.ly == pytest.approx(4 * 4.9019 * np.sqrt(3.0))

    def test_minimum_one_cell(self):
        slab = QuartzBuilder(lx_target=2.0, ly_target=3.0,
                             thickness=12.0).build()
        assert slab.lx == pytest.approx(4.9019)
        assert slab.ly == pytest.approx(4.9019 * np.sqrt(3.0))


class TestStructure:
    def test_every_si_four_coordinated_and_neutral(self):
        slab = _build()
        n_si = sum(a.role == "Si" for a in slab.atoms)
        n_ob = sum(a.role == "OB" for a in slab.atoms)
        n_oh = sum(a.role == "OH" for a in slab.atoms)
        n_ho = sum(a.role == "HO" for a in slab.atoms)
        # full Si coordination <=> 2*N_OB + N_OH = 4*N_Si
        assert 2 * n_ob + n_oh == 4 * n_si
        # hydroxylation: one H per silanol O
        assert n_ho == n_oh
        # charge neutrality
        q = sum(a.charge for a in slab.atoms)
        assert abs(q) < 1e-6

    def test_slab_fits_envelope_and_centers(self):
        slab = _build(thickness=12.0)
        pos = np.array([a.position for a in slab.atoms])
        assert slab.z_extent <= 12.0
        # centered laterally and at mid-plane in z
        assert pos[:, 0].min() >= -slab.lx / 2 - 1e-6
        assert pos[:, 0].max() <= slab.lx / 2 + 1e-6
        assert pos[:, 2].min() == pytest.approx(-slab.z_extent / 2)
        assert pos[:, 2].max() == pytest.approx(slab.z_extent / 2)

    def test_intrinsic_q2_density(self):
        """alpha-quartz(0001) full hydroxylation: 2 silanols per surface
        Si plane (Q2 geminal), i.e. ~9.6/nm^2 on this cut."""
        slab = _build(oh_density=4.6)  # below intrinsic: stays Q2
        area_nm2 = slab.lx * slab.ly / 100.0
        top = slab.n_silanol_top / area_nm2
        assert top == pytest.approx(9.6, abs=0.6)

    def test_reproducible_with_seed(self):
        a = _build(seed=1)
        b = _build(seed=1)
        pa = np.array([x.position for x in a.atoms])
        pb = np.array([x.position for x in b.atoms])
        assert len(a.atoms) == len(b.atoms)
        np.testing.assert_allclose(pa, pb)


class TestForceFields:
    @pytest.mark.parametrize("ff", ["interface", "clayff"])
    def test_builtin_tables_used(self, ff):
        slab = _build(slab_ff=ff)
        table = {"interface": INTERFACE_FF, "clayff": CLAYFF}[ff]
        by_role = {}
        for a in slab.atoms:
            by_role.setdefault(a.role, a)
        for role, atom in by_role.items():
            assert atom.type_name == table[role]["type"]
            assert atom.charge == pytest.approx(table[role]["charge"])

    def test_custom_table(self):
        types = {"Si": "xsi", "OB": "xob", "OH": "xoh", "HO": "xho"}
        charges = {"Si": 1.5, "OB": -0.75, "OH": -0.85, "HO": 0.35}
        lj = {r: (0.1, 3.0) for r in SLAB_ROLES}
        slab = _build(slab_ff="custom", slab_types=types,
                      slab_charges=charges, slab_lj=lj)
        by_role = {}
        for a in slab.atoms:
            by_role.setdefault(a.role, a)
        assert by_role["Si"].type_name == "xsi"
        assert by_role["Si"].charge == pytest.approx(1.5)
        assert by_role["HO"].type_name == "xho"

    def test_custom_table_must_cover_all_roles(self):
        with pytest.raises(ValidationError, match="custom"):
            _build(slab_ff="custom",
                   slab_types={"Si": "xsi"},
                   slab_charges={"Si": 1.5},
                   slab_lj={"Si": (0.1, 3.0)})

    def test_unknown_ff_rejected(self):
        with pytest.raises(ValidationError, match="slab_ff"):
            _build(slab_ff="nosuchff")


class TestLtOutput:
    def test_write_lt(self, tmp_path):
        slab = _build()
        path = slab.write_lt(tmp_path)
        text = path.read_text()
        assert f"{slab.class_name} {{" in text
        assert 'write_once("Data Masses")' in text
        assert 'write_once("In Settings")' in text
        assert 'write("Data Atoms")' in text
        # one line per atom, charges included
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
        for a in slab.atoms:
            assert f"@atom:{a.type_name}" in text

    def test_type_names_do_not_collide_with_gaff(self):
        """Built-in slab type names are namespaced (GAFF uses lowercase
        like c3/oh; OPLS uses numeric opls_xxx)."""
        for table in (INTERFACE_FF, CLAYFF):
            for role, entry in table.items():
                t = entry["type"]
                assert t.startswith(("i15_", "cff_"))


class TestSpecValidation:
    def test_builder_source(self):
        spec = SubstrateSpec(builder="alpha_quartz", thickness=12.0)
        assert spec.is_builder
        assert not spec.is_external

    def test_exactly_one_source(self):
        from AutoPoly import Molecule
        with pytest.raises(ValidationError, match="exactly one source"):
            SubstrateSpec(builder="alpha_quartz",
                          model=Molecule(Count=1, Smiles="CCO",
                                         Name="x"))
        with pytest.raises(ValidationError, match="exactly one source"):
            SubstrateSpec(thickness=10.0)

    def test_unknown_builder(self):
        with pytest.raises(ValidationError, match="Unknown substrate builder"):
            SubstrateSpec(builder="cristobalite")

    def test_thickness_floor(self):
        with pytest.raises(ValidationError, match="thickness"):
            _build(thickness=6.0)

    def test_negative_density(self):
        with pytest.raises(ValidationError):
            SubstrateSpec(builder="alpha_quartz", oh_density=-1.0)
