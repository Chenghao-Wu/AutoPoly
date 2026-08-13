#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Force-field parameter tables for built-in substrate surface builders.

Two fixed parameter sets are shipped for silica (SiO2) slabs.  Both are
applied *as published* to the slab only; the film is typed by the regular
AutoPoly pipeline (GAFF/OPLS), and cross-interactions follow LAMMPS mixing
rules.

Tables are keyed by *role*:

- "Si": tetrahedral silicon (4 O neighbors)
- "OB": bridging oxygen (2 Si neighbors)
- "OH": hydroxyl (silanol) oxygen (1 Si + 1 H neighbor)
- "HO": hydroxyl hydrogen

Each entry gives the moltemplate @atom type name, mass (amu), charge (e),
and 12-6 Lennard-Jones parameters in LAMMPS ``lj/cut/coul/long`` form
(epsilon in kcal/mol, sigma in Angstrom, converted from the published
well-minimum distance Rmin via ``sigma = Rmin / 2**(1/6)``).

Since the built slabs are intended to be held rigid/frozen in MD, only
nonbonded parameters are tabulated.  Published bonded terms (quoted in the
references below) can be added by the user via ``slab_types``/custom .lt
files if a flexible slab is wanted.

References
----------
INTERFACE FF v1.5 silica parameters (types SC4/OC23/OC24/HOY):
    Emami, F. S.; Puddu, V.; Berry, R. J.; Varshney, V.; Patwardhan, S. V.;
    Perry, C. C.; Heinz, H. "Force Field and a Surface Model Database for
    Silica to Simulate Interfacial Properties in Atomic Resolution",
    Chem. Mater. 2014, 26, 2647-2658.  DOI: 10.1021/cm500365c.
    Parameter file: charmm27_interface_v1_5.prm (H. Heinz lab distribution).

CLAYFF (types st/ob/oh/ho):
    Cygan, R. T.; Liang, J.-J.; Kalinichev, A. G. "Molecular Models of
    Hydroxide, Oxyhydroxide, and Clay Phases and the Development of a
    General Force Field", J. Phys. Chem. B 2004, 108, 1255-1266.
    DOI: 10.1021/jp0363287.

Created on 2026-08-12
@author: zwu
"""

#: 2**(1/6): converts a well-minimum distance Rmin to LJ sigma.
_TWO_POW_1_6 = 2.0 ** (1.0 / 6.0)


def _entry(type_name, mass, charge, epsilon, rmin):
    """Build a table entry; sigma derived from the well minimum Rmin."""
    return {
        "type": type_name,
        "mass": mass,
        "charge": charge,
        "epsilon": epsilon,
        "sigma": rmin / _TWO_POW_1_6 if rmin > 0.0 else 0.0,
    }


#: INTERFACE FF v1.5, silica set (CHARMM form: epsilon, Rmin/2 in .prm;
#: Rmin = 2 * Rmin/2).  Charges: neutral slab requires every Si to stay
#: 4-coordinated (the quartz builder guarantees this).
INTERFACE_FF = {
    "Si": _entry("i15_sc4",  28.0855, +1.100, 0.093, 2.0 * 2.075),
    "OB": _entry("i15_oc23", 15.9994, -0.550, 0.054, 2.0 * 1.735),
    "OH": _entry("i15_oc24", 15.9994, -0.675, 0.122, 2.0 * 1.735),
    "HO": _entry("i15_hoy",   1.0080, +0.400, 0.015, 2.0 * 0.5425),
}

#: CLAYFF (E = D0[(R0/r)^12 - 2(R0/r)^6]; epsilon = D0, sigma = R0/2^(1/6)).
#: D0 converted from kJ/mol to kcal/mol (divide by 4.184).
CLAYFF = {
    "Si": _entry("cff_st", 28.0855, +2.100, 7.70065e-6 / 4.184, 3.30203),
    "OB": _entry("cff_ob", 15.9994, -1.050, 0.650194 / 4.184,   3.16554),
    "OH": _entry("cff_oh", 15.9994, -0.950, 0.650194 / 4.184,   3.16554),
    "HO": _entry("cff_ho",  1.0080, +0.425, 0.0,                0.0),
}

#: Registry of built-in slab force fields.
SLAB_FF_REGISTRY = {
    "interface": INTERFACE_FF,
    "clayff": CLAYFF,
}

#: Equilibrium O-H bond length used when capping dangling surface oxygens
#: (INTERFACE FF: OC23-HOY 0.945 A; CLAYFF commonly uses ~1.0 A).
OH_BOND_LENGTH = {
    "interface": 0.945,
    "clayff": 1.0,
}

#: All roles a slab force field must cover.
SLAB_ROLES = ("Si", "OB", "OH", "HO")
