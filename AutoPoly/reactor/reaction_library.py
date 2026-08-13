#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reaction Library for the AutoPoly Reactor

Reaction SMARTS patterns for step-growth polymerization, adapted from the
AutoREACTER project (NanoCIPHER-Lab, MIT license). Each reaction entry
pairs reactant functional-group names (see
:mod:`AutoPoly.reactor.functional_groups`) with an RDKit reaction SMARTS
whose atom-map numbers drive template construction:

- Map numbers 1 and 2 mark the two *initiator* atoms (the atoms between
  which ``fix bond/react`` forms/breaks the connecting bond).
- All other map numbers < 999 mark the *first shell* (reaction center
  atoms tracked into the product).
- ``delete_atom=True`` reactions split off a small byproduct fragment
  (water, HCl, ...); the byproduct atoms are listed under ``deleteIDs``
  in the generated map file so LAMMPS removes them.

Created on 2026-08-04
@author: zwu
"""
from typing import Any, Dict, Optional

# Each entry:
#   same_reactants: True for A+A homo-reactions (monomer reacts with itself)
#   reactant_1 / reactant_2: functional-group names (reactant_2 may equal
#       reactant_1; same_reactants=False with reactant_2 set also matches
#       two different monomers carrying the respective groups)
#   product: descriptive product label
#   delete_atom: whether a byproduct fragment is split off
#   reaction: RDKit reaction SMARTS (atom-mapped)
#   reference: literature pointers for the SMARTS / mechanism

REACTIONS: Dict[str, Dict[str, Any]] = {
    # ============================================================
    # Polyesterification: hydroxy acids / acid halides
    # ============================================================
    "Hydroxy Carboxylic Acid Polycondensation(Polyesterification)": {
        "same_reactants": True,
        "reactant_1": "hydroxy_carboxylic_acid",
        "product": "polyester_chain",
        "delete_atom": True,
        "reaction": "[OX2H1;!$([O][C,S]=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]",
        "reference": {
            "smarts": "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329",
            "reaction_and_mechanism": [
                "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                "https://pubs.acs.org/doi/10.1021/ed073pA312",
            ],
        },
        "comments": None,
    },
    "Hydroxy Carboxylic and Hydroxy Carboxylic Polycondensation(Polyesterification)": {
        "same_reactants": False,
        "reactant_1": "hydroxy_carboxylic_acid",
        "reactant_2": "hydroxy_carboxylic_acid",
        "product": "polyester_chain",
        "delete_atom": True,
        "reaction": "[OX2H1;!$([O][C,S]=*):1]-[H:3].[CX3:2](=[O:5])[OX2H1:4]>>[OX2:1]-[CX3:2](=[O:5]).[O:4]-[H:3]",
        "reference": {
            "smarts": "https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329",
            "reaction_and_mechanism": [
                "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                "https://pubs.acs.org/doi/10.1021/ed073pA312",
            ],
        },
        "comments": None,
    },
    "Hydroxy Acid Halides Polycondensation(Polyesterification)": {
        "same_reactants": True,
        "reactant_1": "hydroxy_acid_halide",
        "product": "polyester_chain",
        "delete_atom": True,
        "reaction": "[OX2H1;!$([O][C,S]=*):1]-[H:3].[CX3:2](=[O:5])[Cl,Br,I:4]>>[OX2:1]-[CX3:2](=[O:5]).[Cl,Br,I:4]-[H:3]",
        "reference": {
            "smarts": None,
            "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed073pA312"],
        },
        "comments": None,
    },
    "Hydroxy Acid Halides Hydroxy Acid Halides Polycondensation(Polyesterification)": {
        "same_reactants": False,
        "reactant_1": "hydroxy_acid_halide",
        "reactant_2": "hydroxy_acid_halide",
        "product": "polyester_chain",
        "delete_atom": True,
        "reaction": "[OX2H1;!$([O][C,S]=*):1]-[H:3].[CX3:2](=[O:5])[Cl,Br,I:4]>>[OX2:1]-[CX3:2](=[O:5]).[Cl,Br,I:4]-[H:3]",
        "reference": {
            "smarts": None,
            "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed073pA312"],
        },
        "comments": None,
    },
    # ============================================================
    # Polyesterification: diols + diacids / diacid halides / esters
    # ============================================================
    "Diol and Di-Carboxylic Acid Polycondensation(Polyesterification)": {
        "same_reactants": False,
        "reactant_1": "diol",
        "reactant_2": "di_carboxylic_acid",
        "product": "polyester_chain",
        "delete_atom": True,
        "reaction": "[CX3:1](=[O:3])[OX2H1:4].[OX2H1;!$([O][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[OX2:2].[O:4]-[H:5]",
        "reference": {
            "smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
            "reaction_and_mechanism": [
                "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                "https://pubs.acs.org/doi/10.1021/ed073pA312",
            ],
        },
        "comments": None,
    },
    "Diol and Di-Acid Halide Polycondensation(Polyesterification)": {
        "same_reactants": False,
        "reactant_1": "diol",
        "reactant_2": "di_carboxylic_acid_halide",
        "product": "polyester_chain",
        "delete_atom": True,
        "reaction": "[CX3:1](=[O:3])[Cl,Br,I:4].[OX2H1;!$([O][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[OX2:2].[Cl,Br,I:4]-[H:5]",
        "reference": {
            "smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
            "reaction_and_mechanism": [
                "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                "https://pubs.acs.org/doi/10.1021/ed073pA312",
            ],
        },
        "comments": None,
    },
    "Diol and Di-Carboxylic Ester Polycondensation(Transesterification)": {
        "same_reactants": False,
        "reactant_1": "diol",
        "reactant_2": "di_carboxylic_ester",
        "product": "polyester_chain",
        "delete_atom": True,
        "reaction": "[OX2H1;!$([O][C,S]=*):1]-[H:3].[CX3:2](=[O:5])[OX2H0:4][#6:6]>>[OX2:1]-[CX3:2](=[O:5]).[OX2:4](-[H:3])-[#6:6]",
        "reference": {"smarts": None, "reaction_and_mechanism": None},
        "comments": None,
    },
    # ============================================================
    # Polyanhydride formation
    # ============================================================
    "Carboxylic Acid and Acid Halide Polycondensation(Polyanhydride Formation)": {
        "same_reactants": True,
        "reactant_1": "carboxylic_acid_acid_halide",
        "product": "polyanhydride_chain",
        "delete_atom": True,
        "reaction": "[CX3:1](=[O:3])[Cl,Br,I:4].[CX3:2](=[O:5])[OX2H1:6]-[H:7]>>[CX3:1](=[O:3])-[OX2:6]-[CX3:2](=[O:5]).[Cl,Br,I:4]-[H:7]",
        "reference": {"smarts": None, "reaction_and_mechanism": None},
        "comments": None,
    },
    # ============================================================
    # Polythioesterification
    # ============================================================
    "Dithiol and Di-Carboxylic Acid Halide Polycondensation(Polythioesterification)": {
        "same_reactants": False,
        "reactant_1": "dithiol",
        "reactant_2": "di_carboxylic_acid_halide",
        "product": "polythioester_chain",
        "delete_atom": True,
        "reaction": "[CX3:1](=[O:3])[Cl,Br,I:4].[SX2H1;!$([S][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[SX2:2].[Cl,Br,I:4]-[H:5]",
        "reference": {"smarts": None, "reaction_and_mechanism": None},
        "comments": None,
    },
    "Dithiol and Di-Carboxylic Acid Polycondensation(Polythioesterification)": {
        "same_reactants": False,
        "reactant_1": "dithiol",
        "reactant_2": "di_carboxylic_acid",
        "product": "polythioester_chain",
        "delete_atom": True,
        "reaction": "[CX3:1](=[O:3])[OX2H1:4].[SX2H1;!$([S][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[SX2:2].[O:4]-[H:5]",
        "reference": {"smarts": None, "reaction_and_mechanism": None},
        "comments": "Thioesterification with water elimination; generally less "
                    "straightforward than the acid-halide route.",
    },
    # ============================================================
    # Polyamidation
    # ============================================================
    "Amino Acid Polycondensation (Polyamidation)": {
        "same_reactants": True,
        "reactant_1": "amino_acid",
        "product": "polyamide_chain",
        "delete_atom": True,
        "reaction": "[NX3;H2,H1;!$([N][C,S]=*):1]-[H:3].[CX3:2](=[O:4])[OX2H1:5]>>[NX3:1]-[CX3:2](=[O:4]).[O:5]-[H:3]",
        "reference": {
            "smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
            "reaction_and_mechanism": [
                "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                "https://pubs.acs.org/doi/10.1021/ed073pA312",
            ],
        },
        "comments": None,
    },
    "Amino Acid and Amino Acid Polycondensation (Polyamidation)": {
        "same_reactants": False,
        "reactant_1": "amino_acid",
        "reactant_2": "amino_acid",
        "product": "polyamide_chain",
        "delete_atom": True,
        "reaction": "[NX3;H2,H1;!$([N][C,S]=*):1]-[H:3].[CX3:2](=[O:4])[OX2H1:5]>>[NX3:1]-[CX3:2](=[O:4]).[O:5]-[H:3]",
        "reference": {
            "smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
            "reaction_and_mechanism": [
                "https://pubs.acs.org/doi/10.1021/ed048pA734.1",
                "https://pubs.acs.org/doi/10.1021/ed073pA312",
            ],
        },
        "comments": None,
    },
    "Di-Amine and Di-Carboxylic Acid Polycondensation (Polyamidation)": {
        "same_reactants": False,
        "reactant_1": "di_amine",
        "reactant_2": "di_carboxylic_acid",
        "product": "polyamide_chain",
        "delete_atom": True,
        "reaction": "[NX3;H2,H1;!$([N][C,S]=*):1]-[H:3].[CX3:2](=[O:4])[OX2H1:5]>>[NX3:1]-[CX3:2](=[O:4]).[O:5]-[H:3]",
        "reference": {
            "smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
            "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed048pA734"],
        },
        "comments": None,
    },
    "Di-Amine and Di-Carboxylic Acid Halide Polycondensation (Polyamidation)": {
        "same_reactants": False,
        "reactant_1": "di_amine",
        "reactant_2": "di_carboxylic_acid_halide",
        "product": "polyamide_chain",
        "delete_atom": True,
        "reaction": "[NX3;H2,H1;!$([N][C,S]=*):1]-[H:3].[CX3:2](=[O:4])[Cl,Br,I:5]>>[NX3:1]-[CX3:2](=[O:4]).[Cl,Br,I:5]-[H:3]",
        "reference": {
            "smarts": ["https://pubs.acs.org/doi/10.1021/acs.jcim.3c00329"],
            "reaction_and_mechanism": ["https://pubs.acs.org/doi/10.1021/ed048pA734"],
        },
        "comments": None,
    },
    # ============================================================
    # Mixed polyester / polythioester formation
    # ============================================================
    "Hydroxy-Thiol and Di-Carboxylic Acid Halide Polycondensation through Hydroxy Group": {
        "same_reactants": False,
        "reactant_1": "hydroxy_thiol",
        "reactant_2": "di_carboxylic_acid_halide",
        "product": "mixed_polyester_polythioester_chain",
        "delete_atom": True,
        "reaction": "[CX3:1](=[O:3])[Cl,Br,I:4].[OX2H1;!$([O][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[OX2:2].[Cl,Br,I:4]-[H:5]",
        "reference": {"smarts": None, "reaction_and_mechanism": None},
        "comments": None,
    },
    "Hydroxy-Thiol and Di-Carboxylic Acid Halide Polycondensation through Thiol Group": {
        "same_reactants": False,
        "reactant_1": "hydroxy_thiol",
        "reactant_2": "di_carboxylic_acid_halide",
        "product": "mixed_polyester_polythioester_chain",
        "delete_atom": True,
        "reaction": "[CX3:1](=[O:3])[Cl,Br,I:4].[SX2H1;!$([S][C,S]=*):2]-[H:5]>>[CX3:1](=[O:3])-[SX2:2].[Cl,Br,I:4]-[H:5]",
        "reference": {"smarts": None, "reaction_and_mechanism": None},
        "comments": None,
    },
    # ============================================================
    # Polyurethane formation (no byproduct)
    # ============================================================
    "Diol and Di-Isocyanate Polyaddition(Polyurethane Formation)": {
        "same_reactants": False,
        "reactant_1": "diol",
        "reactant_2": "di_isocyanate",
        "product": "polyurethane_chain",
        "delete_atom": False,
        "reaction": "[OX2H1;!$([O][C,S]=*):1]-[H:3].[NX2:4]=[CX2:2]=[OX1:5]>>[OX2:1]-[CX3:2](=[OX1:5])-[NX3:4]-[H:3]",
        "reference": {"smarts": None, "reaction_and_mechanism": None},
        "comments": None,
    },
}


def get_reaction_library(
    names: Optional[list] = None,
) -> Dict[str, Dict[str, Any]]:
    """
    Return the reaction library, optionally restricted to a name subset.

    Args:
        names: Reaction names to include; None returns the full library.

    Returns:
        Dict of reaction name -> reaction info.
    """
    if names is None:
        return dict(REACTIONS)
    return {n: REACTIONS[n] for n in names if n in REACTIONS}
