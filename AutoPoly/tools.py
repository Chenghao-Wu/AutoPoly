# -*- coding: utf-8 -*-
"""
Pre-built LangChain tool wrappers for AutoPoly.

Provides Pydantic schemas so LLMs see every field's name, type, and description.
LangChain is an optional dependency - imported lazily.

Usage:
    from AutoPoly.tools import get_autopoly_tools
    tools = get_autopoly_tools()  # Returns [info, generate_atomistic, generate_bead_spring, describe]
"""
import json
from typing import List, Optional

from . import agent


# ---------------------------------------------------------------------------
# Pydantic schemas (no LangChain dependency needed for these)
# ---------------------------------------------------------------------------

try:
    from pydantic import BaseModel, Field
except ImportError:
    BaseModel = None
    Field = None


def _define_schemas():
    """Define Pydantic schemas. Returns None if pydantic unavailable."""
    if BaseModel is None:
        return None

    class PolymerInput(BaseModel):
        """A polymer chain definition."""
        chain_num: int = Field(description="Number of polymer chains to generate")
        sequence: List[str] = Field(
            description="Complement SMILES with [*] wildcards. "
            "First monomer: 1 wildcard right (e.g. 'CC[*]'). "
            "Middle: 2 wildcards (e.g. '[*]CC[*]'). "
            "Last: 1 wildcard left (e.g. '[*]CC')."
        )
        topology: str = Field(default="linear", description="'linear' or 'ring'")
        tacticity: str = Field(default="atactic", description="'atactic', 'isotactic', or 'syndiotactic'")

    class MoleculeInput(BaseModel):
        """A small molecule definition."""
        count: int = Field(description="Number of molecules")
        smiles: str = Field(description="SMILES string WITHOUT wildcards (e.g. 'O' for water)")
        name: Optional[str] = Field(default=None, description="Molecule name (auto-generated if omitted)")

    class BeadTypeInput(BaseModel):
        """Bead type for coarse-grained model."""
        name: str = Field(description="Bead type name (e.g. 'A', 'B')")
        mass: float = Field(default=1.0, description="Bead mass")
        epsilon: float = Field(default=1.0, description="LJ epsilon parameter")
        sigma: float = Field(default=1.0, description="LJ sigma parameter")

    class AtomisticGenerateInput(BaseModel):
        """Config for atomistic polymer system generation."""
        name: str = Field(description="System name for output files")
        output_dir: str = Field(default="./output", description="Output directory")
        force_field: str = Field(
            default="oplsaa",
            description="Force field: 'oplsaa', 'lopls', 'gaff', 'gaff2', 'dreiding', or 'compass'"
        )
        polymers: Optional[List[PolymerInput]] = Field(
            default=None, description="Polymer definitions (need at least one of polymers/molecules)"
        )
        molecules: Optional[List[MoleculeInput]] = Field(
            default=None, description="Small molecule definitions"
        )
        placement_method: str = Field(default="mc_random", description="'grid' or 'mc_random'")

    class BeadSpringGenerateInput(BaseModel):
        """Config for bead-spring coarse-grained system."""
        name: str = Field(description="System name")
        output_dir: str = Field(default="./output", description="Output directory")
        n_chains: int = Field(description="Number of polymer chains")
        bead_types: List[BeadTypeInput] = Field(description="Bead type definitions")
        sequence: list = Field(
            description="Chain structure: list of [name, count] pairs (e.g. [['A', 50]]), "
            "or string (e.g. 'AABB'), or explicit list (e.g. ['A','A','B','B'])"
        )
        topology: str = Field(default="linear", description="'linear' or 'ring'")
        bond_style: str = Field(default="harmonic", description="'harmonic' or 'fene'")
        pair_style: str = Field(default="lj", description="'lj' (full LJ) or 'wca' (repulsive only)")
        generation_method: str = Field(default="saw", description="'geometric', 'saw', or 'mc'")
        density: Optional[float] = Field(default=None, description="Bead density (beads/sigma^3)")

    return {
        "AtomisticGenerateInput": AtomisticGenerateInput,
        "BeadSpringGenerateInput": BeadSpringGenerateInput,
    }


def get_autopoly_tools():
    """Return list of LangChain tools for DeepAgents integration.

    Requires: pip install langchain-core
    """
    try:
        from langchain_core.tools import tool
    except ImportError:
        raise ImportError(
            "LangChain tools require langchain-core. "
            "Install with: pip install 'AutoPoly[agent]' or pip install langchain-core"
        )

    schemas = _define_schemas()

    @tool
    def autopoly_info() -> str:
        """Discover all AutoPoly options: force fields, topologies, bond styles, limits, and
        complete working example configs. CALL THIS FIRST before generating."""
        return json.dumps(agent.info(), indent=2)

    atomistic_schema = schemas["AtomisticGenerateInput"] if schemas else None
    bead_spring_schema = schemas["BeadSpringGenerateInput"] if schemas else None

    if atomistic_schema:
        @tool(args_schema=atomistic_schema)
        def autopoly_generate_atomistic(
            name: str, output_dir: str = "./output", force_field: str = "oplsaa",
            polymers: list = None, molecules: list = None, placement_method: str = "mc_random",
        ) -> str:
            """Generate an atomistic polymer system with full chemistry for LAMMPS simulation.

            Example - polyethylene 10-mer, 2 chains:
              name="pe_system", polymers=[{"chain_num": 2,
                "sequence": ["CC[*]"] + ["[*]CC[*]"]*8 + ["[*]CC"]}]

            Example - polymer + solvent:
              name="pe_water", polymers=[{"chain_num": 1, "sequence": ["[*]CC[*]"]*5}],
              molecules=[{"count": 100, "smiles": "O", "name": "water"}]
            """
            config = {"type": "atomistic", "name": name, "output_dir": output_dir,
                      "force_field": force_field, "placement_method": placement_method}
            if polymers:
                config["polymers"] = [p.model_dump() if hasattr(p, "model_dump") else p for p in polymers]
            if molecules:
                config["molecules"] = [m.model_dump() if hasattr(m, "model_dump") else m for m in molecules]
            return str(agent.generate(config))
    else:
        @tool
        def autopoly_generate_atomistic(config: dict) -> str:
            """Generate an atomistic polymer system. Pass a dict with keys:
            name, output_dir, force_field, polymers, molecules, placement_method."""
            config["type"] = "atomistic"
            return str(agent.generate(config))

    if bead_spring_schema:
        @tool(args_schema=bead_spring_schema)
        def autopoly_generate_bead_spring(
            name: str, n_chains: int, bead_types: list, sequence: list,
            output_dir: str = "./output", topology: str = "linear",
            bond_style: str = "harmonic", pair_style: str = "lj",
            generation_method: str = "saw", density: float = None,
        ) -> str:
            """Generate a coarse-grained bead-spring polymer system for LAMMPS.

            Example - homopolymer melt:
              name="melt", n_chains=100, bead_types=[{"name": "A"}],
              sequence=[["A", 50]], bond_style="fene", pair_style="wca"

            Example - diblock copolymer:
              name="diblock", n_chains=50, bead_types=[{"name":"A"}, {"name":"B","epsilon":1.2}],
              sequence=[["A", 25], ["B", 25]]
            """
            config = {"type": "bead_spring", "name": name, "output_dir": output_dir,
                      "n_chains": n_chains, "sequence": sequence,
                      "topology": topology, "bond_style": bond_style, "pair_style": pair_style,
                      "generation_method": generation_method}
            config["bead_types"] = [bt.model_dump() if hasattr(bt, "model_dump") else bt for bt in bead_types]
            if density is not None:
                config["density"] = density
            return str(agent.generate(config))
    else:
        @tool
        def autopoly_generate_bead_spring(config: dict) -> str:
            """Generate a bead-spring polymer system. Pass a dict with keys:
            name, n_chains, bead_types, sequence, topology, bond_style, pair_style, etc."""
            config["type"] = "bead_spring"
            return str(agent.generate(config))

    @tool
    def autopoly_describe_smiles(smiles: str) -> str:
        """Get structural info about a SMILES string: atom count, elements, molecular weight,
        wildcard count. Use to verify SMILES before generating."""
        return json.dumps(agent.describe_smiles(smiles), indent=2)

    return [autopoly_info, autopoly_generate_atomistic,
            autopoly_generate_bead_spring, autopoly_describe_smiles]
