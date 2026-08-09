"""
Polymer Architecture Graph Core

This module provides the graph representation that underlies all bead-spring
polymer architectures. A chain is described as an explicit graph:

- nodes = beads (each carrying a bead-type name)
- edges = bonds (explicit index pairs)

Linear chains, rings, stars, combs/grafts, tadpoles, dendrimers and arbitrary
user-defined topologies are all the same kind of object; they differ only in
their edge list. Everything downstream (bond/angle enumeration, SAW growth
order, branched MC moves, data-file writing, mixtures) is derived from this
single representation.

Two levels of construction are provided:

1. ``MonomerTemplate`` — a reusable multi-bead monomer unit (e.g. one backbone
   bead + one side-group bead) with named connection points ("head", "tail",
   "side"...). Architectures are assembled by connecting monomer instances
   through their named connection points.
2. Factory functions (``linear``, ``ring``, ``star``, ``comb``, ``graft``,
   ``tadpole``, ``dendrimer``, ``custom``) — convenient constructors returning
   bead-level ``BeadArchitecture`` graphs directly.

Copolymer sequence generators (``block_sequence``, ``alternating_sequence``,
``random_sequence``, ``gradient_sequence``) produce bead-type sequences that
can be used in any sequence slot (backbone, arm, side chain, ...).
"""

from collections import deque
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from ..core.logger import setup_logger

logger = setup_logger()


# Type alias for sequence input formats:
#   "AABB"                  -> ["A", "A", "B", "B"]
#   ["A", "A", "B"]         -> explicit per-bead list
#   [("A", 20), ("B", 10)]  -> block format
SequenceInput = Union[str, List[str], List[Tuple[str, int]], Tuple]


def normalize_sequence(seq: SequenceInput) -> List[str]:
    """
    Convert any supported sequence input to an explicit bead-type list.

    Args:
        seq: Sequence as string ("AABB"), explicit list (["A", "B"]),
            block format ([("A", 20), ("B", 10)]), or a single bare block
            tuple (("A", 20)).

    Returns:
        Explicit list of bead type names.

    Raises:
        ValueError: If the sequence is empty or malformed.
    """
    # Bare block tuple: ("A", 20) -> ["A"] * 20
    if (isinstance(seq, tuple) and len(seq) == 2
            and isinstance(seq[0], str) and isinstance(seq[1], int)):
        seq = [seq]
    if isinstance(seq, str):
        result = list(seq)
    else:
        result = []
        for item in seq:
            if isinstance(item, tuple):
                bead_name, count = item
                if count < 1:
                    raise ValueError(f"Block count must be >= 1, got {count}")
                result.extend([bead_name] * count)
            else:
                result.append(item)
    if not result:
        raise ValueError("Sequence must contain at least one bead")
    return result


# =============================================================================
# Monomer Templates (explicit monomer-level structure)
# =============================================================================

@dataclass
class MonomerTemplate:
    """
    A reusable monomer unit: one or more beads with internal bonds and named
    connection points.

    Example — a graft monomer with one backbone bead and one side-group bead::

        MonomerTemplate(
            name="G",
            beads=["A", "B"],                 # bead 0 = backbone, bead 1 = side
            internal_bonds=[(0, 1)],
            connections={"head": 0, "tail": 0, "side": 1},
        )

    Attributes:
        name: Template name (used to reference it in instance lists).
        beads: Bead type names, one per internal bead.
        internal_bonds: Bonds between internal beads (local indices).
        connections: Named connection points mapping a name ("head", "tail",
            "side", ...) to a local bead index. Inter-monomer bonds are
            defined between connection points.
    """
    name: str
    beads: List[str]
    internal_bonds: List[Tuple[int, int]] = field(default_factory=list)
    connections: Dict[str, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        n = len(self.beads)
        if n == 0:
            raise ValueError(f"MonomerTemplate '{self.name}' has no beads")
        for i, j in self.internal_bonds:
            if not (0 <= i < n and 0 <= j < n):
                raise ValueError(
                    f"MonomerTemplate '{self.name}': internal bond ({i}, {j}) "
                    f"out of range for {n} beads"
                )
            if i == j:
                raise ValueError(
                    f"MonomerTemplate '{self.name}': self bond ({i}, {i})"
                )
        for conn_name, idx in self.connections.items():
            if not 0 <= idx < n:
                raise ValueError(
                    f"MonomerTemplate '{self.name}': connection '{conn_name}' "
                    f"points to bead {idx}, out of range for {n} beads"
                )

    @classmethod
    def single_bead(cls, bead_type: str, name: Optional[str] = None,
                    connections: Tuple[str, ...] = ("head", "tail")) -> "MonomerTemplate":
        """Convenience: a single-bead monomer with connection points on the bead."""
        return cls(
            name=name or bead_type,
            beads=[bead_type],
            internal_bonds=[],
            connections={c: 0 for c in connections},
        )


# =============================================================================
# Bead Architecture (the graph core)
# =============================================================================

@dataclass
class BeadArchitecture:
    """
    Explicit bead-level graph of one polymer chain.

    Attributes:
        bead_types: Bead type name per bead (node labels), length N.
        bonds: Explicit bond list as (i, j) bead index pairs, i != j.
        name: Architecture name ("linear", "ring", "star", "comb", ...).

    Everything else — angle triplets, SAW growth order, MC subtree structure —
    is derived from these two lists.
    """
    bead_types: List[str]
    bonds: List[Tuple[int, int]]
    name: str = "custom"

    # ------------------------------------------------------------------ #
    # Basic properties
    # ------------------------------------------------------------------ #

    @property
    def n_beads(self) -> int:
        return len(self.bead_types)

    @property
    def n_bonds(self) -> int:
        return len(self.bonds)

    def adjacency(self) -> List[List[int]]:
        """Adjacency list (sorted neighbor indices per bead)."""
        adj: List[List[int]] = [[] for _ in range(self.n_beads)]
        for i, j in self.bonds:
            adj[i].append(j)
            adj[j].append(i)
        for nbrs in adj:
            nbrs.sort()
        return adj

    def degrees(self) -> List[int]:
        return [len(nbrs) for nbrs in self.adjacency()]

    @property
    def max_degree(self) -> int:
        return max(self.degrees()) if self.n_beads else 0

    @property
    def is_branched(self) -> bool:
        """True if any bead has degree > 2 (star center, graft point, ...)."""
        return self.max_degree > 2

    @property
    def is_cyclic(self) -> bool:
        """True if the (connected) graph contains at least one cycle."""
        return self.n_bonds >= self.n_beads

    # ------------------------------------------------------------------ #
    # Validation
    # ------------------------------------------------------------------ #

    def validate(self, known_bead_types: Optional[Sequence[str]] = None) -> None:
        """
        Validate the graph.

        Checks: non-empty, bond indices in range, no duplicate/self bonds,
        connectedness, and (optionally) that all bead types are known.

        Raises:
            ValueError: On any violation.
        """
        n = self.n_beads
        if n == 0:
            raise ValueError("Architecture has no beads")

        seen = set()
        for i, j in self.bonds:
            if not (0 <= i < n and 0 <= j < n):
                raise ValueError(
                    f"Bond ({i}, {j}) out of range for {n} beads"
                )
            if i == j:
                raise ValueError(f"Self bond ({i}, {i}) is not allowed")
            key = (min(i, j), max(i, j))
            if key in seen:
                raise ValueError(f"Duplicate bond ({i}, {j})")
            seen.add(key)

        if known_bead_types is not None:
            known = set(known_bead_types)
            unknown = [t for t in self.bead_types if t not in known]
            if unknown:
                raise ValueError(
                    f"Unknown bead type(s) in architecture: {sorted(set(unknown))}. "
                    f"Known types: {sorted(known)}"
                )

        # Connectedness via BFS from bead 0
        adj = self.adjacency()
        visited = {0}
        queue = deque([0])
        while queue:
            u = queue.popleft()
            for v in adj[u]:
                if v not in visited:
                    visited.add(v)
                    queue.append(v)
        if len(visited) != n:
            missing = sorted(set(range(n)) - visited)
            raise ValueError(
                f"Architecture graph is not connected; unreachable beads: {missing}"
            )

    # ------------------------------------------------------------------ #
    # Angle triplets (branch-aware)
    # ------------------------------------------------------------------ #

    def angle_triplets(self, include_branch: bool = True) -> List[Tuple[int, int, int]]:
        """
        All angle triplets as (i, j, k) with j the center bead.

        Enumerated by center bead in node order; for each center, all
        neighbor pairs (i, k) with i < k. Triplets centered on a branch
        point (degree > 2) are included only if ``include_branch`` is True.

        For linear chains this yields exactly (i, i+1, i+2); for rings it
        additionally yields the wrap-around triplets — matching the legacy
        path-based enumeration.
        """
        return [t for t, _ in self.angle_triplets_classified(include_branch)]

    def angle_triplets_classified(
        self, include_branch: bool = True
    ) -> List[Tuple[Tuple[int, int, int], bool]]:
        """
        Like :meth:`angle_triplets` but returns ((i, j, k), is_branch) pairs,
        where is_branch marks triplets centered on a branch point (degree > 2).
        """
        adj = self.adjacency()
        result: List[Tuple[Tuple[int, int, int], bool]] = []
        for center in range(self.n_beads):
            nbrs = adj[center]
            is_branch = len(nbrs) > 2
            if is_branch and not include_branch:
                continue
            for a in range(len(nbrs)):
                for b in range(a + 1, len(nbrs)):
                    result.append(((nbrs[a], center, nbrs[b]), is_branch))
        return result

    # ------------------------------------------------------------------ #
    # Growth order (spanning tree) for SAW / geometric generation
    # ------------------------------------------------------------------ #

    def growth_order(
        self, root: int = 0
    ) -> Tuple[List[int], Dict[int, Optional[int]], List[Tuple[int, int]]]:
        """
        DFS spanning-tree growth order from ``root``.

        DFS (rather than BFS) keeps cycles as single long arcs with one
        closing edge at the far end — the same growth pattern as the
        legacy path-based ring SAW — which makes ring-closure constraints
        far easier to satisfy during placement. For trees (stars, combs,
        dendrimers) the spanning tree coincides with the full edge set
        either way.

        Returns:
            order: Bead indices in placement order (each bead's parent
                appears before it).
            parents: Map bead -> parent bead (None for root).
            closing_edges: Bonds not in the spanning tree (cycle-closing
                edges, e.g. the ring-closure bond). These are constraints
                to satisfy during placement rather than growth directions.
        """
        adj = self.adjacency()
        parents: Dict[int, Optional[int]] = {root: None}
        order: List[int] = []
        tree_edges = set()
        # Iterative DFS (preorder); parents are assigned at push time so
        # each node is pushed exactly once
        stack = [root]
        while stack:
            u = stack.pop()
            order.append(u)
            # Push neighbors in reverse so lower indices are visited first
            for v in reversed(adj[u]):
                if v not in parents:
                    parents[v] = u
                    tree_edges.add((min(u, v), max(u, v)))
                    stack.append(v)
        closing_edges = [
            (i, j) for i, j in self.bonds
            if (min(i, j), max(i, j)) not in tree_edges
        ]
        return order, parents, closing_edges

    # ------------------------------------------------------------------ #
    # Subtree structure for branched MC moves
    # ------------------------------------------------------------------ #

    def subtree_beyond_edge(self, u: int, v: int) -> List[int]:
        """
        Beads reachable from ``v`` without crossing the edge (u, v).

        Cutting edge (u, v) splits the molecule; this returns the component
        containing ``v``. Used by the tree-pivot MC move: rotating this set
        around any axis through bead ``u`` preserves the (u, v) bond length
        and all internal bonds of the subtree.
        """
        adj = self.adjacency()
        visited = {u, v}
        queue = deque([v])
        while queue:
            x = queue.popleft()
            for y in adj[x]:
                if y not in visited:
                    visited.add(y)
                    queue.append(y)
        visited.discard(u)
        return sorted(visited)

    def linear_segments(self, min_length: int = 3) -> List[List[int]]:
        """
        Maximal paths of degree-2 beads (plus their endpoints) suitable for
        crankshaft MC moves on branched architectures.

        Returns:
            List of bead-index paths, each of length >= min_length. For a
            linear chain this returns the whole chain as one segment.
        """
        adj = self.adjacency()
        deg = [len(nbrs) for nbrs in adj]
        n = self.n_beads
        segments: List[List[int]] = []
        visited_edges = set()

        def mark(a: int, b: int) -> None:
            visited_edges.add((min(a, b), max(a, b)))

        # Segments starting at non-degree-2 nodes (branch points, ends)
        starts = [i for i in range(n) if deg[i] != 2]
        for s in starts:
            for nbr in adj[s]:
                if (min(s, nbr), max(s, nbr)) in visited_edges:
                    continue
                path = [s, nbr]
                mark(s, nbr)
                while deg[path[-1]] == 2:
                    a, b = path[-2], path[-1]
                    nxt = adj[b][0] if adj[b][0] != a else adj[b][1]
                    if (min(b, nxt), max(b, nxt)) in visited_edges:
                        # cycle returning to an already-used edge
                        path.append(nxt)
                        mark(b, nxt)
                        break
                    path.append(nxt)
                    mark(b, nxt)
                    if nxt == s:
                        break
                if len(path) >= min_length:
                    segments.append(path)

        # Pure cycles with no branch/end beads (e.g. rings): pick them up
        if not starts and n >= min_length:
            # Whole graph is a single cycle
            path = [0]
            prev, cur = -1, 0
            while True:
                nxt_candidates = [x for x in adj[cur] if x != prev]
                if not nxt_candidates:
                    break
                nxt = nxt_candidates[0]
                if nxt == 0:
                    path.append(nxt)
                    break
                path.append(nxt)
                prev, cur = cur, nxt
            segments.append(path)

        return segments

    # ------------------------------------------------------------------ #
    # Assembly from monomer templates
    # ------------------------------------------------------------------ #

    @classmethod
    def from_monomers(
        cls,
        templates: Dict[str, MonomerTemplate],
        instances: List[str],
        inter_bonds: List[Tuple[int, str, int, str]],
        name: str = "custom",
    ) -> "BeadArchitecture":
        """
        Expand monomer instances into a bead-level architecture.

        Args:
            templates: Map template name -> MonomerTemplate.
            instances: Ordered list of template names (monomer instances).
            inter_bonds: Inter-monomer bonds as
                (instance_i, connection_i, instance_j, connection_j).
            name: Architecture name.

        Returns:
            BeadArchitecture with all internal + inter-monomer bonds.

        Example — three graft monomers polymerized head-to-tail::

            BeadArchitecture.from_monomers(
                templates={"G": graft_monomer},
                instances=["G", "G", "G"],
                inter_bonds=[(0, "tail", 1, "head"), (1, "tail", 2, "head")],
            )
        """
        bead_types: List[str] = []
        bonds: List[Tuple[int, int]] = []
        offsets: List[int] = []

        for inst_name in instances:
            if inst_name not in templates:
                raise ValueError(f"Unknown monomer template '{inst_name}'")
            tpl = templates[inst_name]
            offsets.append(len(bead_types))
            bead_types.extend(tpl.beads)
            for i, j in tpl.internal_bonds:
                bonds.append((offsets[-1] + i, offsets[-1] + j))

        for inst_i, conn_i, inst_j, conn_j in inter_bonds:
            if not (0 <= inst_i < len(instances) and 0 <= inst_j < len(instances)):
                raise ValueError(
                    f"Inter-monomer bond references instance "
                    f"({inst_i}, {inst_j}) but only {len(instances)} instances exist"
                )
            tpl_i = templates[instances[inst_i]]
            tpl_j = templates[instances[inst_j]]
            if conn_i not in tpl_i.connections:
                raise ValueError(
                    f"Template '{tpl_i.name}' has no connection '{conn_i}'"
                )
            if conn_j not in tpl_j.connections:
                raise ValueError(
                    f"Template '{tpl_j.name}' has no connection '{conn_j}'"
                )
            bonds.append((
                offsets[inst_i] + tpl_i.connections[conn_i],
                offsets[inst_j] + tpl_j.connections[conn_j],
            ))

        arch = cls(bead_types=bead_types, bonds=bonds, name=name)
        arch.validate()
        return arch


# =============================================================================
# Factory functions
# =============================================================================

def linear(sequence: SequenceInput, name: str = "linear") -> BeadArchitecture:
    """Linear chain: 0-1-2-...-(N-1)."""
    beads = normalize_sequence(sequence)
    bonds = [(i, i + 1) for i in range(len(beads) - 1)]
    return BeadArchitecture(bead_types=beads, bonds=bonds, name=name)


def ring(sequence: SequenceInput, name: str = "ring") -> BeadArchitecture:
    """Ring (cyclic) chain: linear bonds + closure bond (N-1, 0)."""
    beads = normalize_sequence(sequence)
    if len(beads) < 3:
        raise ValueError(f"Ring topology requires at least 3 beads, got {len(beads)}")
    bonds = [(i, i + 1) for i in range(len(beads) - 1)]
    bonds.append((len(beads) - 1, 0))
    return BeadArchitecture(bead_types=beads, bonds=bonds, name=name)


def star(
    center: str,
    arms: List[SequenceInput],
    name: str = "star",
) -> BeadArchitecture:
    """
    Star polymer: one center bead with f arms.

    Args:
        center: Bead type of the center bead.
        arms: List of arm sequences; arms may differ in length and
            composition (miktoarm / asymmetric stars).

    Bead numbering: 0 = center; arm a occupies consecutive indices in the
    order given.
    """
    if len(arms) < 2:
        raise ValueError(f"Star requires at least 2 arms, got {len(arms)}")
    bead_types = [center]
    bonds: List[Tuple[int, int]] = []
    for arm in arms:
        arm_beads = normalize_sequence(arm)
        start = len(bead_types)
        bead_types.extend(arm_beads)
        bonds.append((0, start))
        bonds.extend((start + k, start + k + 1) for k in range(len(arm_beads) - 1))
    arch = BeadArchitecture(bead_types=bead_types, bonds=bonds, name=name)
    arch.validate()
    return arch


def graft(
    backbone: SequenceInput,
    grafts: Dict[int, SequenceInput],
    name: str = "graft",
) -> BeadArchitecture:
    """
    Graft/comb polymer with explicitly placed side chains.

    Args:
        backbone: Backbone sequence.
        grafts: Map backbone bead index -> side-chain sequence. A side chain
            may be a single bead (a literal "side group", e.g. "B") or an
            oligomer (e.g. [("B", 5)]).

    Bead numbering: backbone beads 0..N-1, then side chains in order of
    sorted graft index.
    """
    bb = normalize_sequence(backbone)
    bead_types = list(bb)
    bonds: List[Tuple[int, int]] = [(i, i + 1) for i in range(len(bb) - 1)]

    for idx in sorted(grafts.keys()):
        if not 0 <= idx < len(bb):
            raise ValueError(
                f"Graft point {idx} out of range for backbone of {len(bb)} beads"
            )
        side = normalize_sequence(grafts[idx])
        start = len(bead_types)
        bead_types.extend(side)
        bonds.append((idx, start))
        bonds.extend((start + k, start + k + 1) for k in range(len(side) - 1))

    arch = BeadArchitecture(bead_types=bead_types, bonds=bonds, name=name)
    arch.validate()
    return arch


def comb(
    backbone: SequenceInput,
    side: SequenceInput,
    every: int,
    offset: int = 0,
    name: str = "comb",
) -> BeadArchitecture:
    """
    Regular comb polymer: identical side chains grafted every ``every``
    backbone beads, starting at ``offset``.

    Args:
        backbone: Backbone sequence.
        side: Side-chain sequence (single bead or oligomer).
        every: Grafting spacing along the backbone.
        offset: Backbone index of the first graft point.
    """
    bb = normalize_sequence(backbone)
    if every < 1:
        raise ValueError(f"'every' must be >= 1, got {every}")
    grafts = {i: side for i in range(offset, len(bb), every)}
    if not grafts:
        raise ValueError(
            f"No graft points: offset={offset}, every={every}, "
            f"backbone length={len(bb)}"
        )
    return graft(bb, grafts, name=name)


def tadpole(
    ring_seq: SequenceInput,
    tail: SequenceInput,
    attach: int = 0,
    name: str = "tadpole",
) -> BeadArchitecture:
    """
    Tadpole (lariat) polymer: a ring with a linear tail attached.

    Args:
        ring_seq: Ring sequence (beads 0..N-1, cyclically bonded).
        tail: Tail sequence, attached by its first bead to ring bead
            ``attach``.
        attach: Ring bead index where the tail attaches.
    """
    ring_beads = normalize_sequence(ring_seq)
    if len(ring_beads) < 3:
        raise ValueError(
            f"Tadpole ring requires at least 3 beads, got {len(ring_beads)}"
        )
    if not 0 <= attach < len(ring_beads):
        raise ValueError(
            f"Attachment index {attach} out of range for ring of {len(ring_beads)}"
        )
    tail_beads = normalize_sequence(tail)

    bead_types = ring_beads + tail_beads
    n_ring = len(ring_beads)
    bonds = [(i, i + 1) for i in range(n_ring - 1)]
    bonds.append((n_ring - 1, 0))  # ring closure
    tail_start = n_ring
    bonds.append((attach, tail_start))
    bonds.extend(
        (tail_start + k, tail_start + k + 1) for k in range(len(tail_beads) - 1)
    )
    arch = BeadArchitecture(bead_types=bead_types, bonds=bonds, name=name)
    arch.validate()
    return arch


def dendrimer(
    core: str,
    branch: str,
    branch_factor: int = 2,
    generations: int = 2,
    name: str = "dendrimer",
) -> BeadArchitecture:
    """
    Dendrimer / regularly branched polymer.

    Args:
        core: Bead type of the core bead.
        branch: Bead type of all branch beads.
        branch_factor: Number of children per node (core degree =
            branch_factor; internal branch beads have degree
            branch_factor + 1).
        generations: Number of branching shells. generations=1 gives a star
            with branch_factor arms of length 1.
    """
    if branch_factor < 2:
        raise ValueError(f"branch_factor must be >= 2, got {branch_factor}")
    if generations < 1:
        raise ValueError(f"generations must be >= 1, got {generations}")

    bead_types = [core]
    bonds: List[Tuple[int, int]] = []
    current_shell = [0]
    for _gen in range(generations):
        next_shell = []
        for parent in current_shell:
            for _ in range(branch_factor):
                child = len(bead_types)
                bead_types.append(branch)
                bonds.append((parent, child))
                next_shell.append(child)
        current_shell = next_shell
    arch = BeadArchitecture(bead_types=bead_types, bonds=bonds, name=name)
    arch.validate()
    return arch


def custom(
    bead_types: List[str],
    bonds: List[Tuple[int, int]],
    name: str = "custom",
) -> BeadArchitecture:
    """
    Arbitrary user-defined architecture from an explicit bead list and
    bond list (the escape hatch for topologies without a factory).
    """
    arch = BeadArchitecture(
        bead_types=list(bead_types),
        bonds=[(int(i), int(j)) for i, j in bonds],
        name=name,
    )
    arch.validate()
    return arch


# =============================================================================
# Copolymer sequence generators
# =============================================================================

def block_sequence(blocks: List[Tuple[str, int]]) -> List[str]:
    """Block copolymer sequence: [("A", 50), ("B", 50)] -> A50-b-B50."""
    return normalize_sequence(blocks)


def alternating_sequence(types: Sequence[str], n_beads: int) -> List[str]:
    """Alternating copolymer: types cycle along the chain ("ABABAB...")."""
    if not types:
        raise ValueError("alternating_sequence requires at least one bead type")
    if n_beads < 1:
        raise ValueError(f"n_beads must be >= 1, got {n_beads}")
    return [types[i % len(types)] for i in range(n_beads)]


def random_sequence(
    types: Sequence[str],
    n_beads: int,
    weights: Optional[Sequence[float]] = None,
    seed: Optional[int] = None,
) -> List[str]:
    """
    Statistical (random) copolymer sequence.

    Args:
        types: Bead type names to draw from.
        n_beads: Chain length.
        weights: Sampling probability per type (uniform if None).
        seed: Random seed for reproducibility.
    """
    if not types:
        raise ValueError("random_sequence requires at least one bead type")
    if n_beads < 1:
        raise ValueError(f"n_beads must be >= 1, got {n_beads}")
    if weights is not None:
        if len(weights) != len(types):
            raise ValueError("weights must have the same length as types")
        if any(w < 0 for w in weights):
            raise ValueError("weights must be non-negative")
        total = sum(weights)
        if total <= 0:
            raise ValueError("weights must sum to a positive value")
        probs = [w / total for w in weights]
    else:
        probs = None
    rng = np.random.default_rng(seed)
    draws = rng.choice(len(types), size=n_beads, p=probs)
    return [types[i] for i in draws]


def gradient_sequence(
    type_start: str,
    type_end: str,
    n_beads: int,
    seed: Optional[int] = None,
) -> List[str]:
    """
    Gradient copolymer: probability of ``type_end`` grows linearly from 0
    at the first bead to 1 at the last bead.
    """
    if n_beads < 2:
        raise ValueError(f"gradient_sequence requires n_beads >= 2, got {n_beads}")
    rng = np.random.default_rng(seed)
    seq = []
    for i in range(n_beads):
        p_end = i / (n_beads - 1)
        seq.append(type_end if rng.random() < p_end else type_start)
    return seq
