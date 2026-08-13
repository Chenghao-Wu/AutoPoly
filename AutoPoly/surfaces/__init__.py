#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Built-in substrate surface builders.

Generates physical substrate slabs (hydroxylated crystalline silica:
alpha-quartz(0001) and beta-cristobalite(111)) as self-contained
moltemplate .lt classes for the "on_substrate" packing strategy.

Created on 2026-08-12
@author: zwu
"""
from .ff_tables import (
    CLAYFF,
    INTERFACE_FF,
    SLAB_FF_REGISTRY,
    SLAB_ROLES,
)
from .base import SlabBuilder, SurfaceSlab
from .quartz import QuartzBuilder, QuartzSlab
from .cristobalite import CristobaliteBuilder, CristobaliteSlab

#: Built-in slab builders by SubstrateSpec(builder=...) name.
SLAB_BUILDERS = {
    "alpha_quartz": QuartzBuilder,
    "beta_cristobalite": CristobaliteBuilder,
}

__all__ = [
    "CLAYFF",
    "INTERFACE_FF",
    "SLAB_FF_REGISTRY",
    "SLAB_ROLES",
    "SlabBuilder",
    "SurfaceSlab",
    "QuartzBuilder",
    "QuartzSlab",
    "CristobaliteBuilder",
    "CristobaliteSlab",
    "SLAB_BUILDERS",
]
