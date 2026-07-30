# generate

One-shot convenience function composing the three pipeline stages
(`GeometryBuilder` → `UnitTyper` → `BoxPacker`). Use the stage classes
directly when you need control between stages (e.g. typing one geometry
under multiple force fields); use `generate()` for single-force-field
one-shot runs.

::: AutoPoly.pipeline.workflow
