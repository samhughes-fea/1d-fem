"""Adapters from other repositories' data to 1d-fem jobs.

The core (``pre_processing``, ``processing``, ``simulation_runner``, ``post_processing``)
knows beams, sections and laminates only. Code that reads a particular kind of structure,
such as a blade's geometry, lives here and is accepted duck typed, so ``src/onedfem``
imports none of the sibling repositories.
"""
