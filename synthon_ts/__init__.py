"""Synthon-aware Thompson-Sampling screen over onepot CORE.

Two interchangeable docking oracles share one loop (`core.run_loop`):
  - oracle_muni.MuniBatchDockingOracle   — docking via the muni CLI (muni credits)
  - oracle_rowan.RowanAnalogueDockingOracle — docking via the direct Rowan SDK

Retrieval is always the direct onepot REST API.
"""
from .core import LoopConfig, Target, run_loop, diverse_top  # noqa: F401
from .clustering import (  # noqa: F401
    cluster_aware_mmr_top,
    cluster_building_blocks,
    product_building_block_similarity,
    write_building_block_cluster_plot,
    write_building_block_clusters,
)
