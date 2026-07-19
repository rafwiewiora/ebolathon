"""Synthon-aware Thompson-Sampling screen over onepot CORE.

Two interchangeable docking oracles share one loop (`core.run_loop`):
  - oracle_muni.MuniBatchDockingOracle   — docking via the muni CLI (muni credits)
  - oracle_rowan.RowanAnalogueDockingOracle — docking via the direct Rowan SDK

Retrieval is always the direct onepot REST API.
"""
from .core import LoopConfig, Target, run_loop, diverse_top  # noqa: F401
