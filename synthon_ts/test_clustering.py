"""Focused tests for final-hit building-block clustering."""
from __future__ import annotations

import csv
import json
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from .clustering import (
    _HAVE_RDKIT,
    cluster_aware_mmr_top,
    cluster_building_blocks,
    write_building_block_cluster_plot,
    write_building_block_clusters,
)
from .core import Product


@unittest.skipUnless(_HAVE_RDKIT, "RDKit is not installed")
class BuildingBlockClusteringTests(unittest.TestCase):
    def test_cluster_aware_mmr_balances_signature_diversity_and_affinity(self):
        products = [
            Product("CCO", "rxn_a", [(0, "CCO"), (1, "c1ccccc1")], score=-10.0),
            Product("CCN", "rxn_a", [(0, "CCO"), (1, "c1ccccc1")], score=-9.9),
            Product("CCC", "rxn_a", [(0, "CCO"), (1, "C")], score=-9.5),
            Product("CCCl", "rxn_b", [(0, "CCO"), (1, "c1ccccc1")], score=-9.0),
            Product("CCBr", "rxn_a", [(0, "N"), (1, "C")], score=-8.0),
        ]

        selected, metadata = cluster_aware_mmr_top(
            products,
            3,
            similarity_threshold=0.60,
            diversity_weight=0.60,
            affinity_weight=0.40,
        )

        # The exact-signature duplicate CCN is collapsed. The different
        # reaction route and maximally different BB signature outrank CCC even
        # though CCC has the stronger docking score.
        self.assertEqual([p.smiles for p in selected], ["CCO", "CCCl", "CCBr"])
        self.assertEqual(metadata["unique_signature_count"], 4)
        self.assertEqual(metadata["diversity_weight"], 0.60)
        self.assertEqual(metadata["affinity_weight"], 0.40)
        self.assertEqual(metadata["selected"]["CCCl"]["diversity_score"], 1.0)

    def test_cluster_aware_mmr_fills_k_if_signatures_are_identical(self):
        products = [
            Product("CCO", "rxn", [(0, "CC")], score=-9.0),
            Product("CCN", "rxn", [(0, "CC")], score=-8.0),
        ]
        selected, metadata = cluster_aware_mmr_top(
            products, 2, similarity_threshold=0.60)
        self.assertEqual([p.smiles for p in selected], ["CCO", "CCN"])
        self.assertFalse(metadata["selected"]["CCO"]["fallback"])
        self.assertTrue(metadata["selected"]["CCN"]["fallback"])

    def test_top_n_position_groups_and_weighted_centroid(self):
        products = [
            Product("compound-worst", "rxn_a", [(0, "CCO"), (1, "N")], score=-5.0),
            Product("compound-best", "rxn_a", [(0, "CCCO"), (1, "N")], score=-10.0),
            Product("compound-second", "rxn_a", [(0, "CCCO"), (1, "N")], score=-9.0),
            Product("compound-other-rxn", "rxn_b", [(0, "CCO")], score=-8.0),
        ]

        report = cluster_building_blocks(
            products, top_n=3, similarity_threshold=0.0, morgan_radius=2)

        self.assertEqual(report["compounds_analyzed"], 3)
        self.assertEqual(
            [p["smiles"] for p in report["top_compounds"]],
            ["compound-best", "compound-second", "compound-other-rxn"],
        )
        positions = {(p["reaction_class"], p["bb_index"]): p
                     for p in report["positions"]}
        self.assertEqual(set(positions), {("rxn_a", 0), ("rxn_a", 1), ("rxn_b", 0)})
        # rxn_a@0 sees CCCO twice and excludes CCO from the worse fourth-ranked
        # rxn_a compound; its sole unique member is necessarily the centroid.
        cluster = positions[("rxn_a", 0)]["clusters"][0]
        self.assertEqual(cluster["centroid_smiles"], "CCCO")
        self.assertEqual(cluster["occurrences"], 2)
        self.assertEqual(cluster["supporting_compounds"], 2)

    def test_occurrence_weighted_medoid(self):
        products = [
            Product("a", "rxn", [(0, "CCO")], score=-7.0),
            Product("b", "rxn", [(0, "CCCO")], score=-8.0),
            Product("c", "rxn", [(0, "CCCO")], score=-9.0),
        ]
        report = cluster_building_blocks(products, top_n=100, similarity_threshold=0.0)
        cluster = report["positions"][0]["clusters"][0]
        self.assertEqual(cluster["centroid_smiles"], "CCCO")
        self.assertEqual(cluster["unique_building_blocks"], 2)

    def test_writes_complete_json_and_centroid_csv(self):
        products = [Product("a", "rxn", [(0, "c1ccccc1")], score=-7.0)]
        report = cluster_building_blocks(products)
        with tempfile.TemporaryDirectory() as tmp:
            json_path = Path(tmp) / "clusters.json"
            csv_path = Path(tmp) / "centroids.csv"
            svg_path = Path(tmp) / "clusters.svg"
            write_building_block_clusters(report, str(json_path), str(csv_path))
            write_building_block_cluster_plot(report, str(svg_path))

            loaded = json.loads(json_path.read_text())
            self.assertEqual(loaded["positions"][0]["cluster_count"], 1)
            with csv_path.open(newline="") as fh:
                rows = list(csv.DictReader(fh))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["centroid_smiles"], "c1ccccc1")
            svg = svg_path.read_text()
            self.assertIn("Building-block similarity clusters", svg)
            self.assertIn("c1ccccc1", svg)
            ET.parse(svg_path)


if __name__ == "__main__":
    unittest.main()
