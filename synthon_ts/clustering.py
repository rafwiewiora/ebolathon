"""Chemical clustering of building blocks from the best-scoring products.

Building blocks are clustered independently for each ``(reaction_class,
bb_index)`` position.  Mixing positions would make the clusters less useful for
the reaction-aware search: two similar fragments are not interchangeable when
they belong to different reactions or coupling positions.

The implementation uses Butina clustering over Tanimoto distances between
Morgan fingerprints.  Because a bit-vector fingerprint has no molecule-valued
arithmetic centroid, the reported "centroid" is a medoid: the member with the
largest occurrence-weighted mean Tanimoto similarity to the cluster.
"""
from __future__ import annotations

import csv
import json
import os
from collections import defaultdict
from html import escape
from math import log2
from statistics import fmean

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.ML.Cluster import Butina

    _HAVE_RDKIT = True
except Exception:  # pragma: no cover - exercised only without the optional dep
    _HAVE_RDKIT = False


def _require_rdkit() -> None:
    if not _HAVE_RDKIT:
        raise RuntimeError(
            "building-block clustering requires RDKit; install it with `pip install rdkit`")


def _score_summary(scores: list[float]) -> dict:
    return {
        "best_score": min(scores),
        "mean_score": fmean(scores),
    }


def _product_building_block_fingerprints(product, generator) -> dict[int, object]:
    """Morgan fingerprints keyed by building-block position for one product."""
    fingerprints = {}
    for bb_index, bb_smiles in (getattr(product, "bbs", None) or []):
        mol = Chem.MolFromSmiles(bb_smiles)
        if mol is not None:
            fingerprints[int(bb_index)] = generator.GetFingerprint(mol)
    return fingerprints


def product_building_block_similarity(left, right, generator=None) -> float:
    """Route-aware similarity between two full products' building blocks.

    Building blocks are only compared when the products share a reaction class,
    and only position-matched blocks are compared.  The returned value is the
    mean Morgan/Tanimoto similarity across those positions.  Products from
    different reaction classes (or without comparable decompositions) have
    similarity zero for this selection purpose.
    """
    _require_rdkit()
    if not (getattr(left, "reaction_class", None)
            and left.reaction_class == getattr(right, "reaction_class", None)):
        return 0.0
    generator = generator or rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=2048)
    left_fps = _product_building_block_fingerprints(left, generator)
    right_fps = _product_building_block_fingerprints(right, generator)
    common_positions = sorted(set(left_fps) & set(right_fps))
    if not common_positions:
        return 0.0
    return fmean(
        DataStructs.TanimotoSimilarity(left_fps[index], right_fps[index])
        for index in common_positions
    )


def cluster_aware_mmr_top(
    ranked,
    k: int,
    similarity_threshold: float = 0.60,
    candidate_limit: int = 100,
    diversity_weight: float = 0.60,
    affinity_weight: float = 0.40,
) -> tuple[list, dict]:
    """Post-generation cluster-aware MMR selection of full products.

    Building blocks from the best ``candidate_limit`` products are first
    clustered per reaction class and position. Every attributed full product is
    assigned a signature ``(reaction_class, ((bb_index, cluster_id), ...))``.
    The best-affinity representative of each exact signature is retained, then
    greedy max-min selection balances signature diversity and normalized docking
    affinity. No part of this function is used by the generation/search loop.
    """
    _require_rdkit()
    if k <= 0:
        return [], {"selected": {}, "cluster_report": None}
    if not 0.0 <= similarity_threshold <= 1.0:
        raise ValueError("similarity_threshold must be between 0 and 1")
    if candidate_limit <= 0:
        raise ValueError("candidate_limit must be greater than zero")
    if diversity_weight < 0 or affinity_weight < 0:
        raise ValueError("selection weights must be non-negative")
    total_weight = diversity_weight + affinity_weight
    if total_weight <= 0:
        raise ValueError("at least one selection weight must be positive")
    diversity_weight /= total_weight
    affinity_weight /= total_weight

    candidates = sorted(
        (p for p in ranked if getattr(p, "score", None) is not None),
        key=lambda p: p.score,
    )[:candidate_limit]
    if not candidates:
        return [], {
            "method": "post-generation cluster-signature MMR",
            "diversity_weight": diversity_weight,
            "affinity_weight": affinity_weight,
            "candidate_count": 0,
            "unique_signature_count": 0,
            "selected": {},
            "cluster_report": None,
        }

    cluster_report = cluster_building_blocks(
        candidates,
        top_n=candidate_limit,
        similarity_threshold=similarity_threshold,
    )

    # product SMILES -> [(reaction_class, bb_index, cluster_id), ...]
    signature_parts: dict = defaultdict(list)
    for position in cluster_report["positions"]:
        reaction_class = position["reaction_class"]
        bb_index = int(position["bb_index"])
        for cluster in position["clusters"]:
            cluster_id = int(cluster["cluster_id"])
            for member in cluster["members"]:
                for product_smiles in member["supporting_compound_smiles"]:
                    signature_parts[product_smiles].append(
                        (reaction_class, bb_index, cluster_id))

    def signature_for(product):
        reaction_class = getattr(product, "reaction_class", None)
        if not reaction_class:
            return None
        positions = sorted(
            (bb_index, cluster_id)
            for rxn, bb_index, cluster_id in signature_parts.get(product.smiles, [])
            if rxn == reaction_class
        )
        return (reaction_class, tuple(positions)) if positions else None

    signature_by_id = {id(p): signature_for(p) for p in candidates}

    # Candidates are score-sorted, so the first product seen for an exact
    # signature is its strongest-affinity full-molecule representative.
    representatives = []
    represented_signatures = set()
    for candidate in candidates:
        signature = signature_by_id[id(candidate)]
        if signature is None or signature in represented_signatures:
            continue
        represented_signatures.add(signature)
        representatives.append(candidate)

    def signature_distance(left_signature, right_signature):
        left_rxn, left_positions = left_signature
        right_rxn, right_positions = right_signature
        if left_rxn != right_rxn:
            return 1.0
        left_map = dict(left_positions)
        right_map = dict(right_positions)
        indices = set(left_map) | set(right_map)
        if not indices:
            return 0.0
        same = sum(left_map.get(index) == right_map.get(index) for index in indices)
        return 1.0 - same / len(indices)

    best_score = min(p.score for p in representatives) if representatives else None
    worst_score = max(p.score for p in representatives) if representatives else None

    def affinity_value(product):
        if best_score is None or worst_score is None or worst_score == best_score:
            return 1.0
        value = (worst_score - product.score) / (worst_score - best_score)
        return max(0.0, min(1.0, value))

    selected = []
    details = {}
    remaining = list(representatives)
    while remaining and len(selected) < k:
        scored_rows = []
        for candidate in remaining:
            signature = signature_by_id[id(candidate)]
            diversity = min(
                (signature_distance(signature, signature_by_id[id(chosen)])
                 for chosen in selected),
                default=1.0,
            )
            affinity = affinity_value(candidate)
            selection_score = diversity_weight * diversity + affinity_weight * affinity
            scored_rows.append((
                -selection_score,
                -diversity,
                candidate.score,
                candidate.smiles,
                candidate,
                diversity,
                affinity,
                selection_score,
            ))
        # Round only for ordering so mathematically tied weighted scores are not
        # decided by binary floating-point noise. Prefer greater diversity,
        # then stronger affinity, for an exact 60/40 tie.
        chosen_row = min(
            scored_rows,
            key=lambda row: (round(row[0], 12), round(row[1], 12), row[2], row[3]),
        )
        chosen, diversity, affinity, selection_score = chosen_row[4:]
        selected.append(chosen)
        remaining.remove(chosen)
        reaction_class, positions = signature_by_id[id(chosen)]
        details[chosen.smiles] = {
            "cluster_signature": {
                "reaction_class": reaction_class,
                "positions": [
                    {"bb_index": bb_index, "cluster_id": cluster_id}
                    for bb_index, cluster_id in positions
                ],
            },
            "diversity_score": round(diversity, 6),
            "affinity_score": round(affinity, 6),
            "selection_score": round(selection_score, 6),
            "fallback": False,
        }

    # If fewer than k cluster signatures exist, fill by full-molecule affinity.
    # These rows are explicitly marked so the diversity guarantee is auditable.
    selected_ids = {id(p) for p in selected}
    for candidate in candidates:
        if len(selected) >= k:
            break
        if id(candidate) in selected_ids:
            continue
        signature = signature_by_id[id(candidate)]
        diversity = min(
            (signature_distance(signature, signature_by_id[id(chosen)])
             for chosen in selected
             if signature is not None and signature_by_id[id(chosen)] is not None),
            default=0.0,
        )
        affinity = affinity_value(candidate) if representatives else 1.0
        selection_score = diversity_weight * diversity + affinity_weight * affinity
        selected.append(candidate)
        selected_ids.add(id(candidate))
        details[candidate.smiles] = {
            "cluster_signature": None if signature is None else {
                "reaction_class": signature[0],
                "positions": [
                    {"bb_index": bb_index, "cluster_id": cluster_id}
                    for bb_index, cluster_id in signature[1]
                ],
            },
            "diversity_score": round(diversity, 6),
            "affinity_score": round(affinity, 6),
            "selection_score": round(selection_score, 6),
            "fallback": True,
        }

    return selected, {
        "method": "post-generation cluster-signature MMR",
        "diversity_weight": diversity_weight,
        "affinity_weight": affinity_weight,
        "bb_cluster_similarity_threshold": similarity_threshold,
        "candidate_limit": candidate_limit,
        "candidate_count": len(candidates),
        "attributed_candidate_count": sum(
            signature_by_id[id(p)] is not None for p in candidates),
        "unique_signature_count": len(represented_signatures),
        "selected": details,
        "cluster_report": cluster_report,
    }


def _weighted_medoid(
    indices: tuple[int, ...],
    fps: list,
    occurrences: list[int],
    best_scores: list[float],
    smiles: list[str],
) -> tuple[int, float]:
    """Return (member index, weighted mean similarity) for a cluster.

    Counts weight the calculation so a building block recurring in many of the
    top compounds has the same influence it would have in an undeduplicated
    list.  Including self-similarity is intentional: it also represents those
    observed occurrences and makes singleton centroids have similarity 1.0.
    """
    ranked = []
    denom = sum(occurrences[j] for j in indices)
    for i in indices:
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], [fps[j] for j in indices])
        mean_sim = sum(sim * occurrences[j] for sim, j in zip(sims, indices)) / denom
        ranked.append((
            -mean_sim,            # most central first
            -occurrences[i],      # then most frequently observed
            best_scores[i],       # then strongest docking support (lower wins)
            smiles[i],            # final deterministic tie breaker
            i,
            mean_sim,
        ))
    _, _, _, _, medoid_i, medoid_mean = min(ranked)
    return medoid_i, medoid_mean


def cluster_building_blocks(
    ranked,
    top_n: int = 100,
    similarity_threshold: float = 0.60,
    morgan_radius: int = 2,
    fingerprint_bits: int = 2048,
) -> dict:
    """Cluster building blocks attributed to the top-scoring compounds.

    Parameters
    ----------
    ranked:
        Iterable of ``Product``-like records with ``score``, ``smiles``,
        ``reaction_class`` and ``bbs`` attributes.  Input order is not trusted;
        compounds are sorted by ascending docking score.
    top_n:
        Maximum number of scored compounds to analyze.
    similarity_threshold:
        Minimum Morgan-fingerprint Tanimoto similarity used by Butina.  A value
        of 1.0 produces only exact-fingerprint groups; lower values merge a
        broader neighborhood.

    Returns a JSON-serializable report.  Invalid building-block SMILES are
    retained in ``invalid_building_blocks`` but excluded from clustering.
    """
    _require_rdkit()
    if top_n <= 0:
        raise ValueError("top_n must be greater than zero")
    if not 0.0 <= similarity_threshold <= 1.0:
        raise ValueError("similarity_threshold must be between 0 and 1")
    if morgan_radius < 0:
        raise ValueError("morgan_radius must be non-negative")
    if fingerprint_bits <= 0:
        raise ValueError("fingerprint_bits must be greater than zero")

    scored = [p for p in ranked if getattr(p, "score", None) is not None]
    selected = sorted(scored, key=lambda p: p.score)[:top_n]
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=morgan_radius, fpSize=fingerprint_bits)

    # position -> canonical BB -> observations.  A product may contribute one
    # observation at each reaction position; retaining observations lets the
    # report show both recurrence and score support.
    grouped: dict = defaultdict(lambda: defaultdict(lambda: {
        "scores": [], "compounds": set(), "mol": None,
    }))
    invalid = []
    attributed_compounds = set()
    for p in selected:
        if not (getattr(p, "reaction_class", None) and getattr(p, "bbs", None)):
            continue
        attributed_compounds.add(p.smiles)
        for bb_index, bb_smiles in p.bbs:
            mol = Chem.MolFromSmiles(bb_smiles)
            if mol is None:
                invalid.append({
                    "reaction_class": p.reaction_class,
                    "bb_index": bb_index,
                    "smiles": bb_smiles,
                    "compound_smiles": p.smiles,
                    "compound_score": p.score,
                })
                continue
            canonical = Chem.MolToSmiles(mol, canonical=True)
            row = grouped[(p.reaction_class, int(bb_index))][canonical]
            row["scores"].append(float(p.score))
            row["compounds"].add(p.smiles)
            row["mol"] = mol

    positions = []
    for (reaction_class, bb_index), bb_rows in sorted(
            grouped.items(), key=lambda item: (item[0][0], item[0][1])):
        # Recurrence-first ordering makes the Butina input and therefore cluster
        # assignment deterministic when several candidates have equal density.
        smiles = sorted(
            bb_rows,
            key=lambda smi: (-len(bb_rows[smi]["scores"]), smi),
        )
        fps = [generator.GetFingerprint(bb_rows[smi]["mol"]) for smi in smiles]
        occurrences = [len(bb_rows[smi]["scores"]) for smi in smiles]
        best_scores = [min(bb_rows[smi]["scores"]) for smi in smiles]

        distances = []
        for i in range(1, len(fps)):
            distances.extend(1.0 - sim for sim in
                             DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i]))
        raw_clusters = Butina.ClusterData(
            distances,
            len(fps),
            1.0 - similarity_threshold,
            isDistData=True,
            reordering=True,
        )

        clusters = []
        for member_indices in raw_clusters:
            centroid_i, centroid_mean = _weighted_medoid(
                tuple(member_indices), fps, occurrences, best_scores, smiles)

            member_rows = []
            cluster_scores = []
            cluster_compounds = set()
            for i in member_indices:
                smi = smiles[i]
                obs = bb_rows[smi]
                cluster_scores.extend(obs["scores"])
                cluster_compounds.update(obs["compounds"])
                member_rows.append({
                    "smiles": smi,
                    "similarity_to_centroid": DataStructs.TanimotoSimilarity(
                        fps[centroid_i], fps[i]),
                    "occurrences": len(obs["scores"]),
                    "supporting_compounds": len(obs["compounds"]),
                    "supporting_compound_smiles": sorted(obs["compounds"]),
                    **_score_summary(obs["scores"]),
                })
            member_rows.sort(key=lambda row: (
                -row["similarity_to_centroid"], -row["occurrences"], row["smiles"]))
            clusters.append({
                "centroid_smiles": smiles[centroid_i],
                "centroid_mean_similarity": centroid_mean,
                "unique_building_blocks": len(member_indices),
                "occurrences": len(cluster_scores),
                "supporting_compounds": len(cluster_compounds),
                "supporting_compound_smiles": sorted(cluster_compounds),
                **_score_summary(cluster_scores),
                "members": member_rows,
            })

        # Largest / most-supported chemical families first.  IDs are assigned
        # only after this sort so they remain stable and meaningful in the CSV.
        clusters.sort(key=lambda row: (
            -row["supporting_compounds"], -row["occurrences"],
            row["best_score"], row["centroid_smiles"]))
        for cluster_id, cluster in enumerate(clusters, start=1):
            cluster["cluster_id"] = cluster_id

        position_compounds = set().union(
            *(obs["compounds"] for obs in bb_rows.values()))
        positions.append({
            "reaction_class": reaction_class,
            "bb_index": bb_index,
            "unique_building_blocks": len(bb_rows),
            "occurrences": sum(occurrences),
            "supporting_compounds": len(position_compounds),
            "cluster_count": len(clusters),
            "clusters": clusters,
        })

    return {
        "method": "Butina clustering over Morgan-fingerprint Tanimoto distance",
        "centroid_definition": (
            "cluster member with maximum occurrence-weighted mean Tanimoto similarity"),
        "top_n_requested": top_n,
        "compounds_analyzed": len(selected),
        "top_compounds": [
            {"rank": rank, "smiles": p.smiles, "score": float(p.score)}
            for rank, p in enumerate(selected, start=1)
        ],
        "compounds_with_attributed_building_blocks": len(attributed_compounds),
        "similarity_threshold": similarity_threshold,
        "morgan_radius": morgan_radius,
        "fingerprint_bits": fingerprint_bits,
        "position_count": len(positions),
        "unique_building_blocks": sum(p["unique_building_blocks"] for p in positions),
        "invalid_building_blocks": invalid,
        "positions": positions,
    }


def write_building_block_clusters(report: dict, json_path: str, csv_path: str) -> None:
    """Write the complete cluster report plus a centroid-only CSV."""
    for path in (json_path, csv_path):
        parent = os.path.dirname(os.path.abspath(path))
        os.makedirs(parent, exist_ok=True)

    with open(json_path, "w") as fh:
        json.dump(report, fh, indent=2)

    fields = [
        "reaction_class", "bb_index", "cluster_id", "centroid_smiles",
        "centroid_mean_similarity", "unique_building_blocks", "occurrences",
        "supporting_compounds", "best_score", "mean_score",
    ]
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for position in report["positions"]:
            for cluster in position["clusters"]:
                writer.writerow({
                    "reaction_class": position["reaction_class"],
                    "bb_index": position["bb_index"],
                    **{field: cluster[field] for field in fields[2:]},
                })


def _mds_coordinates(
    smiles: list[str],
    radius: int,
    fingerprint_bits: int,
) -> list[tuple[float, float]]:
    """Classical 2D MDS coordinates from Morgan/Tanimoto distances."""
    try:
        import numpy as np
    except ImportError as exc:  # RDKit distributions normally include NumPy
        raise RuntimeError("cluster plotting requires NumPy") from exc

    if not smiles:
        return []
    if len(smiles) == 1:
        return [(0.0, 0.0)]

    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius, fpSize=fingerprint_bits)
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    if any(mol is None for mol in mols):
        raise ValueError("cluster report contains an invalid building-block SMILES")
    fps = [generator.GetFingerprint(mol) for mol in mols]

    n = len(fps)
    distances = np.zeros((n, n), dtype=float)
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        for j, similarity in enumerate(sims):
            distances[i, j] = distances[j, i] = 1.0 - similarity

    centering = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * centering @ (distances ** 2) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    dimensions = [i for i in np.argsort(eigenvalues)[::-1] if eigenvalues[i] > 1e-12][:2]
    coords = np.zeros((n, 2), dtype=float)
    for dimension, eigen_index in enumerate(dimensions):
        coords[:, dimension] = (
            eigenvectors[:, eigen_index] * np.sqrt(eigenvalues[eigen_index]))
        # Eigenvector signs are arbitrary. Orient each dimension deterministically
        # so the same report produces the same plot across runs.
        anchor = int(np.argmax(np.abs(coords[:, dimension])))
        if coords[anchor, dimension] < 0:
            coords[:, dimension] *= -1
    return [(float(x), float(y)) for x, y in coords]


def _plot_color(cluster_id: int) -> str:
    """Colorblind-friendly categorical colors; IDs remain visible as labels."""
    palette = (
        "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
        "#56B4E9", "#6A3D9A", "#B15928", "#1B9E77", "#7570B3",
    )
    return palette[(cluster_id - 1) % len(palette)]


def _short_smiles(smiles: str, max_chars: int = 52) -> str:
    if len(smiles) <= max_chars:
        return smiles
    return smiles[:max_chars - 1] + "…"


def write_building_block_cluster_plot(report: dict, svg_path: str) -> None:
    """Write a similarity-map SVG with cluster medoids highlighted.

    Each reaction position gets its own panel. Points are unique building blocks
    embedded by classical multidimensional scaling of Morgan/Tanimoto distances;
    circles are ordinary members and labeled diamonds are cluster medoids.
    """
    _require_rdkit()
    positions = report.get("positions") or []
    panel_width = 1120
    plot_left, plot_top, plot_width, plot_height = 72, 66, 650, 330
    legend_left = 760
    panel_heights = [max(440, 100 + 23 * p["cluster_count"]) for p in positions]
    total_height = sum(panel_heights) + 40 if positions else 220

    parent = os.path.dirname(os.path.abspath(svg_path))
    os.makedirs(parent, exist_ok=True)
    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (f'<svg xmlns="http://www.w3.org/2000/svg" width="{panel_width}" '
         f'height="{total_height}" viewBox="0 0 {panel_width} {total_height}" '
         'role="img" aria-labelledby="plot-title plot-desc">'),
        '<title id="plot-title">Building-block similarity clusters</title>',
        ('<desc id="plot-desc">Morgan fingerprint Tanimoto distance projected '
         'into two dimensions. Centroid medoids are shown as labeled diamonds.</desc>'),
        '<style>',
        '.background{fill:#FFFFFF}',
        'text{font-family:sans-serif;fill:#1F2937}',
        '.panel-title{font-size:17px;font-weight:500}.sub{font-size:12px;fill:#4B5563}',
        '.axis{stroke:#9CA3AF;stroke-width:1}.axis-label{font-size:11px;fill:#4B5563}',
        '.cluster-label{font-size:11px;font-weight:500}.legend{font-size:12px}',
        '.centroid{stroke:#1F2937}.member-point{stroke:#FFFFFF}',
        '@media(prefers-color-scheme:dark){.background{fill:#111827}text{fill:#F3F4F6}'
        '.sub,.axis-label{fill:#D1D5DB}.axis{stroke:#6B7280}.centroid{stroke:#F3F4F6}'
        '.member-point{stroke:#111827}}',
        '</style>',
        '<rect width="100%" height="100%" class="background"/>',
    ]

    if not positions:
        parts.extend([
            '<text x="560" y="100" text-anchor="middle" class="panel-title">'
            'No attributed building blocks to plot</text>',
            '<text x="560" y="130" text-anchor="middle" class="sub">'
            'The cluster report contains no reaction-position groups.</text>',
        ])

    panel_y = 20
    radius = int(report.get("morgan_radius", 2))
    fingerprint_bits = int(report.get("fingerprint_bits", 2048))
    for position, panel_height in zip(positions, panel_heights):
        members = []
        centroids = {}
        for cluster in position["clusters"]:
            cluster_id = int(cluster["cluster_id"])
            centroids[cluster_id] = cluster["centroid_smiles"]
            for member in cluster["members"]:
                members.append({**member, "cluster_id": cluster_id})
        members.sort(key=lambda row: (row["cluster_id"], row["smiles"]))
        smiles = [member["smiles"] for member in members]
        coords = _mds_coordinates(smiles, radius, fingerprint_bits)

        x_values = [xy[0] for xy in coords] or [0.0]
        y_values = [xy[1] for xy in coords] or [0.0]
        x_min, x_max = min(x_values), max(x_values)
        y_min, y_max = min(y_values), max(y_values)

        def scale(value, low, high, start, length):
            if abs(high - low) < 1e-12:
                return start + length / 2
            return start + (value - low) / (high - low) * length

        title = (f'{escape(str(position["reaction_class"]))} · building block '
                 f'{position["bb_index"]}')
        parts.extend([
            f'<g transform="translate(0,{panel_y})">',
            f'<text x="32" y="26" class="panel-title">{title}</text>',
            (f'<text x="32" y="47" class="sub">'
             f'{position["unique_building_blocks"]} unique building blocks · '
             f'{position["cluster_count"]} clusters · similarity threshold '
             f'{report["similarity_threshold"]:.2f}</text>'),
            (f'<line x1="{plot_left}" y1="{plot_top + plot_height}" '
             f'x2="{plot_left + plot_width}" y2="{plot_top + plot_height}" class="axis"/>'),
            (f'<line x1="{plot_left}" y1="{plot_top}" x2="{plot_left}" '
             f'y2="{plot_top + plot_height}" class="axis"/>'),
            (f'<text x="{plot_left + plot_width / 2}" y="{plot_top + plot_height + 30}" '
             'text-anchor="middle" class="axis-label">MDS 1 · Tanimoto distance</text>'),
            (f'<text x="20" y="{plot_top + plot_height / 2}" text-anchor="middle" '
             'transform="rotate(-90 20 '
             f'{plot_top + plot_height / 2})" class="axis-label">'
             'MDS 2 · Tanimoto distance</text>'),
        ])

        for member, (mds_x, mds_y) in zip(members, coords):
            x = scale(mds_x, x_min, x_max, plot_left + 18, plot_width - 36)
            y = scale(mds_y, y_min, y_max, plot_top + plot_height - 18, -(plot_height - 36))
            cluster_id = member["cluster_id"]
            color = _plot_color(cluster_id)
            tooltip = escape(
                f'C{cluster_id} | {member["smiles"]} | occurrences '
                f'{member["occurrences"]} | best score {member["best_score"]:.3f}')
            if member["smiles"] == centroids[cluster_id]:
                size = 8
                points = (f'{x:.1f},{y-size:.1f} {x+size:.1f},{y:.1f} '
                          f'{x:.1f},{y+size:.1f} {x-size:.1f},{y:.1f}')
                parts.append(
                    f'<polygon points="{points}" fill="{color}" class="centroid" '
                    f'stroke-width="1.5"><title>{tooltip} | centroid</title></polygon>')
                parts.append(
                    f'<text x="{x+10:.1f}" y="{y-8:.1f}" class="cluster-label">'
                    f'C{cluster_id}</text>')
            else:
                point_radius = min(8.0, 3.5 + log2(max(1, member["occurrences"])))
                parts.append(
                    f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{point_radius:.1f}" '
                    f'fill="{color}" fill-opacity="0.72" class="member-point" '
                    f'stroke-width="1"><title>{tooltip}</title></circle>')

        parts.append(
            f'<text x="{legend_left}" y="{plot_top - 18}" class="sub">'
            'Centroids (diamond)</text>')
        for row, cluster in enumerate(position["clusters"]):
            cluster_id = int(cluster["cluster_id"])
            color = _plot_color(cluster_id)
            y = plot_top + 8 + row * 23
            centroid = escape(_short_smiles(cluster["centroid_smiles"]))
            full_centroid = escape(cluster["centroid_smiles"])
            points = (f'{legend_left},{y-6} {legend_left+6},{y} '
                      f'{legend_left},{y+6} {legend_left-6},{y}')
            parts.extend([
                (f'<polygon points="{points}" fill="{color}" class="centroid" '
                 'stroke-width="1"/>'),
                (f'<text x="{legend_left+14}" y="{y+4}" class="legend">'
                 f'C{cluster_id} · n={cluster["unique_building_blocks"]} · '
                 f'{centroid}<title>{full_centroid}</title></text>'),
            ])
        parts.extend([
            (f'<text x="{plot_left}" y="{plot_top + plot_height + 52}" class="sub">'
             'Point size reflects recurrence; distance is approximate after 2D projection.</text>'),
            '</g>',
        ])
        panel_y += panel_height

    parts.append('</svg>')
    with open(svg_path, "w") as fh:
        fh.write("\n".join(parts) + "\n")
