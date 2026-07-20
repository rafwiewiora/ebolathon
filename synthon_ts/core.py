"""Shared synthon-aware Thompson-Sampling loop over onepot CORE.

This is the oracle-agnostic engine used by BOTH the muni-backed runner
(`run_muni.py`) and the direct-API runner (`run_direct.py`). The only thing
that differs between the two versions is the *docking oracle* object you pass
in; everything below (onepot retrieval, the synthon bandit, loop-until-dry,
scaffold diversity) is identical.

Oracle contract
---------------
Pass any object exposing::

    oracle.score(smiles_list: list[str], template_smiles: str | None = None)
        -> dict[str, float]

where the returned value is a docking score and **lower is better** (Vina/
Vinardo convention). Missing/failed ligands may be omitted from the dict.

Retrieval is always the *direct* onepot REST API (`pip install onepot`), because
the muni `onepot` tool is broken (the `credits_remaining` gateway bug). Docking
is what varies per version.

onepot response shape (important — this drove the attribution redesign)
-----------------------------------------------------------------------
`Client.search(smiles_list, max_results=N, decompose=True)` returns::

    resp["queries"][i] = {
        "query_smiles": ..., "query_inchikey": ...,
        "results": [ {smiles, inchikey, similarity, price_usd, supplier_risk}, ... ],
        "decompositions": [ {reaction_class: "rxn_...",
                             bbs:[{bb_index:0, smiles:...}, {bb_index:1, smiles:...}]}, ... ],
    }

Crucially, onepot decomposes only the **query** you searched with — the
`results` analogs carry NO synthons. So to credit a docked analog to its own
building blocks we must decompose the analog itself. Two attribution modes:

* **precise** (default, `precise_attribution=True`): after docking a batch of
  analogs, issue ONE `search(analog_list, decompose=True, max_results=1)` call
  and read each query's `decompositions` — i.e. every docked analog's OWN
  synthons — then credit that analog's reward to its own
  `(reaction_class, bb_index, bb_smiles)` arms. Costs ~1 onepot credit per
  decomposed analog (the approved tradeoff). Analogs that don't decompose keep
  their docking score in the pool but are not attributed.
* **coarse** (`precise_attribution=False`): the cheap anchor-based
  approximation. Credit each analog to the ANCHOR's synthons (the anchor is the
  query that retrieved it), weighting full credit to the exploited position (the
  `bb_filters` window guarantees the analog shares that synthon with the anchor)
  and little/none to the varying explored position. The seed round (no window,
  all positions free) gives uniform small credit to each query synthon. No extra
  onepot calls.
"""
from __future__ import annotations

import json
import math
import os
import random
from collections import defaultdict
from dataclasses import dataclass, field

from onepot import Client

try:  # rdkit is strongly recommended but the loop degrades gracefully without it
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    _HAVE_RDKIT = True
except Exception:  # pragma: no cover
    _HAVE_RDKIT = False


# --------------------------------------------------------------------------- #
# Config
# --------------------------------------------------------------------------- #
@dataclass
class Target:
    """The docking target + box. `protein` is a PDB id (e.g. '1HCK') or a
    backend-specific prepared-protein handle (a muni/Rowan protein UUID)."""
    protein: str
    pocket: list  # [[cx, cy, cz], [sx, sy, sz]]
    executable: str = "qvina2"
    scoring_function: str = "vina"
    exhaustiveness: int = 8


@dataclass
class LoopConfig:
    query_smiles: str
    # Optional list of seed SMILES for the seed round. When set, EVERY seed is
    # searched + docked (used for sample_space multi-seed starts). When None the
    # single `query_smiles` is used, so existing single-query behaviour is intact.
    seed_smiles: list | None = None
    seed_max_results: int = 60      # onepot hits for the seed round
    round_max_results: int = 40     # onepot hits per elaboration query
    exploit_min_sim: float = 0.75   # tight window on the exploited position
    explore_min_sim: float = 0.10   # loose floor on the explored position
    explore_max_sim: float = 0.60   # capped ceiling on the explored position
    n_anchors_per_round: int = 3    # elaboration queries submitted per round
    max_rounds: int = 6
    patience: int = 2               # loop-until-dry: rounds w/o improvement
    improve_eps: float = 0.15       # min top-k mean reward gain to count as "improved"
    max_docks: int = 600            # hard budget on ligand evaluations
    top_k: int = 20
    random_seed: int = 0
    # Attribution: True = decompose each docked analog and credit its OWN
    # synthons (exact, +~1 onepot credit/analog). False = coarse anchor-based
    # credit (no extra onepot calls). See module docstring.
    precise_attribution: bool = True
    seed_credit: float = 0.25       # coarse mode: weak uniform credit per seed synthon
    explore_credit: float = 0.0     # coarse mode: credit to the explored (varying) synthon
    # Cost / supplier constraints applied NATIVELY on every onepot retrieval
    # search (seed + elaboration). None = unconstrained. onepot prices are two
    # tiers ($125 / $295), so max_price=200 keeps only the $125 tier.
    max_price: int | None = None
    max_supplier_risk: str | None = None  # "low" | "medium" | "high"


# --------------------------------------------------------------------------- #
# Small chem helpers (rdkit optional)
# --------------------------------------------------------------------------- #
def canon(smiles: str) -> str:
    if not _HAVE_RDKIT:
        return smiles
    m = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(m) if m else smiles


def murcko(smiles: str) -> str:
    if not _HAVE_RDKIT:
        return smiles
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(smiles=smiles)
    except Exception:
        return smiles


# --------------------------------------------------------------------------- #
# Data records
# --------------------------------------------------------------------------- #
@dataclass
class Decomp:
    """One retrosynthetic route of a molecule: a reaction class + the building
    blocks at each indexed position."""
    reaction_class: str
    bbs: list  # [(bb_index, bb_smiles), ...]


@dataclass
class Product:
    """A docked analog.

    `reaction_class` / `bbs` are the **attribution source** used to credit this
    analog's reward to synthon arms and to elaborate it:

    * precise mode -> this analog's OWN decomposition (recovered by decomposing
      the analog itself), so credit is exact.
    * coarse mode  -> the ANCHOR's decomposition (the query that retrieved this
      analog); a deliberate approximation (onepot never returns an analog's own
      synthons in a similarity search).

    `bb_weights` maps bb_index -> credit weight applied when the bandit observes
    this product (precise: 1.0 on every own-synthon; coarse: 1.0 on the exploited
    position, `explore_credit`/0 elsewhere; seed coarse: `seed_credit` uniform).
    """
    smiles: str
    reaction_class: str | None = None
    bbs: list | None = None          # [(bb_index, bb_smiles), ...] attribution source
    bb_weights: dict | None = None   # bb_index -> credit weight
    similarity: float | None = None  # analog<->query similarity from onepot
    score: float | None = None       # docking score (lower better); None until scored

    @property
    def reward(self) -> float | None:
        return None if self.score is None else -self.score


def parse_query(q: dict):
    """Split one onepot query block into (analog dicts, list[Decomp]).

    Analogs carry only {smiles, similarity, price_usd}; the decompositions
    belong to the QUERY, not to the analogs."""
    analogs = [{"smiles": h["smiles"], "similarity": h.get("similarity"),
                "price_usd": h.get("price_usd")}
               for h in q.get("results", [])]
    decomps = []
    for d in q.get("decompositions", []) or []:
        bbs = [(b["bb_index"], b["smiles"]) for b in d.get("bbs", [])]
        if bbs:
            decomps.append(Decomp(d["reaction_class"], bbs))
    return analogs, decomps


def search_one(client, smiles: str, max_results: int, bb_filters=None,
               max_price=None, max_supplier_risk=None):
    """Run a single-query onepot search; return (analogs, decomps) for it.

    `max_price` / `max_supplier_risk` are onepot-native filters, so the API only
    returns cheap, low-risk analogs (applied to every retrieval search)."""
    kw = {}
    if max_price is not None:
        kw["max_price"] = max_price
    if max_supplier_risk is not None:
        kw["max_supplier_risk"] = max_supplier_risk
    resp = client.search([smiles], max_results=max_results, decompose=True,
                         bb_filters=bb_filters, **kw)
    q = (resp.get("queries") or [{}])[0]
    return parse_query(q)


# --------------------------------------------------------------------------- #
# Synthon bandit
# --------------------------------------------------------------------------- #
class SynthonBandit:
    """Per-(reaction_class, bb_index, bb_smiles) reward posteriors + Thompson
    sampling over which anchor to elaborate and which position to exploit.

    Reward = -docking_score (higher is better). Observations are *weighted* (the
    attribution scheme decides how much of an analog's reward a given synthon
    earns). Docking is ~deterministic, so a synthon's *uncertainty* is really
    about its neighbourhood: little effective evidence => sample around it more
    aggressively. We model each synthon's mean reward as Gaussian with
    posterior-of-the-mean std ~ sigma0 / sqrt(n_eff), n_eff = sum of weights."""

    def __init__(self, rng: random.Random):
        self.rng = rng
        self.obs: dict = defaultdict(list)   # (rxn, idx, bb) -> [(reward, weight), ...]
        self._all_rewards: list = []

    # -- updates -----------------------------------------------------------
    def observe(self, p: Product):
        if p.reward is None:
            return
        self._all_rewards.append(p.reward)
        if p.reaction_class and p.bbs:
            weights = p.bb_weights or {}
            for idx, bb in p.bbs:
                w = weights.get(idx, 1.0)
                if w > 0:
                    self.obs[(p.reaction_class, idx, bb)].append((p.reward, float(w)))

    # -- cross-run persistence --------------------------------------------
    def load(self, path: str):
        """Seed this bandit with synthon posteriors accumulated by PRIOR runs, so
        a new campaign starts already knowing which synthons dock well (priors),
        instead of relearning from scratch. Safe if the file is missing/corrupt."""
        if not path or not os.path.exists(path):
            return 0
        try:
            data = json.load(open(path))
        except Exception:
            return 0
        n = 0
        for rxn, idx, bb, rows in data.get("obs", []):
            self.obs[(rxn, int(idx), bb)].extend((float(r), float(w)) for r, w in rows)
            n += 1
        self._all_rewards.extend(float(r) for r in data.get("all_rewards", []))
        return n

    def save(self, path: str):
        """Persist the accumulated synthon posteriors (this run's obs merged with
        whatever was loaded) so the NEXT run compounds on them."""
        data = {"obs": [[k[0], k[1], k[2], [[r, w] for r, w in rows]]
                        for k, rows in self.obs.items()],
                "all_rewards": self._all_rewards}
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(data, fh)
        os.replace(tmp, path)   # atomic-ish so a concurrent reader never sees half

    def _wstats(self, key):
        """(weighted_mean, effective_n) for a synthon arm, or None if unseen."""
        rows = self.obs.get(key)
        if not rows:
            return None
        sw = sum(w for _, w in rows)
        if sw <= 0:
            return None
        m = sum(r * w for r, w in rows) / sw
        return m, sw

    @property
    def sigma0(self) -> float:
        if len(self._all_rewards) < 2:
            return 1.0
        m = sum(self._all_rewards) / len(self._all_rewards)
        var = sum((r - m) ** 2 for r in self._all_rewards) / (len(self._all_rewards) - 1)
        return max(math.sqrt(var), 1e-6)

    # -- Thompson samples --------------------------------------------------
    def ts_synthon(self, key) -> float:
        """Thompson sample of a synthon's mean reward. Unseen synthons get an
        optimistic prior (global mean + one sigma) so exploration is favoured."""
        st = self._wstats(key)
        if st is None:
            base = (sum(self._all_rewards) / len(self._all_rewards)) if self._all_rewards else 0.0
            return base + self.sigma0
        m, n_eff = st
        return m + self.rng.gauss(0.0, self.sigma0 / math.sqrt(max(n_eff, 1e-6)))

    def anchor_value(self, p: Product) -> float:
        """TS value of elaborating around a scored, decomposed product: the max
        Thompson sample over its synthons (best position drives the neighbourhood)."""
        if not (p.reaction_class and p.bbs):
            return p.reward if p.reward is not None else -1e9
        return max(self.ts_synthon((p.reaction_class, idx, bb)) for idx, bb in p.bbs)

    def positions(self, p: Product):
        """(exploit_idx, explore_idx) for an anchor: exploit = highest posterior
        mean synthon (keep it, tighten); explore = least-observed (rotate)."""
        assert p.reaction_class and p.bbs
        def mean(idx, bb):
            st = self._wstats((p.reaction_class, idx, bb))
            return st[0] if st else -1e9
        def count(idx, bb):
            st = self._wstats((p.reaction_class, idx, bb))
            return st[1] if st else 0.0
        exploit = max(p.bbs, key=lambda t: mean(*t))[0]
        others = [t for t in p.bbs if t[0] != exploit] or p.bbs
        explore = min(others, key=lambda t: count(*t))[0]
        return exploit, explore


def build_bb_filters(reaction_class: str, exploit_idx: int, explore_idx: int,
                     cfg: LoopConfig):
    """Windows anchored to the anchor's OWN bbs: tighten the exploited position,
    open (but cap) the explored position. Other positions left unconstrained.

    onepot computes the Tanimoto of each candidate BB against the *query's* BB at
    that (reaction_class, bb_index), so we only need the class + index + bounds —
    not the BB SMILES."""
    filters = [{"reaction_class": reaction_class, "bb_index": exploit_idx,
                "min_similarity": cfg.exploit_min_sim}]
    if explore_idx != exploit_idx:
        filters.append({"reaction_class": reaction_class, "bb_index": explore_idx,
                        "min_similarity": cfg.explore_min_sim,
                        "max_similarity": cfg.explore_max_sim})
    return filters


# --------------------------------------------------------------------------- #
# Attribution
# --------------------------------------------------------------------------- #
def _attribute_coarse(products, decomp: Decomp | None, weights: dict):
    """Coarse mode: stamp the anchor's decomposition + credit weights onto each
    analog (no onepot call)."""
    for p in products:
        if decomp is not None:
            p.reaction_class = decomp.reaction_class
            p.bbs = list(decomp.bbs)
            p.bb_weights = dict(weights)


def _attribute_precise(client, products, log):
    """Precise mode: ONE batched `search(analogs, decompose=True, max_results=1)`
    recovers each docked analog's OWN synthons; credit its reward to them
    (uniform full weight). Order-aligned with the input list. Analogs that don't
    decompose are left unattributed (score still kept in the pool)."""
    if not products:
        return
    smis = [p.smiles for p in products]
    resp = client.search(smis, max_results=1, decompose=True)
    log(f"[attrib] decomposed {len(smis)} analogs "
        f"(onepot credits_used={resp.get('credits_used')})")
    queries = resp.get("queries") or []
    n_ok = 0
    for p, q in zip(products, queries):
        _analogs, decomps = parse_query(q)
        if not decomps:
            continue
        d = decomps[0]   # primary route: this analog's own synthons
        p.reaction_class = d.reaction_class
        p.bbs = list(d.bbs)
        p.bb_weights = {idx: 1.0 for idx, _ in d.bbs}
        n_ok += 1
    log(f"[attrib] {n_ok}/{len(products)} analogs attributed to their own synthons")


# --------------------------------------------------------------------------- #
# The loop
# --------------------------------------------------------------------------- #
def run_loop(oracle, target: Target, cfg: LoopConfig, onepot_key: str | None = None,
             onepot_base_url: str = "https://api.onepot.ai", client=None,
             log=print, round_callback=None, mol_filter=None,
             belief_store: str | None = None) -> dict:
    """Drive the synthon-TS screen. `oracle.score(smiles, template)` does the
    docking (lower better). Returns the scored pool + ranked hits.

    `client` may be any object with a `.search(smiles_list, max_results,
    decompose, bb_filters=...)` method returning the onepot response dict; when
    None a real onepot `Client` is built from `onepot_key`. Injecting a client
    lets you dry-run the whole loop offline.

    `round_callback`, if given, is called after each elaboration round as
    `round_callback(ranked_products, round_index)` where `ranked_products` is the
    current scored pool sorted best-first (lower score first). It lets a caller
    export docked poses live as the screen progresses. Backward compatible:
    when None (default) behaviour is unchanged.

    `mol_filter`, if given, is a `smiles -> bool` predicate applied to every
    onepot analog BEFORE docking (both the seed round and elaboration rounds), so
    off-profile molecules are dropped before spending a dock. The number dropped
    is logged per search. When None (default) no analogs are filtered."""
    rng = random.Random(cfg.random_seed)
    if client is None:
        client = Client(api_key=onepot_key, base_url=onepot_base_url)
    bandit = SynthonBandit(rng)
    if belief_store:          # start smart: load synthon posteriors from prior runs
        n = bandit.load(belief_store)
        log(f"[belief-store] loaded {n} synthon arms + {len(bandit._all_rewards)} "
            f"prior observations from {belief_store}")
    pool: dict = {}           # canonical smiles -> Product (scored)
    seen: set = set()         # canonical smiles ever retrieved
    n_docks = 0

    def dock(products: list[Product], template: str | None) -> list[Product]:
        """Dock the fresh (unseen, in-budget) subset; return scored products
        (NOT yet attributed to the bandit)."""
        nonlocal n_docks
        fresh = []
        for p in products:
            c = canon(p.smiles)
            if c in seen:
                continue
            seen.add(c)
            fresh.append(p)
        if not fresh:
            return []
        budget = max(0, cfg.max_docks - n_docks)
        if budget == 0:
            return []
        fresh = fresh[:budget]
        scores = oracle.score([p.smiles for p in fresh], template_smiles=template)  # lower better
        n_docks += len(fresh)
        scored = []
        for p in fresh:
            s = scores.get(p.smiles)
            if s is None:
                continue
            p.score = float(s)
            pool[canon(p.smiles)] = p
            scored.append(p)
        return scored

    def filter_analogs(analogs, tag):
        """Drop analogs failing `mol_filter` (drug-like window) before docking."""
        if mol_filter is None:
            return analogs
        kept = [a for a in analogs if mol_filter(a["smiles"])]
        dropped = len(analogs) - len(kept)
        if dropped:
            log(f"[{tag}] druglike filter dropped {dropped}/{len(analogs)} analogs")
        return kept

    def attribute_and_observe(scored, coarse_decomp=None, coarse_weights=None):
        if not scored:
            return
        if cfg.precise_attribution:
            _attribute_precise(client, scored, log)
        else:
            _attribute_coarse(scored, coarse_decomp, coarse_weights or {})
        for p in scored:
            bandit.observe(p)

    # -- seed round (one or many seed SMILES) ------------------------------
    seeds = [s for s in (cfg.seed_smiles or [cfg.query_smiles]) if s]
    log(f"[seed] {len(seeds)} seed(s); max_results={cfg.seed_max_results} each")
    for si, seed_smi in enumerate(seeds):
        log(f"[seed {si}] onepot search: {seed_smi!r}")
        analogs, seed_decomps = search_one(
            client, seed_smi, cfg.seed_max_results,
            max_price=cfg.max_price, max_supplier_risk=cfg.max_supplier_risk)
        analogs = filter_analogs(analogs, f"seed {si}")
        seed_decomp = seed_decomps[0] if seed_decomps else None
        log(f"[seed {si}] {len(analogs)} analogs (post-filter); decomposes to "
            f"{seed_decomp.reaction_class if seed_decomp else 'None'}; docking...")
        products = [Product(smiles=a["smiles"], similarity=a["similarity"]) for a in analogs]
        scored = dock(products, template=seed_smi)
        # coarse: uniform weak credit to each of this seed's synthons
        seed_weights = ({idx: cfg.seed_credit for idx, _ in seed_decomp.bbs}
                        if seed_decomp else {})
        attribute_and_observe(scored, coarse_decomp=seed_decomp, coarse_weights=seed_weights)
    log(f"[seed] scored pool={len(pool)}  docks={n_docks}  synthon arms={len(bandit.obs)}")

    # -- elaboration rounds (loop-until-dry) -------------------------------
    def topk_mean():
        best = sorted((p.score for p in pool.values()))[:cfg.top_k]
        return (sum(best) / len(best)) if best else 0.0

    best_metric = topk_mean()
    stale = 0
    for rnd in range(1, cfg.max_rounds + 1):
        if n_docks >= cfg.max_docks:
            log(f"[round {rnd}] budget exhausted ({n_docks}/{cfg.max_docks})"); break

        # choose anchors: decomposed, scored, Thompson-ranked, scaffold-diverse
        cands = [p for p in pool.values() if p.reaction_class and p.bbs]
        cands.sort(key=bandit.anchor_value, reverse=True)
        anchors, used_scaffolds = [], set()
        for p in cands:
            sc = murcko(p.smiles)
            if sc in used_scaffolds:
                continue
            used_scaffolds.add(sc)
            anchors.append(p)
            if len(anchors) >= cfg.n_anchors_per_round:
                break
        if not anchors:
            log(f"[round {rnd}] no decomposed anchors; stop"); break

        got = 0
        for a in anchors:
            ex, xp = bandit.positions(a)
            filters = build_bb_filters(a.reaction_class, ex, xp, cfg)
            log(f"[round {rnd}] anchor score={a.score:.2f} rxn={a.reaction_class} "
                f"exploit@{ex} explore@{xp}")
            try:
                analogs, _ = search_one(client, a.smiles, cfg.round_max_results,
                                        bb_filters=filters, max_price=cfg.max_price,
                                        max_supplier_risk=cfg.max_supplier_risk)
            except Exception as e:
                log(f"[round {rnd}]   onepot error: {str(e)[:120]}"); continue
            analogs = filter_analogs(analogs, f"round {rnd}")
            products = [Product(smiles=al["smiles"], similarity=al["similarity"])
                        for al in analogs]
            scored = dock(products, template=a.smiles)
            # coarse: full credit to the exploited synthon, little/none to explored
            cw = {ex: 1.0, xp: cfg.explore_credit} if a.bbs else {}
            attribute_and_observe(scored, coarse_decomp=a, coarse_weights=cw)
            got += len(scored)
            log(f"[round {rnd}]   +{len(scored)} scored (pool={len(pool)}, docks={n_docks}, "
                f"arms={len(bandit.obs)})")

        if round_callback is not None:
            try:
                ranked_now = sorted(pool.values(), key=lambda p: p.score)
                round_callback(ranked_now, rnd)
            except Exception as e:  # never let pose export break the screen
                log(f"[round {rnd}] round_callback error: {str(e)[:160]}")

        metric = topk_mean()
        improved = (metric - best_metric) < -cfg.improve_eps  # scores: more negative = better
        log(f"[round {rnd}] top{cfg.top_k} mean score={metric:.3f} "
            f"(prev {best_metric:.3f}); fresh={got}")
        if improved:
            best_metric = metric; stale = 0
        else:
            stale += 1
            if stale >= cfg.patience:
                log(f"[round {rnd}] loop-until-dry: no improvement x{stale}; stop"); break

    if belief_store:          # compound: persist this run's synthon learning
        bandit.save(belief_store)
        log(f"[belief-store] saved {len(bandit.obs)} synthon arms + "
            f"{len(bandit._all_rewards)} observations to {belief_store}")

    ranked = sorted(pool.values(), key=lambda p: p.score)
    return {
        "pool": pool, "ranked": ranked, "n_docks": n_docks,
        "top": [{"smiles": p.smiles, "score": p.score,
                 "reaction_class": p.reaction_class,
                 "scaffold": murcko(p.smiles)} for p in ranked[:cfg.top_k]],
    }


def diverse_top(ranked: list[Product], k: int) -> list[Product]:
    """Scaffold-diverse skim of the ranked hits (one per Murcko scaffold)."""
    out, seen = [], set()
    for p in ranked:
        s = murcko(p.smiles)
        if s in seen:
            continue
        seen.add(s)
        out.append(p)
        if len(out) >= k:
            break
    return out
