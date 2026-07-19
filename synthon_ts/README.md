# synthon_ts — synthon-aware Thompson-Sampling screen over onepot CORE

Get a hit → decompose it to building blocks (synthons) → learn which synthons
dock well → have onepot fetch new molecules carrying *similar* synthons at the
same reaction positions → dock those → repeat. An offline-seeded bandit over
synthons, with a **docking oracle** as the fitness function.

There are **two interchangeable docking backends** behind one shared loop
(`core.run_loop`). Retrieval is always the **direct onepot REST API** (the muni
`onepot` tool is broken — the `credits_remaining` gateway bug).

```
core.py          shared loop: onepot retrieval, synthon bandit, loop-until-dry, diversity
oracle_muni.py   backend A — docking via muni CLI `rowan_batch_docking` (QVina2)
oracle_rowan.py  backend B — docking via the direct Rowan SDK (no muni)
run_muni.py      entrypoint A
run_direct.py    entrypoint B
mocks.py         offline mocks for a credit-free dry run
```

## One command (recommended): `run.py`

Point it at a **target PDB** + a **ligand that already sits in the pocket**, and
everything runs — box, protein prep, seeding, the screen, and **docked poses on
disk for PyMOL**:

```bash
export ONEPOT_API_KEY=...          # retrieval + sample_space
export ROWAN_API_KEY=...           # docking backend AND pose export (always)

python -m synthon_ts.run \
  --pdb 1HCK --pocket-ligand ATP \
  --query "CCC(CO)Nc1nc(NCc2ccccc2)c2ncn(C(C)C)c2n1" \
  --backend direct --out-dir runs/cdk2 --pymol
```

What happens:
1. **Box from the pocket ligand** — `pocket.box_from_ligand` parses the named
   ligand (`ATP`) out of the **raw** PDB (id auto-downloaded from RCSB, or a local
   `.pdb`/`.cif`) and builds the docking box: center = ligand centroid, size =
   ligand extent + `2*--padding`. The receptor is prepped/stripped separately but
   the frame is preserved, so the box stays valid. (Verified: Rowan `prepare()`
   keeps the raw coordinate frame — CA centroids match within 0.3 Å.)
2. **Protein prep** — direct: `create_protein_from_pdb_id` + `prepare()` → a
   prepared UUID; muni: the PDB id/UUID + box are passed to the batch job.
3. **Seed** — from `--query` (a known binder) **or**, if you omit `--query`, by
   sampling drug-like CORE-native molecules from onepot `sample_space` (below).
4. **Screen** — `core.run_loop` (synthon Thompson Sampling), with a **local RDKit
   drug-like filter** dropping off-profile analogs *before* docking.
5. **Poses for PyMOL** — after **every round** and at the end,
   `poses.export_top_poses` writes the current top-`--top-k` poses so PyMOL can
   refresh live. On the **direct** backend it **reuses the pose from the screening
   dock** (via its pose uuid) — no re-dock, so the pose you view is the exact one
   that was scored/ranked. Only the **muni** backend (batch docking = scores-only,
   no poses) re-docks the top hits through direct Rowan to get geometry.

### Seeding without a known binder (`sample_space`)
Omit `--query` and the loop is seeded from onepot `sample_space` — property-
filtered, **makeable** molecules straight out of CORE (they decompose cleanly, so
they make better seeds than an arbitrary drug):

```bash
python -m synthon_ts.run --pdb 1HCK --pocket-ligand ATP \
  --backend direct --out-dir runs/cdk2_denovo \
  --sample-seeds 3 --sample-strategy diverse --qed-min 0.5
```

The drawn seeds (SMILES) and the returned `seed` int are logged; replay an exact
draw with `--sample-seed <int>`.

### Drug-like window (two places)
The same thresholds drive **both** the `sample_space` seed draw (native property
filters) **and** a **local RDKit filter** on retrieved analogs (onepot's analog
`search` has *no* property filters, so analogs are filtered locally before
docking). Defaults: MW ≤ `--mw-max` 550 · cLogP `--logp-min/-max` 1–5 · HBD ≤
`--hbd-max` 5 · HBA ≤ `--hba-max` 10 · TPSA ≤ `--tpsa-max` 140 · rot. bonds ≤
`--rot-max` 8 (plus `--mw-min` 300 / `--qed-min` 0.5 for seeds only). **Caveat:**
the med-chem target is cLogD(pH 7.4) 1–5, but neither RDKit nor onepot expose
cLogD — RDKit `Crippen.MolLogP` (matching onepot's `clogp`) is used as a proxy.
Disable the analog filter with `--no-druglike-filter`.

### Cost & supplier constraints
`--max-price` (default **200**) and `--max-supplier-risk` (default **"low"**) keep
only cheap, low-risk compounds. onepot prices are two tiers (**$125 / $295**), so
`--max-price 200` keeps just the $125 tier. Both are **native** onepot `search`
filters, applied to every retrieval (seed + elaboration) so the API returns only
qualifying analogs. `sample_space` has no price/supplier param, so seeds are
**over-sampled and price-cut locally** on the returned `price_usd` (supplier risk
isn't returned at draw time, but the downstream analog searches enforce it, so
docked hits stay low-risk). Drops are logged.

### Convergence metrics (`convergence.json`)
Each round the callback also writes `<out_dir>/convergence.json` (full history,
overwritten every round) for a dashboard: per round it records `best_score`,
`topk_mean`, the accumulated `scores` distribution, the synthon-space `positions`
(per `reaction_class` → per `bb_index`, the top ~12 building blocks ranked by mean
reward = −score, with count + std), and the pruning funnel `n_synthons_seen` /
`n_synthons_confident` (confident = seen ≥ 2×). It is a serialization of the
bandit state — no extra docking.

### Building-block clusters and centroids

After **all screening rounds have completed**, the best `--cluster-top-n`
compounds from the final accumulated docking-score ranking (default **100**) are
analyzed using their attributed building blocks. Those building blocks are
clustered separately per
`(reaction_class, bb_index)` using Butina clustering over radius-2, 2048-bit
Morgan-fingerprint Tanimoto similarity. The default minimum similarity is 0.60;
change it with `--bb-cluster-similarity` or set `--cluster-top-n 0` to disable
the analysis.

`building_block_clusters.json` contains every cluster and member. Because a
binary fingerprint has no molecule-valued arithmetic centroid, each reported
centroid is a **medoid**: the observed building block with the highest
occurrence-weighted mean Tanimoto similarity to its cluster. A compact list of
those representatives is also written to `building_block_centroids.csv`, with
cluster size, compound support, and docking-score summaries. Precise attribution
(the default) analyzes each hit's own decomposition; `--coarse` can only report
the anchor building blocks used by its approximate attribution.

`building_block_clusters.svg` plots each reaction position separately. Unique
building blocks are projected into two dimensions by classical MDS over their
Morgan/Tanimoto distance matrix, colored by cluster, sized by recurrence, and
shown with the centroid medoids as labeled diamonds.

### Full-molecule top hits with building-block diversity

The final `--top-k` output contains complete product molecules, not isolated
building blocks. This ranking runs **only after `run_loop` has finished every
generation round**; it does not affect OnePot generation, Thompson-sampling
anchors, pruning, or live per-round pose displays.

The post-generation ranking takes the best `--cluster-top-n` scored products
(default 100), maps each product to its reaction-aware building-block cluster
signature, and retains the strongest-affinity full molecule for each exact
signature. It then greedily selects full molecules using:

```text
selection score = 0.60 × minimum signature diversity
                + 0.40 × normalized docking affinity
```

Signature diversity is the fraction of building-block positions assigned to a
different cluster from the nearest already selected molecule. Different
reaction classes are treated as maximally different because their building-block
positions are not interchangeable. The first selection is the best-affinity
signature representative. If fewer than `--top-k` unique signatures exist, the
remaining slots are filled by affinity and explicitly marked as fallbacks.
Set `--cluster-top-n 0` to disable both this post-generation ranking and its
cluster artifacts; the existing scaffold-diverse affinity ranking is then used.

`hits.json` records each selected full molecule, its reaction class, building
blocks, cluster signature, diversity component, affinity component, combined
selection score, and fallback status. `hits.csv`, `top_hits.sdf`, and the rank
PDB files contain the same selected full molecules.

### Viewing the poses in PyMOL
The out-dir gets: `receptor.pdb` (grey cartoon), `rank{N}_{score}.pdb` (one docked
ligand per hit, green-carbon sticks), `top_hits.sdf` (best-effort combined SDF),
`hits.csv` (rank, smiles, score, mmgbsa, posebusters_valid), `view.pml`,
`hits.json` (ranked summary), and `convergence.json` (per-round metrics, below).

```bash
/Users/bb/.local/share/mamba/envs/cadd-pymol/bin/pymol runs/cdk2/view.pml
```

`view.pml` `reinitialize`s, `cd`s into the out-dir, loads everything, colours by
element (green ligand carbons), labels each pose with its rank+score, and zooms
the pocket. **Live refresh while a screen is still running:** in the PyMOL command
line re-run `@view.pml` (or *File → Reload All*) to pull in the newest poses.
`--pymol` auto-launches PyMOL on the output.

> Pose export **always** needs `ROWAN_API_KEY`, even with `--backend muni` — muni
> batch docking returns scores only (no poses), so poses are re-docked through
> direct Rowan against a separately-prepared receptor. `mmgbsa` is written when the
> Rowan result includes it (absent on some accounts → left blank).

## The two backends (lower-level runners)

| | **A — muni** (`run_muni.py`) | **B — direct** (`run_direct.py`) |
|---|---|---|
| Docking | muni `rowan_batch_docking` (QVina2, **batch**) | Rowan SDK single `docking` (box, **per-ligand**) |
| Credits | muni credits (~**0.05 cr/ligand**) | your Rowan credits (~**1 cr/ligand**) |
| Keys | `ONEPOT_API_KEY` only (muni holds Rowan) | `ONEPOT_API_KEY` + `ROWAN_API_KEY` |
| Runs on | the shared muni workspace | any machine, no muni |
| Speed | fast (one batch job/round) | slower (one workflow/ligand, submitted concurrently) |

**Why not Rowan `batch_docking` in version B?** It is feature-gated on a standard
subscription — `400 You do not have access to this feature` even with credits
(verified 2026-07-19). Single `docking` and `analogue_docking` *are* accessible;
`analogue_docking` needs a bound template pose (fails on arbitrary anchors), so
version B uses box-based single `docking`, which mirrors the muni batch semantics.

## Setup

```bash
pip install onepot rdkit          # both versions
pip install rowan-python          # version B only
export ONEPOT_API_KEY=...          # both
export ROWAN_API_KEY=...           # version B only
```

## Run

```bash
# A — through muni (bind a page first so trials land there)
muni page use "Synthon-TS oracle calibration (CDK2 1HCK)"
python -m synthon_ts.run_muni   --query "<seed SMILES>" --max-docks 200

# B — direct APIs only
python -m synthon_ts.run_direct --query "<seed SMILES>" --max-docks 120
```

Both default to the **CDK2 / PDB 1HCK** ATP-site target (the calibrated box +,
for version A, the prepared protein UUID). Override with `--protein` / `--pocket`
for another target. A credit-free dry run:

```bash
python -m synthon_ts.selftest      # mocks onepot + docking, exercises the loop
```

## The loop (what `core.run_loop` does)

1. **Seed** — `onepot.search(query, decompose=True)` → a batch of purchasable
   **analogs** (which carry NO synthons) plus the **query's** own decomposition.
   onepot only ever decomposes the *query*, never the returned analogs.
2. **Dock** every new analog (oracle; lower score = better).
3. **Attribute** each analog's reward (= −score) to synthon arms
   (`(reaction_class, bb_index, bb_smiles)`). Two modes:
   - **precise** (default, `precise_attribution=True`) — after docking a batch,
     ONE `search(analogs, decompose=True, max_results=1)` recovers **each docked
     analog's OWN synthons**; the reward is credited to them (true per-product
     credit). Costs ~1 onepot credit per decomposed analog. Analogs onepot can't
     decompose keep their score but aren't attributed.
   - **coarse** (`precise_attribution=False`, `--coarse`) — credit each analog to
     the **anchor's** synthons instead (the query that retrieved it), full credit
     to the exploited position (the `bb_filters` window guarantees the analog
     shares that synthon), little/none to the varying explored one; seed round
     gives uniform weak credit. No extra onepot calls, so much cheaper.
4. **Thompson-sample**: rank candidate anchors by the best TS draw over their
   synthons (favours high mean *and* high uncertainty); take the top few,
   scaffold-diverse.
5. For each anchor, **exploit** its best position (tight Tanimoto window) and
   **explore** its least-observed position (loose, capped window) via onepot
   `bb_filters`; fetch new products.
6. Go to 2. **Loop-until-dry**: stop when the top-k mean score stops improving
   for `patience` rounds (or `--max-docks` / `--max-rounds`).

### Known approximations (documented, not bugs)
- Credit assignment to synthons is lossy — used only to *prioritise*; every
  product is always docked in full.
- Precise attribution credits an analog's **primary** decomposition (onepot may
  offer several retrosynthetic routes); coarse attribution credits the anchor's
  synthons, not the analog's own.
- Per-position independence is an approximation.
- onepot caps: `max_results ≤ 100`, `max_depth = 1`.

## Calibration (CDK2 / 1HCK, QVina2 via muni)
Known actives vs decoys, docking-ready box from Rowan Pocketeer (rank-1 pocket,
druggability 8.562). Scores (lower = better): roscovitine −7.2 (best), ibuprofen
−6.5, caffeine −5.7, adenine −5.2, benzene −3.8, acetic acid −3.2. The oracle
ranks the real CDK2 inhibitor top and tiny fragments bottom — box + oracle
discriminate. Cost: 0.31 muni credits for 6 ligands.
