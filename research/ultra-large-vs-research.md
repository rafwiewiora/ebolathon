# Ultra-Large Library Virtual Screening — Research Dossier

Compiled from four parallel research agents (2019–2026 literature). Focus: strategies for
billion+ compound screening, with a deep dive on **building a GPU-powered, Pharmit-style
shape + pharmacophore engine** (GPU conformer generation + GPU Gaussian shape overlay).

Sections:
1. Strategic landscape & taxonomy
2. Hierarchical / synthon-based docking
3. Shape & pharmacophore methods (the Pharmit-clone deep dive)
4. GPU cheminformatics stack (build-your-own-tool)
5. Pitfalls, benchmarks, and the modern layered pipeline
6. Consolidated bibliography

---

## 1. Strategic landscape & taxonomy

Enumerable make-on-demand space grew ~10⁷ (2019) → 10¹⁰–10¹¹ tangible (2026): Enamine REAL
~48–95B, xREAL ~4.4T, WuXi GalaXi, Mcule Ultimate, Freedom, OTAVA CHEMriya. Two coupled
problems: (i) how to *search* a space too big to enumerate-and-dock; (ii) how to *triage*
hit lists numerically dominated by scoring artifacts.

Complementary layers (nearly every serious campaign combines 3+):
- **Brute-force docking** — DOCK3.7/3.8 (Shoichet/Irwin + ZINC22), VirtualFlow (Gorgulla).
- **Hierarchical / synthon docking** — dock building blocks, prune, elaborate (V-SYNTHES,
  Chemical Space Docking). Cost scales with reagents, not products.
- **Active-learning / ML surrogate docking** — Deep Docking, MolPAL, HASTEN, AL-Glide.
- **Shape / pharmacophore / 2D prefilter** — cheap gate to carve a dockable subset.
- **ML scoring & co-folding** — Boltz-2, AlphaFold3, DiffDock-L; late-stage rescorers.

Canonical modern stack: cheap ligand-based/synthon reduction → docking (brute or AL) →
orthogonal artifact-stripping rescore → MM/GBSA or FEP → human inspection → synthesis.
Organizing principle: **consensus across orthogonal error modes**, diversity over top-score.

---

## 2. Hierarchical / synthon-based docking

**Premise:** a REAL-type library is reactions × reagents. Attach a minimal "cap" to a
reagent's reactive group → a **synthon** carrying an extension vector. Dock synthons first,
keep those whose pose points the growth vector into open pocket volume, enumerate full
products only for surviving combinations. Combinatorial cost O(reagents^n) → ~O(reagents×n).

### V-SYNTHES (Katritch lab, USC) — reference method
Sadybekov et al., *Nature* 601:452 (2022). Four stages:
1. **Minimal Enumeration Library (MEL)** ~600K: one R-group fully enumerated with real
   synthons, other attachment points capped (methyl/phenyl) to keep fragments realistic.
2. **Dock MEL** (published with ICM-Pro/Molsoft; engine-agnostic).
3. **Productive-pose filter**: keep high-scoring seeds whose growth vector points to
   unoccupied volume.
4. **Hierarchical elaboration**: replace caps with full synthon set one position at a time,
   re-dock focused library, prune, repeat. Only final enumeration is docked in full.

Numbers: ~0.5–1M docked in last step vs 11B → **>5,000× reduction**; computational EF100 ≈
250 (2-cmpt) / 460 (3-cmpt). Prospective CB2: **33% hit rate** (vs 15% for standard VLS of
115M), best Ki ≈ 0.28 nM, optimized to 1 nM with 50–200× CB2/CB1 selectivity.

**V-SYNTHES2** (Nazarova et al., *npj Drug Discovery* 2026): automated fragment selection via
geometry-based **CapSelect** (removes manual pose inspection), scaled to ~36B, validated on
harder targets (shallow pockets, RNA site, GPCRs).

### Chemical Space Docking (Beroza/Lemmen — Genentech/BioSolveIT)
Beroza et al., *Nat. Commun.* 13:6447 (2022). Incremental-construction docking (FlexX/SeeSAR
lineage) on reaction-based space: dock building block with extension vector, grow products
*in situ* by compatible reagents — never enumerate the full product space. ROCK1: 69
purchased, **27 (39%) Ki < 10 µM**. Commercialized as **Chemical Space Docking®** (BioSolveIT).
Distinction from V-SYNTHES: CSD grows atom-by-fragment inside pocket via FlexX; V-SYNTHES
caps/enumerates and re-docks discrete focused libraries.

### Rarey/ZBH non-enumerative substrate (powers BioSolveIT)
- **FTrees** — fuzzy feature-tree similarity over combinatorial spaces.
- **SpaceLight** (Bellmann et al., *JCIM* 2021) — 2D fingerprint Tanimoto over billion–trillion
  spaces in seconds on a laptop.
- **SpaceMACS** — max-common-substructure over fragment spaces.
- **SpaceGrow** (Hönig et al., *JCAMD* 2024) — **shape-based (ROCS-like) 3D** screening of
  billion-sized combinatorial fragment spaces in hours on **a single CPU**; resources scale
  with synthons, not molecules.
Packaged in **infiniSee** (Scaffold Hopper / Analog Hunter / Motif Matcher) and infiniSee
xREAL (trillions).

### Fragment-context caveat
Fragments docked in isolation often DON'T reproduce their pose in the full ligand (the rest
of the molecule provides anchoring contacts). Every synthon method therefore re-docks full
products in the final stage. This motivated Thompson Sampling (below), which scores full
products but chooses which via reagent-space statistics.

### ML / active-learning surrogate docking
Dock small seed → train cheap surrogate (fingerprint→dock-score) → acquisition picks next
batch → dock, retrain, iterate. 20–100× fewer docking calls.
- **Deep Docking** (Gentile/Cherkasov, *ACS Cent. Sci.* 2020; *Nat. Protoc.* 2021) — QSAR
  classifier prunes; up to 100× fewer dockings; 1.36B ZINC15, 1.3B SARS-CoV-2 Mpro.
- **Deep Docking Ultra (DDU)** (Pandey et al., *Chem. Sci.* 2026) — >10.1B REAL in ~10 days,
  ~28,500× reduction; current high-water mark for full-REAL surrogate docking.
- **MolPAL** (Graff/Coley, *Chem. Sci.* 2021) — batched Bayesian optimization; recovered
  94.8% of true top-50k (UCB) docking ~2.5% of a 100M pool (~40× speedup). Open-source.
- **Regression-based AL** (*JCIM* 2024, 64:2955) — simple **linear regression** retrieves ~90%
  of top-1% after docking ~10%; deep models unnecessary at shallow sampling depth.
- **Active Learning Glide** (Schrödinger, commercial); **HASTEN**/**SpaceHASTEN** (Kalliokoski,
  hybrid ML + SpaceLight non-enumerative expansion).

### Thompson Sampling (Klarich/Walters) — the key bridge
Klarich et al., *JCIM* 64:1158 (2024). Scores **full products** but chooses which via a
multi-armed bandit over **reagent space**: each reagent = an "arm" with a Beta/Gaussian
belief over product scores; sample posteriors, assemble & score one product (docking / ROCS /
FP / ML), update posteriors. Concentrates on promising reagent combos **without enumerating**.
Evaluates <0.1–1% of a 10⁹–10¹¹ library; modality-agnostic. Open-source (github.com/PatWalters/TS).
Advantage over Deep Docking/MolPAL: natively combinatorial (reasons in reagent space).
**Enhanced TS** (Yu et al., *J. Cheminform.* 2025) adds roulette-wheel selection + thermal
cycling for better diversity/exploration.

**Note:** "CHINN" could not be verified in the docking/synthon literature — likely a
mis-remembered name (nearest concepts: CNN pose scoring, or TS reagent scoring).

---

## 3. Shape & pharmacophore methods — the Pharmit-clone deep dive

### Gaussian shape overlay (Grant/Pickup formalism, the ROCS foundation)
Grant & Pickup, *J. Phys. Chem.* 99:3503 (1995): replace hard-sphere atomic volume with a
**sum of continuous Gaussians**, making molecular volume — and *overlap* — an analytically
differentiable function of relative orientation.

Each atom i: ρᵢ(r) = p·exp(−αᵢ‖r−Rᵢ‖²), prefactor **p = 2.7**, width αᵢ set so the Gaussian
integral reproduces the atom's vdW sphere volume. The **Gaussian product theorem** makes each
pairwise overlap integral closed-form:

  O_AB = Σ_{i∈A} Σ_{j∈B} π^{3/2} p² (1/(αᵢ+αⱼ))^{3/2} exp(−αᵢαⱼ/(αᵢ+αⱼ) · r_ij²)

Self-overlaps O_AA, O_BB are orientation-independent constants. OpenEye carries
inclusion-exclusion to ~6th order; all GPU implementations truncate at **first-order
(atom-atom pairwise)** for speed.

**Scoring:**
- Shape Tanimoto = O_AB / (O_AA + O_BB − O_AB) ∈ [0,1]
- Shape Tversky (α=0.95, β=0.05 default) — asymmetric, for "query fragment contained in DB
  molecule" (can exceed 1.0)
- Distance form S² = O_AA + O_BB − 2·O_AB (minimize = maximize overlap)

**Alignment optimization (your GPU engine must nail this):**
1. Rigid conformers (flexibility = pre-enumerated multi-conformer libraries).
2. Multiple starts: typically **4 inertial/"star" starts** (align principal moments, flip axes).
3. Parameterize by **unit quaternion + translation**; compute O_AB and its analytic gradient;
   refine to local overlap max by **quasi-Newton (BFGS)**.
4. Smooth surface → converges in a handful of iterations. Each conformer × start is an
   independent optimization → **massive data parallelism** (one thread-block per alignment).

**Color / pharmacophore ("color force field"):** feature types (donor, acceptor, hydrophobe,
anion, cation, ring) from SMARTS, each placed as a harder feature-Gaussian; Color Tanimoto
scored the same way. **ComboScore = Shape + Color Tanimoto ∈ [0,2]** — consistently best for
enrichment. Force fields: ImplicitMillsDean (pH-7 pKa model, robust) / ExplicitMillsDean.

### Pharmit architecture (Koes lab) — exactly what you're cloning
Sunseri & Koes, *Nucleic Acids Res.* 44:W442 (2016). Web front-end over three engines:
**Pharmer** (pharmacophore), **shapedb** (shape), **3Dmol.js** (viz). The performance secret is
the backend (Koes, *IBM J. Res. Dev.* 62, 2018):

**Pharmacophore search (Pharmer):**
- Decompose an n-feature pharmacophore into **triangles**: every feature triple → a 3-tuple of
  pairwise distances = a point in "triangle space."
- Index triplets in a **KDB-tree** (disk-friendly k-d/B-tree hybrid) → sub-linear range queries.
- Each leaf (TripleData) stores triangle vertex coords, conformer ID, metadata, and a **Bloom
  filter** of the conformer's *other* features (fast multi-feature pruning).
- Query: range-search triplet index → candidate poses → verify remaining features + exclusion
  spheres → RMSD align. Result is **exact**, not approximate.

**Shape search (shapedb):**
- Molecules → **molecular octrees**, stored in a **GSS tree** with **MSV/MIV** (max-surrounding /
  min-included volume) bounds → **branch-and-bound pruning** against a Tanimoto/Tversky
  threshold. Supports inclusive/exclusive shape constraints.
- Repos: github.com/dkoes/shapedb, github.com/dkoes/pharmer.

**Systems engineering (the "sub-second over 10⁸ conformers" magic):**
- Libraries **striped at molecule level** across many drives (24×4TB JBOD), searched
  independently, communication-free.
- Trees in depth-first order; leaves separate; molecule fetches **sorted by on-disk position**
  for sequential I/O; **mmap** (MAP_PRIVATE|MAP_NORESERVE + FADV_SEQUENTIAL); tcmalloc, bump
  allocators, pipeline parallelism.
- Scale: MolPort ~104M conformers = 2.5TB (72% is TripleData); >28TB hosted total.
- **Limitation:** bulk-load only (no incremental updates; PubChem rebuild takes days);
  semi-custom server hard to cloud-scale. **Pharmit is disk-I/O bound, NOT compute bound —
  this is exactly where a GPU rewrite of the shape-overlay stage differentiates.**

### Other open-source shape/pharmacophore tools
- **Shape-it** (github.com/rdkit/shape-it) — open Gaussian shape alignment, CPU/C++. Correctness
  baseline.
- **Align-it / Pharao** (Taminau et al. 2008) — open pharmacophore alignment as Gaussian volumes.
- **PheSA** (Rufener/von Korff, *JCIM* 2024) — open Pharmacophore-Enhanced Shape Alignment (Java,
  OpenChemLib); clean readable combined shape+color objective; competitive with ROCS.
- **RDKit**: `rdShapeAlign`, `ShapeTanimotoDist`/`ShapeProtrudeDist` (grid), **Open3DAlign (O3A)**
  (field-based flexible overlay).
- **ODDT** — Python; USR-family, pharmacophore, rescoring.

### USR / USRCAT / ElectroShape — alignment-FREE moment methods
No alignment → screening = nearest-neighbor on fixed-length vectors → trivially billion-scale.
- **USR** (Ballester & Richards 2007): 4 reference points (centroid, closest, farthest, farthest-
  from-farthest); first 3 moments of atom-distance distributions each → **12-D** descriptor.
- **USRCAT** (Schreyer & Blundell, *J. Cheminform.* 2012): moments per atom subtype (all,
  hydrophobic, aromatic, HBD, HBA) → **60-D**. **In RDKit** (`GetUSRCAT`/`GetUSR`/`GetUSRScore`).
- **ElectroShape** — adds partial-charge/lipophilicity dims.
- Best as an **ultra-fast in-RAM prefilter** (60 floats/conformer) feeding alignment rescoring.

### Pharmacophore matching paradigms
1. **Alignment-based / geometric** (Pharmer, Phase, LigandScout, Catalyst): rigid superposition
   of feature points within tolerance spheres respecting exclusion volumes. Pharmer's triangle-
   KDB indexing makes this exact-and-fast at scale.
2. **Pharmacophore fingerprints**: hash k-point pharmacophores (feature-type + binned distance)
   → bit/count vector; screen by Tanimoto (no alignment). RDKit Gobbi2D, ErG. Fast but lossy.
Commercial: Phase (Schrödinger; shape screening seeds trial alignments from atom/feature-pair
triplets — same triangle idea), LigandScout (best structure-based pharmacophore perception),
Catalyst/DS, MOE.

### Deep-learning pharmacophore (2023–2024)
- **PharmacoNet** (Seo & Kim, NeurIPS 2023) — protein-based pharmacophore model for fast large-
  scale scoring.
- **PharmacoMatch** (arXiv 2409.06316) — neural subgraph matching for 3D pharmacophore screening
  at scale.

---

## 4. GPU cheminformatics stack — build-your-own-tool

### "NVIDIA version of RDKit" = **nvMolKit**
github.com/NVIDIA-BioNeMo/nvMolKit (mirror NVIDIA-Digital-Bio). The de facto GPU RDKit — a
**core-functions accelerator**, NOT a full port. GPU-accelerates: Morgan fingerprints,
Tanimoto/cosine similarity, **MMFF minimization**, **ETKDG conformer generation**. 1–4 orders of
magnitude over CPU RDKit. v0.5.1 (Jun 2026), CUDA (≥cc7.0/V100+, CUDA ≥12.5), needs PyTorch-CUDA,
RDKit 2025.03–2026.03, Python 3.11–3.14. `conda install -c conda-forge nvmolkit`. One benchmark:
ETKDG 63.8 conf/s GPU vs 3.9 conf/s CPU. (Caveat: `AceConfgen` reported ~40× faster than nvMolKit
on one large benchmark — nvMolKit isn't necessarily the fastest GPU conformer route.)

Broader NVIDIA stack (mostly not "GPU RDKit"): **BioNeMo** (ML framework, DiffDock/AF2 NIMs),
**cuEquivariance** (E(3)-equivariant NN kernels — relevant only for ML conformer/scoring),
**RAPIDS/cuGraph/cuDF** (GPU dataframes/fingerprints/clustering — data wrangling, not overlay).
There is **no** official GPU fork of RDKit core.

### GPU Gaussian shape overlay — the core engine
- **roshambo** (Atwi et al., *JCIM* 64:8945, 2024) — open-source **FastROCS analog** from
  **Biogen** (NOT the Volkamer lab). RDKit I/O + conformers, **PAPER** (CUDA) as the Gaussian-
  overlap backend via Cython, shape + color (RDKit features) + ComboTanimoto, aligned-SDF
  output. **GPL-3.0.** Friction: needs a **from-source-compiled RDKit** (2023.03.1) — main
  install pain / maturity caveat.
- **roshambo2** (Atwi et al., *JCIM* 2025, 5c01322) — **USE THIS.** Near-complete rewrite,
  **>200×** over roshambo via GPU memory/batching engineering, multi-GPU **server mode**,
  **HDF5** large-library mode, **MIT license**, and crucially **does NOT need custom-compiled
  RDKit**. Targets billion-molecule libraries. Repo github.com/molecularinformatics/roshambo2.
  API:
  ```python
  from roshambo2 import Roshambo2
  calc = Roshambo2('query.sdf', 'library.sdf', color=True)
  scores = calc.compute(optim_mode='combination')
  calc.write_best_fit_structures(hits_sdf_prefix='hits')
  ```
- **PAPER** (Haque & Pande, *JCC* 2010) — the open CUDA Gaussian-overlap **optimizer** underneath
  roshambo. Maps each molecule-orientation to an independent CUDA thread. Historically >1 order
  of magnitude over commercial ROCS. simtk.org group_id 339. **This is the kernel to fork rather
  than reimplement the overlap optimizer.** (Related: GWEGA — GPU WEGA.)
- **FastROCS / ROCS X** (OpenEye/Cadence, commercial) — the **performance bar**: >2M
  conformers/s/GPU; 1B compounds in <1 day on 4 GPUs. Your open stack will be slower but free
  and hackable.

### GPU / high-throughput conformer generation
| Tool | GPU? | Method | Notes |
|---|---|---|---|
| RDKit ETKDGv3 | CPU, embarrassingly parallel | distance-geometry + torsion knowledge | Default; parallelize across cores/nodes. |
| **nvMolKit ETKDG** | **True GPU** | GPU ETKDG + GPU MMFF | RDKit-native GPU option. |
| **Auto3D** (Isayev, *JCIM* 2022) | **GPU** | isomer enum + embed + **ANI/AIMNet NNP** optimization + rank | Energy-ranked quality, slower. github.com/isayevlab/Auto3D_pkg |
| Torsional Diffusion (Jing/MIT) | GPU (ML) | diffusion on torsion hypertorus, RDKit-seeded | Fast sample but needs seed; research-grade. |
| GeoMol / GeoDiff | GPU (ML) | GNN over torsions / equivariant diffusion | Earlier ML; research-grade. |
| Flow-matching (2024–25) | GPU (ML) | SO(3) flow matching, ET-Flow, FlexiFlow | Active frontier; not production-hardened. |
| OMEGA (OpenEye) | mostly CPU | systematic torsion search | Commercial; NVIDIA collab claims 30× throughput. |

For a shape screen you need many conformers fast, "good enough for shape": **nvMolKit GPU-ETKDG**
(fastest to integrate) or massively parallel **RDKit ETKDGv3**. Use **Auto3D** only if you need
energy-ranked/NNP quality. roshambo/roshambo2 also have optional built-in conformer generation.

### GPU building blocks for a custom overlay
- **CuPy** — NumPy-on-CUDA for Gaussian-volume math, Tanimoto matrices, batched linear algebra.
- **Numba CUDA** — write a bespoke overlap kernel in Python if not forking PAPER.
- **PyTorch / JAX** — **differentiable shape alignment**: shape = sum of atom-centered Gaussians,
  volume overlap as a differentiable function of pose (quaternion + translation), optimize with
  autograd/Adam over batched poses. JAX `vmap`+`jit` natural for thousands of pose starts. Modern
  "reimplement ROCS as a differentiable objective" route; unifies shape + color in one autograd
  objective.
- **GNINA 1.3** (*J. Cheminform.* 2025) — GPU CNN docking (PyTorch scoring, knowledge-distilled
  CNN for HTVS); optional structure-based rescoring on top of ligand-based shape hits.

### Recommended open-source architecture for a GPU Pharmit-clone
1. **Chemistry core: RDKit** (BSD) — I/O, standardization, feature/pharmacophore definitions
   (`FeatureFactory` gives donor/acceptor/aromatic/hydrophobe = your "color" layer).
2. **GPU conformers: nvMolKit GPU-ETKDG + GPU-MMFF** (or parallel CPU ETKDGv3); Auto3D for
   curated high-quality subsets.
3. **GPU shape+color overlay core: reuse roshambo2** (MIT, multi-GPU, HDF5). Gives Shape + Color
   + Combo Tanimoto + best-fit pose export. **Don't rewrite the Gaussian-overlap optimizer** —
   roshambo2 (modern) or PAPER (classic) already solve the hard part.
4. **Pharmacophore constraint layer (main net-new code):** the thing that makes it "Pharmit," not
   just "ROCS" — geometric feature-point matching with tolerance spheres + exclusion volumes.
   Build on RDKit features + a GPU-batched RMSD/constraint check (CuPy/Numba/PyTorch). Consider
   Pharmit's triplet-KDB-tree + Bloom-filter indexing as the blueprint (but the I/O-bound part;
   your GPU win is overlay + RMSD verification).
5. **Prefilter + orchestration:** GPU fingerprints/Tanimoto (nvMolKit/RAPIDS) and USRCAT (in-RAM,
   60-D) for cheap prefilters; Dask/Ray for multi-GPU fan-out; HDF5/Parquet conformer store.
6. **Optional tail: GNINA 1.3** CNN rescoring; shape-then-dock funnel.

**Reuse vs build:** reuse RDKit, nvMolKit, roshambo2/PAPER, GNINA. Build the pharmacophore-
constraint query engine, library indexing/sharding, orchestration/UI. Optionally a JAX/PyTorch
differentiable-overlap module for a unified batched shape+pharmacophore objective.

**License watch:** RDKit BSD, nvMolKit open (NVIDIA), roshambo2 **MIT**, roshambo v1 **GPL-3.0**,
GNINA Apache-2.0, PAPER open. Use roshambo2, not v1.

**Corrections to common assumptions:** roshambo is from **Biogen**, not Volkamer; "NVIDIA RDKit"
= **nvMolKit** (core accelerator, not a fork); roshambo's CUDA backend is **PAPER** (Haque & Pande).

---

## 5. Pitfalls, benchmarks, and the modern layered pipeline

### Landmark campaigns & lessons
- **Lyu et al., *Nature* 2019** — 170M vs AmpC/D4: bigger libraries → higher hit rates (11–24%),
  better potency (180 pM D4 agonist), novel chemotypes. Founded ZINC/DOCK3.7 workflow.
- **Alon et al., *Nature* 2021** — σ2 receptor: large-library docking generalizes to membrane
  receptors; selectivity engineerable from geometry.
- **V-SYNTHES *Nature* 2022** — you don't have to enumerate to win.
- **Lyu et al., *Nat. Chem. Biol.* 2025** — 1.7B vs 99M AmpC, tested **1,521** compounds: larger
  screen ~**doubled hit rate**, ~50× more inhibitors. BUT hit-rate estimates **only converged
  after testing several hundred** — a rebuke to reporting "hit rate" from 20–50 compounds.

**Cumulative lessons:** (1) bigger genuinely helps for *discovery* — but only with enough testing
and artifact control; (2) dock score **enriches, does not rank** individual affinity;
(3) **chemical diversity beats raw top-score** — cluster and pick sane, well-posed reps across a
score band, not the literal top-N.

### The "top of the list is noise" problem
As libraries grow, the extreme top is enriched in **rare artifacts exploiting scoring
pathologies** (over-counted electrostatics/H-bonds, unphysical buried nonpolar surface, strained
poses).
- **"Identifying Artifacts from Large Library Docking"** (Shoichet, *J. Med. Chem.* 2024): from the
  top of a 1.7B AmpC list, synthesized 128 — **0/39 "cheaters"** active, **57% of 89 "plausible"**
  active. Mitigate with implicit-solvent rescoring + AB-FEP.
- **Sindt/Rognan, *JCIM* 2025** — confirms orthogonal rescoring separates cheaters (0%) from
  plausible (57%) but stresses reliable rescoring of *millions* of poses is genuinely hard; the
  bottleneck is **diversity-aware post-processing**, not more docking.

### Benchmarks
- **DUD-E** (Mysinger 2012) — analog/decoy bias, learnable; treat skeptically for ML.
- **DUDE-Z** / property-matched decoys — de-biased.
- **LIT-PCBA** (Tran-Nguyen/Rognan, *JCIM* 2020) — 7,844 actives / 407,381 experimentally
  confirmed inactives, 15 targets; EF1%/EF10%/BEDROC. **Current preferred rigorous benchmark**;
  explicitly benchmarks shape vs docking.
- **Shape vs docking** (Hawkins/Skillman/Nicholls, *J. Med. Chem.* 2007) — ligand-based ROCS often
  matches or beats docking when a good query ligand exists, far cheaper. Ligand-centric shape/
  pharmacophore methods are more *consistent* across targets than docking.
- Benchmark your GPU engine for **correctness against Shape-it / Align-it / PheSA / ROCS**, and
  for enrichment on **LIT-PCBA (BEDROC/EF1%)**, not DUD-E.

### Co-folding / ML scoring frontier (2024–2026)
- **Boltz-2** (Passaro/Corso et al., bioRxiv 2025) — AF3-style structure module + **trained
  affinity head** approaching **FEP accuracy at ~1000× lower cost** (~20 s/compound/GPU). Best on
  CASP16 affinity; ~doubles average precision vs ML/docking on MF-PCBA. **MIT, open-source.**
  Position: **late-stage rescorer / triage oracle** on top 10³–10⁵, not a primary screen (too slow
  for 10⁹). Efficiency wrappers: Boltzina, FlashAffinity. Reliability uneven off-distribution.
- **AlphaFold3** (Abramson, *Nature* 2024) — complex prediction; useful for consensus but costly,
  and "have co-folding methods moved beyond memorisation?" is an open concern.
- **DiffDock-L** (Corso 2024) — better pose *hypothesis generation* than *scoring*; use diffusion/
  co-folding for pose + consensus, trained affinity heads / FEP for ranking.
- **Generative synthon-aware design** — SyntheMol (MCTS over BBs+reactions, *Nat. Mach. Intell.*
  2024), SynFlowNet (GFlowNet, ICLR 2025): steer generation through 10¹⁵ synthesizable space
  toward the pocket; inherit docking's scoring noise as reward.

### Modern layered pipeline (reference)
- **L0 Target/structure prep** — validate site with retrospective DUDE-Z/property-matched
  enrichment before screening; protonation/tautomers/waters/ensembles.
- **L1 Cheap ligand-based reduction (10¹⁰–10¹² → 10⁸–10⁹)** — property windows, PAINS/reactive
  filters, FTrees/SpaceLight/shape/pharmacophore on combinatorial space. Cap MW/logP to counter
  property inflation.
- **L2 Accelerated docking (→ 10⁴–10⁶)** — synthon-hierarchical (V-SYNTHES2) or active-learning
  (Deep Docking/MolPAL/HASTEN; simple regression suffices). Keep a score band + scaffold diversity.
- **L3 Orthogonal artifact stripping (→ 10²–10³)** — the most-skipped, most-critical step:
  implicit-solvent rescoring, strain/buried-nonpolar checks, cross-engine pose consistency,
  optionally Boltz-2 affinity / shape-pharmacophore consensus. Discard cheaters.
- **L4 Physics ranking (→ 10¹–10²)** — MM/GBSA then FEP+/OpenFE on survivors; Boltz-2 as pre-FEP
  triage.
- **L5 Human visual inspection (mandatory)** — pose sanity, interaction plausibility, make-on-
  demand feasibility (Rognan, *J. Med. Chem.* 2021). Diversity & novelty over score.
- **L6 Experimental validation & iteration** — test **hundreds** to estimate hit rate honestly;
  SAR-by-catalog analog search (SpaceLight/Analog Hunter) around confirmed hits for potency jumps.

**Guiding principles:** score enriches not ranks; consensus across orthogonal error modes;
diversity over top-scores; test enough to know your hit rate; the assay is ground truth.

---

## 6. Consolidated bibliography

**Landmark campaigns / library size**
- Lyu et al. *Ultra-large library docking for discovering new chemotypes.* Nature 566:224 (2019). 10.1038/s41586-019-0917-9
- Gorgulla et al. *VirtualFlow.* Nature 580:663 (2020). 10.1038/s41586-020-2117-z; VirtualFlow 2.0 bioRxiv 2023.04.25.537981
- Alon et al. *σ2 receptor structures enable docking.* Nature 600:759 (2021). 10.1038/s41586-021-04175-x
- Sadybekov et al. *V-SYNTHES.* Nature 601:452 (2022). 10.1038/s41586-021-04220-9
- Nazarova et al. *V-SYNTHES2.* npj Drug Discovery (2026). 10.1038/s44386-026-00053-6
- Lyu, Irwin, Shoichet. *Impact of library size and scale of testing.* Nat. Chem. Biol. 21 (2025). 10.1038/s41589-024-01797-w
- Tingle et al. *ZINC-22.* JCIM 63:1166 (2023). 10.1021/acs.jcim.2c01253
- Warr, Nicklaus, Nicolaou, Rarey. *Exploration of Ultralarge Compound Collections.* JCIM 62:2021 (2022).

**Hierarchical / synthon / combinatorial**
- Beroza et al. *Chemical space docking (ROCK1).* Nat. Commun. 13:6447 (2022). 10.1038/s41467-022-33981-8
- Bellmann, Penner, Rarey. *SpaceLight.* JCIM 61:238 (2021). 10.1021/acs.jcim.0c00850
- Hönig et al. *SpaceGrow.* JCAMD 38:14 (2024). 10.1007/s10822-024-00551-7
- Klarich et al. *Thompson Sampling.* JCIM 64:1158 (2024). 10.1021/acs.jcim.3c01790
- Yu et al. *Enhanced Thompson sampling.* J. Cheminform. 17:105 (2025). 10.1186/s13321-025-01105-1
- Gentile et al. *Deep Docking.* ACS Cent. Sci. 6:939 (2020). 10.1021/acscentsci.0c00229; Nat. Protoc. 16:5761 (2021)
- Pandey et al. *Deep Docking part 2 (DDU).* Chem. Sci. (2026). 10.1039/D5SC09599A
- Graff, Shakhnovich, Coley. *MolPAL.* Chem. Sci. 12:7866 (2021). 10.1039/D0SC06805E
- *Regression-based active learning.* JCIM 64:2955 (2024). 10.1021/acs.jcim.3c01661
- Kalliokoski et al. *SpaceHASTEN.* JCIM 65 (2025). 10.1021/acs.jcim.4c01790

**Shape / pharmacophore**
- Grant & Pickup. *A Gaussian Description of Molecular Shape.* J. Phys. Chem. 99:3503 (1995). 10.1021/j100011a016
- Grant, Gallardo, Pickup. *Fast molecular shape comparison.* JCC 17:1653 (1996).
- Haque & Pande. *PAPER — Accelerating Parallel Evaluations of ROCS.* JCC 31:117 (2010). 10.1002/jcc.21307
- Atwi et al. *ROSHAMBO.* JCIM 64:8945 (2024). 10.1021/acs.jcim.4c01225
- Atwi et al. *ROSHAMBO2.* JCIM (2025). 10.1021/acs.jcim.5c01322
- Koes & Camacho. *Pharmer.* JCIM 51:1307 (2011). 10.1021/ci200097m
- Sunseri & Koes. *Pharmit.* Nucleic Acids Res. 44:W442 (2016). 10.1093/nar/gkw287
- Koes. *The Pharmit Backend.* IBM J. Res. Dev. 62 (2018). 10.1147/JRD.2018.2883977
- Taminau, Thijs, De Winter. *Pharao.* J. Mol. Graph. Model. 27:161 (2008). 10.1016/j.jmgm.2008.04.003
- Rufener, von Korff et al. *PheSA.* JCIM 64:5443 (2024). 10.1021/acs.jcim.4c00516
- Ballester & Richards. *USR.* JCC 28:1711 (2007). 10.1002/jcc.20681
- Schreyer & Blundell. *USRCAT.* J. Cheminform. 4:27 (2012). 10.1186/1758-2946-4-27
- Tosco, Balle, Shiri. *Open3DAlign (O3A).* JCAMD 25:777 (2011). 10.1007/s10822-011-9462-9
- Dixon et al. *PHASE.* JCAMD 20:647 (2006); Sastry et al. *Phase Shape.* JCIM 51:2455 (2011). 10.1021/ci2002704
- Wolber & Langer. *LigandScout.* JCIM 45:160 (2005). 10.1021/ci049885e
- Seo & Kim. *PharmacoNet.* arXiv 2310.00681 (NeurIPS 2023). *PharmacoMatch* arXiv 2409.06316 (2024).

**GPU stack**
- nvMolKit: github.com/NVIDIA-BioNeMo/nvMolKit
- roshambo2: github.com/molecularinformatics/roshambo2 · roshambo: github.com/molecularinformatics/roshambo
- Auto3D: Isayev lab, JCIM 2022. 10.1021/acs.jcim.2c00817 · github.com/isayevlab/Auto3D_pkg
- McNutt et al. *GNINA 1.0.* J. Cheminform. 13:43 (2021). 10.1186/s13321-021-00522-2; GNINA 1.3 (2025) 10.1186/s13321-025-00973-x
- FastROCS: eyesopen.com/fastrocs · ROCS X: eyesopen.com/rocsx
- Torsional Diffusion: github.com/gcorso/torsional-diffusion · GeoMol arXiv 2106.07802

**Pitfalls / benchmarks / co-folding**
- Shoichet et al. *Identifying Artifacts from Large Library Docking.* J. Med. Chem. 67:16292 (2024). 10.1021/acs.jmedchem.4c01632
- Sindt, Bret, Rognan. *On the Difficulty to Rescore Hits.* JCIM 65:5553 (2025). 10.1021/acs.jcim.5c00730
- Sindt & Rognan. *SBVS of ultra-large chemical spaces: Advances and pitfalls.* Eur. J. Med. Chem. (2026). S0223-5234(26)00021-8
- Corrêa Veríssimo et al. *Ultra-Large Virtual Screening: Definition, Advances, Challenges.* Mol. Inform. 44:e202400305 (2025).
- Mysinger et al. *DUD-E.* J. Med. Chem. 55:6582 (2012). 10.1021/jm300687e
- Tran-Nguyen, Jacquemard, Rognan. *LIT-PCBA.* JCIM 60:4263 (2020). 10.1021/acs.jcim.0c00155
- Hawkins, Skillman, Nicholls. *Shape-Matching vs Docking.* J. Med. Chem. 50:74 (2007). 10.1021/jm0603365
- Bender, Gahbauer et al. *A practical guide to large-scale docking.* Nat. Protoc. 16:4799 (2021). 10.1038/s41596-021-00597-z
- Rognan et al. *Decision Making in SBDD: Visual Inspection.* J. Med. Chem. 64:2489 (2021). 10.1021/acs.jmedchem.0c02227
- Passaro, Corso et al. *Boltz-2.* bioRxiv 2025.06.14.659707
- Abramson et al. *AlphaFold3.* Nature 630:493 (2024). 10.1038/s41586-024-07487-w
- Corso et al. *DiffDock-L.* arXiv 2402.18396 (ICLR 2024).
- Swanson, Zou et al. *SyntheMol.* Nat. Mach. Intell. (2024). Cretu et al. *SynFlowNet.* ICLR 2025 / arXiv 2405.01155.

**Unverified:** "CHINN" (synthon method) — not found; likely mis-remembered.
