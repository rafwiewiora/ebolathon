# Runbook — spend ~40k Rowan credits on a value-adding deep screen

**Audience:** the agent executing this. **Goal:** consume ~40,000 Rowan credits on
one target, launched within ~30 min, spent on *valuable* deep computation (not
waste). **Reality check up front:** QVina2 docking is ~0.1–2 cr/ligand — you
cannot burn 40k on docking in 30 min. The budget is consumed by the **expensive
per-hit tiers (MD pose-validation, SQM binding affinity, optional FEP)**, which
**submit in minutes but run over hours**. So "30 minutes" = time to LAUNCH; the
credits charge as jobs complete over the next several hours.

## Prereqs
- Env: `/Users/bb/.local/share/mamba/envs/researchia/bin/python` (has `rowan`, `onepot`, `rdkit`).
- Keys (env): `ROWAN_API_KEY=$(cat ~/.config/ebolathon_rowan_key)`, `ONEPOT_API_KEY` (see repo).
- `rowan.api_key = os.environ["ROWAN_API_KEY"]`.
- Target: default **1HCK / ATP** (runnable now). **Swap in the real target** by changing `--pdb`/`--pocket-ligand` and re-preparing the protein UUID.
- HARD GUARDRAIL: put `max_credits=<cap>` on EVERY workflow submission. Track cumulative `credits_charged`; stop submitting at ~38k (leave margin). There is no un-charge.

## Credit budget (calibrate, then scale to hit 40k)
| Phase | Workflow | ~cost/job (VERIFY) | count to ~40k | value |
|---|---|---|---|---|
| 0 Screen | `batch_docking` (qvina2) | ~0.1 cr | thousands | find + rank hits (cheap) |
| 1 Careful dock | single `docking`, `do_csearch=True`, exh 32, max_poses 8 | ~6–15 cr | ~1–2k | robust poses/scores |
| 2 SQM affinity | `binding_affinity` | ~tens (MEASURE) | ~300–600 | consensus rescore |
| 3 MD validation | `pose_analysis_md` | ~hundreds (MEASURE) | ~80–150 | **main burn**; pose stability |
| 4 (opt) FEP | `rbfe_graph` + `rbfe` | ~1000s/edge (MEASURE) | tens of edges | ΔΔG on top analogs |

Priors are rough — **Phase A below measures the real per-job cost, then you scale the counts so Phases 2–4 sum to ~40k.**

## Phase A — get a ranked hit pool (~10 min, cheap)
Two options; pick one:
- **Focused (preferred):** run the synthon screen on the target (see repo `synthon_ts`), e.g.
  ```bash
  python -m synthon_ts.run --pdb 1HCK --pocket-ligand ATP --sample-seeds 3 \
    --backend muni --out-dir runs/burn_screen --max-docks 400 --top-k 200
  ```
  (muni batch = fast/cheap ranking; save Rowan credits for the deep tiers.)
- **Broad pool:** draw a big drug-like set straight from onepot and batch-dock it:
  ```python
  from onepot import Client
  mols = Client(api_key=OP_KEY).sample_space(count=5000, strategy="diverse",
      properties={"molecular_weight":{"min":300,"max":550},"clogp":{"min":1,"max":5},
                  "tpsa":{"max":140},"hbd":{"max":5},"hba":{"max":10},
                  "rotatable_bonds":{"max":8},"qed":{"min":0.5}})["molecules"]
  # dock in chunks of ~50 via rowan.submit_batch_docking_workflow(protein=<prepared uuid>, pocket=<box>, ...)
  ```
Produce a ranked list `hits = [(smiles, score), ...]` (best/lowest first). Prepare the protein ONCE to a UUID: `rowan.create_protein_from_pdb_id("1HCK").prepare().uuid`.

## Phase B — CALIBRATE the expensive tiers (~10 min, ~a few hundred cr)
Submit **one** of each on a single top hit and read `credits_charged`:
```python
import rowan, os; rowan.api_key=os.environ["ROWAN_API_KEY"]
smi = hits[0][0]; PROT = "<prepared uuid>"
ba  = rowan.submit_binding_affinity_workflow(protein=PROT, ligand_structures=[rowan.Molecule.from_smiles(smi)], name="cal-ba",  max_credits=500)
md  = rowan.submit_pose_analysis_md_workflow(protein=PROT, initial_smiles=smi, num_trajectories=4, simulation_time_ns=10, name="cal-md", max_credits=3000)
for wf in (ba, md):
    wf.wait_for_result(poll_interval=10); wf.fetch_latest(in_place=True)
    print(wf.model_dump().get("name"), "->", wf.model_dump().get("credits_charged"))
```
Let `c_ba`, `c_md` be the measured costs. (Optionally calibrate one RBFE edge if you want Phase 4.)

## Phase C — LAUNCH the burn (~10 min to submit; runs for hours)
Allocate the remaining budget and compute counts from the measured costs, e.g. split
the ~40k as **~40% MD, ~35% SQM affinity, ~25% careful docking** (tune to taste):
```python
BUDGET = 40000; spent = <credits used in A+B>
n_md = int(0.40*(BUDGET-spent)/c_md)          # top n_md hits -> pose_analysis_md
n_ba = int(0.35*(BUDGET-spent)/c_ba)          # top n_ba hits -> binding_affinity
# submit in concurrent waves (submit all, don't block per-job); cap each with max_credits.
for smi,_ in hits[:n_md]:
    rowan.submit_pose_analysis_md_workflow(protein=PROT, initial_smiles=smi,
        num_trajectories=4, simulation_time_ns=10, name="burn-md", max_credits=int(c_md*1.5))
for smi,_ in hits[:n_ba]:
    rowan.submit_binding_affinity_workflow(protein=PROT,
        ligand_structures=[rowan.Molecule.from_smiles(smi)], name="burn-ba", max_credits=int(c_ba*1.5))
# careful docking on the next tranche to soak the remaining ~25%:
from stjames.workflows.docking import VinaSettings
for smi,_ in hits[n_md:n_md+n_dock]:
    rowan.submit_docking_workflow(protein=PROT, pocket=BOX,
        initial_molecule=rowan.Molecule.from_smiles(smi), do_csearch=True,
        docking_settings=VinaSettings(executable="qvina2", scoring_function="vina",
            exhaustiveness=32, max_poses=8), name="burn-dock", max_credits=30)
```
Submit fast (these are POSTs). They queue/run on Rowan over the next hours.

## Monitor + stop
- Sum `credits_charged` across your submitted workflows (or watch the Rowan web console) and **stop submitting once cumulative ≈ 38k**.
- `max_credits` on every job is the per-job kill-switch; the sum + your stop-submitting rule is the campaign-level cap.
- Nothing here is un-chargeable — err on the side of under-submitting and topping up, not over-submitting.

## Deliverables
- `runs/burn_screen/` ranked hits + poses.
- A results table: per top hit → QVina2 score, careful-dock score, SQM affinity, MD pose stability (RMSD/occupancy). This is the *value* the 40k bought: a deeply-validated shortlist, not just a docking rank.

## Honest caveats
- 40k in 30 min of *finished* compute is not physical; MD/FEP run for hours. The 30 min buys the LAUNCH + calibration; charging trails.
- Cost priors above are guesses — Phase B is mandatory; scale from measured numbers.
- If the goal were scientific value per credit (not spending the budget), you'd stop after Phase 2 — Phases 3–4 are where the big credits go and are worth it only if MD/FEP-grade validation matters to you.
