# Homology-overlay pipeline (simplest starting version)

Take a **protein–inhibitor complex PDB** → build a **variant/homolog model** by
threading a target sequence onto the reference backbone → **overlay** the two.

Zero setup beyond the `cadd-pymol` env (PyMOL + Biopython). No MODELLER, no license.

## Run

```bash
micromamba run -n cadd-pymol pymol -cq pipeline/homology_overlay.py -- \
    --complex 3ptb \                    # PDB file OR 4-letter id to fetch
    --target  target.fasta \            # FASTA file OR raw one-letter string
    --receptor-chain A \                # optional; default = chain with most CAs
    --ligand-resn BEN \                 # optional; default = all non-polymer HETATM
    --out-dir out
```

Outputs in `--out-dir`:
- `template.pdb`     — receptor pulled from the reference complex
- `variant_model.pdb`— the threaded variant model
- `inhibitor.pdb`    — the kept ligand
- `overlay.pse`      — ready-to-open PyMOL session (template gray, variant blue, ligand yellow)

## How it works

1. Load complex, split receptor chain + inhibitor.
2. Read the template sequence off the CA atoms.
3. Pairwise-align target→template (Biopython, BLOSUM62).
4. **Thread:** per aligned column, mutate the template residue to the target
   identity (PyMOL mutagenesis, best rotamer). Deletions remove the residue.
5. Overlay + write session/PDBs.

## Known limitations (what to upgrade next)

- **Insertions not built.** Target residues with no template backbone are
  skipped and reported. For indel-heavy homologs, swap step 4 for MODELLER —
  everything else stays the same.
- **Fixed backbone, no refinement.** Sidechains are placed on the template
  backbone; no loop remodeling or minimization. Add a quick OpenMM/`pdbfixer`
  relax (in the `dock` env) if you need cleaner geometry.

## Smoke test

`_test_target.fasta` is 3PTB's own sequence with ~8% substitutions and a 3-residue
deletion — a synthetic homolog that exercises every branch. Re-run the command
above with `--target pipeline/_test_target.fasta`; expect ~11 mutated, 3 deleted.
