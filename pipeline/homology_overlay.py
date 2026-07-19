"""
Simplest starting pipeline: template complex (PDB) -> variant homology model -> overlay.

Strategy (zero-setup, PyMOL-only "threading"):
  1. Load a reference protein-inhibitor complex.
  2. Read the template protein sequence straight off the structure.
  3. Pairwise-align the target (variant/homolog) sequence to the template.
  4. Thread: for every aligned position, mutate the template residue to the
     target identity (PyMOL mutagenesis wizard, best rotamer). Deletions in the
     target remove the template residue; insertions are skipped (see NOTE).
  5. Overlay the variant model back onto the reference and report RMSD. The
     inhibitor is kept from the reference so you can see it in the variant pocket.

NOTE / known limitations of this "simplest" version:
  - Insertions in the target (residues with no template backbone) are NOT built.
    You get substitutions + deletions only. For indel-heavy homologs, graduate
    this step to MODELLER; the rest of the pipeline is unchanged.
  - No loop refinement / minimization. Backbone stays fixed at the template.

Run:
  micromamba run -n cadd-pymol pymol -cq pipeline/homology_overlay.py -- \
      --complex ref.pdb --target target.fasta --receptor-chain A \
      --ligand-resn LIG --out-dir out
"""

import argparse
import os
import sys

from pymol import cmd

# Biopython for the one thing PyMOL can't do: align two raw sequences.
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.SeqUtils import seq1


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--complex", required=True,
                   help="Reference protein-inhibitor complex PDB (path or 4-letter PDB id to fetch).")
    p.add_argument("--target", required=True,
                   help="Target/variant sequence: FASTA file or a raw one-letter string.")
    p.add_argument("--receptor-chain", default=None,
                   help="Chain to model. Default: the chain with the most CA atoms.")
    p.add_argument("--ligand-resn", default=None,
                   help="Residue name of the inhibitor to keep (e.g. LIG). "
                        "Default: keep all non-polymer HETATM in the receptor chain.")
    p.add_argument("--out-dir", default="out", help="Output directory.")
    # argparse sees args after PyMOL's '--'
    return p.parse_args(argv)


def load_complex(spec, name):
    if os.path.isfile(spec):
        cmd.load(spec, name)
    elif len(spec) == 4 and spec.isalnum():
        cmd.fetch(spec, name, type="pdb", async_=0)
    else:
        sys.exit(f"--complex '{spec}' is neither a file nor a 4-letter PDB id.")
    if cmd.count_atoms(name) == 0:
        sys.exit(f"Loaded nothing from '{spec}'.")


def pick_receptor_chain(obj, chain):
    if chain:
        return chain
    counts = {}
    for ch in cmd.get_chains(obj):
        counts[ch] = cmd.count_atoms(f"{obj} and chain {ch} and name CA and polymer.protein")
    if not counts or max(counts.values()) == 0:
        sys.exit("No protein chains found.")
    best = max(counts, key=counts.get)
    print(f"[receptor] auto-picked chain {best} ({counts[best]} CA atoms)")
    return best


def read_target_sequence(spec):
    if os.path.isfile(spec):
        seq = []
        with open(spec) as fh:
            for line in fh:
                if line.startswith(">"):
                    continue
                seq.append(line.strip())
        s = "".join(seq)
    else:
        s = spec.strip()
    s = "".join(c for c in s.upper() if c.isalpha())
    if not s:
        sys.exit("Empty target sequence.")
    return s


def template_residues(obj, chain):
    """Ordered (resi, one-letter) for the receptor chain, plus the sequence string."""
    residues = []
    cmd.iterate(
        f"{obj} and chain {chain} and polymer.protein and name CA",
        "residues.append((resi, resn))",
        space={"residues": residues},
    )
    ordered, seq = [], []
    for resi, resn in residues:
        one = seq1(resn.capitalize()) or "X"
        ordered.append(resi)
        seq.append(one)
    return ordered, "".join(seq)


def align(template_seq, target_seq):
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aln = aligner.align(template_seq, target_seq)[0]
    # aln.aligned -> ((tpl_blocks), (tgt_blocks)); rebuild column-wise mapping.
    tpl_idx = tgt_idx = 0
    columns = []  # (template_pos or None, target_char or None)
    t_str, q_str = str(aln[0]), str(aln[1])
    for tc, qc in zip(t_str, q_str):
        columns.append(
            (None if tc == "-" else tpl_idx,
             None if qc == "-" else qc),
        )
        if tc != "-":
            tpl_idx += 1
        if qc != "-":
            tgt_idx += 1
    return aln, columns


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else sys.argv[1:]
    args = parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)

    # 1. reference complex ---------------------------------------------------
    load_complex(args.complex, "ref")
    cmd.remove("solvent")
    chain = pick_receptor_chain("ref", args.receptor_chain)

    # isolate receptor (protein of the chosen chain) and the inhibitor
    cmd.create("template", f"ref and chain {chain} and polymer.protein")
    lig_sel = (f"ref and chain {chain} and resn {args.ligand_resn}"
               if args.ligand_resn
               else f"ref and chain {chain} and not polymer and not solvent")
    if cmd.count_atoms(lig_sel):
        cmd.create("inhibitor", lig_sel)
        print(f"[inhibitor] kept {cmd.count_atoms('inhibitor')} atoms")
    else:
        print("[inhibitor] none found (continuing without a ligand)")

    ordered, tpl_seq = template_residues("template", chain)
    print(f"[template] chain {chain}: {len(tpl_seq)} residues")

    # 2/3. target + alignment -----------------------------------------------
    tgt_seq = read_target_sequence(args.target)
    print(f"[target]   {len(tgt_seq)} residues")
    aln, columns = align(tpl_seq, tgt_seq)
    identity = sum(1 for a, b in zip(str(aln[0]), str(aln[1])) if a == b and a != "-")
    aligned_len = sum(1 for a, b in zip(str(aln[0]), str(aln[1])) if a != "-" and b != "-")
    print(f"[align]    {identity}/{aligned_len} identical over aligned "
          f"({100*identity/max(aligned_len,1):.0f}% id)")

    # 4. thread: build the variant onto the template backbone ----------------
    cmd.create("variant", "template")
    cmd.set("retain_order", 0)

    to_mutate = []   # (resi, target_three_letter)
    to_delete = []   # resi with no target residue (deletion)
    n_insert = 0
    for tpl_pos, tgt_char in columns:
        if tpl_pos is None:            # insertion in target -> can't build backbone
            n_insert += 1
            continue
        resi = ordered[tpl_pos]
        if tgt_char is None:           # deletion in target
            to_delete.append(resi)
            continue
        if tgt_char != tpl_seq[tpl_pos]:
            to_mutate.append((resi, seq1_to_three(tgt_char)))

    # mutate
    cmd.wizard("mutagenesis")
    cmd.refresh_wizard()
    n_mut = 0
    for resi, three in to_mutate:
        if three is None:
            continue
        cmd.get_wizard().set_mode(three)
        # NB: the mutagenesis wizard only commits with an object-qualified
        # selector of the form /object//chain/resi/ (headless PyMOL).
        cmd.get_wizard().do_select(f"/variant//{chain}/{resi}/")
        cmd.get_wizard().apply()
        n_mut += 1
    cmd.set_wizard()

    # delete
    for resi in to_delete:
        cmd.remove(f"variant and chain {chain} and resi {resi}")

    print(f"[thread]   {n_mut} mutated, {len(to_delete)} deleted, "
          f"{n_insert} target insertions skipped (not built)")

    # 5. overlay + outputs ---------------------------------------------------
    # variant was built from template, so it already shares the frame; align to
    # confirm and to place it if you later relax the backbone.
    rms = cmd.align("variant", "template")[0]
    print(f"[overlay]  variant vs template CA RMSD after align: {rms:.3f} A")

    tpl_out = os.path.join(args.out_dir, "template.pdb")
    var_out = os.path.join(args.out_dir, "variant_model.pdb")
    ses_out = os.path.join(args.out_dir, "overlay.pse")
    cmd.save(tpl_out, "template")
    cmd.save(var_out, "variant")
    if cmd.count_atoms("inhibitor"):
        cmd.save(os.path.join(args.out_dir, "inhibitor.pdb"), "inhibitor")

    # a ready-to-open overlay session
    cmd.hide("everything")
    cmd.show("cartoon", "template or variant")
    cmd.color("gray70", "template")
    cmd.color("marine", "variant")
    if cmd.count_atoms("inhibitor"):
        cmd.show("sticks", "inhibitor")
        cmd.color("yellow", "inhibitor")
    cmd.save(ses_out)
    print(f"[done]     wrote:\n  {tpl_out}\n  {var_out}\n  {ses_out}")


def seq1_to_three(one):
    table = {
        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
        "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
        "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
        "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
    }
    return table.get(one.upper())


if __name__ == "__main__" or "pymol" in sys.modules:
    main()
