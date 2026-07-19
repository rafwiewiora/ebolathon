"""Derive a docking box from a *pocket ligand* already bound in a structure.

The one-command entrypoint (`run.py`) takes a target PDB + the residue name of a
ligand sitting in the pocket of interest, and turns that into the Vina/onepot
`[[cx,cy,cz],[sx,sy,sz]]` box that both docking backends consume.

Key idea: the box is derived from the **raw** structure (which still contains the
ligand). The receptor is prepared/stripped separately (Rowan `prepare()` removes
heteroatoms), but coordinates are preserved, so a box computed here stays valid
against the prepared receptor.

`pdb` may be a 4-char PDB id (downloaded from RCSB) or a local `.pdb`/`.cif`
path. Only the standard coordinate columns are parsed — no external deps.
"""
from __future__ import annotations

import os
import urllib.request


# --------------------------------------------------------------------------- #
# Fetching
# --------------------------------------------------------------------------- #
def _looks_like_pdb_id(s: str) -> bool:
    return len(s) == 4 and s[0].isdigit() and s.isalnum()


def _fetch_structure(pdb: str) -> tuple[str, str]:
    """Return (text, fmt) where fmt is 'pdb' or 'cif'. Accepts a PDB id (fetched
    from RCSB) or a local .pdb/.cif file path."""
    if os.path.exists(pdb):
        fmt = "cif" if pdb.lower().endswith((".cif", ".mmcif")) else "pdb"
        with open(pdb) as fh:
            return fh.read(), fmt
    if _looks_like_pdb_id(pdb):
        url = f"https://files.rcsb.org/download/{pdb.upper()}.pdb"
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode("utf-8", "replace"), "pdb"
        except Exception as e:  # noqa: BLE001
            raise RuntimeError(
                f"could not download PDB {pdb!r} from RCSB ({url}): {e}") from e
    raise FileNotFoundError(
        f"{pdb!r} is neither an existing file nor a 4-char PDB id")


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #
def _iter_pdb_atoms(text: str):
    """Yield (resname, chain, x, y, z) for every ATOM/HETATM record in a PDB."""
    for line in text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            resname = line[17:20].strip()
            chain = line[21:22].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            continue
        yield resname, chain, x, y, z


def _iter_cif_atoms(text: str):
    """Yield (resname, chain, x, y, z) for every _atom_site record in an mmCIF.

    Parses the `_atom_site` loop generically by mapping the declared column
    order, so it tolerates PDBx column layouts."""
    lines = text.splitlines()
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        if line.strip() == "loop_":
            # collect the tag block that follows
            tags = []
            j = i + 1
            while j < n and lines[j].lstrip().startswith("_"):
                tags.append(lines[j].strip())
                j += 1
            if any(t.startswith("_atom_site.") for t in tags):
                col = {t: k for k, t in enumerate(tags)}
                def _idx(*cands):
                    for c in cands:
                        if c in col:
                            return col[c]
                    return None
                ci_res = _idx("_atom_site.label_comp_id", "_atom_site.auth_comp_id")
                ci_ch = _idx("_atom_site.auth_asym_id", "_atom_site.label_asym_id")
                ci_x = _idx("_atom_site.Cartn_x")
                ci_y = _idx("_atom_site.Cartn_y")
                ci_z = _idx("_atom_site.Cartn_z")
                k = j
                while k < n:
                    row = lines[k]
                    s = row.strip()
                    if s == "" or s.startswith("#") or s.startswith("loop_") \
                            or s.startswith("_") or s.startswith("data_"):
                        break
                    parts = row.split()
                    try:
                        resname = parts[ci_res]
                        chain = parts[ci_ch] if ci_ch is not None else ""
                        x = float(parts[ci_x]); y = float(parts[ci_y]); z = float(parts[ci_z])
                    except (ValueError, IndexError, TypeError):
                        k += 1
                        continue
                    yield resname, chain, x, y, z
                    k += 1
                i = k
                continue
        i += 1


def _hetatm_resnames(text: str, fmt: str) -> list[str]:
    it = _iter_cif_atoms(text) if fmt == "cif" else _iter_pdb_atoms(text)
    common_solvent = {"HOH", "WAT", "DOD"}
    names = {}
    for resname, _chain, _x, _y, _z in it:
        if resname in common_solvent:
            continue
        names[resname] = names.get(resname, 0) + 1
    # order by descending atom count; hetero ligands tend to be mid-sized
    return sorted(names, key=lambda k: -names[k])


# --------------------------------------------------------------------------- #
# Public API
# --------------------------------------------------------------------------- #
def box_from_ligand(pdb: str, resname: str, chain: str | None = None,
                    padding: float = 8.0) -> list:
    """Compute a docking box centred on the bound ligand `resname`.

    Parameters
    ----------
    pdb : 4-char PDB id (downloaded from RCSB) or a local .pdb/.cif path.
    resname : residue name of the pocket ligand (e.g. ``"ATP"``).
    chain : optional chain id to disambiguate (e.g. ``"A"``).
    padding : Angstroms added on every side beyond the ligand extent.

    Returns
    -------
    ``[[cx, cy, cz], [sx, sy, sz]]`` — center (mean of ligand coords) and size
    (per-axis extent + 2*padding), the format both docking backends consume.
    """
    text, fmt = _fetch_structure(pdb)
    resname = resname.strip().upper()
    it = _iter_cif_atoms(text) if fmt == "cif" else _iter_pdb_atoms(text)
    xs, ys, zs = [], [], []
    for rn, ch, x, y, z in it:
        if rn.upper() != resname:
            continue
        if chain and ch and ch != chain:
            continue
        xs.append(x); ys.append(y); zs.append(z)
    if not xs:
        avail = _hetatm_resnames(text, fmt)
        hint = (", ".join(avail[:25]) + (" ..." if len(avail) > 25 else "")) \
            if avail else "(none found)"
        raise ValueError(
            f"no atoms with residue name {resname!r}"
            + (f" on chain {chain!r}" if chain else "")
            + f" in {pdb!r}. Available non-water HET/residue names: {hint}")
    cx = sum(xs) / len(xs); cy = sum(ys) / len(ys); cz = sum(zs) / len(zs)
    sx = (max(xs) - min(xs)) + 2 * padding
    sy = (max(ys) - min(ys)) + 2 * padding
    sz = (max(zs) - min(zs)) + 2 * padding
    return [[cx, cy, cz], [sx, sy, sz]]
