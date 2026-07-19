"""Docking oracle backend #1 — muni `rowan_batch_docking` (QVina2).

Shells out to the muni CLI. Runs on muni's Rowan integration and bills muni
credits (~0.05 cr/ligand at exhaustiveness 8, measured on the CDK2/1HCK
calibration). No Rowan API key needed — muni holds it. This is the version that
"runs through muni": another agent on the shared muni workspace can drive it.

Requires: the `muni` CLI installed, logged in, and (ideally) a page bound via
`muni page use` so trial jobs land on the right page.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
import time

from .core import Target

_JOB_RE = re.compile(r"(job_[0-9a-f-]{8,})")


class MuniBatchDockingOracle:
    def __init__(self, target: Target, muni_bin: str | None = None,
                 timeout_s: int = 7200, chunk: int = 50, page_id: str | None = None,
                 poll_interval: int = 30, name: str = "synthon-TS dock", log=print):
        self.t = target
        self.muni = muni_bin or os.environ.get("MUNI_BIN", "muni")
        self.timeout_s = timeout_s      # muni batch docks take ~1h; wait long enough
        self.chunk = chunk
        self.page_id = page_id
        self.poll_interval = poll_interval
        self.name = name
        self.log = log

    # oracle contract: lower score = better; template unused (box defines site)
    _TERMINAL = ("completed", "completed_ok", "success", "failed", "error", "cancelled")

    def score(self, smiles_list, template_smiles=None) -> dict:
        chunks = [smiles_list[i:i + self.chunk]
                  for i in range(0, len(smiles_list), self.chunk)]
        # 1) submit ALL chunks up front so muni runs them CONCURRENTLY (muni
        #    parallelises jobs; submitting serially wasted ~1h per chunk).
        jobs = []  # (chunk, job_id)
        for ch in chunks:
            jid = self._submit(ch)
            if jid:
                jobs.append((ch, jid))
        self.log(f"[muni-oracle] submitted {len(jobs)}/{len(chunks)} jobs concurrently "
                 f"({sum(len(c) for c, _ in jobs)} ligands); polling every {self.poll_interval}s")
        # 2) poll them all together until each is terminal
        pending = {jid for _, jid in jobs}
        deadline = time.time() + self.timeout_s
        while pending and time.time() < deadline:
            time.sleep(self.poll_interval)
            for jid in list(pending):
                s = subprocess.run([self.muni, "job", "status", jid, "--json"],
                                   capture_output=True, text=True)
                try:
                    status = str(json.loads(s.stdout).get("status", "")).lower()
                except Exception:
                    continue
                if status in self._TERMINAL:
                    pending.discard(jid)
            self.log(f"[muni-oracle] {len(jobs) - len(pending)}/{len(jobs)} jobs done")
        # 3) read scores from every job
        out: dict = {}
        for ch, jid in jobs:
            out.update(self._read_scores(jid, ch))
        return out

    def _submit(self, batch) -> str | None:
        """Submit one batch-docking job (no --follow); return its job id."""
        params = {
            "smiles_list": batch, "pocket": self.t.pocket,
            "executable": self.t.executable,
            "scoring_function": self.t.scoring_function,
            "exhaustiveness": self.t.exhaustiveness, "name": self.name,
        }
        if _looks_like_uuid(self.t.protein):
            params["protein_uuid"] = self.t.protein
        else:
            params["protein"] = self.t.protein
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
            json.dump(params, f)
            pfile = f.name
        cmd = [self.muni, "run", "rowan_batch_docking", "--params-file", pfile,
               "--json", "--title", self.name]
        if self.page_id:
            cmd += ["--page-id", self.page_id]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        jid = _extract_job_id(proc.stdout) or _extract_job_id(proc.stderr)
        if not jid:
            self.log(f"[muni-oracle] submit failed: {proc.stdout[-200:]}")
        return jid

    def _read_scores(self, job_id: str, batch) -> dict:
        q = subprocess.run([self.muni, "job", "query", job_id, "--json"],
                           capture_output=True, text=True)
        try:
            row = json.loads(q.stdout)["rows"][0]
        except Exception as e:
            self.log(f"[muni-oracle] parse fail for {job_id}: {e}")
            return {}
        scores = row.get("best_scores") or []
        sent = row.get("initial_smiles_list") or batch
        out = {}
        for smi, sc in zip(sent, scores):
            if sc is not None:
                out[smi] = float(sc)
        return out


def _looks_like_uuid(s: str) -> bool:
    return isinstance(s, str) and s.count("-") == 4 and len(s) >= 32


def _extract_job_id(text: str) -> str | None:
    # muni --follow --json interleaves "Status: ..." lines with a pretty-printed
    # JSON block; a regex for the job_ token is the most robust extraction.
    m = _JOB_RE.search(text or "")
    return m.group(1) if m else None
