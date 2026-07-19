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
    def score(self, smiles_list, template_smiles=None) -> dict:
        out: dict = {}
        for i in range(0, len(smiles_list), self.chunk):
            batch = smiles_list[i:i + self.chunk]
            out.update(self._dock_batch(batch))
        return out

    def _dock_batch(self, batch) -> dict:
        params = {
            "smiles_list": batch,
            "pocket": self.t.pocket,
            "executable": self.t.executable,
            "scoring_function": self.t.scoring_function,
            "exhaustiveness": self.t.exhaustiveness,
            "name": self.name,
        }
        # thread a prepared protein UUID if given, else pass a PDB id
        if _looks_like_uuid(self.t.protein):
            params["protein_uuid"] = self.t.protein
        else:
            params["protein"] = self.t.protein

        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
            json.dump(params, f)
            pfile = f.name

        # Submit WITHOUT --follow: muni's --follow drops the connection on the
        # ~1-hour batch jobs (that's what silently lost the earlier run's scores).
        # Instead grab the job id and poll `muni job status` until it's terminal.
        cmd = [self.muni, "run", "rowan_batch_docking", "--params-file", pfile,
               "--json", "--title", self.name]
        if self.page_id:
            cmd += ["--page-id", self.page_id]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        job_id = _extract_job_id(proc.stdout) or _extract_job_id(proc.stderr)
        if not job_id:
            self.log(f"[muni-oracle] no job_id; stdout tail: {proc.stdout[-300:]}")
            return {}
        self.log(f"[muni-oracle] submitted {job_id} ({len(batch)} ligands); polling...")
        deadline = time.time() + self.timeout_s
        while time.time() < deadline:
            time.sleep(self.poll_interval)
            st = subprocess.run([self.muni, "job", "status", job_id, "--json"],
                                capture_output=True, text=True)
            try:
                status = str(json.loads(st.stdout).get("status", "")).lower()
            except Exception:
                continue
            if status in ("completed", "completed_ok", "success"):
                break
            if status in ("failed", "error", "cancelled"):
                self.log(f"[muni-oracle] job {job_id} {status}")
                return {}
        else:
            self.log(f"[muni-oracle] {job_id} not terminal after {self.timeout_s}s; reading anyway")
        return self._read_scores(job_id, batch)

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
