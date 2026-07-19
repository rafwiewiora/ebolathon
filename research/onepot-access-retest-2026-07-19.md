# onepot access — retest report

**Date:** 2026-07-19 (UTC; jobs timestamped 00:02–00:03 UTC)
**Tester:** rafwiewiora@gmail.com · muni CLI 0.1.61 · space Ebolathon
**Results page:** *Synthon-TS oracle calibration (CDK2 1HCK)* (`page_bd72c0e2-0799-44c3-b45b-1d4ffddaad17`)

## TL;DR

onepot is **still down** for us, with the exact same failure seen on 2026-07-18:

```
error: "OnePot response missing required credits_remaining"
```

It is **not** a balance problem (59.9 credits available), **not** query-specific, and **not**
fixable by us client-side — muni CLI 0.1.61 is already the latest release on PyPI. The onepot
search itself runs (each attempt gets a provider job id and burns 2–5 s of provider time) but the
response envelope muni receives back from onepot is missing the `credits_remaining` field that
muni's gateway treats as required, so muni marks the job failed on response validation. This is a
**muni↔onepot integration contract break on the provider/gateway side.**

## What was tested

Three independent live calls, varying query and params to isolate the cause:

| Job ID | Query | Params | Provider job id | Runtime | Result |
|---|---|---|---|---|---|
| `job_33fba5bf` | ibuprofen `CC(C)Cc1ccc(cc1)C(C)C(=O)O` | `max_results=10, decompose=true` | `op-b499e470` | 4.58 s | **failed** — missing `credits_remaining` |
| `job_83cebfa9` | benzene `c1ccccc1` | `max_results=5` (minimal, no decompose) | `op-83a07b80` | 2.77 s | **failed** — missing `credits_remaining` |

(Plus the two pre-existing failed onepot jobs from the 2026-07-18 session, same error.)

## Diagnosis

1. **Reaches the provider.** Every attempt is assigned an `op-…` provider job id and consumes
   real provider runtime (2.8–4.6 s), so the request is authenticated and the search executes on
   onepot's side. The failure is on the *response*, not the request.
2. **Universal, not query- or param-specific.** The minimal benzene / `max_results=5` /
   no-decompose call fails identically to the richer decompose call. Rules out our query, our
   `bb_filters`/`decompose` usage, and result-size limits.
3. **Not a credit/balance issue.** Balance is 59.9 credits; the error is a *missing field named*
   `credits_remaining` in onepot's reply, i.e. a response-schema mismatch — not "insufficient
   credits."
4. **Not client-fixable.** Installed muni is 0.1.61 and PyPI latest is also 0.1.61 — no newer CLI
   to upgrade into. The validation that rejects the response lives in the muni gateway wrapping
   `onepot==0.2.0`, which we don't control.

**Most likely root cause:** onepot's API response format changed (or this account's plan returns a
response without `credits_remaining`), while muni's gateway still requires that field. Either
onepot must restore the field, or muni must relax it to optional.

## Impact on the synthon-TS project

onepot is the **retriever leg** of the loop (seed pool + `decompose`/`bb_filters` per-position BB
search). With it down, the loop's chemistry-generation half cannot run live. This does **not** block:

- the **oracle leg** (Rowan pocket detection + `rowan_batch_docking`/QVina2), which is working — see
  the CDK2/1HCK calibration on the same page; or
- **scaffolding the full loop with onepot mocked**, so it runs the moment the endpoint returns.

## Action / ask for onepot + muni

Please flag to the onepot/muni team:

> Every `muni run onepot …` (any SMILES, with or without `decompose`) fails with
> `OnePot response missing required credits_remaining`. The provider assigns an `op-…` job id and
> the search runs (~2–5 s), so the request path is fine — the returned envelope is missing the
> `credits_remaining` field muni's gateway requires. Balance is healthy (59.9 credits) and muni CLI
> is at latest (0.1.61). Fix on the provider/gateway side: restore `credits_remaining` in onepot's
> response, or make it optional in muni's onepot client. Reproducible provider job ids:
> `op-b499e470-8bfa-4a8d-95f6-555c6e60fe74`, `op-83a07b80-63eb-4f5b-8a82-526c617bae3a`.
