# onepot direct API

Direct client for `https://api.onepot.ai` — bypasses the muni wrapper (and its
`credits_remaining` integration bug). The `onepot/` package here is vendored from
`pip install onepot` (v0.2.1); `pip install onepot` also works. The only runtime
dep is `httpx==0.28.1` — install it (or just `pip install onepot`) before running.

## Auth
Header `X-API-Key: <key>`. Set the key in the environment (do not hardcode):

```bash
export ONEPOT_API_KEY=<your-onepot-api-key>
```

## Endpoints (all verified live 2026-07-18)
- `client.search(smiles_list=[...], max_results, decompose, bb_filters, substructure_search, exact_lookup, max_price, ...)` → `POST /v1/search`. 1 credit/query. Hits carry `similarity`, `price_usd`, and (when enumerated) `reaction_class` + `bbs`. `decompose=True` adds a `decompositions` block per query.
- `client.search_stream(smiles=...)` → `POST /v1/search/stream`, SSE, single molecule.
- `client.sample_space(count, strategy="diverse"|"random", properties={...}, include_properties, exclude_generic_scaffolds, seed)` → `POST /v1/space/sample`. Samples makeable CORE molecules (up to 10k) with property filters. Direct Scenario-B sampler.
- `client.order(smiles=[...], email=..., notes=...)` → `POST /v1/order`.

Notes: `max_depth` fixed at 1. `max_results` default 100. Price tiers $125 / $295.
`bb_filters` = per-`(reaction_class, bb_index)` Tanimoto `min/max` windows applied
before enumeration — the knob the synthon-TS loop uses. See `example.py`.
