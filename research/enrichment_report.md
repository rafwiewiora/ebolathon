# TS vs random enrichment (statistical control)

**Question:** does the synthon-TS strategy enrich for good binders vs random draws
from the same makeable, drug-like CORE space, docked in the same 5JQ7/T0R box?

**Control:** 399 random drug-like molecules (`sample_space(strategy="random")`, same
MW/cLogP/TPSA/HBD/HBA/rotB/QED window + $≤200 price), docked identically (QVina2,
exh 8, 5JQ7 frame). **Test:** all 4,513 TS-campaign docked scores.

| metric | random (n=399) | TS full (n=4,513) |
|---|---|---|
| best | −9.6 | **−10.3** |
| median | **−7.30** | −7.00 |
| mean | **−7.31** | −6.96 |
| frac < −8.5 | 6.0% | 9.1% (1.5×) |
| frac < −9.0 | 0.8% | 3.0% (**4.0×**) |
| frac < −9.5 | 0.3% | 0.5% (2.0×) |

**Stats:** Mann-Whitney U (TS<random) p≈1.0 · Cliff's δ = −0.18 · KS = 0.26, p ≈ 7e-23.

## Verdict (honest)
- TS does **not** beat random on the **average** — median/mean are slightly *worse*,
  the cost of exploration (TS deliberately docks off-target molecules to learn).
- TS **does** enrich the **good tail** — ~4× more sub-−9.0 hits and a better champion
  (−10.3 vs −9.6); the distributions differ significantly (KS p≈7e-23).
- For a screen (where only the top hits matter) that tail enrichment is the win —
  real and significant, but modest. The earlier "43×" figure was a biased artifact
  (TS top-20 vs random full) and is retracted.

## Caveats / next
- The TS "full" pool is diluted by seed rounds (~random) + exploration. A cleaner
  test of the *steering* is **elaboration-round scores only** vs random.
- More exploit-weighted windows (raise `exploit_min_sim`, lower explore ceiling)
  would trade tail-breadth for a better median.
- A sharper scorer (GNINA/MM-GBSA on the tail) would tighten the ranking.
