# Run notes — revision on RStudio Server

Working document for the revision run. Status doc for this project: record what
was run, what came out, and what it decides.

Started: 2026-10-29

---

## 1. How to run

```bash
cd global_coastal_bc_analysis
git pull                              # brings data/bc_analysis.rds
bash scripts/get_external_data.sh     # skips anything already present
Rscript R/blue_carbon_full.R 2>&1 | tee output/run_$(date +%Y%m%d).log
```

Run from the repository root. The script resolves `data/`, `output/` and
`R/02_reviewer_analyses.R` against the **working directory**, and if it cannot
find them it prints a skip message and carries on rather than stopping — so a
run from the wrong directory looks like it worked.

`data/bc_analysis.rds` and `data/ecklonia_records.csv` both ship with the
repository; `R/01_ecklonia_gbif.R` only needs rerunning to refresh the latter
from GBIF. The log's first lines print `Occurrence table: data/bc_analysis.rds`
— if it names the root copy instead, `data/bc_analysis.rds` is missing (the
pull did not bring it) and the run is using whatever old copy sits in the root.

Budget: the main pipeline as before, plus 10–20 min for PART 8 (≈165 maxnet
fits: 3 ecosystems × 11 model specifications × 5 fits each).

### Checkpoints to watch in the log

| message | meaning |
|---|---|
| `Ecklonia: N recovered coordinates, M with complete environment` | the merge landed; M should be in the hundreds |
| `bc after Ecklonia merge: ... records` | total should exceed 2,040 |
| `=> best distance form: ...` | drives the `DIST_VAR` decision below |
| `Pearson r(log area, log discharge)` | how much area and discharge actually differ |

---

## 2. Run log

All runs on the RTX server (conda env `gcbc`: R 4.5.3, terra 1.9.50, sf 1.1.3,
maxnet 0.1.4, xgboost 3.4.2.1, ranger 0.18.0; EOT20 tides present).

| date | what | status | notes |
|---|---|---|---|
| 2026-10-29 | run 1 — revision as pushed (`eb145fe`), no cap | done (6 min) | exposed the scenario bug (§3e) and *E. cava* swamping kelp (§3f); outputs kept in `output_run1_prefix/` |
| 2026-10-29 | v5 code (`f5593e1`) on the same server | done | reference for "what the revision changed" (§3g); worktree `../gcbc_v5_baseline` |
| 2026-10-29 | run 2 — scenario fix + per-species cap 100 | done | **main results**, `output/` |
| 2026-10-29 | sensitivity — cap 200, no cap | done | `output_cap200/`, `output_capInf/` |
| 2026-10-29 | sensitivity — cap 50 | done | `output_cap50/` |
| 2026-10-29 | `R/03_robustness.R`, 50 bootstrap replicates, caps 50/100/200/none | done | §3h–3i; `TabS_robust_*.csv` |
| 2026-10-29 | v5 code + terra 1.8.86 (side library) | done | reproduces the v5 manuscript (§3g); worktree `../gcbc_v5_terra18` |
| 2026-10-30 | `R/04_structural.R`, co-author's structural model, 50 replicates | done | §3k |
| 2026-10-30 | main pipeline with `BC_FEAT=struct` (S, σk, log D), cap 100 | done | §3l; `output_struct/` |
| 2026-10-30 | 10 background draws × {current, struct}, 20 runs | done | §3m; `output_bgseeds/` |
| 2026-10-30 | struct made default; main run → `output/`; figures compared with v5 | done | report §図の比較 |
| 2026-10-30 | σk definitions STx and T, ten draws each | done | §3n |

Environment notes: conda-forge GDAL needs `libgdal-hdf5 libgdal-netcdf` or terra
cannot open the EOT20 `.nc` files (run halts at the tidal step). Running three
pipelines at once needs `OMP_NUM_THREADS` set — xgboost otherwise takes every
core in each process and the runs thrash. EOT20 downloads without
registration from https://www.seanoe.org/data/00683/79489/data/85762.zip; the
`ocean_tides.zip` the script wants is the inner of the two zips in it.

Checkpoints (run 2):

| message | value |
|---|---|
| `Ecklonia: N recovered coordinates, M with complete environment` | 1,447 → 1,388 |
| `bc after Ecklonia merge` | 3,428 (kelp 1,843); 3,049 after the depth filter |
| `Per-species cap (100)` | *E. cava* 1,273 → 100; *Z. marina* 173 → 100; bc 3,049 → 1,803 |
| `=> best distance form` | `exp_dist`, by 0.002 — within noise, see §3b |
| `Pearson r(log area, log discharge)` | 0.764 |

---

## 3. Results

Main results are run 2 (cap 100). Spatial-block CV AUC throughout.

### 3a. Turbidity decomposition

**Question.** Kd490 at a coastal cell sums ocean optics (phytoplankton,
resuspension — both latitudinal) and terrestrial delivery (CDOM, mineral
sediment). How much is each?

**Reads.** `output/TabS_Kd490_variance_partition.csv`,
`output/TabS_Kd490_residual_model.csv`, `output/FigR1_Turbidity_Decomposition`

| component | deviance explained |
|---|---|
| Total (climate + river) | 40.7% |
| Unique to climate | 19.1% |
| Unique to river | 23.0% |
| Shared | −1.4% (suppression; read as ~0) |

Residual model (`kd_resid ~ log_discharge + log_dist`, n = 2,986), R² = 0.224

| term | estimate | p | ×per 10-fold |
|---|---|---|---|
| log_discharge | 0.151 | 6 × 10⁻¹⁶ | 1.42 |
| log_dist | −0.442 | 1 × 10⁻¹⁵³ | 0.36 |

Stable across caps (river-unique 22.2–23.0%, elasticities 1.39–1.42 and
0.36–0.38).

**Decision.** The river-unique share is as large as the climate-unique share.
The Discussion can state that roughly a quarter of the Kd490 variation is
attributable to river terms alone, rather than hedging that turbidity is hard
to interpret. This is the decomposition the review asked for (Kd ~ SST + |lat|,
then residual ~ river size + distance on the log scale, i.e. multiplicative).

---

### 3b. Distance functional form

**Reads.** `output/TabS_distance_form_comparison.csv`, `FigR2` panel a

| form | seagrass | mangrove | kelp | mean |
|---|---|---|---|---|
| `log_dist` log₁₀(d+1) | 0.675 | 0.689 | 0.697 | 0.687 |
| `inv_dist` 1/(d+1) | 0.669 | 0.684 | 0.694 | 0.682 |
| `exp_dist` exp(−d/50) | 0.683 | 0.677 | 0.706 | 0.689 |

Best: `exp_dist` by 0.002 over `log_dist`; per-ecosystem SD 0.05–0.09. The
ranking flips between runs (no cap: `log_dist` first, cap 200: `exp_dist`
first by 0.0005).

**Decision.** Keep `DIST_VAR = "log_dist"`. The three forms are
indistinguishable (spread ≤ 0.007 in mean AUC), so switching on a 0.002
difference would be fitting noise. One line in Methods: the result does not
depend on the transform; 1/d, the form suggested in review for radial
dilution, performs the same.

---

### 3c. River size: catchment area vs discharge

**Reads.** `output/TabS_riversize_comparison.csv`

r(log area, log discharge) = 0.764

| predictor set | seagrass | mangrove | kelp |
|---|---|---|---|
| area only | 0.651 | 0.700 | 0.698 |
| discharge only | 0.675 | 0.689 | 0.697 |
| both | 0.652 | 0.744 | 0.688 |
| neither | 0.667 | 0.692 | 0.712 |

Discharge importance: 1.3–11.6% in the ecosystem MaxEnt models, 1.6–2.5%
median across thinned species models.

**Decision.** "Proximity matters, river size does not" stands, now on the
better measure: no row beats `neither` beyond the CV noise, and discharge
importance is as small as area's was.

One thing to mention but not claim: mangrove `both` (0.744) is 0.05 above
`neither`, about one SD, and identical in every run (mangrove data are
untouched by the cap). Area conditional on discharge is in effect inverse
runoff — an aridity signal. That is the nearest the data come to the review's
point about catchment composition (vegetated vs barren).

---

### 3d. One stage or two

**Reads.** `output/TabS_onestage_vs_twostage.csv`, `FigR2` panel b

Run 2 changed the importance to 10 permutations per variable and kept the raw
AUC loss (`dAUC`) beside the share. The share is not comparable across
designs: the one-stage model has six predictors (SST and lat included), stage 2
has four, so distance's *share* is larger in stage 2 partly by arithmetic. The
raw ΔAUC has no such dependence.

AUC lost when river distance is permuted (raw ΔAUC):

| design | seagrass | mangrove | kelp |
|---|---|---|---|
| 1-stage, ocean-wide bg | 0.061 | 0.053 | 0.051 |
| 1-stage, Kd residual | 0.043 | 0.053 | 0.052 |
| 2-stage (stage 2) | 0.143 | 0.101 | 0.168 |
| 2-stage, Kd residual | 0.124 | 0.134 | 0.174 |

Share of total importance, for reference: 12–18% one-stage, 12–15% one-stage
with Kd residual, 50–72% stage 2.

**Decision.** The central methodological claim holds on the raw scale:
distance is worth two to three times more AUC within the climatic range than in
a one-stage model. The review's proposal — remove the temperature signal from
turbidity and fit one stage — does not recover it (ΔAUC 0.04–0.05, no better
than the plain one-stage model). The reason is that a one-stage model is
dominated by the range question (where the climate allows the group at all),
and the within-range river gradient is a small part of that discrimination
whatever is done to turbidity. Same pattern in every cap setting.

> **Superseded by §3h.** These are training-data numbers from one background
> draw. Out of sample with 50 bootstrap replicates the gap is about 2×, not
> 2–3×, and the refit-based contrast includes zero.

---

### 3e. Scenario bug (new)

PART 7 rescales `dist_km` for the "dist↑" pathway but the models read
`log_dist`, which was never recomputed. In v5 (`f5593e1`) and in run 1,
"dist↑ only" therefore returned exactly zero change for all three ecosystems,
and every combined scenario carried the Kd490 pathway alone. Fixed in both the
ecosystem and species scenario loops (`wg[[DIST_VAR]] <- km_to_dist(...)`).

Effect, run 2: "dist↑ only (50%)" now lowers mean suitability by 0.006–0.010,
the same order as "Kd490↓ only" (0.007–0.015). 70% regulation, mean % change
across species: seagrass −10.1, mangrove −12.1, kelp −5.7 (v5 code on this
server: −6.5 / −7.7 / −3.0). The regulation effects reported in v5 were
underestimates by roughly a factor of 1.5–2, and FigSc3 (distance pathway
only) was flat. This needs a sentence in the response letter.

### 3f. Per-species cap (new)

After the merge, *E. cava* was 1,273 of 1,677 kelp presences (76%), at
19.4 ± 1.6 °C and median |lat| 34° against 11.9 ± 4.0 °C and 49.6° for the other
kelps, and the main SDM is not thinned. Without a cap the kelp models describe
*E. cava*: CV AUC 0.649 ± 0.145, Stage 1 AUC up to 0.945 because the range
model learns one narrow thermal band, and an occupancy hump at 15–20 km from
Japanese/Korean observer density (`output_run1_prefix/FigD4`).

Decided 2026-10-29 (option b): cap each species at 100 records, removing the
excess by round-robin over 1° cells so the kept records stay spread. maxnet has
no case weights, hence a cap. Genus-only records are untouched. On these data
the cap touches only *E. cava* (1,273 → 100) and *Z. marina* (173 → 100).

| kelp | no cap | cap 200 | **cap 100** |
|---|---|---|---|
| presences | 1,677 | 604 | 504 |
| *E. cava* share | 76% | 33% | 20% |
| SST (°C), median \|lat\| | 17.6 ± 4.0, 34.3° | 14.2 ± 4.8, 41.5° | 13.3 ± 4.6, 44.9° |
| CV AUC | 0.649 ± 0.145 | 0.708 ± 0.099 | 0.697 ± 0.079 |
| Stage 1 AUC | 0.945 | 0.902 | 0.891 |
| distance / salinity (ecosystem MaxEnt) | 57 / 29% | 68 / 19% | 72 / 18% |
| 70% regulation | −8.3% | −7.1% | −5.7% |

Caps 100 and 200 agree on everything that matters; no cap is the outlier.
Seagrass under cap 100 (*Z. marina* capped): CV AUC 0.675 vs 0.683, Stage 1
0.793 vs 0.804 — small.

The kelp group is now temperate to subpolar rather than
"cool-temperate to subpolar" (13.3 ± 4.6 °C, median |lat| 44.9°); the v5
wording has to change.

### 3g. Stage 1 AUC drift — cause found (new)

Run 1 gave seagrass and mangrove Stage 1 AUC 0.804 / 0.816 against v5's
0.834 / 0.880, although the revision does not touch those groups' Stage 1
inputs. The v5 code itself (`f5593e1`) gives the same 0.804 / 0.816 under
terra 1.9.50. Cause: `terra::spatSample` draws a completely different set of
20,000 background points for the same `set.seed(42)` in terra 1.8.86
(RStudio Server) and 1.9.50 (here) — 1,437 vs 1,183 survive the depth filter.

Confirmed by running the v5 code with terra 1.8.86 built into a side library
(`~/R/terra18_lib`; needs `PKG_CXXFLAGS=-fpermissive` against GDAL 3.13):

| | v5 manuscript | v5 code + terra 1.8.86 | v5 code + terra 1.9.50 |
|---|---|---|---|
| Stage 1 AUC sg / mg / kp | 0.834 / 0.880 / 0.914 | 0.834 / 0.880 / 0.914 | 0.804 / 0.816 / 0.907 |
| CV AUC sg | 0.611 ± 0.104 | 0.612 ± 0.099 | 0.655 ± 0.099 |
| CV AUC mg | 0.687 ± 0.045 | 0.681 ± 0.049 | 0.700 ± 0.036 |
| CV AUC kp | 0.749 ± 0.014 | 0.742 ± 0.017 | 0.719 ± 0.068 |
| 70% regulation sg / mg / kp | −8.3 / −6.0 / −3.1% | −8.3 / −6.0 / −3.8% | −6.5 / −7.7 / −3.0% |

The point that matters for the paper: one background redraw moves Stage 1 AUC
by 0.03–0.06 and CV AUC by 0.02–0.04. Headline numbers should be the median
over several background draws, from one environment.

### 3h. Robustness of the design comparison (new, `R/03_robustness.R`)

50 replicates, each with presences and background resampled with replacement
and the four longitude blocks rotated by a random offset; everything scored on
the held-out block; designs paired within a replicate. Distance contribution
measured two ways: permutation on the held-out block, and drop-variable
(refit without distance).

Paired difference, two-stage minus one-stage, cap 100 (median [95% interval]):

| metric | seagrass | mangrove | kelp |
|---|---|---|---|
| permutation ΔAUC | +0.078 [+0.024, +0.153] | +0.061 [+0.030, +0.092] | +0.080 [−0.039, +0.131] |
| drop-variable ΔAUC | +0.013 [−0.037, +0.095] | +0.020 [−0.010, +0.050] | +0.047 [−0.071, +0.141] |
| 1-stage Kd residual minus 1-stage (perm.) | −0.008 [−0.034, +0.038] | −0.001 [−0.008, +0.004] | −0.003 [−0.008, +0.002] |

Same pattern at caps 50, 200 and none (`output_cap*/TabS_robust_design_contrasts.csv`):
the permutation contrast excludes zero for seagrass and mangrove in all four,
for kelp at cap 200 and none (borderline, p≤0 = 0.06–0.08, at 50 and 100).

**Decision.** The in-sample result of §3d (2–3× on raw ΔAUC) overstated the
case. Out of sample, distance is worth about twice as much in stage 2 by
permutation; the unique (refit) contribution points the same way but its
interval includes zero, and distance is detectable in the one-stage model too.
The manuscript's "the river signal collapses in one stage" becomes "the
two-stage design measures about twice the distance contribution". The
reviewer's Kd-residual one-stage model does not help — robust everywhere.

### 3i. Runoff depth (new)

log10(discharge / catchment area, mm per year + 1). r = −0.79 with log area,
−0.21 with log discharge. CV AUC gain over `none`, 50 replicates, cap 100:

| set | seagrass | mangrove | kelp |
|---|---|---|---|
| area | −0.012 [−0.040, +0.009] | −0.003 [−0.019, +0.012] | +0.011 [−0.032, +0.059] |
| discharge | +0.006 [−0.011, +0.041] | −0.006 [−0.020, +0.000] | −0.004 [−0.046, +0.036] |
| both | +0.004 [−0.018, +0.054] | +0.010 [−0.009, +0.027] | +0.010 [−0.025, +0.057] |
| runoff | −0.005 [−0.039, +0.028] | **+0.025 [+0.007, +0.051]** | +0.021 [−0.010, +0.074] |
| discharge + runoff | +0.003 [−0.033, +0.055] | +0.019 [−0.000, +0.043] | +0.014 [−0.013, +0.063] |

Mangrove presences sit at LOWER runoff than climate-matched background
(median 748 vs 1,090 mm/yr; kelp the reverse, 841 vs 507). So the gain is not
"arid catchments exclude mangroves", and GBIF's lean toward dry coasts
(N Australia, Red Sea, Florida) could produce it. One sentence as a lead, not a
claim. River size (area, discharge) still adds nothing anywhere.

### 3j. Cap 50 (new)

Touches eight species (incl. three mangroves). Kelp: 428 presences, *E. cava*
12%, 12.8 ± 4.5 °C, |lat| 46.5°, CV AUC 0.713, Stage 1 0.890, distance /
salinity 72 / 22%, 70% regulation −8.0%. Caps 50–200 agree; no cap is the
outlier.

---

### 3k. Co-author's structural model (new, `R/04_structural.R`, 2026-10-30)

Proposal: salinity carries land input, S = S_ref − c·B^β·exp(−D/L); turbidity
log K = a0 + a1·S + a2·T + σk with σk read as ocean-derived; SDM on T, S, σk
with D and B acting only through S. S_ref = 90th percentile of Bio-ORACLE
salinity in a 2° block (1° and 4° as checks). Bootstrap CV as in §3h, cap 100.

| step | result |
|---|---|
| (1) salinity equation, coastal grid | ΔS ∝ B^0.03, L ≈ 590 km (log fit; nls 325 km), R² 0.02. Median ΔS 0.31 PSU at 0–10 km vs 0.20 at 200–500 km — Bio-ORACLE salinity (~9 km climatology) does not resolve plumes |
| (2) log K ~ S + T | R² 0.17; σk still tracks distance (r = −0.45; σk ~ log D + log Q R² 0.22) — not ocean-only |
| (3) structural (S, σk) − current, 2-stage | −0.041 / −0.047 / −0.089 (sg / mg / kp) |
| (3) adding exp(−D/L̂) back, 1-stage | +0.040 / +0.024 / +0.045, all intervals > 0 |
| (4) doubling B, via S | < 0.03% change in suitability |
| (5) plume index log10(1 + Q·exp(−D/L)) in place of D + Q | never better than current; L = 10 km closest (−0.018 / −0.009 / −0.028), longer L worse; S, K, log D alone ≈ current |

**Update — b0 by IDW (co-author's suggestion).** Mask ocean cells within M km
of land, fill by IDW from offshore 0.5° means (radius 3°, power 2); M = 50
main, 25 and 100 as checks. With this reference the decay becomes visible:
median ΔS 0.15 PSU at 0–10 km, 0.09 at 25–50, 0.02 at 50–100, ~0 beyond 100
(90th pct 2.4 → 1.6 → 1.1 → 0.07). Nonlinear fit L = 58 km (M = 25: 36, M =
100: 92 — L follows the mask width, partly by construction). But catchment
area still does not explain ΔS (elasticity −0.005, p = 0.60), only 38% of
coastal cells are fresher than the reference, R² ≈ 0.05, B-via-S effect zero,
and every SDM contrast is unchanged (σk does not depend on b0). Written as
S = b0 − b1·B^β·exp(−D/L), b1 > 0 and the log form is well defined.

**Update — literal form and 1/D (2026-10-30).** Fitting the comment's form
literally, log(b0 − S) = log(b1)·log(B) − b2·D with no intercept, gives
log(b1) = −0.16 (it absorbs the mean level, since ΔS is mostly < 1 PSU),
b2 = 0.0042 km⁻¹ and R² = 0.046 — the same poor fit. Distance as 1/(D+1)
added to (S, σk) recovers the same AUC as log or exp (2-stage +0.045 / +0.040 /
+0.053; 1-stage +0.044 / +0.025 / +0.038, all 1-stage intervals > 0): the form
of D does not matter. (S, σk, D) matches the current model (0.664 / 0.741 /
0.698 vs 0.670 / 0.734 / 0.698) — an interpretable variable set at no cost.

**Decision.** The process model is sound but the data cannot carry it: global
salinity does not resolve river plumes, so S cannot stand in for D or B, and
σk keeps the land signal salinity misses. Distance stays in the SDM. If an
exponential form is wanted for distance, its scale is short (~10 km), not the
300–600 km that salinity implies. River size × decay adds nothing over
distance alone, consistent with §3c/§3i. Catchment vegetation (RiverATLAS)
not tested: figshare's bot check blocks a server-side download; the
2.33 GB geodatabase (figshare file 20087321) has to be fetched in a browser.

---

### 3l. Main pipeline with salinity, σk, distance (new, 2026-10-30)

`BC_FEAT=struct` sets stage 2 to Salinity, σk, log_dist; σk = residual of
log Kd490 ~ Salinity + SST (R² 0.166), refitted on presences + background and
reapplied to the grid and inside every scenario (`sigk_of()`), so the Kd490
pathway reaches it. Default stays `current`. Cap 100, `output_struct/`.

| | current (S, Kd490, log D, log Q) | struct (S, σk, log D) |
|---|---|---|
| CV AUC sg / mg / kp | 0.675 ± 0.082 / 0.689 ± 0.050 / 0.697 ± 0.079 | 0.670 ± 0.091 / 0.699 ± 0.048 / 0.699 ± 0.077 |
| Stage 1 AUC | 0.793 / 0.816 / 0.891 | same |
| ecosystem MaxEnt, distance | 56 / 42 / 72% | 70 / 53 / 74% |
| ecosystem MaxEnt, turbidity term | Kd 9 / 56 / 6% | σk 10 / 46 / 5% |
| thinned species median, sal / dist / turb | sg 14 / 14 / 8; mg 6 / 56 / 26; kp 52 / 37 / 4 | sg 12 / 20 / 6; mg 8 / 61 / 31; kp 51 / 27 / 11 |
| 70% regulation, species mean | −10.1 / −12.1 / −5.7% | −15.1 / −14.4 / −12.6% |

Response to σk is unimodal in all three groups, peaking at σk ≈ +0.2 to +0.5
(somewhat more turbid than salinity and SST predict); clearing the water near
mouths pushes cells down the left limb, which is why the regulation effects
grow, most for kelp. Same predictive skill, one variable fewer, cleaner
reading — but the regulation magnitudes depend on how turbidity enters, so
whichever set is adopted, report the other as sensitivity.

Also fixed in passing: the PART 5 partial-dependence curve for distance varied
`dist_km`, which is not a predictor, so it was flat; it now varies DIST_VAR.

---

### 3m. Ten background draws × two variable sets (new, 2026-10-30)

`BC_BG_SEED` (default 42) sets the background draw; seeds 42 and 101–109,
cap 100, PART 8 skipped (`BC_SKIP_PART8=1`), each for `current` and `struct`.
`output_bgseeds/`, summarised by `R/06_bgseed_summary.R`. Background points
after the depth filter: 1,151–1,266.

Median [min, max] over the ten draws:

| | current | struct |
|---|---|---|
| CV AUC sg | 0.651 [0.621, 0.678] | 0.657 [0.628, 0.675] |
| CV AUC mg | 0.725 [0.689, 0.749] | 0.732 [0.699, 0.758] |
| CV AUC kp | 0.737 [0.685, 0.765] | 0.731 [0.691, 0.765] |
| Stage 1 AUC sg / mg / kp | 0.789 / 0.823 / 0.894 (ranges 0.78–0.81, 0.82–0.84, 0.88–0.90) | same |
| distance importance (ecosystem MaxEnt) sg / mg / kp | 65 / 52 / 70% | 75 / 57 / 69% |
| 70% regulation sg | −10.6 [−12.1, −9.1]% | −15.1 [−17.7, −10.1]% |
| 70% regulation mg | −11.2 [−13.1, −9.1]% | −13.5 [−16.7, −11.6]% |
| 70% regulation kp | −7.2 [−9.7, −5.1]% | −11.3 [−12.6, −9.0]% |

Paired by draw (struct − current): CV AUC −0.001 / +0.009 / +0.004 (all
within ±0.03); 70% regulation −5.1 / −2.3 / −3.4 points, more negative in
9/10, 10/10, 10/10 draws.

**Decision input.** The two sets predict equally well in every draw. The
regulation effect is larger under struct consistently, not by chance of one
draw. Seed 42 sat at the low end for kelp (CV AUC 0.697 vs median 0.737):
single-draw numbers should not be quoted. Revised Table S1:
`data/TableS1_taxa_revised.csv` (`R/05_tableS1.R`) — 25 genera, 94 species,
3,428 records, 3,049 after the depth filter, 1,803 used (cap 100); Ecklonia
2 species, 1,439 records.

---

### 3n. σk definition: interaction or temperature only (new, 2026-10-30)

Co-author: land-derived material also decays on the way out at a
temperature-dependent rate, so try log K = a0 + a1·S + a2·T + a3·S·T + σk, or
σk from log K ~ T alone. `BC_SIGK` = ST (default) / STx / T; ten draws each
(`output_bgseeds/structSTx_s*`, `structT_s*`), compared by `R/07_sigk_compare.R`.

Turbidity model: ST R² 0.166; STx R² 0.191, S·T term t = −9.6; T R² 0.058.
σk still tracks distance under all three (r = −0.45 / −0.43 / −0.47); under T
it also overlaps salinity (r = −0.32). The interaction runs the other way from
the decay idea: per PSU fresher, Kd rises ~5% at 10 °C and ~10% at 28 °C, i.e.
the land–turbidity link is stronger in warm water (sediment- and CDOM-rich
tropical rivers outweigh faster decay).

| ten draws, median | ST | STx | T |
|---|---|---|---|
| CV AUC sg / mg / kp | 0.657 / 0.732 / 0.731 | 0.661 / 0.725 / 0.732 | 0.634 / 0.733 / 0.744 |
| σk importance sg / mg / kp | 9 / 40 / 4% | 11 / 36 / 1% | 2 / 41 / 6% |
| distance importance sg / mg / kp | 75 / 57 / 69% | 75 / 59 / 70% | 80 / 58 / 72% |
| 70% regulation sg / mg / kp | −15.1 / −13.5 / −11.3% | −14.0 / −12.9 / −9.7% | −11.4 / −12.4 / −7.5% |
| σk at peak suitability sg / mg / kp | 0.53 / 0.43 / 0.19 | 0.52 / 0.39 / 0.11 | 0.52 / 0.47 / −0.02 |

Paired with ST in the same draw: STx AUC ±0.01, regulation +0.6 to +1.4
points; T AUC −0.021 for seagrass (others ±0.01), regulation +1.0 to +4.0.

**Decision (agreed 2026-10-30).** Keep ST as the main definition; report STx and T as
sensitivity. The interaction is real but changes nothing material, and its
sign does not support the decay rationale; T-only σk overlaps salinity and
loses seagrass skill. Quote the regulation effect as a range across
definitions (70%: sg −11 to −15, mg −12 to −14, kp −8 to −11%).

---

## 4. Departures from v5

| quantity | v5 (manuscript) | v5 code, this server | revision (cap 100) |
|---|---|---|---|
| records after depth filter (sg / mg / kp) | — | 666 / 706 / 404 | 593 / 706 / 504 |
| records in bc_analysis.rds (sg / mg / kp) | 777 / 808 / 455 | same | + 1,388 *E. cava* before cap |
| CV AUC seagrass | 0.611 ± 0.104 | 0.655 | 0.675 ± 0.082 |
| CV AUC mangrove | 0.687 ± 0.045 | 0.700 | 0.689 ± 0.050 |
| CV AUC kelp | 0.749 ± 0.014 | 0.719 | 0.697 ± 0.079 |
| Stage 1 AUC (sg / mg / kp) | 0.834 / 0.880 / 0.914 | 0.804 / 0.816 / 0.907 | 0.793 / 0.816 / 0.891 |
| mangrove distance importance* | 47.4% (30.0–78.2) | 52.3% (12.1–75.7) | 56.4% (11.5–75.4) |
| kelp salinity / distance* | 68.8 / 24.6% | 61.9 / 16.9% | 52.0 / 37.3% |
| seagrass sal / dist / turb* | 15.4 / 14.7 / 13.4% | 11.1 / 13.6 / 7.2% | 13.8 / 14.2 / 7.7% |
| river size* | area 0.3–3.3% | area 1.5–3.5% | discharge 1.6–2.5% |
| median distance mg / kp / sg (km) | 31.8 / 33.0 / 43.7 | 28.5 / 32.0 / 38.6 | 28.5 / 33.2 / 39.7 |
| occupancy at mouth mg / sg / kp | 0.98 / 0.81 / 0.48 | — | ≈ 0.99 / 0.80 / 0.56 |
| 0.50 crossing mg / sg | 21 km / 17 km | — | ≈ 25 km / 16 km |
| 70% regulation sg / mg / kp | −8.3 / −6.0 / −3.1% | −6.5 / −7.7 / −3.0% (bug) | −10.1 / −12.1 / −5.7% |
| kelp climate | 11.9 ± 4.2 °C, \|lat\| 49.7° | (as v5) | 13.3 ± 4.6 °C, \|lat\| 44.9° |

\* median (IQR) across thinned species models (FigSp1). Occupancy values are
read off FigD4.

What the revision itself changed (column 2 → 3): kelp, through *E. cava* (more
distance, less salinity); the regulation effects, through the bug fix; and the
river-size row, through area → discharge. Seagrass and mangrove otherwise
move within noise.

---

## 5. Open decisions

- [x] `DIST_VAR` — keep `log_dist` (§3b)
- [x] River-size claim stands, now on discharge (§3c)
- [x] Kelp group — per-species cap 100, caps 200 and none as sensitivity (§3f)
- [x] Scenario distance pathway — fixed (§3e)
- [ ] **One environment for all manuscript numbers** (§3g — cause now known:
      terra version changes the background draw). Recommend terra 1.9.50 here;
      or rerun the revision under terra 1.8.86 to keep v5's draw. Not a mix.
- [x] **σk definition**: ST (log Kd490 ~ S + SST) main; STx and T as
      sensitivity, regulation quoted as the range across the three (§3n).
- [x] **Stage-2 variable set**: struct (S, σk, D) adopted 2026-10-30 and made
      the default; current (S, Kd490, D, Q) kept as sensitivity (`BC_FEAT=current`).
      Main outputs in `output/` (seed 42); the previous current-set run is in
      `output_current_cap100/`. Figure comparison with v5 in the report.
- [x] **Report headline numbers over several background draws** — done for
      ten draws (§3m); quote medians and ranges.
- [ ] **Soften the one-stage claim** (§3h): "about twice", with both
      permutation and drop-variable intervals.
- [ ] Response letter: disclose the scenario bug and the corrected regulation
      effects (§3e)
- [x] Taxonomy table (Supplementary Table S1) — regenerated,
      `data/TableS1_taxa_revised.csv` (§3m)
- [ ] **RiverATLAS not implemented.** Unchanged judgement: discharge covers
      the substance. The one hint that composition might matter is the mangrove
      `both` row (§3c); revisit only if a reviewer presses it.
- [ ] *Sargassum* not added (decided 2026-10-28): too speciose and too tropical.
      *S. muticum* as a single widespread introduced species remains an option.

---

## 6. Manuscript sections to update once results are in

| section | what changes |
|---|---|
| Methods — predictors | catchment area → discharge, with the rationale |
| Methods — distance | functional form selected rather than assumed; resolves the inconsistency with the scenario model |
| Methods — taxa | *Ecklonia cava* complex, Akita et al. (2020), the GBIF backbone limitation, and the per-species cap |
| Methods — scenarios | the distance pathway now actually enters the prediction |
| Results — importance | every number, from the FEAT change and the cap |
| Results — scenarios | regulation effects ~1.5–2× larger; FigSc3 no longer flat |
| Results — new paragraph | turbidity decomposition |
| Discussion — turbidity | replace "hard to interpret" hedging with the measured split |
| Discussion — design | one-stage comparison as evidence, with intervals (§3h); "about twice", not "collapses" |
| Kelp characterisation | "cool-temperate to subpolar" → temperate to subpolar (13.3 ± 4.6 °C) |
| Supplementary | Table S1 regenerated; TabS_* and FigR1/FigR2 added; cap sensitivity table |
| References | Akita et al. 2020 *Phycologia* 59:330–339; Freckleton 2002 *J. Anim. Ecol.* 71:542–545 |

Citation numbering is first-appearance order (Nature style) — renumber after
the text settles, not before.
