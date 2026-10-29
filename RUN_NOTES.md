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

### 3g. Stage 1 AUC drift (new)

Run 1 gave seagrass and mangrove Stage 1 AUC 0.804 / 0.816 against v5's
0.834 / 0.880, although the revision does not touch those groups' Stage 1
inputs (SST, lat). The v5 code itself (`f5593e1`) run on this server gives the
same 0.804 / 0.816, and its other headline numbers also miss the v5 anchors
(CV AUC 0.655 / 0.700 / 0.719 vs 0.611 / 0.687 / 0.749). So the drift is the
environment, not the revision. Leading suspect: the 20,000 background points
from `terra::spatSample` under `set.seed(42)` differ between terra 1.8.86
(RStudio Server) and 1.9.50 (here); the depth filter keeps 1,183 of them here.

Consequence: the v5 manuscript numbers cannot serve as the "before" for the
revision. The fair comparison is v5 code vs revised code in one environment
(column 2 vs 3 in §4), and every number in the revised manuscript should come
from one environment.

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
- [ ] **One environment for all manuscript numbers** (§3g). Recommend: take
      every number from this server and record its session info; or rerun
      v5 + revision on RStudio Server and use that. Not a mix.
- [ ] Response letter: disclose the scenario bug and the corrected regulation
      effects (§3e)
- [ ] Taxonomy table (Supplementary Table S1) — regenerate with the merged
      *E. cava*; kelp genus counts and the species total both change
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
| Discussion — design | one-stage comparison as evidence (raw ΔAUC), not assertion |
| Kelp characterisation | "cool-temperate to subpolar" → temperate to subpolar (13.3 ± 4.6 °C) |
| Supplementary | Table S1 regenerated; TabS_* and FigR1/FigR2 added; cap sensitivity table |
| References | Akita et al. 2020 *Phycologia* 59:330–339; Freckleton 2002 *J. Anim. Ecol.* 71:542–545 |

Citation numbering is first-appearance order (Nature style) — renumber after
the text settles, not before.
