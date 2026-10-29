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

| date | what | status | notes |
|---|---|---|---|
| | first full run | | |

---

## 3. Results

### 3a. Turbidity decomposition

**Question.** Kd490 at a coastal cell sums ocean optics (phytoplankton,
resuspension — both latitudinal) and terrestrial delivery (CDOM, mineral
sediment). How much is each?

**Reads.** `output/TabS_Kd490_variance_partition.csv`,
`output/TabS_Kd490_residual_model.csv`, `output/FigR1_Turbidity_Decomposition`

| component | deviance explained |
|---|---|
| Total (climate + river) | |
| Unique to climate | |
| Unique to river | |
| Shared | |

Residual model (`kd_resid ~ log_discharge + log_dist`), R² = ______

| term | estimate | p | ×per 10-fold |
|---|---|---|---|
| log_discharge | | | |
| log_dist | | | |

**What it decides.** If the river-unique share is non-trivial, the turbidity
predictor can be described in the Discussion as carrying a real terrestrial
signal rather than being an uninterpretable mixture — which is the criticism
this answers. If it is near zero, turbidity is essentially an ocean variable
and the paper should say so plainly and lean on distance and salinity instead.

The coefficients read as elasticities because the residual is on the log scale:
the additive model is already multiplicative, `Kd_land ∝ discharge^b1 ×
(dist+1)^b2`.

---

### 3b. Distance functional form

**Question.** v1 used `log10(d+1)` in the SDM but `exp(−d/L)` in the scenarios —
two different claims about one process. Radial dilution from a point source
argues for `1/d`.

**Reads.** `output/TabS_distance_form_comparison.csv`, `FigR2` panel a

| form | seagrass | mangrove | kelp | mean |
|---|---|---|---|---|
| `log_dist` log₁₀(d+1) | | | | |
| `inv_dist` 1/(d+1) | | | | |
| `exp_dist` exp(−d/50) | | | | |

Best: ______

**What it decides.** If the winner is not `log_dist`, set `DIST_VAR` in
`R/blue_carbon_full.R` and rerun — the response-curve figures follow it
automatically. Either way the paper gains a sentence saying the form was
selected rather than assumed, and the inconsistency with PART 7 is resolved.

If the three are within noise of each other, say so: that is itself the answer
(the result does not hinge on the transform) and is worth one line in Methods.

---

### 3c. River size: catchment area vs discharge

**Question.** Area counts an arid basin the same as a wet one. Discharge
integrates area with runoff.

**Reads.** `output/TabS_riversize_comparison.csv`

r(log area, log discharge) = ______

| predictor set | seagrass | mangrove | kelp |
|---|---|---|---|
| area only | | | |
| discharge only | | | |
| both | | | |
| neither | | | |

**What it decides.** v1 reported watershed area at 0.3–3.3% importance and
concluded "proximity matters, river size does not". If discharge scores
materially higher, that conclusion was partly an artefact of a bad proxy and
the claim must be softened. If discharge is also near zero, the original claim
stands and is now much better supported — a null on the better measure is worth
more than a null on the worse one.

Watch the `both` row for collinearity: if r is high, the two will split
importance and neither will look strong.

---

### 3d. One stage or two

**Question.** Does the two-stage design earn its complexity? Would decomposing
turbidity let one stage do the job?

**Reads.** `output/TabS_onestage_vs_twostage.csv`, `FigR2` panel b

River-distance importance share (%):

| design | seagrass | mangrove | kelp |
|---|---|---|---|
| 1-stage, ocean-wide bg | | | |
| 1-stage, Kd residual | | | |
| 2-stage (stage 2) | | | |
| 2-stage, Kd residual | | | |

**What it decides.** The expected pattern is that distance importance collapses
in the one-stage rows because climate absorbs it — that is the paper's central
methodological claim, and this turns it from an assertion into a measurement.

The second row is the specific hypothesis raised in review. If it recovers the
distance signal without the second stage, the two-stage design is harder to
justify and the framing needs rethinking. Report it either way.

---

## 4. Expected departures from v5

The kelp group and every importance number will move. Anchors from v5
(over-water distance, pre-revision) for comparison:

| quantity | v5 | new |
|---|---|---|
| records / species / genera | 2,040 / 93 / 25 | |
| seagrass / mangrove / kelp records | 777 / 808 / 455 | |
| CV AUC seagrass | 0.611 ± 0.104 | |
| CV AUC mangrove | 0.687 ± 0.045 | |
| CV AUC kelp | 0.749 ± 0.014 | |
| Stage 1 AUC (sg/mg/kp) | 0.834 / 0.880 / 0.914 | |
| mangrove distance importance | 47.4% (30.0–78.2) | |
| kelp salinity / distance | 68.8% / 24.6% | |
| seagrass sal / dist / turb | 15.4 / 14.7 / 13.4% | |
| watershed area (all) | 0.3–3.3% | |
| median distance mg / kp / sg (km) | 31.8 / 33.0 / 43.7 | |
| occupancy at mouth mg / sg / kp | 0.98 / 0.81 / 0.48 | |
| 0.50 crossing mg / sg | 21 km / 17 km | |
| 70% regulation sg / mg / kp | −8.3 / −6.0 / −3.1% | |

**Why kelp moves most.** 1,447 East Asian coordinates enter a group that had
455 records in total. That widens the kelp thermal and latitudinal range, which
shifts both the Stage 1 envelope and the climate-matched background window for
Stage 2. Expect the kelp numbers to depart from v5 substantially; seagrass and
mangrove should move only through the area→discharge swap.

**Watch for.** The kelp group was "cool-temperate to subpolar (11.9 ± 4.2 °C,
median |lat| 49.7°)" in v5. Japanese and Korean *E. cava* sits warmer and lower
than that. If the kelp thermal envelope widens a lot, the Stage 1 range model
weakens and the "climate sets the range" half of the story needs checking
rather than assuming it still holds.

---

## 5. Open decisions

- [ ] `DIST_VAR` — keep `log_dist` or switch, per §3b
- [ ] Whether the area→discharge result changes the "river size does not
      matter" claim, per §3c
- [ ] Taxonomy table (Supplementary Table S1) — regenerate with the merged
      *E. cava*; kelp genus counts and the species total both change
- [ ] **RiverATLAS not implemented.** Review suggested weighting catchment area
      by vegetation / cropland / urban fraction. That needs RiverATLAS (3–4 GB
      additional download). Current judgement: discharge already addresses the
      substance of the point, because an arid barren catchment has low discharge
      despite large area. Land-use *composition* is a question about nutrient
      and sediment quality, which is a different claim from "proximity matters".
      Revisit only if a reviewer presses it.
- [ ] *Sargassum* not added (decided 2026-10-28): too speciose and too tropical,
      would break the kelp group's cool-temperate characterisation. *S. muticum*
      as a single widespread introduced species remains an option.

---

## 6. Manuscript sections to update once results are in

| section | what changes |
|---|---|
| Methods — predictors | catchment area → discharge, with the rationale |
| Methods — distance | functional form selected rather than assumed; resolves the inconsistency with the scenario model |
| Methods — taxa | *Ecklonia cava* complex, Akita et al. (2020), and the GBIF backbone limitation |
| Results — importance | every number, from the FEAT change |
| Results — new paragraph | turbidity decomposition |
| Discussion — turbidity | replace "hard to interpret" hedging with the measured split |
| Discussion — design | one-stage comparison as evidence, not assertion |
| Supplementary | Table S1 regenerated; TabS_* and FigR1/FigR2 added |
| References | Akita et al. 2020 *Phycologia* 59:330–339; Freckleton 2002 *J. Anim. Ecol.* 71:542–545 |

Citation numbering is first-appearance order (Nature style) — renumber after
the text settles, not before.
