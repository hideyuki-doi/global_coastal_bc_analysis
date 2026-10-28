# Blue carbon ecosystems and river proximity

Code and derived data for a global analysis of whether riverine connectivity
constrains where seagrasses, mangroves and kelp forests occur, once
biogeography is accounted for.

Two-stage MaxEnt species distribution models separate the climatic envelope
(stage 1: sea surface temperature, latitude) from within-range habitat
selection (stage 2: salinity, turbidity, distance to the nearest river mouth,
river discharge), with stage-2 background sampling matched to each taxon's own
climate so that climate is controlled rather than discarded.

---

## Running the pipeline

```bash
Rscript R/01_ecklonia_gbif.R     # once: recover the Ecklonia cava complex
Rscript R/blue_carbon_full.R     # main pipeline; sources R/02_reviewer_analyses.R
```

`blue_carbon_full.R` expects `bc_analysis.rds` (the screened GBIF occurrence
table with Bio-ORACLE predictors attached) in the working directory, and the
external datasets listed below. It writes figures and tables to `output/`.

Runtime is dominated by the over-water distance surface and the species-level
models; budget a few hours on a first run, less afterwards because the distance
surface and the prediction grid are cached in `output/`.

---

## What changed in this revision

The revision responds to four points raised in review.

### 1. River size: catchment area replaced by discharge

Catchment area treats an arid, sparsely vegetated basin and a wet one of the
same size alike, so it is a poor proxy for the freshwater and dissolved load a
river actually delivers — which is what the predictor is meant to represent.
Area was also the weakest predictor in the first version (0.3–3.3% importance),
and some of that is measurement rather than ecology.

HydroRIVERS carries `DIS_AV_CMS`, long-term average discharge at each reach, so
no additional download is needed. Discharge integrates area with runoff and
replaces `log_upstream` in `FEAT`. Area is retained in the data frames and
compared against discharge head to head in `R/02_reviewer_analyses.R` (§8c).

### 2. Turbidity decomposed into ocean and land components

Kd490 at a coastal cell sums local ocean optics (phytoplankton, resuspension —
both strongly latitudinal) and terrestrial delivery (CDOM, mineral sediment).
The first version could not say which half carried the signal, which is why the
predictor was hard to interpret.

`R/02_reviewer_analyses.R` §8a does this two ways:

* **Variance partition.** `log Kd490 ~ s(SST) + s(|lat|) + s(distance) +
  s(discharge)` fitted over presences and background together, with deviance
  explained split into climate-unique, river-unique and shared.
* **Residual model.** Kd490 residual after removing the ocean-climate signal,
  regressed on the two river terms — the decomposition suggested in review.

One clarification on that second model: because the residual is on the **log**
scale, the additive regression is **already the multiplicative model**
(`Kd_land ∝ discharge^b1 × (distance + 1)^b2`). Its coefficients read as
elasticities, so a tenfold change in discharge multiplies land-attributable
turbidity by `10^b1`. No log link is needed on top.

The SDM is then refitted with the land-attributable residual in place of raw
Kd490, as a sensitivity check rather than as the headline model: using a
residual as a predictor is known to bias coefficients relative to fitting both
terms jointly (Freckleton 2002, *J. Anim. Ecol.* 71, 542–545), and the joint
variance partition above is the cleaner statement of the same thing.

### 3. Distance functional form

The first version was internally inconsistent: the SDM used `log10(d + 1)`
while the flow-regulation scenarios perturbed distance through `exp(−d / L)` —
two different claims about one process. A plume diluting radially from a point
source argues for `1/d` instead.

All three forms are now declared together (`add_dist_forms()`), carried by
presences, background and the prediction grid alike, and compared under the
same spatial-block cross-validation (§8b). `DIST_VAR` selects the active form
in one place, and the response-curve figures follow it rather than hard-coding
log10.

### 4. One stage or two?

Answered by experiment rather than argument (§8d). Four designs are fitted and
their river-distance importance compared:

| design | background | turbidity |
|---|---|---|
| 1-stage | ocean-wide | raw Kd490 |
| 1-stage | ocean-wide | land-attributable residual |
| 2-stage (stage 2) | climate-matched | raw Kd490 |
| 2-stage (stage 2) | climate-matched | land-attributable residual |

The second row tests the specific hypothesis raised in review — that removing
the temperature signal from turbidity might make a single stage sufficient.

### 5. *Ecklonia cava* complex recovered

Akita et al. (2020, *Phycologia* 59, 330–339) showed *E. cava*, *E. kurome* and
*E. stolonifera* to be one species, the latter two becoming subspecies of
*E. cava*; WoRMS and AlgaeBase have adopted this.

The GBIF backbone carries only **two** species in the genus — *E. radiata* and
*E. maxima*. *E. cava* and its synonyms match nothing, so every record bearing
those names is indexed at genus rank with an empty `species` field and is
invisible to any species-level query:

```r
name_backbone("Ecklonia cava")                 # matchType = NONE
occ_search(scientificName = "Ecklonia cava")   # 0 records
```

That is why the genus contributed a single species (*E. radiata*, 51 thinned
records) plus 13 unusable genus-only records in the first version, leaving the
kelp group effectively a Southern-Hemisphere and North-Atlantic sample.

`R/01_ecklonia_gbif.R` recovers them through `verbatimScientificName`, which
matches the name the publisher supplied:

| queried name | georeferenced |
|---|---|
| *Ecklonia cava* | 1,335 |
| *Ecklonia kurome* | 735 |
| *Ecklonia stolonifera* | 176 |
| *E. cava* subsp. *kurome* | 34 |
| *E. cava* subsp. *stolonifera* | 81 |
| **merged, unique coordinates** | **1,447** (JP 1,247 · KR 197 · CN 3) |

*Eckloniopsis radicosa* was moved into *Ecklonia* by the same paper but kept as
its own species, so it is retrieved separately (190 coordinates) and not
merged. *E. radiata* and *E. maxima* resolve through the backbone already and
are deliberately **not** re-downloaded here.

---

## Layout

```
R/01_ecklonia_gbif.R        recover the Ecklonia cava complex from GBIF
R/blue_carbon_full.R        main pipeline (PART 1-7, then sources PART 8)
R/02_reviewer_analyses.R    PART 8: the four analyses described above
R/blue_carbon_figures.R     standalone figure rebuilds
R/peninsula_check.R         over-water vs great-circle distance diagnostic
scripts/get_external_data.sh  download the third-party datasets
data/                       derived data small enough to version
output/                     figures and tables (written by the pipeline)
```

---

## External data (not in this repository)

These are large third-party datasets, redistributed by their own providers
under their own licences. `scripts/get_external_data.sh` fetches them.

| dataset | used for | size | licence |
|---|---|---|---|
| [HydroRIVERS v1.0](https://www.hydrosheds.org/products/hydrorivers) | river mouths, catchment area, discharge | 2.1 GB | CC BY 4.0 |
| [WorldClim 2.1 elevation 30″](https://worldclim.org/data/worldclim21.html) | land mask for the over-water distance surface | 326 MB | CC BY-SA 4.0 |
| [Bio-ORACLE](https://bio-oracle.org/) via `sdmpredictors` | SST, salinity, Kd490 | ~1 GB | CC BY 4.0 |
| [MARSPEC MS_bathy_5m](http://www.marspec.org/) via `sdmpredictors` | bathymetry | — | CC BY-NC-SA 3.0 |
| [EOT20 ocean tides](https://doi.org/10.17882/79489) | tidal range (optional) | ~2 GB | CC BY 4.0 |

`bc_analysis.rds` — the screened occurrence table — is produced upstream and is
not in this repository.

---

## Reproducibility

Random seeds are fixed for background sampling, model fitting and spatial
thinning, so the pipeline is deterministic. Package versions used for the
submitted analysis are recorded in `doc/session_info.txt`.

## Licence

Code: MIT (`LICENSE`). Derived data in `data/`: CC BY 4.0. Third-party datasets
keep their own licences, listed above.
