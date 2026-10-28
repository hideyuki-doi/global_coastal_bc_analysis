#!/usr/bin/env bash
# =============================================================================
# Fetch the third-party datasets the pipeline needs.
#
# These are large and belong to their own providers, so they are not kept in
# this repository. Run this once from the repository root; everything lands
# where blue_carbon_full.R looks for it.
#
#   bash scripts/get_external_data.sh
#
# Bio-ORACLE and MARSPEC are NOT fetched here — sdmpredictors downloads them on
# first use into sdm_data/, which the pipeline creates itself.
# =============================================================================
set -euo pipefail

say() { printf '\n\033[1m==> %s\033[0m\n' "$*"; }

need() {
  command -v "$1" >/dev/null 2>&1 || { echo "missing required tool: $1" >&2; exit 1; }
}
need curl
need unzip

# ── HydroRIVERS v1.0 (global, shapefile) ─────────────────────────────────────
# River mouths, catchment area (UPLAND_SKM) and discharge (DIS_AV_CMS).
# CC BY 4.0 — Lehner, B. & Grill, G. (2013) Hydrol. Process. 27, 2171-2186.
if [ -d HydroRIVERS_v10_shp ]; then
  say "HydroRIVERS already present — skipping"
else
  say "HydroRIVERS v1.0 (~2.1 GB extracted)"
  curl -L --fail --retry 5 --retry-delay 5 -C - \
    -o HydroRIVERS_v10_shp.zip \
    "https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_shp.zip"
  unzip -q HydroRIVERS_v10_shp.zip -d HydroRIVERS_v10_shp
  rm -f HydroRIVERS_v10_shp.zip
fi

# ── WorldClim 2.1 elevation, 30 arc-seconds ──────────────────────────────────
# Land mask for the over-water (gridDist) distance surface.
# CC BY-SA 4.0 — Fick, S.E. & Hijmans, R.J. (2017) Int. J. Climatol. 37, 4302-4315.
if [ -f wc2.1_30s_elev.tif ]; then
  say "WorldClim elevation already present — skipping"
else
  say "WorldClim 2.1 elevation 30s (~326 MB)"
  curl -L --fail --retry 5 --retry-delay 5 -C - \
    -o wc2.1_30s_elev.zip \
    "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_elev.zip"
  unzip -q -j wc2.1_30s_elev.zip "*.tif" -d .
  rm -f wc2.1_30s_elev.zip
fi

say "Done."
cat <<'EOF'

Still needed, and not fetched by this script:

  bc_analysis.rds   the screened GBIF occurrence table with Bio-ORACLE
                    predictors attached, produced upstream of this repository.
                    Place it in the repository root.

  EOT20 ocean tides (optional, ~2 GB) — only used for the tidal-range
  predictor, which the models drop when it is absent. Register and download
  from https://doi.org/10.17882/79489, then place ocean_tides.zip in the root.

Bio-ORACLE and MARSPEC are downloaded automatically by sdmpredictors into
sdm_data/ the first time the pipeline runs.
EOF
