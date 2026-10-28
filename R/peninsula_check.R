# =============================================================================
# peninsula_check.R — how wrong is dist_km where land blocks the straight line?
#
# THE QUESTION
#   dist_km is a great-circle distance to the nearest river mouth. Where a spit,
#   peninsula or isthmus separates two water bodies, the geometrically nearest
#   mouth may be one the cell cannot receive water from. Reviewer 2 asked for the
#   size of that residual error.
#
# WHAT DOES NOT WORK (all three were tried and discarded)
#   (a) "does the straight line cross land?" -> ~90% of records, meaningless.
#       HydroRIVERS places an outlet at the downstream end of a mapped channel,
#       which sits inland of any independent coastline product, so almost every
#       path clips land within the first kilometre.
#   (b) "how long is the over-land run?" -> median ~11 km against Natural Earth
#       medium-scale coastlines. That is the coastline mismatch, not a peninsula.
#   (c) Snapping the endpoints to water does not rescue (a) or (b): at 1 km
#       resolution 49% of these records, and 63% of the mangrove records, fall on
#       a land cell, because the taxa live exactly at the land-water boundary.
#   (d) A "strict" water mask (a cell is water only if every 1 km subcell is
#       water) is unusable here: it deletes the coastal strip, so only 60 of the
#       6,202 river mouths survive and every distance becomes meaningless.
#
# WHAT THIS SCRIPT DOES
#   It never classifies an endpoint. Using one terra::gridDist() call it computes,
#   for every water cell on Earth, the shortest path to ANY river mouth that
#   stays in water, then compares that with the great-circle distance the
#   analysis actually used:
#         inflation = over-water distance / great-circle distance
#   Inflation is 1 where the straight line is a fair approximation and large
#   where land blocks it. No endpoint decision, no coastline-mismatch artefact.
#
#   A cell counts as water if ANY of its 1 km subcells is water, which keeps the
#   coastal strip continuous. The consequence is that a land barrier thinner than
#   one cell stays passable, so the reported shares are LOWER BOUNDS. The
#   calculation is therefore repeated at three grid resolutions: refining the
#   grid closes more narrow barriers and raises the estimate, and the trend shows
#   how far from converged the answer is.
#
# RESULT ON THIS DATASET (2,040 records, 6,202 mouths)
#   cell size   median inflation   >25% inflated   >50%    >2x
#     2.8 km          1.11             31.1%       21.1%   12.8%
#     5.6 km          1.09             24.3%       14.9%    8.1%
#    11.1 km          1.07             21.4%       10.4%    4.7%
#
#   Reading: for the typical record the straight-line distance is within ~10% of
#   the over-water distance. For roughly a fifth of records it understates the
#   over-water distance by more than half, and for about one in eight by more
#   than a factor of two.
#
#   Direction of the bias: a straight line is never longer than a water path, so
#   affected records are always recorded as CLOSER to a river than they
#   hydrologically are. Presences and background are measured identically, so the
#   error adds noise to the predictor and attenuates its estimated importance
#   rather than manufacturing it.
#
# INPUTS   bc_analysis.rds, HydroRIVERS_v10.shp, and a 1 km land mask whose
#          cells are non-NA on land and NA on water (WorldClim elevation works).
# RUNTIME  seconds per resolution; the aggregation dominates.
# =============================================================================

suppressMessages({library(sf); library(terra); library(dplyr)})

BC_RDS    <- "bc_analysis.rds"
HYDRO     <- list.files(".", "HydroRIVERS_v10\\.shp$", recursive = TRUE,
                        full.names = TRUE)[1]
LAND_MASK <- "wc2.1_30s_elev.tif"        # land = elevation, water = NA
OUT_DIR   <- "output"; dir.create(OUT_DIR, showWarnings = FALSE)

AGG_FACTS <- c(3, 6, 12)                 # 30 arc-sec x fact -> working grids
LAT_LIM   <- c(-70, 80)                  # trim to the latitudes the data occupy
SEARCH_DEG <- 0.4                        # fallback search radius for land records

lastcol <- function(e) if (is.null(dim(e))) e else e[, ncol(e)]

# ── inputs ──────────────────────────────────────────────────────────────────
bc  <- readRDS(BC_RDS)
pts <- as.matrix(bc[, c("lon", "lat")])
message(sprintf("occurrences: %d", nrow(pts)))

riv <- sf::st_read(HYDRO, quiet = TRUE, query = paste(
  "SELECT HYRIV_ID FROM HydroRIVERS_v10",
  "WHERE ORD_FLOW <= 5 AND NEXT_DOWN = 0 AND ENDORHEIC = 0"))
mouths <- do.call(rbind, lapply(sf::st_geometry(riv), \(g) {
  m <- sf::st_coordinates(g); m[nrow(m), 1:2]
}))
colnames(mouths) <- c("lon", "lat")
message(sprintf("river mouths: %d", nrow(mouths)))

# the distance the analysis actually used
psf  <- sf::st_as_sf(as.data.frame(pts),    coords = c("lon", "lat"), crs = 4326)
msf  <- sf::st_as_sf(as.data.frame(mouths), coords = c("lon", "lat"), crs = 4326)
nn   <- sf::st_nearest_feature(psf, msf)
d_gc <- as.numeric(sf::st_distance(psf, msf[nn, ], by_element = TRUE)) / 1000

message("Reading the 1 km land mask ...")
elev0 <- terra::crop(terra::rast(LAND_MASK),
                     terra::ext(-180, 180, LAT_LIM[1], LAT_LIM[2]))

overwater_distance <- function(fct) {
  # mean with na.rm = FALSE is NA if ANY subcell is water -> water mask
  agg <- terra::aggregate(elev0, fact = fct, fun = "mean", na.rm = FALSE)
  cellkm <- terra::res(agg)[1] * 111
  water <- terra::ifel(is.na(agg), 1, NA)
  mc <- terra::cellFromXY(water, mouths)
  ok <- !is.na(mc) & !is.na(lastcol(terra::extract(water, mouths)))
  water[mc[ok]] <- 2
  dw <- terra::gridDist(water, target = 2)

  v <- lastcol(terra::extract(dw, pts)) / 1000
  miss <- which(is.na(v))
  step <- terra::res(water)[1]
  for (i in miss) {
    for (rr in seq(step, SEARCH_DEG, by = step)) {
      a <- seq(0, 2 * pi, length.out = 17)[-17]
      cand <- cbind(pts[i, 1] + rr * cos(a) / cos(pts[i, 2] * pi / 180),
                    pts[i, 2] + rr * sin(a))
      vv <- lastcol(terra::extract(dw, cand)) / 1000
      if (any(!is.na(vv))) { v[i] <- min(vv, na.rm = TRUE); break }
    }
  }
  list(v = v, cellkm = cellkm, n_mouth = sum(ok), n_fallback = length(miss))
}

res <- data.frame(Type = factor(bc$Type, levels = c("Seagrass", "Mangrove", "Kelp")),
                  d_greatcircle_km = d_gc)
summary_rows <- list()

for (fct in AGG_FACTS) {
  o <- overwater_distance(fct)
  r <- o$v / pmax(d_gc, 0.5)
  tag <- sprintf("%.1fkm", o$cellkm)
  res[[paste0("d_overwater_", tag)]] <- o$v
  res[[paste0("inflation_",   tag)]] <- r
  message(sprintf("  %-7s cells | mouths on water %d/%d | land-cell fallbacks %d",
                  tag, o$n_mouth, nrow(mouths), o$n_fallback))
  summary_rows[[tag]] <- data.frame(
    cell_km      = round(o$cellkm, 1),
    median_infl  = round(median(r, na.rm = TRUE), 2),
    pct_gt_1.25  = round(100 * mean(r > 1.25, na.rm = TRUE), 1),
    pct_gt_1.5   = round(100 * mean(r > 1.5,  na.rm = TRUE), 1),
    pct_gt_2     = round(100 * mean(r > 2,    na.rm = TRUE), 1))
}

cat("\n=== over-water vs great-circle distance to the nearest river mouth ===\n")
cat("(finer grids close more narrow barriers, so the shares are lower bounds)\n\n")
print(do.call(rbind, summary_rows), row.names = FALSE)

fine <- paste0("inflation_", names(summary_rows)[1])
cat(sprintf("\n=== per ecosystem, finest grid (%s) ===\n", names(summary_rows)[1]))
print(res |> dplyr::mutate(.i = .data[[fine]]) |> dplyr::group_by(Type) |>
  dplyr::summarise(n = dplyr::n(),
                   median_infl = round(median(.i, na.rm = TRUE), 2),
                   pct_gt_1.5  = round(100 * mean(.i > 1.5, na.rm = TRUE), 1),
                   pct_gt_2    = round(100 * mean(.i > 2,   na.rm = TRUE), 1),
                   .groups = "drop") |> as.data.frame())

write.csv(res, file.path(OUT_DIR, "peninsula_check.csv"), row.names = FALSE)
saveRDS(res, file.path(OUT_DIR, "peninsula_check.rds"))
message("\nSaved output/peninsula_check.{csv,rds}")

# =============================================================================
# OPTIONAL follow-up, if a reviewer wants the analysis re-run without the
# affected records. Drop the worst cases and refit the species-level models:
#
#   keep <- res[[paste0("inflation_", names(summary_rows)[1])]] <= 2
#   bc   <- bc[keep, ]          # then rerun blue_carbon_full.R from PART 5c
#
# This is a genuine robustness check rather than a correction: the right fix is
# to define dist_km itself on the over-water surface computed above.
# =============================================================================
