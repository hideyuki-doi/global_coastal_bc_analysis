# =============================================================================
# Blue Carbon × River Proximity — Full Analysis Pipeline
# Models: GAM (spatial smooth) + MaxEnt + XGBoost ensemble + Flow scenarios
# Usage: place bc_analysis.rds in getwd(), then source this file.
# =============================================================================

# ── 0. Packages ───────────────────────────────────────────────────────────────
REQ_PKGS <- c(
  "dplyr","tidyr","stringr","tibble","purrr",
  "ggplot2","patchwork","ggridges","scales","ggh4x","sf","terra",
  "sdmpredictors","mgcv","gratia","spdep",
  "maxnet","xgboost","pdp","pROC","ranger",
  "rnaturalearth","rnaturalearthdata"
)
new_pkgs <- REQ_PKGS[!sapply(REQ_PKGS, requireNamespace, quietly = TRUE)]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)

library(dplyr);  library(tidyr);  library(stringr); library(tibble); library(purrr)
library(ggplot2); library(patchwork); library(ggridges); library(scales); library(ggh4x)
library(sf);     library(terra);  library(sdmpredictors)
library(mgcv);   library(gratia); library(spdep)
library(maxnet); library(xgboost); library(pdp); library(pROC)
library(ranger)
library(rnaturalearth); library(rnaturalearthdata)

# fastshap was archived on CRAN (last release 0.1.1), so it is treated as
# OPTIONAL: it is used only for the TreeSHAP attribution behind the discarded
# tree-based robustness check (FigS2b), which no figure in the manuscript uses.
# Without it, PART 5b falls back to RF permutation importance automatically.
HAS_SHAP <- requireNamespace("fastshap", quietly = TRUE) &&
            requireNamespace("shapviz",  quietly = TRUE)
if (HAS_SHAP) { library(fastshap); library(shapviz) } else {
  message("NOTE: fastshap/shapviz not available — SHAP attribution skipped, ",
          "RF permutation importance used instead (no manuscript figure affected).")
}

# ── 1. Paths — edit only if your files are elsewhere ─────────────────────────
BC_RDS  <- "bc_analysis.rds"    # upload to working directory (getwd())
OUT_DIR <- "output"             # figures saved here
dir.create(OUT_DIR,    showWarnings = FALSE)
dir.create("sdm_data", showWarnings = FALSE)

if (!file.exists(BC_RDS))
  stop("bc_analysis.rds not found. Current dir: ", getwd(),
       "\nUpload the file via RStudio Files pane → Upload.")

# ── 2. Shared theme & save helper ─────────────────────────────────────────────
ECO_PAL  <- c(Seagrass = "#009E73", Mangrove = "#E69F00", Kelp = "#0072B2")
ECO_LVLS <- c("Seagrass","Mangrove","Kelp")
eco_strip <- function() {
  ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(fill = alpha(unname(ECO_PAL[ECO_LVLS]), 0.22)),
    text_x       = ggh4x::elem_list_text(colour = unname(ECO_PAL[ECO_LVLS]), face = "bold")
  )
}

theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey94", linewidth = 0.25),
    strip.background  = element_rect(fill = "grey90", colour = NA),
    strip.text        = element_text(colour = "grey15", face = "bold", size = rel(0.9)),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.key        = element_rect(fill = NA, colour = NA),
    legend.key.size   = unit(0.4, "cm"),
    plot.title        = element_text(face = "bold", size = rel(1.05), hjust = 0,
                                     margin = margin(b = 4)),
    plot.subtitle     = element_text(colour = "grey40", size = rel(0.85), hjust = 0,
                                     margin = margin(b = 8)),
    plot.caption      = element_text(colour = "grey55", size = rel(0.8), hjust = 1),
    axis.text         = element_text(colour = "grey25", size = rel(0.85)),
    axis.ticks        = element_line(colour = "grey70", linewidth = 0.3),
    axis.title        = element_text(size = rel(0.9))
  )
}

# ── Distance to the river mouth: functional form ─────────────────────────────
# Three candidate transforms, all decreasing in river influence:
#   log_dist  log10(d + 1)       used throughout v1; steepest near the mouth
#   inv_dist  1 / (d + 1)        concentration under 2-D radial dilution of a
#                                point source, which is what a plume is
#   exp_dist  exp(-d / L_PLUME)  plume decay — the form the flow-regulation
#                                scenarios in PART 7 already assume
# v1 was internally inconsistent here: the SDM used log_dist while the scenarios
# perturbed distance through exp(-d / L). Declaring the set makes the choice
# explicit, lets the data decide (R/02_reviewer_analyses.R compares all three
# under spatial-block CV), and lets the SDM and the scenarios share one form.
L_PLUME   <- 50          # km — plume decay length, as used by PART 7
DIST_VAR  <- "log_dist"  # active transform; see R/02_reviewer_analyses.R

add_dist_forms <- function(d) {
  d$log_dist <- log10(d$dist_km + 1)
  d$inv_dist <- 1 / (d$dist_km + 1)
  d$exp_dist <- exp(-d$dist_km / L_PLUME)
  d
}

# Back-transform for whichever form is active, so the response-curve figures
# follow DIST_VAR instead of hard-coding log10.
dist_to_km <- function(v, form = DIST_VAR) switch(form,
  log_dist = 10^v - 1,
  inv_dist = 1 / v - 1,
  exp_dist = -L_PLUME * log(v))
km_to_dist <- function(km, form = DIST_VAR) switch(form,
  log_dist = log10(km + 1),
  inv_dist = 1 / (km + 1),
  exp_dist = exp(-km / L_PLUME))

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p,
         width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p,
         width = w, height = h, units = "mm", dpi = 300, bg = "white")
  message("Saved: ", name)
}

# =============================================================================
# PART 1 — LOAD DATA
# =============================================================================

message("Loading bc_analysis.rds ...")
bc <- readRDS(BC_RDS)

if (!"Type_f"   %in% names(bc)) bc$Type_f   <- factor(bc$Type, levels = ECO_LVLS)
if (!"abs_lat"  %in% names(bc)) bc$abs_lat  <- abs(bc$lat)
bc$log_dist  <- log10(bc$dist_km + 1)  # log10 scale: emphasises near-river gradient
if (!"log_Kd"   %in% names(bc)) bc$log_Kd   <- log(bc$Kd490)

message(sprintf("Loaded: %s records", format(nrow(bc), big.mark = ",")))
print(table(bc$Type))

# ── Rivers (for SDM background point distances) ───────────────────────────────
message("Loading rivers...")
hydro_path <- list.files(".", "HydroRIVERS_v10\\.shp", recursive = TRUE,
                          full.names = TRUE)[1]
if (!is.na(hydro_path)) {
  rivers <- sf::st_read(hydro_path, quiet = TRUE) |>
    dplyr::filter(ORD_FLOW <= 5) |>
    (\(x) if (is.na(sf::st_crs(x))) sf::st_set_crs(x, 4326)
           else sf::st_transform(x, 4326))()
  message("  HydroRIVERS: ", nrow(rivers), " segments")
} else {
  message("  HydroRIVERS not found — using Natural Earth (lower resolution)")
  rivers <- rnaturalearth::ne_download(
    scale = 10, type = "rivers_lake_centerlines",
    category = "physical", returnclass = "sf"
  ) |> sf::st_transform(4326)
}

# ── River mouths (outlet termini) ─────────────────────────────────────────────
# All distance predictors are measured to the nearest river MOUTH, not to the
# nearest reach anywhere in the network.  `rivers` spans headwaters to outlet, so
# st_nearest_feature() against it can match a marine cell to an inland reach
# whose own mouth is far away: on the occurrence records, 33% of matches were
# non-outlet reaches and 1.4-1.8% sat >200 km inland along the network (max
# 1,143 km).  Plume dilution, the salinity ecocline and sediment delivery all
# propagate from the mouth, so the mouth is the correct anchor.
#   outlet = NEXT_DOWN == 0 (no downstream reach); ENDORHEIC == 1 excluded
#   because those basins terminate in inland sinks, not the sea.
# The mouth location is the last vertex of the outlet reach (its downstream end).
message("Building river-mouth point layer from network outlets...")
if ("NEXT_DOWN" %in% names(rivers)) {
  .outlets <- rivers |>
    dplyr::filter(NEXT_DOWN == 0,
                  if ("ENDORHEIC" %in% names(rivers)) ENDORHEIC == 0 else TRUE)
  .mouth_xy <- do.call(rbind, lapply(sf::st_geometry(.outlets), \(g) {
    m <- sf::st_coordinates(g); m[nrow(m), 1:2]
  }))
  river_mouths <- sf::st_as_sf(
    as.data.frame(.mouth_xy) |> setNames(c("lon", "lat")),
    coords = c("lon", "lat"), crs = 4326)
  river_mouths$UPLAND_SKM <- if ("UPLAND_SKM" %in% names(.outlets))
    .outlets$UPLAND_SKM else 0
  river_mouths$ORD_FLOW   <- if ("ORD_FLOW" %in% names(.outlets))
    .outlets$ORD_FLOW else NA_integer_
  # DIS_AV_CMS — HydroRIVERS long-term average discharge at the outlet (m³/s).
  # Discharge measures freshwater and dissolved-load delivery better than
  # catchment area does, because it already integrates area with the
  # catchment's runoff: a large arid or sparsely vegetated basin is correctly
  # down-weighted against a small wet one, which area alone cannot do.
  # UPLAND_SKM is kept alongside so the two can be compared directly.
  river_mouths$DIS_AV_CMS <- if ("DIS_AV_CMS" %in% names(.outlets))
    .outlets$DIS_AV_CMS else NA_real_
  rm(.outlets, .mouth_xy)
} else {
  # Natural Earth fallback has no topology: use reach endpoints as mouths
  message("  no NEXT_DOWN field — falling back to all reach endpoints")
  .mouth_xy <- do.call(rbind, lapply(sf::st_geometry(rivers), \(g) {
    m <- sf::st_coordinates(g); m[nrow(m), 1:2]
  }))
  river_mouths <- sf::st_as_sf(
    as.data.frame(.mouth_xy) |> setNames(c("lon", "lat")),
    coords = c("lon", "lat"), crs = 4326)
  river_mouths$UPLAND_SKM <- 0
  river_mouths$ORD_FLOW   <- NA_integer_
  river_mouths$DIS_AV_CMS <- NA_real_
  rm(.mouth_xy)
}
message(sprintf("  river mouths: %d", nrow(river_mouths)))

# ── Cities (for dist_city predictor) ──────────────────────────────────────────
CITY_CACHE <- file.path(OUT_DIR, "cities_sf.rds")
if (file.exists(CITY_CACHE)) {
  cities_sf <- readRDS(CITY_CACHE)
} else {
  message("Downloading Natural Earth populated places (scale=10)...")
  cities_sf <- tryCatch(
    rnaturalearth::ne_download(
      scale = 10, type = "populated_places",
      category = "cultural", returnclass = "sf"
    ) |> sf::st_transform(4326),
    error = \(e) {
      message("  Cities download failed: ", e$message, " — dist_city_km set to NA")
      NULL
    }
  )
  if (!is.null(cities_sf)) saveRDS(cities_sf, CITY_CACHE)
}
message("  Cities: ", if (!is.null(cities_sf)) nrow(cities_sf) else "unavailable")

# ── Bio-ORACLE rasters (for SDM background env values + global grid) ──────────
message("Loading Bio-ORACLE layers...")
options(sdmpredictors_datadir = "sdm_data")

# Layer codes confirmed from list_layers(datasets = "Bio-ORACLE"):
#   BO_sstmean      = Sea surface temperature (mean)
#   BO_salinitymean = Sea surface salinity (mean)
#   BO_damean       = Diffuse attenuation coefficient at 490 nm (mean) = Kd490
env_codes <- c(SST = "BO_sstmean", Salinity = "BO_salinity", Kd490 = "BO_damean")

rstack   <- sdmpredictors::load_layers(unname(env_codes), equalarea = FALSE)
env_rast <- terra::rast(rstack)
message("  Rasters loaded: ", paste(names(env_rast), collapse = ", "))

# ── Ecklonia cava complex — Akita et al. (2020) taxonomy ─────────────────────
# Akita, S., Hashimoto, K., Hanyuda, T. & Kawai, H. (2020) Phycologia 59, 330-339.
# E. cava, E. kurome and E. stolonifera are one species; the latter two become
# subspecies of E. cava. WoRMS and AlgaeBase have adopted this.
#
# The GBIF backbone carries only TWO species in the genus — E. radiata and
# E. maxima — so E. cava and its synonyms match nothing and every record
# bearing those names is indexed at GENUS rank with an empty species field.
# Any species-level query therefore drops the entire East Asian Ecklonia
# material: in v1 the genus contributed one species (E. radiata) and 13
# unusable genus-only records, which is why the kelp group was effectively a
# Southern-Hemisphere and North-Atlantic sample.
#
# R/01_ecklonia_gbif.R recovers the records through the verbatimScientificName
# parameter of the REST API and writes data/ecklonia_records.csv. Here they are
# merged in: synonyms already present in bc are relabelled, and the recovered
# coordinates are appended with environment extracted from the SAME env_rast as
# every other record, so there is one extraction path. Everything downstream —
# depth filter, over-water distance, discharge, 1-degree thinning — then treats
# them like any other occurrence.
ECK_CSV <- "data/ecklonia_records.csv"
if (file.exists(ECK_CSV)) {
  .spc <- intersect(c("species","Species","scientificName","taxon",
                      "verbatimScientificName"), names(bc))
  if (length(.spc) == 0) {
    message("  Ecklonia merge skipped: bc has no species column")
  } else {
    .spc <- .spc[1]
    .syn <- c("Ecklonia kurome", "Ecklonia stolonifera",
              "Ecklonia cava subsp. kurome", "Ecklonia cava subsp. stolonifera",
              "Ecklonia cava subsp. cava", "Ecklonia cava f. stolonifera")
    .hit <- !is.na(bc[[.spc]]) & bc[[.spc]] %in% .syn
    if (any(.hit)) {
      message(sprintf("  Ecklonia: relabelled %d existing record(s) to E. cava",
                      sum(.hit)))
      bc[[.spc]][.hit] <- "Ecklonia cava"
    }

    .eck <- utils::read.csv(ECK_CSV, stringsAsFactors = FALSE)
    .eck <- .eck[.eck$species_accepted == "Ecklonia cava", c("lon","lat")]
    # drop coordinates bc already holds for this species
    .have <- paste(round(bc$lon, 4), round(bc$lat, 4))[
               !is.na(bc[[.spc]]) & bc[[.spc]] == "Ecklonia cava"]
    .eck <- .eck[!paste(round(.eck$lon, 4), round(.eck$lat, 4)) %in% .have, ,
                 drop = FALSE]

    if (nrow(.eck) > 0) {
      .ev <- terra::extract(env_rast,
               terra::vect(as.matrix(.eck[, c("lon","lat")]), crs = "EPSG:4326"))[, -1,
               drop = FALSE]
      names(.ev) <- names(env_codes)[match(names(.ev), unname(env_codes))]
      if (any(is.na(names(.ev)))) names(.ev) <- names(env_codes)
      .new <- data.frame(lon = .eck$lon, lat = .eck$lat, Type = "Kelp",
                         .ev, stringsAsFactors = FALSE)
      .new[[.spc]] <- "Ecklonia cava"
      .keep <- is.finite(.new$SST) & is.finite(.new$Salinity) &
               is.finite(.new$Kd490) & .new$Kd490 > 0
      .new <- .new[.keep, , drop = FALSE]
      message(sprintf("  Ecklonia: %d recovered coordinates, %d with complete environment",
                      nrow(.eck), nrow(.new)))
      if (nrow(.new) > 0) bc <- dplyr::bind_rows(bc, .new)
    } else {
      message("  Ecklonia: no new coordinates to add")
    }

    # derived columns set before the merge have to be rebuilt for the new rows
    bc$Type_f  <- factor(bc$Type, levels = ECO_LVLS)
    bc$abs_lat <- abs(bc$lat)
    bc$log_Kd  <- log(bc$Kd490)
    message(sprintf("  bc after Ecklonia merge: %s records", format(nrow(bc), big.mark = ",")))
    print(table(bc$Type))
  }
} else {
  message("  ", ECK_CSV, " not found — run R/01_ecklonia_gbif.R to add the ",
          "East Asian Ecklonia records (see Methods)")
}

# ── Over-water distance to the nearest river mouth ───────────────────────────
# A great-circle distance can assign a cell to a mouth on the far side of a
# peninsula. Measured against a 2.8 km water mask, the great-circle distance
# understates the over-water distance by >50% for about a fifth of the
# occurrence records and by >2x for about one in eight (peninsula_check.R).
#
# The fix is a distance that has to stay in water. terra::gridDist() gives, for
# every water cell on Earth, the shortest in-water path to ANY river mouth, in a
# single pass. Two details matter:
#
#   1. Mouths and points that fall on a land cell are snapped to the nearest
#      water cell rather than dropped. At 2.8 km, 26% of mouths and half of the
#      occurrence records sit on a land cell, because both lie on the coastline
#      itself; dropping them would inflate distances everywhere.
#
#   2. The grid distance is quantised at the cell size, so near a mouth - where
#      the ecological signal is strongest and occupancy changes fastest - it is
#      far coarser than the great-circle distance. Replacing dist_km outright
#      would therefore destroy resolution in exactly the interval that matters.
#      We instead override the great-circle distance only where the water path
#      is longer by more than two cells, i.e. only where a genuine land barrier
#      exists. Everywhere else the metrically finer great-circle value is kept.
OW_AGG_FACT <- 3                       # 30 arc-sec x 3 -> ~2.8 km working grid
OW_CACHE    <- file.path(OUT_DIR, "overwater_dist.tif")
LAND_MASK   <- "wc2.1_30s_elev.tif"    # land = elevation, water = NA

.lastcol <- function(e) if (is.null(dim(e))) e else e[, ncol(e)]

if (file.exists(OW_CACHE)) {
  message("Loading cached over-water distance surface...")
  ow_dist <- terra::rast(OW_CACHE)
  OW_CELL_KM <- terra::res(ow_dist)[1] * 111
} else {
  message("Building over-water distance surface (one-time, a few minutes)...")
  if (file.exists(LAND_MASK)) {
    .elev <- terra::crop(terra::rast(LAND_MASK), terra::ext(-180, 180, -70, 80))
    # mean with na.rm = FALSE is NA if ANY subcell is water -> water mask that
    # keeps the coastal strip continuous
    .agg <- terra::aggregate(.elev, fact = OW_AGG_FACT, fun = "mean", na.rm = FALSE)
    .water <- terra::ifel(is.na(.agg), 1, NA)
    rm(.elev, .agg)
  } else {
    message("  ", LAND_MASK, " not found — falling back to the Bio-ORACLE domain\n",
            "  (coarser, ~9 km, but it is the ocean domain the models already use)")
    .water <- terra::ifel(is.na(env_rast[[1]]), NA, 1)
  }
  OW_CELL_KM <- terra::res(.water)[1] * 111
  message(sprintf("  water mask cell size: %.1f km", OW_CELL_KM))

  # snap every mouth to the nearest water cell so none is lost
  .mxy <- sf::st_coordinates(river_mouths)
  .onland <- is.na(.lastcol(terra::extract(.water, .mxy)))
  message(sprintf("  mouths on a land cell: %d of %d — snapping to nearest water",
                  sum(.onland), nrow(.mxy)))
  .step <- terra::res(.water)[1]
  .ang  <- seq(0, 2 * pi, length.out = 17)[-17]
  for (i in which(.onland)) {
    for (rr in seq(.step, 0.5, by = .step)) {
      .cand <- cbind(.mxy[i, 1] + rr * cos(.ang) / cos(.mxy[i, 2] * pi / 180),
                     .mxy[i, 2] + rr * sin(.ang))
      .ok <- which(!is.na(.lastcol(terra::extract(.water, .cand))))
      if (length(.ok)) { .mxy[i, ] <- .cand[.ok[1], ]; break }
    }
  }
  .mc <- terra::cellFromXY(.water, .mxy)
  .water[.mc[!is.na(.mc)]] <- 2
  message("  computing in-water distance to the nearest mouth...")
  ow_dist <- terra::gridDist(.water, target = 2)
  terra::writeRaster(ow_dist, OW_CACHE, overwrite = TRUE)
  rm(.water, .mxy, .mc); gc()
}
message(sprintf("  over-water surface ready (cell %.1f km)", OW_CELL_KM))

# Distance to the nearest river mouth for arbitrary points: great-circle,
# overridden by the in-water distance only where land genuinely blocks the path.
dist_to_mouth <- function(lon, lat) {
  psf  <- sf::st_as_sf(data.frame(lon = lon, lat = lat),
                       coords = c("lon", "lat"), crs = 4326)
  nnm  <- sf::st_nearest_feature(psf, river_mouths)
  d_gc <- as.numeric(sf::st_distance(psf, river_mouths[nnm, ],
                                     by_element = TRUE)) / 1000
  xy   <- cbind(lon, lat)
  d_ow <- .lastcol(terra::extract(ow_dist, xy)) / 1000
  # Points on a land cell take the nearest water cell's value. Done ring by ring
  # across all outstanding points at once: the prediction grid has tens of
  # thousands of cells and a per-point loop is far too slow.
  step <- terra::res(ow_dist)[1]
  ang  <- seq(0, 2 * pi, length.out = 17)[-17]
  miss <- which(is.na(d_ow))
  rr   <- step
  while (length(miss) && rr <= 0.5) {
    best <- rep(NA_real_, length(miss))
    for (a in ang) {
      cand <- cbind(xy[miss, 1] + rr * cos(a) / cos(xy[miss, 2] * pi / 180),
                    xy[miss, 2] + rr * sin(a))
      vv <- .lastcol(terra::extract(ow_dist, cand)) / 1000
      best <- pmin(best, vv, na.rm = TRUE)
    }
    d_ow[miss] <- best
    miss <- miss[is.na(best)]
    rr <- rr + step
  }
  blocked <- !is.na(d_ow) & d_ow > d_gc + 2 * OW_CELL_KM
  out <- ifelse(blocked, d_ow, d_gc)
  attr(out, "n_blocked")   <- sum(blocked)
  # nearest_idx is the great-circle nearest mouth and is used only to attach
  # UPLAND_SKM. gridDist returns a distance, not the identity of the mouth it
  # came from, so for a land-blocked cell the watershed area is still taken from
  # the geometrically nearest mouth rather than the reachable one. That is a
  # known inconsistency; it is immaterial because watershed area carries only
  # 0.4-2.8% of within-range importance and the paper's conclusion is that it is
  # uninformative. Resolving it would need a per-mouth basin allocation.
  attr(out, "nearest_idx") <- nnm
  out
}


# ── Bathymetry (MARSPEC MS_bathy_5m via sdmpredictors) ───────────────────────
# MS_bathy_5m: positive values = water depth (m) below sea surface; land = NA
BATHY_CACHE_RAST <- file.path(OUT_DIR, "depth_rast.tif")
if (!file.exists(BATHY_CACHE_RAST)) {
  message("Downloading bathymetry layer (MS_bathy_5m, ~1 min) ...")
  bathy_raw  <- sdmpredictors::load_layers("MS_bathy_5m", equalarea = FALSE)
  depth_rast <- terra::rast(bathy_raw); names(depth_rast) <- "depth"
  rm(bathy_raw)
  terra::writeRaster(depth_rast, BATHY_CACHE_RAST, overwrite = TRUE)
} else {
  message("Loading cached bathymetry raster...")
  depth_rast <- terra::rast(BATHY_CACHE_RAST)
}

# ── Tidal range from EOT20 (ocean_tides.zip) ─────────────────────────────────
# Source: Hart-Davis et al. (2021) EOT20, SEANOE doi:10.17882/79489
# Upload ocean_tides.zip to the project directory, then this block builds the
# cached tidal_rast.tif (amplitude in metres) from 8 main tidal constituents.
# Tidal range ≈ 2 × Σ amplitude  (M2+S2+K1+O1+N2+K2+P1+Q1), units: metres.
tidal_rast  <- NULL
TIDAL_SCALE <- 1
TIDAL_CACHE_RAST <- file.path(OUT_DIR, "tidal_rast.tif")
EOT20_ZIP    <- "ocean_tides.zip"           # place in working directory
EOT20_DIR    <- file.path(OUT_DIR, "eot20") # default extraction target

# Auto-detect pre-extracted EOT20 directory (searches common locations)
for (.cand in c(EOT20_DIR, "ocean_tides", file.path(OUT_DIR, "ocean_tides"),
                file.path(dirname(getwd()), "ocean_tides"))) {
  if (file.exists(file.path(.cand, "M2_ocean_eot20.nc"))) {
    EOT20_DIR <- .cand; message("  EOT20 found at: ", EOT20_DIR); break
  }
}

if (file.exists(TIDAL_CACHE_RAST)) {
  message("Loading cached tidal raster...")
  tidal_rast <- terra::rast(TIDAL_CACHE_RAST)
} else if (file.exists(EOT20_ZIP) || file.exists(file.path(EOT20_DIR, "M2_ocean_eot20.nc"))) {
  message("Building tidal range raster from EOT20...")
  if (!dir.exists(EOT20_DIR) || !file.exists(file.path(EOT20_DIR, "M2_ocean_eot20.nc"))) {
    message("  Extracting ocean_tides.zip (~2 GB, one-time) ...")
    unzip(EOT20_ZIP, exdir = EOT20_DIR)
    # Files land in EOT20_DIR/ocean_tides/ — flatten one level
    sub_dir <- file.path(EOT20_DIR, "ocean_tides")
    if (dir.exists(sub_dir))
      file.rename(list.files(sub_dir, full.names = TRUE),
                  file.path(EOT20_DIR, basename(list.files(sub_dir))))
  }
  # 8 main constituents: semidiurnal (M2,S2,N2,K2) + diurnal (K1,O1,P1,Q1)
  TIDAL_CONST <- c("M2","S2","K1","O1","N2","K2","P1","Q1")
  amp_sum <- NULL
  for (tc in TIDAL_CONST) {
    nc_f <- file.path(EOT20_DIR, sprintf("%s_ocean_eot20.nc", tc))
    if (!file.exists(nc_f)) { message("  Missing: ", basename(nc_f)); next }
    amp  <- terra::rast(nc_f, lyrs = "amplitude")   # units: cm
    amp_sum <- if (is.null(amp_sum)) amp else amp_sum + amp
    message(sprintf("  Added %s (max amplitude %.1f cm)", tc,
                    as.numeric(terra::global(amp, "max", na.rm = TRUE))))
  }
  # tidal range = 2 × Σ amp; convert cm → m
  tidal_rast <- terra::rotate(amp_sum * 2 / 100)   # longitude 0-360 → -180-180
  names(tidal_rast) <- "tidal_range_m"
  # Resample to match env_rast resolution (Bio-ORACLE ~0.083°)
  tidal_rast <- terra::resample(tidal_rast, env_rast[[1]], method = "bilinear")
  # Land cells (NA in env_rast) → 0
  tidal_rast[is.na(env_rast[[1]])] <- NA
  terra::writeRaster(tidal_rast, TIDAL_CACHE_RAST, overwrite = TRUE)
  message(sprintf("  Tidal range: max=%.2f m  — saved to %s",
                  as.numeric(terra::global(tidal_rast, "max", na.rm=TRUE)), TIDAL_CACHE_RAST))
} else {
  message("  ocean_tides.zip not found — tidal_m=0 (log_tidal=0, excluded by MaxEnt)")
}
if (!is.null(tidal_rast)) {
  raw_max     <- as.numeric(terra::global(tidal_rast, "max", na.rm = TRUE))
  TIDAL_SCALE <- if (raw_max > 100) 100 else 1   # already metres; scale=1 expected
  message(sprintf("  Tidal raster ready: max=%.2f m (TIDAL_SCALE=%d)", raw_max, TIDAL_SCALE))
}
# Extract depth at occurrence points
# MS_bathy_5m uses NEGATIVE values for ocean (ETOPO convention: surface=0, deep=-50).
# abs() converts to non-negative depth (metres below surface); land/NA → 0.
bc_vect  <- terra::vect(as.matrix(bc[, c("lon","lat")]), crs = "EPSG:4326")
depth_raw <- terra::extract(depth_rast, bc_vect)[, 2]
bc$depth     <- pmax(0, replace(abs(depth_raw), !is.finite(depth_raw), 0))
# Restrict to blue carbon depth range: Kelp max ~40 m; use 50 m as conservative ceiling
MAX_DEPTH_M <- 100
n_bc_before <- nrow(bc)
bc <- bc |> dplyr::filter(depth <= MAX_DEPTH_M)
bc$log_depth <- log10(bc$depth + 1)
# River-mouth distance, upstream catchment area and city distance for occurrences.
# dist_km is recomputed here (overriding whatever bc_analysis.rds supplied) so
# that presences, background and the prediction grid all use the same definition:
# distance to the nearest river MOUTH, measured through water where land blocks
# the straight line (see dist_to_mouth()).
bc_sf_up         <- sf::st_as_sf(bc[, c("lon","lat")], coords = c("lon","lat"), crs = 4326)
bc_vect_dep      <- terra::vect(bc_sf_up)
bc$dist_km       <- dist_to_mouth(bc$lon, bc$lat)
nn_bc_m          <- attr(bc$dist_km, "nearest_idx")
message(sprintf("  occurrences with a land-blocked straight line: %d of %d (%.1f%%)",
                attr(bc$dist_km, "n_blocked"), nrow(bc),
                100 * attr(bc$dist_km, "n_blocked") / nrow(bc)))
bc$dist_km       <- as.numeric(bc$dist_km)
bc               <- add_dist_forms(bc)
bc$upstream_km2  <- river_mouths$UPLAND_SKM[nn_bc_m]
bc$log_upstream  <- log10(bc$upstream_km2 + 1)
bc$discharge_cms <- river_mouths$DIS_AV_CMS[nn_bc_m]
bc$log_discharge <- log10(pmax(bc$discharge_cms, 0) + 1)
message(sprintf("  dist to river mouth (bc): median %.1f km | IQR %.1f-%.1f km",
                median(bc$dist_km), quantile(bc$dist_km, .25), quantile(bc$dist_km, .75)))
if (!is.null(tidal_rast)) {
  tidal_raw_bc <- terra::extract(tidal_rast, bc_vect_dep)[, 2]
  bc$tidal_m   <- pmax(0, replace(tidal_raw_bc, !is.finite(tidal_raw_bc), 0)) / TIDAL_SCALE
  message(sprintf("  Tidal range bc: median %.2f m | max %.2f m",
                  median(bc$tidal_m, na.rm=TRUE), max(bc$tidal_m, na.rm=TRUE)))
} else {
  bc$tidal_m <- 0
}
bc$log_tidal <- log10(bc$tidal_m + 1)
if (!is.null(cities_sf)) {
  nn_bc_city      <- sf::st_nearest_feature(bc_sf_up, cities_sf)
  bc$dist_city_km <- as.numeric(
    sf::st_distance(bc_sf_up, cities_sf[nn_bc_city, ], by_element = TRUE)) / 1000
} else {
  bc$dist_city_km <- 0  # no cities_sf: constant 0 → zero-var → MaxEnt drops it
}
bc$log_city <- log10(bc$dist_city_km + 1)
message(sprintf("  Depth filter (<=%.0f m): bc %d → %d records (%.0f%% kept)",
                MAX_DEPTH_M, n_bc_before, nrow(bc), 100 * nrow(bc) / n_bc_before))
message(sprintf("  Depth for bc: median %.0f m | max %.0f m | n_zero %d",
                median(bc$depth), max(bc$depth), sum(bc$depth == 0)))

# =============================================================================
# PART 2 — GAM WITH SPATIAL SMOOTH  [DISABLED — not used in final analysis]
# =============================================================================
if (FALSE) {  # skip GAM block
message("Fitting GAM...")

bc_coords <- cbind(bc$lon, bc$lat)
knn5      <- spdep::knearneigh(bc_coords, k = 5)
listw5    <- spdep::nb2listw(spdep::knn2nb(knn5), style = "W")
moran_res <- spdep::moran.test(bc$log_dist, listw5)
message(sprintf("Moran's I = %.4f, p = %.3e",
                moran_res$estimate["Moran I statistic"], moran_res$p.value))

gam_full <- mgcv::gam(
  log_dist ~
    Type_f +
    s(SST,      by = Type_f, k = 8, bs = "tp") +
    s(Salinity, by = Type_f, k = 8, bs = "tp") +
    s(log_Kd,   k = 8, bs = "tp") +
    s(abs_lat,  k = 8, bs = "tp") +
    s(lon, lat, bs = "sos", k = 150),
  data   = bc,
  method = "REML"
)
message(summary(gam_full))

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Fig 1: Map — n per type embedded in legend labels
n_per_type <- bc |> dplyr::count(Type_f) |>
  dplyr::mutate(label = sprintf("%s (n = %s)", Type_f, format(n, big.mark = ",")))
type_labels <- setNames(n_per_type$label, as.character(n_per_type$Type_f))

fig1 <- ggplot() +
  geom_sf(data = world, fill = "grey85", colour = "white", linewidth = 0.15) +
  geom_point(data = bc[sample(nrow(bc)), ],   # shuffle so no type always on top
             aes(lon, lat, colour = Type_f),
             size = 0.35, alpha = 0.55, shape = 16) +
  scale_colour_manual(values = ECO_PAL, labels = type_labels, name = NULL) +
  coord_sf(ylim = c(-65, 75), expand = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 2.8, alpha = 1))) +
  labs(caption = "GBIF occurrences after 50-km spatial thinning | HydroRIVERS ORD_FLOW ≤ 5") +
  theme_pub() +
  theme(axis.title = element_blank(),
        legend.position    = c(0.12, 0.32),
        legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
        legend.key.height  = unit(0.5, "lines"))
save_fig(fig1, "Fig1_Map", 170, 90)

# Fig Sup TypeMap: distribution per ecosystem type (faceted, supplement)
fig_stype_map <- ggplot() +
  geom_sf(data = world, fill = "grey85", colour = "white", linewidth = 0.12) +
  geom_point(data = bc[sample(nrow(bc)), ],
             aes(lon, lat, colour = Type_f),
             size = 0.3, alpha = 0.6, shape = 16) +
  scale_colour_manual(values = ECO_PAL, guide = "none") +
  facet_wrap(~Type_f, ncol = 3) +
  coord_sf(ylim = c(-65, 75), expand = FALSE) +
  theme_pub() +
  theme(axis.title  = element_blank(),
        axis.text   = element_text(size = 6),
        strip.text  = element_text(face = "bold", size = 9))
save_fig(fig_stype_map, "FigSup_TypeMap", 170, 72)

# Fig 2: ECDF with DKW bands
ecdf_df <- bc |>
  dplyr::group_by(Type) |>
  dplyr::mutate(n = dplyr::n()) |>
  dplyr::arrange(dist_km) |>
  dplyr::mutate(
    F    = seq_len(dplyr::n()) / dplyr::n(),
    eps  = sqrt(log(2 / 0.05) / (2 * dplyr::n())),
    F_lo = pmax(0, F - eps),
    F_hi = pmin(1, F + eps)
  ) |> dplyr::ungroup()

fig2 <- ggplot(ecdf_df, aes(dist_km, F, colour = Type, fill = Type)) +
  geom_ribbon(aes(ymin = F_lo, ymax = F_hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = ECO_PAL, name = NULL) +
  scale_fill_manual(  values = ECO_PAL, name = NULL) +
  scale_x_log10(labels = label_comma(suffix = " km"),
                breaks  = c(1, 5, 10, 50, 100, 500)) +
  annotation_logticks(sides = "b", linewidth = 0.2) +
  labs(x = "Distance to nearest major river (km)", y = "Cumulative probability") +
  theme_pub() + theme(legend.position = c(0.18, 0.8))
save_fig(fig2, "Fig2_ECDF", 120, 90)

# CI helper: robust across gratia versions that may lack .lower_ci/.upper_ci
add_ci <- function(df) {
  if (!".lower_ci" %in% names(df))
    df <- dplyr::mutate(df,
      .lower_ci = .estimate - 1.96 * .se,
      .upper_ci = .estimate + 1.96 * .se)
  df
}

# Fig 3: GAM partials — SST and Salinity
smth_sst <- gratia::smooth_estimates(gam_full,
  smooth = paste0("s(SST):Type_f", ECO_LVLS)) |>
  dplyr::mutate(Type = stringr::str_remove(.smooth, "s\\(SST\\):Type_f")) |>
  dplyr::filter(Type %in% ECO_LVLS) |> add_ci()

smth_sal <- gratia::smooth_estimates(gam_full,
  smooth = paste0("s(Salinity):Type_f", ECO_LVLS)) |>
  dplyr::mutate(Type = stringr::str_remove(.smooth, "s\\(Salinity\\):Type_f")) |>
  dplyr::filter(Type %in% ECO_LVLS) |> add_ci()

p3a <- ggplot(smth_sst, aes(SST, .estimate, colour = Type, fill = Type)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = ECO_PAL, name = NULL) +
  scale_fill_manual(  values = ECO_PAL, name = NULL) +
  labs(x = "SST (°C)", y = "Partial effect on log(dist+1)",
       title = "Sea surface temperature") +
  theme_pub() + theme(legend.position = "none")

p3b <- ggplot(smth_sal, aes(Salinity, .estimate, colour = Type, fill = Type)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = ECO_PAL, name = NULL) +
  scale_fill_manual(  values = ECO_PAL, name = NULL) +
  labs(x = "Salinity (PSU)", y = NULL, title = "Salinity") +
  theme_pub() + theme(legend.position = "right")

fig3 <- p3a + p3b + plot_annotation(theme = theme_pub())
save_fig(fig3, "Fig3_GAM_SST_Salinity", 170, 80)

# Fig 4: GAM partials — Kd490 and latitude
smth_kd  <- gratia::smooth_estimates(gam_full, smooth = "s(log_Kd)")  |> add_ci()
smth_lat <- gratia::smooth_estimates(gam_full, smooth = "s(abs_lat)") |> add_ci()

fig4 <- (
  ggplot(smth_kd, aes(log_Kd, .estimate)) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "grey60") +
    geom_line(linewidth = 1, colour = "grey20") +
    labs(x = expression(log(Kd490)~(m^{-1})), y = "Partial effect",
         title = "Turbidity (Kd490, log-scale)") + theme_pub()
) + (
  ggplot(smth_lat, aes(abs_lat, .estimate)) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "grey60") +
    geom_line(linewidth = 1, colour = "grey20") +
    labs(x = "Absolute latitude (°)", y = NULL, title = "Latitude") + theme_pub()
) + plot_annotation(theme = theme_pub())
save_fig(fig4, "Fig4_GAM_Kd_Lat", 170, 80)

# Fig 5: Spatial residuals
bc$spatial_resid <- residuals(gam_full, type = "deviance")
fig5 <- ggplot() +
  geom_sf(data = world, fill = "grey88", colour = "white", linewidth = 0.1) +
  geom_point(data = bc, aes(lon, lat, colour = spatial_resid),
             size = 0.5, alpha = 0.7, shape = 16) +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                          midpoint = 0, name = "Deviance\nresidual",
                          limits = c(-3, 3), oob = scales::squish) +
  coord_sf(ylim = c(-65, 75), expand = FALSE) +
  labs(caption = "Blue = closer to river than predicted | Red = farther") +
  theme_pub() + theme(axis.title = element_blank())
save_fig(fig5, "Fig5_SpatialResiduals", 170, 90)

# Fig 6: Ridge plots
fig6 <- ggplot(bc, aes(x = dist_km, y = Type_f, fill = Type_f, colour = Type_f)) +
  ggridges::geom_density_ridges(scale = 1.4, alpha = 0.6, linewidth = 0.6,
    quantile_lines = TRUE, quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(  values = ECO_PAL, guide = "none") +
  scale_colour_manual(values = ECO_PAL, guide = "none") +
  scale_x_log10(labels = label_comma(suffix = " km"),
                breaks  = c(1, 5, 10, 50, 100, 500)) +
  annotation_logticks(sides = "b", linewidth = 0.2) +
  labs(x = "Distance to nearest major river (km)", y = NULL) +
  theme_pub()
save_fig(fig6, "Fig6_Ridge", 120, 80)

# Fig 7: Observed vs predicted
bc$pred_dist <- expm1(predict(gam_full, bc,
                               exclude = "s(lon,lat)",
                               newdata.guaranteed = TRUE))
fig7 <- ggplot(bc, aes(pred_dist, dist_km, colour = Type_f)) +
  geom_point(size = 0.35, alpha = 0.3, shape = 16) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.5, colour = "grey40") +
  facet_wrap(~Type_f, ncol = 3) +
  scale_colour_manual(values = ECO_PAL, guide = "none") +
  scale_x_log10(labels = label_comma(suffix = " km")) +
  scale_y_log10(labels = label_comma(suffix = " km")) +
  annotation_logticks(sides = "bl", linewidth = 0.2) +
  labs(x = "Predicted (km)", y = "Observed (km)") +
  theme_pub()
save_fig(fig7, "Fig7_Obs_vs_Pred", 170, 75)

} # end if(FALSE) — GAM block

# =============================================================================
# PART 3 — TURBIDITY DEEP-DIVE
# =============================================================================

# Fig T1: dist_km × Kd490 scatter + LOESS
fig_t1 <- ggplot(bc, aes(dist_km, Kd490, colour = Type, fill = Type)) +
  geom_point(size = 0.3, alpha = 0.2, shape = 16) +
  geom_smooth(method = "loess", span = 0.6, linewidth = 1, se = TRUE, alpha = 0.18) +
  scale_colour_manual(values = ECO_PAL, name = NULL) +
  scale_fill_manual(  values = ECO_PAL, name = NULL) +
  scale_x_log10(labels = label_comma(suffix = " km"), breaks = c(1, 10, 100, 500)) +
  scale_y_log10(labels = label_scientific()) +
  annotation_logticks(sides = "bl", linewidth = 0.2) +
  labs(x = "Distance to nearest major river (km)",
       y = expression(Kd490~(m^{-1}))) +
  theme_pub() + theme(legend.position = c(0.84, 0.85))
save_fig(fig_t1, "FigT1_Dist_vs_Kd490", 120, 90)

# Fig T2: Kd490 by distance zone
bc_zones <- bc |>
  dplyr::mutate(DistZone = cut(dist_km,
    breaks = c(0, 10, 25, 50, 100, 200, Inf),
    labels = c("<10","10–25","25–50","50–100","100–200",">200"), right = FALSE))

fig_t2 <- ggplot(bc_zones, aes(DistZone, Kd490, fill = Type, colour = Type)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.25, outlier.alpha = 0.2,
               linewidth = 0.35) +
  stat_summary(fun = median, geom = "point", shape = 18, size = 2,
               colour = "grey10", position = position_dodge(0.75)) +
  ggh4x::facet_wrap2(~Type, ncol = 1, strip = eco_strip()) +
  scale_fill_manual(  values = ECO_PAL, guide = "none") +
  scale_colour_manual(values = ECO_PAL, guide = "none") +
  scale_y_log10(labels = label_scientific()) +
  labs(x = "Distance zone (km)", y = expression(Kd490~(m^{-1}))) +
  theme_pub()
save_fig(fig_t2, "FigT2_Kd490_byZone", 120, 120)

# Fig T3: 2D joint density
fig_t3 <- ggplot(bc, aes(dist_km, Kd490)) +
  stat_density_2d_filled(aes(fill = after_stat(level)),
    contour_var = "ndensity", breaks = seq(0, 1, by = 0.1)) +
  geom_density_2d(colour = "white", linewidth = 0.2, alpha = 0.5) +
  facet_wrap(~Type_f, ncol = 3) +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_scientific()) +
  scale_fill_brewer(palette = "YlOrRd", direction = 1, name = "Relative\ndensity") +
  labs(x = "Distance to nearest river (km)", y = expression(Kd490~(m^{-1}))) +
  theme_pub()
save_fig(fig_t3, "FigT3_2D_Density", 170, 80)

cor_tbl <- bc |>
  dplyr::group_by(Type) |>
  dplyr::summarise(n = dplyr::n(),
                   rho = cor(dist_km, Kd490, method = "spearman"),
                   p   = cor.test(dist_km, Kd490, method = "spearman")$p.value,
                   .groups = "drop")
write.csv(cor_tbl, file.path(OUT_DIR, "Spearman_Kd490_dist.csv"), row.names = FALSE)
message("Spearman rho:"); print(cor_tbl)

# Full predictor correlation matrix (collinearity diagnostic)
cor_matrix <- bc |>
  dplyr::select(SST, Salinity, Kd490, dist_km, abs_lat) |>
  dplyr::mutate(log_Kd490 = log(Kd490), log_dist = log1p(dist_km)) |>
  dplyr::select(SST, Salinity, log_Kd490, log_dist, abs_lat) |>
  cor(method = "spearman", use = "complete.obs")
write.csv(round(as.data.frame(cor_matrix), 3),
          file.path(OUT_DIR, "Spearman_CorMatrix.csv"))
message("Predictor correlation matrix (Spearman):")
print(round(cor_matrix, 3))

# =============================================================================
# PART 4 — DISTANCE ANALYSIS EXTENSIONS
# =============================================================================

# Fig D1: Violin + box + jitter
fig_d1 <- ggplot(bc, aes(Type_f, dist_km, fill = Type_f, colour = Type_f)) +
  geom_violin(alpha = 0.45, scale = "width", trim = TRUE, linewidth = 0.3) +
  geom_boxplot(width = 0.1, fill = "white", colour = "grey20",
               outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(width = 0.2, size = 0.25, alpha = 0.12, shape = 16) +
  stat_summary(fun = median, geom = "point", size = 3.5, colour = "grey10", shape = 18) +
  scale_fill_manual(  values = ECO_PAL, guide = "none") +
  scale_colour_manual(values = ECO_PAL, guide = "none") +
  scale_y_log10(labels = label_comma(suffix = " km"), breaks = c(1, 5, 10, 50, 100, 500)) +
  annotation_logticks(sides = "l", linewidth = 0.2) +
  labs(x = NULL, y = "Distance to nearest major river (km)") +
  theme_pub()
save_fig(fig_d1, "FigD1_Violin", 89, 90)

# Fig D2: Percentile profile
pct_df <- bc |>
  dplyr::group_by(Type) |>
  dplyr::summarise(p10 = quantile(dist_km,.10), p25 = quantile(dist_km,.25),
                   p50 = quantile(dist_km,.50), p75 = quantile(dist_km,.75),
                   p90 = quantile(dist_km,.90), .groups = "drop") |>
  tidyr::pivot_longer(-Type, names_to = "Percentile", values_to = "dist_km") |>
  dplyr::mutate(Percentile = factor(Percentile, levels = c("p10","p25","p50","p75","p90"),
                                    labels = c("10th","25th","50th","75th","90th")))

fig_d2 <- ggplot(pct_df, aes(Percentile, dist_km, colour = Type, group = Type)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_text(aes(label = round(dist_km)), vjust = -0.9, size = 2.7, fontface = "bold") +
  scale_colour_manual(values = ECO_PAL, name = NULL) +
  scale_y_log10(labels = label_comma(suffix = " km")) +
  labs(x = "Percentile", y = "Distance (km)") +
  theme_pub() + theme(legend.position = c(0.18, 0.8))
save_fig(fig_d2, "FigD2_Percentiles", 120, 90)

# Fig D3: Kernel density
fig_d3 <- ggplot(bc, aes(dist_km, fill = Type, colour = Type)) +
  geom_density(alpha = 0.4, linewidth = 0.6, bw = "SJ") +
  scale_fill_manual(  values = ECO_PAL, name = NULL) +
  scale_colour_manual(values = ECO_PAL, name = NULL) +
  scale_x_log10(labels = label_comma(suffix = " km"), breaks = c(1, 5, 10, 50, 100, 500)) +
  annotation_logticks(sides = "b", linewidth = 0.2) +
  labs(x = "Distance from nearest river mouth (km)", y = "Density") +
  theme_pub() + theme(legend.position = c(0.85, 0.8))
save_fig(fig_d3, "FigD3_Density", 120, 80)

# [FigD4 moved to after PART 5b switch — needs bg_df which is created in PART 5]

# =============================================================================
# PART 5 — SDM: MaxEnt + XGBoost Ensemble
# =============================================================================
# MaxEnt (maxnet): ecological gold standard for GBIF presence-background data
# XGBoost: gradient boosting, complements MaxEnt
# Ensemble = mean probability of both models

message("\n=== Building SDM ===")
# abs_lat removed (redundant with SST for tropical/polar gradient)
# Two-stage SDM:
#   Stage 1 — range model: FEAT_range = SST + lat  (full coastal bg, no matching)
#   Stage 2 — habitat model: FEAT = local factors  (SST+lat matched background)
#   Prediction = suit_range × suit_habitat
FEAT_range <- c("SST","lat")
# depth: domain filter only (<=MAX_DEPTH_M). log_depth excluded from FEAT:
# depth discriminates trivially (bg 0-100 m vs presences in shallow zone) and
# absorbs importance from the ecologically meaningful log_dist signal.
# River size enters as DISCHARGE, not catchment area. Area was the weakest
# predictor in v1 (0.3-3.3% importance), and part of that is measurement rather
# than ecology: area counts an arid, sparsely vegetated catchment the same as a
# wet one, so it is a poor proxy for the freshwater and dissolved load actually
# delivered. Discharge integrates area with runoff and is carried by
# HydroRIVERS itself (DIS_AV_CMS). log_upstream stays in the data frames and is
# compared against log_discharge in R/02_reviewer_analyses.R.
FEAT       <- c("Salinity","Kd490", DIST_VAR, "log_discharge")
FEAT_SP    <- FEAT  # same for species level

VAR_LABELS <- c(
  SST          = "SST (°C)",
  lat          = "Latitude (°)",
  Salinity     = "Salinity (PSU)",
  Kd490        = "Turbidity / Kd490 (m⁻¹)",
  log_dist     = "log₁₀(Dist. from river mouth + 1) (km)",
  dist_km      = "Dist. from river mouth (km)",
  log_depth    = "log₁₀(Water depth + 1) (m)",
  depth        = "Water depth (m)",
  log_upstream = "log₁₀(Watershed area + 1) (km²)",
  upstream_km2 = "Watershed area (km²)",
  log_discharge= "log₁₀(River discharge + 1) (m³ s⁻¹)",
  discharge_cms= "River discharge (m³ s⁻¹)",
  inv_dist     = "1 / (Dist. from river mouth + 1) (km⁻¹)",
  exp_dist     = "exp(−Dist. from river mouth / 50 km)",
  log_city     = "log₁₀(Dist. to city + 1) (km)",
  dist_city_km = "Dist. to nearest city (km)",
  log_tidal    = "log₁₀(Tidal range + 1) (m)",
  tidal_m      = "Tidal range (m)"
)

# Background: random ocean sample from env raster
set.seed(42)
bg_pts  <- terra::spatSample(env_rast[[1]], size = 20000,
                              method = "random", na.rm = TRUE, xy = TRUE)
bg_vect <- terra::vect(bg_pts[, c("x","y")], crs = "EPSG:4326")
bg_env  <- terra::extract(env_rast, bg_vect)[, -1]
# Rename columns to match env_codes names
names(bg_env) <- names(env_codes)[match(names(bg_env), unname(env_codes))]
if (any(is.na(names(bg_env)))) names(bg_env) <- names(env_codes)

message("Computing river-mouth distances for background points...")
bg_sf   <- sf::st_as_sf(as.data.frame(bg_pts[, c("x","y")]),
                          coords = c("x","y"), crs = 4326)
.bgd    <- dist_to_mouth(bg_pts$x, bg_pts$y)
nn_bg   <- attr(.bgd, "nearest_idx")
message(sprintf("  background points with a land-blocked straight line: %d of %d (%.1f%%)",
                attr(.bgd, "n_blocked"), length(.bgd),
                100 * attr(.bgd, "n_blocked") / length(.bgd)))
bg_dist <- as.numeric(.bgd)
bg_upstream  <- river_mouths$UPLAND_SKM[nn_bg]
bg_discharge <- river_mouths$DIS_AV_CMS[nn_bg]
if (!is.null(cities_sf)) {
  nn_bg_city   <- sf::st_nearest_feature(bg_sf, cities_sf)
  bg_city_dist <- as.numeric(
    sf::st_distance(bg_sf, cities_sf[nn_bg_city, ], by_element = TRUE)) / 1000
} else {
  bg_city_dist <- rep(0, nrow(bg_pts))  # no cities_sf: constant 0
}

bg_df <- data.frame(lon = bg_pts$x, lat = bg_pts$y, dist_km = bg_dist,
                     abs_lat = abs(bg_pts$y), bg_env,
                     upstream_km2 = bg_upstream, discharge_cms = bg_discharge,
                     dist_city_km = bg_city_dist,
                     presence = 0L) |>
  dplyr::filter(!is.na(SST), !is.na(Salinity), !is.na(Kd490),
                is.finite(dist_km), Kd490 > 0)

bg_vect2    <- terra::vect(as.matrix(bg_df[, c("lon","lat")]), crs = "EPSG:4326")
depth_raw2  <- terra::extract(depth_rast, bg_vect2)[, 2]
bg_df$depth <- pmax(0, replace(abs(depth_raw2), !is.finite(depth_raw2), 0))
if (!is.null(tidal_rast)) {
  tidal_raw2    <- terra::extract(tidal_rast, bg_vect2)[, 2]
  bg_df$tidal_m <- pmax(0, replace(tidal_raw2, !is.finite(tidal_raw2), 0)) / TIDAL_SCALE
} else {
  bg_df$tidal_m <- 0
}
bg_df <- bg_df |> dplyr::filter(depth <= MAX_DEPTH_M)
bg_df               <- add_dist_forms(bg_df)
bg_df$log_depth     <- log10(bg_df$depth        + 1)
bg_df$log_upstream  <- log10(bg_df$upstream_km2 + 1)
bg_df$log_discharge <- log10(pmax(bg_df$discharge_cms, 0) + 1)
bg_df$log_city     <- log10(bg_df$dist_city_km + 1)
bg_df$log_tidal    <- log10(bg_df$tidal_m      + 1)
message(sprintf("  bg_df after depth filter (<=%.0f m): %d points", MAX_DEPTH_M, nrow(bg_df)))

lon_cuts   <- seq(-180, 180, length.out = 5)
assign_blk <- function(lon) as.integer(cut(lon, breaks = lon_cuts, include.lowest = TRUE))

all_mx <- list(); all_xgb <- list(); all_auc <- list()
all_vimp <- list(); all_pdp <- list()

for (tp in ECO_LVLS) {
  message(sprintf("\n--- %s ---", tp))

  sel_cols <- unique(c("lon","lat","abs_lat","SST", FEAT))
  pres <- bc |> dplyr::filter(Type == tp) |>
    dplyr::select(all_of(sel_cols)) |>
    dplyr::mutate(presence = 1L)

  # SST + lat matched background: control for climate zone and biogeographic range
  sst_q <- quantile(pres$SST, c(0.05, 0.95), na.rm = TRUE)
  lat_q <- quantile(pres$lat, c(0.05, 0.95), na.rm = TRUE)
  bg_pool <- bg_df |> dplyr::filter(SST >= sst_q[1], SST <= sst_q[2],
                                     lat >= lat_q[1], lat <= lat_q[2])
  if (nrow(bg_pool) < 100) bg_pool <- bg_df
  n_bg  <- min(nrow(pres) * 10, nrow(bg_pool), 6000)
  bg_sub <- bg_pool |>
    dplyr::select(all_of(unique(c(sel_cols, "presence")))) |>
    dplyr::slice_sample(n = n_bg)

  df <- dplyr::bind_rows(pres, bg_sub) |>
    dplyr::filter(if_all(all_of(FEAT), is.finite), Kd490 > 0) |>
    dplyr::mutate(block = assign_blk(lon))

  X <- as.matrix(df[, FEAT]); colnames(X) <- FEAT
  y <- df$presence

  # Remove zero-variance features from MaxEnt only (depth=0 everywhere crashes maxnet)
  feat_var  <- apply(X, 2, var)
  FEAT_mx   <- names(feat_var[feat_var > 1e-10])
  if (length(FEAT_mx) < length(FEAT))
    message("  Zero-var removed from MaxEnt: ", paste(setdiff(FEAT, FEAT_mx), collapse=", "))
  X_mx <- X[, FEAT_mx, drop = FALSE]

  # MaxEnt
  set.seed(42)
  mx_mod <- maxnet::maxnet(
    p    = y, data = as.data.frame(X_mx),
    f    = maxnet::maxnet.formula(y, as.data.frame(X_mx), classes = "lqph"),
    regmult = 2.0
  )

  # XGBoost
  set.seed(42)
  xgb_mod <- xgboost::xgb.train(
    params = list(objective = "binary:logistic", eval_metric = "auc",
                  max_depth = 5, eta = 0.05, subsample = 0.8,
                  colsample_bytree = 0.8, min_child_weight = 5),
    data    = xgboost::xgb.DMatrix(X, label = y),
    nrounds = 400, verbose = 0
  )

  # Spatial CV AUC
  cv_mx <- cv_xgb <- cv_ens <- numeric(4)
  for (b in 1:4) {
    tr <- df[df$block != b, ]; te <- df[df$block == b, ]
    if (length(unique(te$presence)) < 2 || nrow(te) < 10) next
    Xtr <- as.matrix(tr[, FEAT]); Xte <- as.matrix(te[, FEAT])
    colnames(Xtr) <- colnames(Xte) <- FEAT
    # Use same zero-variance-safe feature set as full model
    Xtr_mx <- Xtr[, FEAT_mx, drop = FALSE]; Xte_mx <- Xte[, FEAT_mx, drop = FALSE]

    mx_cv  <- maxnet::maxnet(tr$presence, as.data.frame(Xtr_mx),
                              maxnet::maxnet.formula(tr$presence, as.data.frame(Xtr_mx),
                                                     classes = "lqph"), regmult = 2)
    px_mx  <- predict(mx_cv, as.data.frame(Xte_mx), type = "cloglog") |> as.numeric()

    xgb_cv <- xgboost::xgb.train(
      params  = list(objective = "binary:logistic", max_depth = 5, eta = 0.05,
                     subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 5),
      data    = xgboost::xgb.DMatrix(Xtr, label = tr$presence),
      nrounds = 300, verbose = 0
    )
    px_xgb <- predict(xgb_cv, xgboost::xgb.DMatrix(Xte))
    px_ens <- (px_mx + px_xgb) / 2

    cv_mx[b]  <- pROC::auc(te$presence, px_mx,  quiet = TRUE) |> as.numeric()
    cv_xgb[b] <- pROC::auc(te$presence, px_xgb, quiet = TRUE) |> as.numeric()
    cv_ens[b] <- pROC::auc(te$presence, px_ens, quiet = TRUE) |> as.numeric()
  }
  v <- cv_mx > 0
  message(sprintf("  AUC — MaxEnt: %.3f | XGBoost: %.3f | Ensemble: %.3f",
                  mean(cv_mx[v]), mean(cv_xgb[v]), mean(cv_ens[v])))

  # Variable importance
  xgb_vimp <- xgboost::xgb.importance(model = xgb_mod)
  # MaxEnt permutation importance — only over features the model actually used
  base_pred <- predict(mx_mod, as.data.frame(X_mx), type = "cloglog") |> as.numeric()
  base_auc  <- pROC::auc(y, base_pred, quiet = TRUE) |> as.numeric()
  mx_vimp_used <- sapply(FEAT_mx, \(fv) {
    Xp <- as.data.frame(X_mx); Xp[[fv]] <- sample(Xp[[fv]])
    pa <- pROC::auc(y, predict(mx_mod, Xp, type="cloglog") |> as.numeric(),
                    quiet=TRUE) |> as.numeric()
    max(0, base_auc - pa)
  })
  # Map back to full FEAT vector (zero-var features get 0 importance)
  mx_vimp <- setNames(rep(0, length(FEAT)), FEAT)
  mx_vimp[FEAT_mx] <- mx_vimp_used
  mx_vimp <- mx_vimp / max(sum(mx_vimp), 1e-9)

  # Partial dependence
  kd_seq   <- exp(seq(log(min(bc$Kd490)), log(max(bc$Kd490)), length.out = 80))
  dist_seq <- seq(0, log10(1 + quantile(bc$dist_km, 0.90)), length.out = 80)
  mean_row_mx <- as.data.frame(t(colMeans(X_mx)))
  mean_row    <- as.data.frame(t(colMeans(X)))

  make_pd <- function(var, vals) {
    purrr::map_dfr(vals, \(v) {
      # MaxEnt uses X_mx; XGBoost uses full X
      r_mx <- mean_row_mx; if (var %in% FEAT_mx) r_mx[[var]] <- v
      r    <- mean_row;    r[[var]] <- v
      Xr   <- as.matrix(r[, FEAT]); colnames(Xr) <- FEAT
      data.frame(x = v,
                 MaxEnt  = predict(mx_mod, r_mx, type = "cloglog") |> as.numeric(),
                 XGBoost = predict(xgb_mod, xgboost::xgb.DMatrix(Xr)))
    }) |> dplyr::mutate(Type = tp, Variable = var)
  }

  all_mx[[tp]]   <- mx_mod
  all_xgb[[tp]]  <- xgb_mod
  all_auc[[tp]]  <- tibble::tibble(Type = tp,
    MaxEnt = mean(cv_mx[v]), XGBoost = mean(cv_xgb[v]), Ensemble = mean(cv_ens[v]))
  xgb_raw <- xgb_vimp$Gain[match(FEAT, xgb_vimp$Feature)]
  xgb_raw[is.na(xgb_raw)] <- 0
  all_vimp[[tp]] <- tibble::tibble(Type = tp, Variable = FEAT,
    MaxEnt  = mx_vimp,
    XGBoost = xgb_raw / sum(xgb_raw, na.rm = TRUE))
  all_pdp[[tp]]  <- dplyr::bind_rows(make_pd("Kd490", kd_seq),
                                      make_pd("dist_km", dist_seq))
}

# FigS3 moved to after PART 5c — now uses species-level models for CI (±SD across species)

# =============================================================================
# PART 5b — OBSERVATION-DENSITY WEIGHTED COASTAL BACKGROUND
# =============================================================================
# Strategy: restrict background to coastal zone (dist_km ≤ 95th pctile of
# all presence records) AND weight each background point by its proximity to
# the nearest occurrence record across all taxa.
#   weight_i = exp(−dist_to_nearest_obs_i / 500 km)
# This approximates survey effort bias: areas visited by marine biologists
# (where any of the three taxa were recorded) get over-represented in the
# background, mimicking the spatial observation bias in the presence data.
# Result: MaxEnt sees background that reflects the same accessible-coast
# sampling pattern as the presences → dist_km and lon/lat explain less.

message("\n=== Weighted coastal background (observation-density corrected) ===")

coast_thresh <- quantile(bc$dist_km, 0.95)
# No river-distance cap on background: wider dist_km range in background
# increases contrast with presences (near rivers) → more informative for log_dist importance.
bg_coastal   <- bg_df  # use full depth-filtered background (no dist_km cap)
message(sprintf("  Coastal threshold for maps: %.0f km | bg_coastal = full bg_df (%d pts, no dist cap)",
                coast_thresh, nrow(bg_coastal)))

# Uniform background: satellite-derived predictors (SST, Salinity, Kd490) are
# geographically unbiased, so geographic sampling bias in GBIF occurrences is
# addressed by spatial thinning of presences rather than background weighting.
bg_coastal$obs_weight <- 1L
message("  obs_weight = 1 (uniform random background)")

bias_auc  <- list(); bias_vimp <- list(); bias_mx <- list(); bias_shap <- list()
bias_roc  <- list()  # pooled CV predictions for ROC curves
range_mx  <- list()  # Stage 1: biogeographic range models (SST + lat)
range_auc <- list()  # training AUC for Stage 1 range models

for (tp in ECO_LVLS) {
  message(sprintf("  %s ...", tp))
  sel_cols2 <- unique(c("lon","lat","abs_lat","SST", FEAT))
  pres   <- bc |> dplyr::filter(Type == tp) |>
    dplyr::select(all_of(sel_cols2)) |>
    dplyr::mutate(presence = 1L)
  # SST + lat matched background: control for climate zone and biogeographic range
  sst_q2 <- quantile(pres$SST, c(0.05, 0.95), na.rm = TRUE)
  lat_q2 <- quantile(pres$lat, c(0.05, 0.95), na.rm = TRUE)
  bg_pool2 <- bg_coastal |> dplyr::filter(SST >= sst_q2[1], SST <= sst_q2[2],
                                           lat >= lat_q2[1], lat <= lat_q2[2])
  if (nrow(bg_pool2) < 100) {
    message("  SST+lat-match fallback to full bg_coastal (", tp, ")")
    bg_pool2 <- bg_coastal
  }
  n_bg   <- min(nrow(pres) * 10, nrow(bg_pool2), 6000)
  bg_sub <- bg_pool2 |>
    dplyr::select(all_of(unique(c(sel_cols2, "presence")))) |>
    dplyr::slice_sample(n = n_bg)

  df2 <- dplyr::bind_rows(pres, bg_sub) |>
    dplyr::filter(if_all(all_of(FEAT), is.finite), Kd490 > 0) |>
    dplyr::mutate(block = assign_blk(lon))
  X2 <- as.matrix(df2[, FEAT]); colnames(X2) <- FEAT; y2 <- df2$presence

  # Zero-variance guard for PART 5b maxnet
  feat_var2 <- apply(X2, 2, var)
  FEAT_mx2  <- names(feat_var2[feat_var2 > 1e-10])
  X2_mx     <- X2[, FEAT_mx2, drop = FALSE]

  set.seed(42)
  mx2 <- maxnet::maxnet(y2, as.data.frame(X2_mx),
    maxnet::maxnet.formula(y2, as.data.frame(X2_mx), classes = "lqph"),
    regmult = 2.0)

  cv2 <- numeric(4); roc_pool_mx <- list()
  for (b in 1:4) {
    tr2 <- df2[df2$block != b, ]; te2 <- df2[df2$block == b, ]
    if (length(unique(te2$presence)) < 2 || nrow(te2) < 10) next
    Xtr2 <- as.matrix(tr2[, FEAT_mx2]); Xte2 <- as.matrix(te2[, FEAT_mx2])
    colnames(Xtr2) <- colnames(Xte2) <- FEAT_mx2
    mx2cv <- maxnet::maxnet(tr2$presence, as.data.frame(Xtr2),
      maxnet::maxnet.formula(tr2$presence, as.data.frame(Xtr2), classes = "lqph"),
      regmult = 2)
    px2    <- predict(mx2cv, as.data.frame(Xte2), type = "cloglog") |> as.numeric()
    cv2[b] <- pROC::auc(te2$presence, px2, quiet = TRUE) |> as.numeric()
    roc_pool_mx[[b]] <- data.frame(y = te2$presence, pred = px2)
  }
  v2 <- cv2 > 0

  bp2 <- predict(mx2, as.data.frame(X2_mx), type = "cloglog") |> as.numeric()
  ba2 <- pROC::auc(y2, bp2, quiet = TRUE) |> as.numeric()
  vi2_used <- sapply(FEAT_mx2, \(fv) {
    Xp <- as.data.frame(X2_mx); Xp[[fv]] <- sample(Xp[[fv]])
    pa <- pROC::auc(y2, predict(mx2, Xp, type = "cloglog") |> as.numeric(),
                    quiet = TRUE) |> as.numeric()
    max(0, ba2 - pa)
  })
  vi2 <- setNames(rep(0, length(FEAT)), FEAT)
  vi2[FEAT_mx2] <- vi2_used
  vi2 <- vi2 / max(sum(vi2), 1e-9)

  # ── Balanced RF + TreeSHAP (same training data as MaxEnt above) ──────────
  # Case weights: upweight presence to balance 1:1 effective ratio.
  # Tree SHAP (Lundberg et al. 2020): exact, model-agnostic directional importance.
  n_pos  <- sum(y2 == 1L); n_neg <- sum(y2 == 0L)
  cw_rf  <- ifelse(y2 == 1L, n_neg / n_pos, 1.0)
  df2_rf <- df2[, c("presence", FEAT)] |>
    dplyr::mutate(presence = factor(presence, levels = c("0","1")))

  rf2 <- ranger::ranger(
    presence ~ .,
    data         = df2_rf,
    num.trees    = 500,
    case.weights = cw_rf,
    probability  = TRUE,
    num.threads  = 2,
    seed         = 42
  )

  # SHAP via fastshap (Monte Carlo permutation SHAP — reliable with ranger)
  # Falls back to RF permutation importance if fastshap is unavailable.
  sv <- NULL; rf2_vimp_sh <- NULL
  X_rf <- as.data.frame(df2[, FEAT])
  # fastshap is archived on CRAN; do not attempt to install it. If it is absent,
  # go straight to the permutation-importance fallback below.
  tryCatch({
    if (!HAS_SHAP) stop("fastshap/shapviz not installed")
    pred_fun <- function(model, newdata) {
      predict(model, data = newdata)$predictions[, "1"]
    }
    shap_mat <- fastshap::explain(
      rf2, X = X_rf, pred_wrapper = pred_fun,
      nsim = 50
    )
    shap_mat <- as.matrix(shap_mat)
    colnames(shap_mat) <- FEAT
    shap_mat_lab <- shap_mat
    colnames(shap_mat_lab) <- VAR_LABELS[FEAT]
    X_rf_lab <- X_rf; colnames(X_rf_lab) <- VAR_LABELS[FEAT]
    sv <- shapviz::shapviz(shap_mat_lab, X = X_rf_lab)
    shap_raw <- colMeans(abs(shap_mat))
    names(shap_raw) <- FEAT
    rf2_vimp_sh <- shap_raw / sum(shap_raw)
    message("    SHAP (fastshap): OK")
  }, error = function(e) {
    message("    WARNING: fastshap failed (", conditionMessage(e),
            ") — falling back to RF permutation importance")
    rf2_pi <- ranger::ranger(
      presence ~ ., data = df2_rf, num.trees = 500,
      importance = "permutation", case.weights = cw_rf,
      probability = TRUE, num.threads = 2, seed = 42
    )
    raw_pi <- pmax(0, rf2_pi$variable.importance)
    rf2_vimp_sh <<- raw_pi / sum(raw_pi)
  })

  # RF spatial CV AUC
  cv_rf2 <- numeric(4); roc_pool_rf <- list()
  for (b in 1:4) {
    tr_b <- df2_rf[df2$block != b, ]; te_b <- df2_rf[df2$block == b, ]
    te_y <- as.integer(as.character(te_b$presence))
    if (length(unique(te_y)) < 2 || nrow(te_b) < 10) next
    tr_y <- as.integer(as.character(tr_b$presence))
    cw_b <- ifelse(tr_y == 1L, sum(tr_y == 0L) / sum(tr_y == 1L), 1.0)
    rf_cv <- ranger::ranger(
      presence ~ ., data = tr_b, num.trees = 300,
      probability = TRUE, case.weights = cw_b, num.threads = 2, seed = 42
    )
    pred_rf   <- predict(rf_cv, te_b[, FEAT])$predictions[, "1"]
    cv_rf2[b] <- pROC::auc(te_y, pred_rf, quiet = TRUE) |> as.numeric()
    roc_pool_rf[[b]] <- data.frame(y = te_y, pred = pred_rf)
  }
  v_rf2 <- cv_rf2 > 0

  # Pool CV folds into one data frame for ROC curves
  pool_mx_df <- do.call(rbind, Filter(Negate(is.null), roc_pool_mx))
  pool_rf_df <- do.call(rbind, Filter(Negate(is.null), roc_pool_rf))

  bias_auc[[tp]]  <- tibble::tibble(Type = tp,
    AUC_random   = all_auc[[tp]]$MaxEnt,
    AUC_weighted = mean(cv2[v2]),
    AUC_rf       = mean(cv_rf2[v_rf2]))
  bias_vimp[[tp]] <- tibble::tibble(Type = tp, Variable = FEAT,
    Random   = all_vimp[[tp]]$MaxEnt,
    Weighted = vi2,
    RF_SHAP  = rf2_vimp_sh)
  bias_mx[[tp]]   <- mx2
  bias_shap[[tp]] <- sv   # shapviz object for beeswarm figure
  bias_roc[[tp]]  <- list(
    mx_pool   = pool_mx_df,
    mx_folds  = Filter(Negate(is.null), roc_pool_mx),
    mx_cv_auc = cv2[v2],
    rf_pool   = pool_rf_df
  )
  message(sprintf("    AUC: random=%.3f  MaxEnt-coastal=%.3f  RF-coastal=%.3f",
                  all_auc[[tp]]$MaxEnt, mean(cv2[v2]), mean(cv_rf2[v_rf2])))

  # Stage 1: biogeographic range model — full coastal bg (no SST/lat matching intentional)
  n_bg_r <- min(nrow(pres) * 10, nrow(bg_coastal), 6000)
  bg_r   <- bg_coastal |>
    dplyr::slice_sample(n = n_bg_r) |>
    dplyr::select(all_of(FEAT_range)) |>
    dplyr::mutate(presence = 0L)
  df_r <- dplyr::bind_rows(
    pres |> dplyr::select(all_of(FEAT_range)) |> dplyr::mutate(presence = 1L),
    bg_r
  ) |>
    dplyr::mutate(dplyr::across(all_of(FEAT_range),
                                \(x) dplyr::if_else(is.na(x), median(x, na.rm = TRUE), x)))
  y_r <- df_r$presence
  X_r <- as.data.frame(df_r[, FEAT_range])
  range_mx[[tp]] <- tryCatch(
    maxnet::maxnet(y_r, X_r,
                   maxnet::maxnet.formula(y_r, X_r, classes = "lqph"),
                   regmult = 1.0),
    error = \(e) { message("  range model FAIL (", tp, "): ", e$message); NULL }
  )
  if (!is.null(range_mx[[tp]])) {
    p_r <- predict(range_mx[[tp]], X_r, type = "cloglog") |> as.numeric()
    range_auc[[tp]] <- tryCatch(pROC::auc(y_r, p_r, quiet = TRUE) |> as.numeric(),
                                error = \(e) NA_real_)
    message(sprintf("  Stage 1 range model (%s): AUC=%.3f (n_pres=%d n_bg=%d)",
                    tp, range_auc[[tp]], sum(y_r==1L), n_bg_r))
  }
}

# Switch primary model to weighted coastal for all downstream predictions
all_shap <- list(); all_roc <- list()
for (tp in ECO_LVLS) {
  all_mx[[tp]]           <- bias_mx[[tp]]
  all_auc[[tp]]$MaxEnt   <- bias_auc[[tp]]$AUC_weighted
  all_auc[[tp]]$RF       <- bias_auc[[tp]]$AUC_rf
  all_vimp[[tp]]$MaxEnt  <- bias_vimp[[tp]]$Weighted
  all_vimp[[tp]]$RF_SHAP <- bias_vimp[[tp]]$RF_SHAP
  all_shap[[tp]]         <- bias_shap[[tp]]
  all_roc[[tp]]          <- bias_roc[[tp]]
}
message("Primary SDM: observation-density weighted coastal background")

message("  Primary SDM updated: all_mx, all_auc, all_vimp switched to weighted coastal")

# ── Checkpoint 1: save SDM results ──────────────────────────────────────────
RDS_DIR <- file.path(OUT_DIR, "rds"); dir.create(RDS_DIR, showWarnings = FALSE)
saveRDS(bc,           file.path(RDS_DIR, "bc.rds"))
saveRDS(bg_df,        file.path(RDS_DIR, "bg_df.rds"))
saveRDS(bg_coastal,   file.path(RDS_DIR, "bg_coastal.rds"))
saveRDS(coast_thresh, file.path(RDS_DIR, "coast_thresh.rds"))
saveRDS(all_mx,       file.path(RDS_DIR, "all_mx.rds"))
saveRDS(range_mx,     file.path(RDS_DIR, "range_mx.rds"))
saveRDS(range_auc,    file.path(RDS_DIR, "range_auc.rds"))   # Stage 1 AUC, for Fig. S4
saveRDS(all_auc,      file.path(RDS_DIR, "all_auc.rds"))
saveRDS(all_vimp,     file.path(RDS_DIR, "all_vimp.rds"))
saveRDS(all_shap,     file.path(RDS_DIR, "all_shap.rds"))
saveRDS(all_roc,      file.path(RDS_DIR, "all_roc.rds"))
saveRDS(all_pdp,      file.path(RDS_DIR, "all_pdp.rds"))
if (exists("RF_SHAP") && !is.null(RF_SHAP))
  saveRDS(RF_SHAP, file.path(RDS_DIR, "RF_SHAP.rds"))
message("Checkpoint 1 saved → ", RDS_DIR)

# =============================================================================
# FIGURE A1 — Predictor correlation matrix (Appendix)
# =============================================================================
CORR_LABELS <- c(Salinity      = "Salinity",
                 Kd490         = "Kd490",
                 log_dist      = "log₁₀(Dist+1)",
                 inv_dist      = "1/(Dist+1)",
                 exp_dist      = "exp(−Dist/50)",
                 log_upstream  = "log₁₀(Wshed+1)",
                 log_discharge = "log₁₀(Discharge+1)")

cor_mat <- bc |>
  dplyr::select(all_of(FEAT)) |>
  cor(use = "pairwise.complete.obs")

# P-values for each pair
sig_raw <- purrr::map_dfr(FEAT, \(v1) {
  purrr::map_dfr(FEAT, \(v2) {
    pval <- if (v1 == v2) 0 else
      tryCatch(cor.test(bc[[v1]], bc[[v2]])$p.value, error = function(e) NA_real_)
    data.frame(Var1_raw = v1, Var2_raw = v2, pval = pval)
  })
})

cor_long <- as.data.frame(as.table(cor_mat)) |>
  setNames(c("Var1_raw", "Var2_raw", "r")) |>
  dplyr::left_join(sig_raw, by = c("Var1_raw", "Var2_raw")) |>
  dplyr::mutate(
    Var1      = factor(CORR_LABELS[Var1_raw], levels = CORR_LABELS[FEAT]),
    Var2      = factor(CORR_LABELS[Var2_raw], levels = CORR_LABELS[FEAT]),
    sig_level = dplyr::case_when(
      Var1_raw == Var2_raw        ~ "sig",
      is.na(pval) | pval < 0.05  ~ "sig",
      TRUE                         ~ "ns"),
    sig_stars = dplyr::case_when(
      Var1_raw == Var2_raw ~ "",
      is.na(pval) | pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE         ~ ""),
    label     = ifelse(sig_stars == "", sprintf("%.2f", r),
                       paste0(sprintf("%.2f", r), "\n", sig_stars)),
    dark_text = abs(r) > 0.55 & sig_level == "sig"
  )

fig_a1 <- ggplot(cor_long, aes(Var1, Var2, fill = r, alpha = sig_level)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = label, colour = dark_text), size = 2.9, lineheight = 0.85) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  scale_alpha_manual(values = c(sig = 1.0, ns = 0.18), guide = "none") +
  scale_colour_manual(values = c("FALSE" = "grey25", "TRUE" = "white"), guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Environmental predictor correlations",
       caption = "* p<0.05  ** p<0.01  *** p<0.001  |  grey = not significant") +
  theme_pub() +
  theme(panel.grid   = element_blank(), axis.ticks = element_blank(),
        axis.text.x  = element_text(angle = 45, hjust = 1),
        plot.margin  = margin(4, 4, 6, 4, "mm"))
save_fig(fig_a1, "FigA1_CorrMatrix", 130, 115)

# =============================================================================
# FIGURE S2 — Variable importance: MaxEnt permutation ΔAUC (ecosystem level)
# =============================================================================
vimp_mx <- dplyr::bind_rows(all_vimp) |>
  dplyr::select(Type, Variable, MaxEnt) |>
  dplyr::filter(!is.na(MaxEnt)) |>
  dplyr::mutate(
    Type     = factor(Type, levels = ECO_LVLS),
    Variable = factor(VAR_LABELS[Variable], levels = rev(VAR_LABELS[FEAT]))
  )

fig_s2 <- ggplot(vimp_mx, aes(MaxEnt, Variable, fill = Type)) +
  geom_col(colour = "white", linewidth = 0.2, alpha = 0.85) +
  ggh4x::facet_wrap2(~Type, ncol = 3, strip = eco_strip()) +
  scale_fill_manual(values = ECO_PAL, guide = "none") +
  scale_x_continuous(labels = scales::label_percent()) +
  labs(x = "Importance share (%) — MaxEnt permutation ΔAUC", y = NULL) +
  theme_pub()
save_fig(fig_s2, "FigS2_VarImp", 170, 85)

# =============================================================================
# FIGURE S2b — SHAP beeswarm (only if treeshap succeeded for all species)
# =============================================================================
shap_ok <- exists("all_shap") &&
           all(purrr::map_lgl(ECO_LVLS, \(tp) !is.null(all_shap[[tp]])))
if (shap_ok) {
  shap_bee_plots <- purrr::map(ECO_LVLS, \(tp) {
    shapviz::sv_importance(all_shap[[tp]], kind = "beeswarm",
                            max_display = length(FEAT)) +
      ggplot2::labs(title = tp, x = "SHAP value") +
      ggplot2::scale_x_continuous(n.breaks = 4) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(plot.title   = ggplot2::element_text(face = "bold", size = 10),
                     axis.title.x = ggplot2::element_text(size = 8))
  })
  fig_s2b <- patchwork::wrap_plots(shap_bee_plots, ncol = 3) +
    patchwork::plot_annotation(theme = ggplot2::theme_minimal())
  save_fig(fig_s2b, "FigS2b_SHAP_Beeswarm", 200, 100)
  message("FigS2b saved (SHAP beeswarm)")
} else {
  message("FigS2b skipped: TreeSHAP not available for all species")
}

# =============================================================================
# FIGURE D4 — Occupancy rate vs distance from river mouth
# (moved here: needs bg_df from PART 5)
#
# Primary layer is a binomial GAM fitted to EVERY presence and background
# record, so the curve involves no binning choice at all.  The binned points
# are kept only as a raw-data check on the smooth, and their size carries a
# visible legend giving the number of records behind each point.
# =============================================================================
N_BINS_D4 <- 12   # log-spaced bins for the check points only; the GAM is unbinned

coast_thresh_d4 <- quantile(bc$dist_km, 0.95)
bg_d4 <- bg_df |> dplyr::filter(dist_km <= coast_thresh_d4) |>
  dplyr::slice_sample(n = 5000) |>
  dplyr::select(dist_km) |> dplyr::mutate(presence = 0L)

# presence + background, per ecosystem, unbinned
occ_all <- purrr::map_dfr(ECO_LVLS, \(tp) {
  dplyr::bind_rows(
    bc |> dplyr::filter(Type == tp, dist_km <= coast_thresh_d4) |>
      dplyr::select(dist_km) |> dplyr::mutate(presence = 1L),
    bg_d4
  ) |> dplyr::mutate(Type = tp)
}) |>
  dplyr::filter(is.finite(dist_km)) |>
  dplyr::mutate(Type = factor(Type, levels = ECO_LVLS))

# ── unbinned binomial GAM: occupancy as a smooth function of log distance ────
occ_gam <- purrr::map_dfr(ECO_LVLS, \(tp) {
  d <- dplyr::filter(occ_all, Type == tp)
  m <- tryCatch(
    mgcv::gam(presence ~ s(log10(dist_km + 1), k = 6), family = binomial(), data = d),
    error = function(e) NULL)
  if (is.null(m)) return(NULL)
  nd <- data.frame(dist_km = exp(seq(log(max(0.3, min(d$dist_km))),
                                     log(max(d$dist_km)), length.out = 250)))
  p  <- predict(m, nd, type = "link", se.fit = TRUE)
  nd |> dplyr::mutate(Type = tp,
                      fit = plogis(p$fit),
                      lo  = plogis(p$fit - 1.96 * p$se.fit),
                      hi  = plogis(p$fit + 1.96 * p$se.fit))
}) |> dplyr::mutate(Type = factor(Type, levels = ECO_LVLS))

# ── binned points, retained as a check on the smooth ─────────────────────────
emp_occ <- purrr::map_dfr(ECO_LVLS, \(tp) {
  pts  <- dplyr::filter(occ_all, Type == tp)
  brks <- exp(seq(log(pmax(1, min(pts$dist_km))), log(max(pts$dist_km)),
                  length.out = N_BINS_D4))
  pts |>
    dplyr::mutate(bin_id = findInterval(dist_km, brks, all.inside = TRUE)) |>
    dplyr::group_by(bin_id) |>
    dplyr::summarize(
      occ_rate = mean(presence),
      n        = dplyr::n(),
      n_bg     = sum(1L - presence),
      bin_mid  = exp(mean(log(dist_km))),
      .groups  = "drop"
    ) |>
    dplyr::filter(n >= 10, n_bg >= 3) |>
    dplyr::mutate(Type = tp) |>
    dplyr::select(-n_bg)
}) |> dplyr::mutate(Type = factor(Type, levels = ECO_LVLS))

# Evenly log-spaced breaks: 1,5,10,50,100 crowds 50 next to 100, which is what
# made the tick labels collide in the previous version.
X_BREAKS_D4 <- c(1, 3, 10, 30, 100, 300)
X_LIMS_D4   <- range(occ_all$dist_km[occ_all$dist_km > 0], na.rm = TRUE) * c(0.85, 1.15)

fig_d4 <- ggplot() +
  geom_ribbon(data = occ_gam, aes(dist_km, ymin = lo, ymax = hi, fill = Type),
              alpha = 0.15, colour = NA) +
  geom_line(data = occ_gam, aes(dist_km, fit, colour = Type), linewidth = 0.9) +
  geom_point(data = emp_occ, aes(bin_mid, occ_rate, colour = Type, size = n),
             alpha = 0.75, shape = 16) +
  scale_colour_manual(values = ECO_PAL, name = NULL) +
  scale_fill_manual(  values = ECO_PAL, name = NULL, guide = "none") +
  scale_size_area(max_size = 4.5, name = "Records\nper bin") +
  scale_x_log10(breaks = X_BREAKS_D4, labels = X_BREAKS_D4, limits = X_LIMS_D4) +
  annotation_logticks(sides = "b", linewidth = 0.2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Distance from nearest river mouth (km)", y = "Occupancy rate") +
  theme_pub() +
  theme(legend.position = "top", legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "mm"))
save_fig(fig_d4, "FigD4_EmpiricalOccupancy", 140, 90)

# =============================================================================
# FIGURE ROC — Spatial block CV ROC curves (pooled across 4 longitude folds)
# =============================================================================
roc_df <- purrr::map_dfr(ECO_LVLS, \(tp) {
  purrr::map_dfr(c("MaxEnt", "RF balanced"), \(mdl) {
    pool <- if (mdl == "MaxEnt") all_roc[[tp]]$mx_pool else all_roc[[tp]]$rf_pool
    if (is.null(pool) || nrow(pool) < 20) return(NULL)
    roc_obj <- pROC::roc(pool$y, pool$pred, quiet = TRUE)
    data.frame(
      Type  = tp,
      Model = mdl,
      FPR   = 1 - rev(roc_obj$specificities),
      TPR   = rev(roc_obj$sensitivities),
      AUC   = as.numeric(pROC::auc(roc_obj))
    )
  })
}) |> dplyr::mutate(Type = factor(Type, levels = ECO_LVLS))

auc_ann <- roc_df |>
  dplyr::distinct(Type, Model, AUC) |>
  dplyr::mutate(
    y     = dplyr::if_else(Model == "MaxEnt", 0.18, 0.08),
    label = sprintf("AUC = %.3f", AUC)
  )

fig_roc <- ggplot(roc_df, aes(FPR, TPR, colour = Model, linetype = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey70", linewidth = 0.4) +
  geom_line(linewidth = 0.85) +
  geom_text(data = auc_ann,
            aes(x = 0.58, y = y, label = label, colour = Model),
            size = 2.6, fontface = "bold", show.legend = FALSE) +
  ggh4x::facet_wrap2(~Type, ncol = 3, strip = eco_strip()) +
  scale_colour_manual(values = c(MaxEnt = "#4477AA", "RF balanced" = "#117733"),
                      name = NULL) +
  scale_linetype_manual(values = c(MaxEnt = "solid", "RF balanced" = "dashed"),
                        name = NULL) +
  scale_x_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = 0.02), limits = c(0, 1)) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = 0.02), limits = c(0, 1)) +
  coord_fixed() +
  labs(x = "False positive rate (1 − Specificity)",
       y = "True positive rate (Sensitivity)") +
  theme_pub() + theme(legend.position = "bottom")
save_fig(fig_roc, "FigROC_CrossValidation", 170, 80)

# =============================================================================
# PART 5c — SPECIES-LEVEL MaxEnt: variable importance distribution
# =============================================================================
# Fits one MaxEnt per species (≥MIN_PRES records) within each ecosystem type.
# Variable importance distribution across species shows whether type-level
# conclusions are consistent or driven by a few outliers.
# Uses same coastal-weighted background as PART 5b (bg_coastal still in scope).

MIN_PRES <- 30   # minimum presence records per species
MAX_SP   <- 25   # cap per type to limit runtime

sp_col <- intersect(
  c("species","Species","scientificName","taxon","verbatimScientificName"),
  names(bc)
)
if (length(sp_col) == 0) {
  message("PART 5c SKIPPED: no species column in bc. Columns: ",
          paste(names(bc), collapse = ", "))
} else {
  sp_col <- sp_col[1]
  message("\n=== PART 5c: species-level MaxEnt | column: '", sp_col, "' ===")

  sp_info <- bc |>
    dplyr::rename(.sp = dplyr::all_of(sp_col)) |>
    dplyr::filter(!is.na(.sp), .sp != "") |>
    dplyr::count(Type, .sp) |>
    dplyr::filter(n >= MIN_PRES) |>
    dplyr::group_by(Type) |>
    dplyr::slice_max(n, n = MAX_SP, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::arrange(Type, dplyr::desc(n))

  message(sprintf("Species qualifying (n >= %d, max %d/type): %d",
                  MIN_PRES, MAX_SP, nrow(sp_info)))
  print(as.data.frame(sp_info))

  # All qualifying species for Appendix maps (no per-type cap)
  top3_per_type <- sp_info  # all species with n >= MIN_PRES, ordered by n desc within type

  # Panel order for the per-species figures: Seagrass | Mangrove | Kelp
  # (row = rank within type).  FigSup_SpeciesMap itself is drawn after the
  # species-model loop below, because it annotates each panel with that
  # species' AUC and sp_vimp does not exist yet at this point.
  top3_reordered_occ <- top3_per_type |>
    dplyr::group_by(Type) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::arrange(rank, factor(Type, levels = ECO_LVLS))

  # Per-ecosystem-type occurrence map (matches FigMain_GlobalSuitability layout)
  occ_by_type <- bc |>
    dplyr::filter(!is.na(lon), !is.na(lat)) |>
    dplyr::mutate(Type = factor(Type, levels = ECO_LVLS))
  fig_occ_type <- ggplot() +
    geom_sf(data = world, fill = "grey85", colour = "white", linewidth = 0.1) +
    geom_point(data = occ_by_type, aes(lon, lat, colour = Type),
               size = 0.4, alpha = 0.6, shape = 16) +
    ggh4x::facet_wrap2(~Type, ncol = 1, strip = eco_strip()) +
    scale_colour_manual(values = ECO_PAL, guide = "none") +
    coord_sf(ylim = c(-65, 75), expand = FALSE) +
    theme_pub() +
    theme(axis.title = element_blank(),
          axis.text  = element_text(size = 7),
          plot.margin = margin(1, 1, 2, 1, "mm"))
  save_fig(fig_occ_type, "FigSup_OccurrenceByType", 89, 145)

  # FEAT_SP = FEAT (SST/lat controlled via background matching, not excluded from FEAT)

  sp_vimp_list <- list()
  sp_models    <- list()  # store fitted models for prediction maps
  sp_roc_list  <- list()  # MaxEnt per-species training ROC
  for (i in seq_len(nrow(sp_info))) {
    tp  <- sp_info$Type[i]
    spp <- sp_info$.sp[i]
    n_s <- sp_info$n[i]

    pres_s <- bc |>
      dplyr::rename(.sp = dplyr::all_of(sp_col)) |>
      dplyr::filter(Type == tp, .sp == spp)

    # SST + lat matched background: same climate zone and latitude band as this species
    sst_qs <- quantile(pres_s$SST, c(0.05, 0.95), na.rm = TRUE)
    lat_qs <- quantile(pres_s$lat, c(0.05, 0.95), na.rm = TRUE)
    bg_pool_s <- bg_coastal |> dplyr::filter(SST >= sst_qs[1], SST <= sst_qs[2],
                                              lat >= lat_qs[1], lat <= lat_qs[2])
    if (nrow(bg_pool_s) < 50) bg_pool_s <- bg_coastal
    n_bg_s <- min(nrow(pres_s) * 10, nrow(bg_pool_s), 3000)
    bg_s   <- bg_pool_s |> dplyr::slice_sample(n = n_bg_s)

    df_s <- dplyr::bind_rows(
      pres_s |> dplyr::select(dplyr::all_of(FEAT_SP)) |> dplyr::mutate(presence = 1L),
      bg_s   |> dplyr::select(dplyr::all_of(FEAT_SP)) |> dplyr::mutate(presence = 0L)
    ) |>
      dplyr::mutate(dplyr::across(
        dplyr::all_of(FEAT_SP),
        \(x) dplyr::if_else(is.na(x), median(x, na.rm = TRUE), x)
      ))

    X_s <- as.matrix(df_s[, FEAT_SP]); colnames(X_s) <- FEAT_SP
    y_s <- df_s$presence

    # Remove zero-variance features (e.g. log_tidal=0 when ocean_tides.zip absent)
    feat_var_s <- apply(X_s, 2, var)
    FEAT_sp_mx <- names(feat_var_s)[!is.na(feat_var_s) & feat_var_s > 1e-10]
    X_s_mx     <- X_s[, FEAT_sp_mx, drop = FALSE]
    if (length(FEAT_sp_mx) < length(FEAT_SP))
      message("  Zero-var removed (sp): ", paste(setdiff(FEAT_SP, FEAT_sp_mx), collapse = ", "))

    mx_s <- tryCatch(
      maxnet::maxnet(y_s, as.data.frame(X_s_mx),
                     maxnet::maxnet.formula(y_s, as.data.frame(X_s_mx), classes = "lqph"),
                     regmult = 2.0),
      error = function(e) { message("  SKIP ", spp, ": ", e$message); NULL }
    )
    if (is.null(mx_s)) next

    bp_s <- predict(mx_s, as.data.frame(X_s), type = "cloglog") |> as.numeric()
    ba_s <- tryCatch(pROC::auc(y_s, bp_s, quiet = TRUE) |> as.numeric(),
                     error = \(e) NA_real_)
    vi_s <- sapply(FEAT_SP, \(fv) {
      Xp <- as.data.frame(X_s); Xp[[fv]] <- sample(Xp[[fv]])
      pa <- tryCatch(
        pROC::auc(y_s, predict(mx_s, Xp, type = "cloglog") |> as.numeric(),
                  quiet = TRUE) |> as.numeric(),
        error = \(e) ba_s
      )
      max(0, ba_s - pa)
    })
    vi_tot <- sum(vi_s, na.rm = TRUE)
    if (is.na(vi_tot) || vi_tot == 0) next

    key <- paste(tp, spp, sep = "__")
    sp_vimp_list[[key]] <- tibble::tibble(
      Type      = tp,
      Species   = spp,
      n_pres    = n_s,
      AUC       = ba_s,
      Variable  = FEAT_SP,
      Importance = vi_s / vi_tot
    )
    sp_models[[key]]   <- mx_s  # save for prediction maps
    sp_roc_list[[key]] <- data.frame(Type = tp, Species = spp, y = y_s, pred = bp_s)
    message(sprintf("  [%d/%d] %-12s | %-40s AUC=%.3f n=%d",
                    i, nrow(sp_info), tp, spp, ba_s, n_s))
  }

  if (length(sp_vimp_list) == 0) {
    message("PART 5c: no models succeeded — check MIN_PRES or species column")
  } else {
    sp_vimp_raw <- dplyr::bind_rows(sp_vimp_list)  # keep raw variable names for PART 5d
    sp_vimp <- sp_vimp_raw |>
      dplyr::mutate(
        Type     = factor(Type, levels = ECO_LVLS),
        Variable = factor(VAR_LABELS[Variable], levels = rev(VAR_LABELS[FEAT_SP]))
      )

    # --- FigSup_SpeciesMap: occurrence records per species, with training AUC ---
    # Drawn here (not before the loop) so each panel can be annotated with the
    # AUC of that species' own MaxEnt model, matching FigSup_SpeciesPred.
    sp_map_panels <- purrr::map(seq_len(nrow(top3_reordered_occ)), \(i) {
      tp       <- top3_reordered_occ$Type[i]
      spp      <- top3_reordered_occ$.sp[i]
      type_col <- ECO_PAL[[tp]]
      pts <- bc |>
        dplyr::rename(.sp = dplyr::all_of(sp_col)) |>
        dplyr::filter(.sp == spp)
      auc_val   <- dplyr::filter(sp_vimp_raw, Species == spp)$AUC
      auc_label <- if (length(auc_val) > 0) sprintf("AUC = %.3f", auc_val[1]) else ""
      ggplot() +
        geom_sf(data = world, fill = "grey85", colour = "white", linewidth = 0.1) +
        geom_point(data = pts, aes(lon, lat), colour = type_col,
                   size = 0.35, alpha = 0.7, shape = 16) +
        annotate("text", x = Inf, y = -Inf, label = auc_label,
                 hjust = 1.05, vjust = -0.4, size = 2.0, colour = "grey35") +
        coord_sf(ylim = c(-65, 75), expand = FALSE) +
        labs(title = paste("●", spp)) +
        theme_pub() +
        theme(axis.title = element_blank(),
              axis.text  = element_blank(),
              axis.ticks = element_blank(),
              plot.title = element_text(size = 7, face = "italic",
                                        colour = type_col,
                                        margin = margin(b = 1)))
    })
    n_occ_rows  <- ceiling(nrow(top3_reordered_occ) / 3)
    fig_ssp_map <- patchwork::wrap_plots(sp_map_panels, ncol = 3)
    save_fig(fig_ssp_map, "FigSup_SpeciesMap", 170, n_occ_rows * 50 + 10)

    # --- FigS3: partial dependence with ±SD across species, all FEAT_SP variables ---
    if (length(sp_models) > 0) {
      message("Computing species-level partial dependence for FigS3...")

      bg_med <- sapply(FEAT_SP, \(v) median(bc[[v]], na.rm = TRUE))

      var_grids <- list(
        Salinity     = seq(quantile(bc$Salinity, 0.02, na.rm = TRUE),
                           quantile(bc$Salinity, 0.98, na.rm = TRUE), length.out = 50),
        Kd490        = exp(seq(log(max(1e-4, quantile(bc$Kd490[bc$Kd490 > 0],
                                                       0.02, na.rm = TRUE))),
                               log(quantile(bc$Kd490, 0.98, na.rm = TRUE)),
                               length.out = 50)),
        dist_seq      = km_to_dist(seq(0, quantile(bc$dist_km, 0.90, na.rm = TRUE),
                                       length.out = 50)),
        log_discharge = seq(quantile(bc$log_discharge, 0.02, na.rm = TRUE),
                            quantile(bc$log_discharge, 0.98, na.rm = TRUE), length.out = 50)
      )
      names(var_grids)[names(var_grids) == "dist_seq"] <- DIST_VAR

      sp_pd_raw <- purrr::map_dfr(names(sp_models), \(key) {
        parts <- strsplit(key, "__")[[1]]
        tp    <- parts[1]; spp <- parts[2]
        mod   <- sp_models[[key]]
        purrr::map_dfr(FEAT_SP, \(vname) {
          xs   <- var_grids[[vname]]
          X_pd <- matrix(rep(bg_med, length(xs)), nrow = length(xs), byrow = TRUE)
          colnames(X_pd) <- FEAT_SP
          X_pd[, vname]  <- xs
          suit <- predict(mod, as.data.frame(X_pd), type = "cloglog") |> as.numeric()
          tibble::tibble(Type = tp, Species = spp, Variable = vname, x = xs, suit = suit)
        })
      })

      sp_pd_summary <- sp_pd_raw |>
        dplyr::group_by(Type, Variable, x) |>
        dplyr::summarise(suit_mean = mean(suit, na.rm = TRUE),
                         suit_sd   = sd(suit,   na.rm = TRUE), .groups = "drop") |>
        dplyr::mutate(suit_lo = pmax(0, suit_mean - suit_sd),
                      suit_hi = pmin(1, suit_mean + suit_sd),
                      Type    = factor(Type, levels = ECO_LVLS))

      make_pd_panel <- function(vname, x_label, log_x = FALSE, show_y = TRUE,
                                back_log10 = FALSE, nice_pow10 = NULL) {
        df <- sp_pd_summary |> dplyr::filter(Variable == vname)
        p  <- ggplot(df, aes(x, suit_mean, colour = Type, fill = Type)) +
          geom_ribbon(aes(ymin = suit_lo, ymax = suit_hi), alpha = 0.20, colour = NA) +
          geom_line(linewidth = 0.9) +
          scale_colour_manual(values = ECO_PAL, name = NULL) +
          scale_fill_manual(values   = ECO_PAL, name = NULL) +
          scale_y_continuous(labels  = scales::label_percent(), limits = c(0, NA),
                             expand  = expansion(mult = c(0, 0.05))) +
          labs(x = x_label, y = if (show_y) "Predicted suitability" else NULL) +
          theme_pub()
        if (log_x) p <- p + scale_x_log10(labels = scales::label_number())
        if (back_log10 && !is.null(nice_pow10)) {
          p <- p + scale_x_continuous(
            breaks = log10(nice_pow10 + 1),
            labels = scales::label_number(big.mark = " ")(nice_pow10))
        } else if (back_log10) {
          p <- p + scale_x_continuous(
            labels = \(v) scales::label_number()(round(10^v - 1)))
        }
        p
      }

      p_sal  <- make_pd_panel("Salinity",     "Salinity (PSU)",            log_x = FALSE, show_y = TRUE)
      p_kd   <- make_pd_panel("Kd490",        "Kd490 (m⁻¹)",              log_x = TRUE,  show_y = FALSE)
      p_dist <- make_pd_panel(DIST_VAR,        "Distance from river mouth (km)",  log_x = FALSE, show_y = TRUE, back_log10 = (DIST_VAR == "log_dist"))
      p_up    <- make_pd_panel("log_discharge", "River discharge (m³ s⁻¹)", log_x = FALSE, show_y = FALSE, back_log10 = TRUE,
                              nice_pow10 = c(1, 10, 100, 1000, 10000))
      fig_s3 <- ((p_sal | p_kd) / (p_dist | p_up)) +
        patchwork::plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
      save_fig(fig_s3, "FigS3_ResponseCurves_CI", 130, 110)
      message("FigS3_ResponseCurves_CI saved")
    }
  }
}  # end PART 5c

# =============================================================================
# FigROC — MaxEnt: per-species training thin lines + ecosystem CV 95% CI ribbon
# =============================================================================
if (exists("sp_roc_list") && length(sp_roc_list) > 0) {
  sp_roc_df <- purrr::map_dfr(names(sp_roc_list), \(key) {
    d <- sp_roc_list[[key]]
    if (length(unique(d$y)) < 2 || nrow(d) < 20) return(NULL)
    roc_s <- tryCatch(pROC::roc(d$y, d$pred, quiet = TRUE), error = function(e) NULL)
    if (is.null(roc_s)) return(NULL)
    data.frame(Type = d$Type[1], Species = d$Species[1],
               FPR  = 1 - rev(roc_s$specificities),
               TPR  = rev(roc_s$sensitivities))
  }) |> dplyr::mutate(Type = factor(Type, levels = ECO_LVLS))

  # Per-fold ROC CI ribbon (95% CI across 4 spatial block CV folds)
  FPR_GRID <- seq(0, 1, length.out = 300)
  roc_ribbon <- purrr::map_dfr(ECO_LVLS, \(tp) {
    folds <- all_roc[[tp]]$mx_folds
    if (is.null(folds) || length(folds) == 0) return(NULL)
    purrr::map_dfr(seq_along(folds), \(b) {
      d <- folds[[b]]
      if (is.null(d) || length(unique(d$y)) < 2) return(NULL)
      roc_b <- tryCatch(pROC::roc(d$y, d$pred, quiet = TRUE), error = function(e) NULL)
      if (is.null(roc_b)) return(NULL)
      tpr_b <- approx(1 - rev(roc_b$specificities), rev(roc_b$sensitivities),
                      xout = FPR_GRID, yleft = 0, yright = 1)$y
      data.frame(Type = tp, fold = b, FPR = FPR_GRID, TPR = tpr_b)
    }) |>
      dplyr::group_by(FPR) |>
      dplyr::summarise(tpr_mean = mean(TPR, na.rm = TRUE),
                       tpr_sd   = sd(TPR,   na.rm = TRUE),
                       n_fold   = dplyr::n(), .groups = "drop") |>
      dplyr::mutate(Type   = tp,
                    tpr_lo = pmax(0, tpr_mean - 1.96 * tpr_sd / sqrt(n_fold)),
                    tpr_hi = pmin(1, tpr_mean + 1.96 * tpr_sd / sqrt(n_fold)))
  }) |> dplyr::mutate(Type = factor(Type, levels = ECO_LVLS))

  # AUC mean ± SD across folds
  auc_cv_summ <- purrr::map_dfr(ECO_LVLS, \(tp) {
    vals <- all_roc[[tp]]$mx_cv_auc
    if (is.null(vals) || length(vals) < 1) {
      fallback <- auc_ann |> dplyr::filter(Type == tp, Model == "MaxEnt") |> dplyr::pull(AUC)
      return(data.frame(Type = tp, label = sprintf("AUC = %.3f", fallback[1])))
    }
    data.frame(Type = tp, label = sprintf("AUC = %.3f ± %.3f", mean(vals), sd(vals)))
  }) |> dplyr::mutate(Type = factor(Type, levels = ECO_LVLS))

  fig_roc2 <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey70", linewidth = 0.4) +
    geom_line(data = sp_roc_df,
              aes(FPR, TPR, group = Species, colour = Type),
              linewidth = 0.28, alpha = 0.32) +
    geom_ribbon(data = roc_ribbon,
                aes(FPR, ymin = tpr_lo, ymax = tpr_hi, fill = Type),
                alpha = 0.35, colour = NA) +
    geom_line(data = roc_ribbon,
              aes(FPR, tpr_mean, colour = Type), linewidth = 0.7) +
    geom_text(data = auc_cv_summ,
              aes(x = 0.30, y = 0.10, label = label, colour = Type),
              hjust = 0, size = 2.4, fontface = "bold") +
    scale_colour_manual(values = ECO_PAL, name = NULL, guide = "none") +
    scale_fill_manual(values   = ECO_PAL, name = NULL, guide = "none") +
    ggh4x::facet_wrap2(~Type, ncol = 3, strip = eco_strip()) +
    scale_x_continuous(labels = scales::label_number(scale = 100),
                       breaks = c(0.25, 0.5, 0.75, 1.0),
                       expand = expansion(mult = 0.02), limits = c(0, 1)) +
    scale_y_continuous(labels = scales::label_number(scale = 100),
                       breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                       expand = expansion(mult = 0.02), limits = c(0, 1)) +
    coord_fixed() +
    labs(x = "False positive rate (%)",
         y = "True positive rate (%)",
         caption = "Thin lines: per-species training ROC  |  Ribbon: ecosystem MaxEnt CV (95% CI, 4-fold spatial block)") +
    theme_pub() + theme(plot.caption = element_text(size = 7, colour = "grey50"))
  save_fig(fig_roc2, "FigROC_CrossValidation", 170, 80)
  message("FigROC saved — MaxEnt CI ribbon with per-species training ROC")
}

# =============================================================================
# PART 5d — DIAGNOSTIC: unthinned vs spatially thinned presence comparison
# =============================================================================
if (exists("sp_vimp_raw") && !is.null(sp_vimp_raw) && exists("sp_col") && exists("bg_coastal")) {
  message("\n=== PART 5d: thinning comparison ===")

  # ── 1. Spatial thinning of presences (1° × 1° grid, 1 record/cell/species) ──
  THIN_DEG <- 1.0
  bc_thin <- bc |>
    dplyr::rename(.sp_tmp = dplyr::all_of(sp_col)) |>
    dplyr::mutate(cell_lon = floor(lon / THIN_DEG),
                  cell_lat = floor(lat / THIN_DEG)) |>
    dplyr::group_by(Type, .sp_tmp, cell_lon, cell_lat) |>
    dplyr::slice_sample(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(-cell_lon, -cell_lat) |>
    dplyr::rename(!!sp_col := .sp_tmp)

  sp_info_thin <- bc_thin |>
    dplyr::rename(.sp = dplyr::all_of(sp_col)) |>
    dplyr::count(Type, .sp) |>
    dplyr::filter(n >= 10) |>
    dplyr::arrange(Type, dplyr::desc(n))

  message(sprintf("  Thinned: %d → %d records | qualifying spp (n>=10): %d",
                  nrow(bc), nrow(bc_thin), nrow(sp_info_thin)))

  # ── 2. Species loop with thinned presences (MaxEnt importance only) ──────────
  sp_vimp_thin_list <- list()

  for (i in seq_len(nrow(sp_info_thin))) {
    tp  <- sp_info_thin$Type[i]
    spp <- sp_info_thin$.sp[i]

    pres_s <- bc_thin |>
      dplyr::rename(.sp = dplyr::all_of(sp_col)) |>
      dplyr::filter(Type == tp, .sp == spp)

    sst_qs <- quantile(pres_s$SST, c(0.05, 0.95), na.rm = TRUE)
    lat_qs <- quantile(pres_s$lat, c(0.05, 0.95), na.rm = TRUE)
    bg_pool_s <- bg_coastal |> dplyr::filter(SST >= sst_qs[1], SST <= sst_qs[2],
                                              lat >= lat_qs[1], lat <= lat_qs[2])
    if (nrow(bg_pool_s) < 50) bg_pool_s <- bg_coastal
    n_bg_s <- min(nrow(pres_s) * 10, nrow(bg_pool_s), 3000)
    bg_s   <- bg_pool_s |> dplyr::slice_sample(n = n_bg_s)

    df_s <- dplyr::bind_rows(
      pres_s |> dplyr::select(dplyr::all_of(FEAT_SP)) |> dplyr::mutate(presence = 1L),
      bg_s   |> dplyr::select(dplyr::all_of(FEAT_SP)) |> dplyr::mutate(presence = 0L)
    ) |>
      dplyr::mutate(dplyr::across(
        dplyr::all_of(FEAT_SP),
        \(x) dplyr::if_else(is.na(x), median(x, na.rm = TRUE), x)
      ))

    X_s <- as.matrix(df_s[, FEAT_SP]); colnames(X_s) <- FEAT_SP
    y_s <- df_s$presence

    feat_var_s2 <- apply(X_s, 2, var)
    FEAT_sp_mx2 <- names(feat_var_s2)[!is.na(feat_var_s2) & feat_var_s2 > 1e-10]
    X_s_mx2     <- X_s[, FEAT_sp_mx2, drop = FALSE]

    mx_s <- tryCatch(
      maxnet::maxnet(y_s, as.data.frame(X_s_mx2),
                     maxnet::maxnet.formula(y_s, as.data.frame(X_s_mx2), classes = "lqph"),
                     regmult = 2.0),
      error = function(e) NULL
    )
    if (is.null(mx_s)) next

    base_auc_s <- tryCatch(
      pROC::auc(y_s, predict(mx_s, as.data.frame(X_s), type = "cloglog") |> as.numeric(),
                quiet = TRUE) |> as.numeric(),
      error = function(e) 0.5
    )

    vi_s <- sapply(FEAT_SP, \(fv) {
      Xp <- as.data.frame(X_s); Xp[[fv]] <- sample(Xp[[fv]])
      pa <- tryCatch(
        pROC::auc(y_s, predict(mx_s, Xp, type = "cloglog") |> as.numeric(),
                  quiet = TRUE) |> as.numeric(),
        error = function(e) base_auc_s
      )
      max(0, base_auc_s - pa)
    })
    vi_s <- vi_s / max(sum(vi_s), 1e-9)

    sp_vimp_thin_list[[paste(tp, spp, sep = "__")]] <- tibble::tibble(
      Type = tp, Species = spp, Variable = FEAT_SP, Importance = vi_s
    )
    message(sprintf("  [thinned %d/%d] %-12s | %-35s n=%d",
                    i, nrow(sp_info_thin), tp, spp, nrow(pres_s)))
  }

  sp_vimp_thinned <- dplyr::bind_rows(sp_vimp_thin_list)

  # ── 3. Comparison figure ──────────────────────────────────────────────────────
  if (nrow(sp_vimp_thinned) > 0) {
    sp_vimp_compare <- dplyr::bind_rows(
      sp_vimp_raw     |> dplyr::mutate(Method = "Unthinned"),
      sp_vimp_thinned |> dplyr::mutate(Method = sprintf("Thinned (%.0f°)", THIN_DEG))
    ) |>
      dplyr::mutate(
        Type     = factor(Type, levels = ECO_LVLS),
        Method   = factor(Method, levels = c("Unthinned", sprintf("Thinned (%.0f°)", THIN_DEG))),
        Variable = factor(VAR_LABELS[Variable], levels = rev(VAR_LABELS[FEAT_SP]))
      )

    fig_compare <- ggplot(sp_vimp_compare,
                          aes(Variable, Importance * 100, fill = Type, colour = Type)) +
      geom_violin(alpha = 0.35, scale = "width", trim = TRUE, linewidth = 0.3) +
      geom_boxplot(width = 0.18, outlier.size = 0.8, linewidth = 0.35,
                   fill = "white", colour = "grey20", alpha = 0.85) +
      geom_jitter(width = 0.07, size = 0.7, alpha = 0.55) +
      scale_fill_manual(values   = ECO_PAL) +
      scale_colour_manual(values = ECO_PAL) +
      scale_y_continuous(labels = scales::label_percent(scale = 1)) +
      facet_grid(Type ~ Method, scales = "free_x") +
      coord_flip() +
      labs(x = NULL, y = "Importance share (%)") +
      theme_pub() +
      theme(legend.position = "none",
            strip.text.x = element_text(face = "bold", size = 8),
            strip.text.y = element_text(face = "bold", size = 8),
            panel.grid.major.y = element_blank())

    save_fig(fig_compare, "FigDiag_ThinningComparison", 170, 120)
    message("FigDiag_ThinningComparison saved")
  }
}  # end PART 5d

# FigSp1 — thinned per-species MaxEnt importance (boxplot + jitter, no violin)
if (exists("sp_vimp_thinned") && !is.null(sp_vimp_thinned) && nrow(sp_vimp_thinned) > 0) {
  sp_vimp_t_plot <- sp_vimp_thinned |>
    dplyr::mutate(
      Type     = factor(Type, levels = ECO_LVLS),
      Variable = factor(VAR_LABELS[Variable], levels = rev(VAR_LABELS[FEAT_SP]))
    )
  fig_sp1 <- ggplot(sp_vimp_t_plot,
                    aes(Variable, Importance, fill = Type, colour = Type)) +
    geom_boxplot(width = 0.45, fill = "white", colour = "grey25",
                 outlier.shape = NA, linewidth = 0.35) +
    geom_jitter(width = 0.20, size = 1.1, alpha = 0.65, shape = 16) +
    ggh4x::facet_wrap2(~Type, ncol = 3, strip = eco_strip()) +
    scale_fill_manual(values   = ECO_PAL, guide = "none") +
    scale_colour_manual(values = ECO_PAL, guide = "none") +
    scale_y_continuous(labels  = scales::label_number(scale = 100),
                       breaks  = c(0.25, 0.5, 0.75, 1.0),
                       expand  = expansion(mult = c(0.02, 0.08))) +
    labs(x = NULL, y = "Importance share (%)") +
    coord_flip() +
    theme_pub()
  save_fig(fig_sp1, "FigSp1_SpeciesImportance", 170, 95)
  message(sprintf("FigSp1 saved (thinned) — %d species across %d ecosystem types",
                  dplyr::n_distinct(sp_vimp_thinned$Species),
                  dplyr::n_distinct(sp_vimp_thinned$Type)))
}

# ── Checkpoint 2: save species-level results ─────────────────────────────────
if (!exists("RDS_DIR")) { RDS_DIR <- file.path(OUT_DIR, "rds"); dir.create(RDS_DIR, showWarnings = FALSE) }
if (exists("sp_models")     && !is.null(sp_models))     saveRDS(sp_models,     file.path(RDS_DIR, "sp_models.rds"))
if (exists("sp_vimp")         && !is.null(sp_vimp))         saveRDS(sp_vimp,         file.path(RDS_DIR, "sp_vimp.rds"))
if (exists("sp_vimp_thinned") && !is.null(sp_vimp_thinned)) saveRDS(sp_vimp_thinned, file.path(RDS_DIR, "sp_vimp_thinned.rds"))
if (exists("sp_roc_list") && !is.null(sp_roc_list)) saveRDS(sp_roc_list, file.path(RDS_DIR, "sp_roc_list.rds"))
if (exists("top3_per_type") && !is.null(top3_per_type)) saveRDS(top3_per_type, file.path(RDS_DIR, "top3_per_type.rds"))
if (exists("sp_col"))                                    saveRDS(sp_col,        file.path(RDS_DIR, "sp_col.rds"))
if (exists("sp_pd_summary") && !is.null(sp_pd_summary)) saveRDS(sp_pd_summary, file.path(RDS_DIR, "sp_pd_summary.rds"))
message("Checkpoint 2 saved → ", RDS_DIR)

# =============================================================================
# PART 6 — GLOBAL PREDICTION MAPS
# =============================================================================

# Cache name is versioned by the distance definition: a world_grid.rds built
# with the old nearest-channel metric must NOT be silently reused now that
# dist_km means distance to the nearest river mouth.
GRID_CACHE <- file.path(OUT_DIR, "world_grid_overwater.rds")
for (.old in c("world_grid.rds", "world_grid_mouth.rds"))
  if (!file.exists(GRID_CACHE) && file.exists(file.path(OUT_DIR, .old)))
    message("NOTE: output/", .old, " uses an earlier distance definition and is\n",
            "      ignored. Rebuilding the grid with over-water distances (~15 min).")

if (file.exists(GRID_CACHE)) {
  message("Loading cached global grid...")
  world_grid <- readRDS(GRID_CACHE)
  if (!"depth" %in% names(world_grid)) {
    message("  Adding depth to cached world_grid...")
    wg_vect <- terra::vect(as.matrix(world_grid[, c("lon","lat")]), crs = "EPSG:4326")
    wg_depth_raw <- terra::extract(depth_rast, wg_vect)[, 2]
    world_grid$depth <- pmax(0, replace(abs(wg_depth_raw), !is.finite(wg_depth_raw), 0))
    saveRDS(world_grid, GRID_CACHE)
  }
  if (!"upstream_km2" %in% names(world_grid)) {
    message("  Adding upstream_km2 to cached world_grid (one-time)...")
    wg_sf_tmp <- sf::st_as_sf(world_grid[, c("lon","lat")], coords = c("lon","lat"), crs = 4326)
    nn_g_up   <- sf::st_nearest_feature(wg_sf_tmp, river_mouths)
    world_grid$upstream_km2  <- river_mouths$UPLAND_SKM[nn_g_up]
    world_grid$discharge_cms <- river_mouths$DIS_AV_CMS[nn_g_up]
    saveRDS(world_grid, GRID_CACHE)
  }
  if (!"discharge_cms" %in% names(world_grid)) {
    message("  Adding discharge_cms to cached world_grid (one-time)...")
    wg_sf_tmp <- sf::st_as_sf(world_grid[, c("lon","lat")], coords = c("lon","lat"), crs = 4326)
    nn_g_dis  <- sf::st_nearest_feature(wg_sf_tmp, river_mouths)
    world_grid$discharge_cms <- river_mouths$DIS_AV_CMS[nn_g_dis]
    saveRDS(world_grid, GRID_CACHE)
  }
  if (!"dist_city_km" %in% names(world_grid) && !is.null(cities_sf)) {
    message("  Adding dist_city_km to cached world_grid (one-time)...")
    wg_sf_tmp <- sf::st_as_sf(world_grid[, c("lon","lat")], coords = c("lon","lat"), crs = 4326)
    nn_g_city <- sf::st_nearest_feature(wg_sf_tmp, cities_sf)
    world_grid$dist_city_km <- as.numeric(
      sf::st_distance(wg_sf_tmp, cities_sf[nn_g_city, ], by_element = TRUE)) / 1000
    saveRDS(world_grid, GRID_CACHE)
  }
  if (!"dist_city_km" %in% names(world_grid)) world_grid$dist_city_km <- 0  # no cities_sf
  if (!"tidal_m" %in% names(world_grid)) {
    message("  Adding tidal_m to cached world_grid (one-time)...")
    wg_vect_tmp  <- terra::vect(as.matrix(world_grid[, c("lon","lat")]), crs = "EPSG:4326")
    if (!is.null(tidal_rast)) {
      wg_tidal_tmp       <- terra::extract(tidal_rast, wg_vect_tmp)[, 2]
      world_grid$tidal_m <- pmax(0, replace(wg_tidal_tmp, !is.finite(wg_tidal_tmp), 0)) / TIDAL_SCALE
    } else {
      world_grid$tidal_m <- 0
    }
    saveRDS(world_grid, GRID_CACHE)
  }
} else {
  message("Building global grid at 0.5° resolution...")
  env_05 <- terra::aggregate(env_rast,
    fact = round(0.5 / terra::res(env_rast)[1]), fun = "mean", na.rm = TRUE)

  world_grid <- as.data.frame(env_05, xy = TRUE, na.rm = TRUE) |>
    setNames(c("lon","lat", names(env_codes))) |>
    dplyr::filter(!is.na(SST), !is.na(Salinity), !is.na(Kd490), Kd490 > 0) |>
    dplyr::mutate(abs_lat = abs(lat))

  message("  Computing river-mouth distances for grid (~15 min)...")
  grid_sf <- sf::st_as_sf(world_grid[, c("lon","lat")],
                            coords = c("lon","lat"), crs = 4326)
  .gd     <- dist_to_mouth(world_grid$lon, world_grid$lat)
  nn_g    <- attr(.gd, "nearest_idx")
  message(sprintf("  grid cells with a land-blocked straight line: %d of %d (%.1f%%)",
                  attr(.gd, "n_blocked"), length(.gd),
                  100 * attr(.gd, "n_blocked") / length(.gd)))
  world_grid$dist_km <- as.numeric(.gd)
  world_grid$upstream_km2  <- river_mouths$UPLAND_SKM[nn_g]
  world_grid$discharge_cms <- river_mouths$DIS_AV_CMS[nn_g]
  if (!is.null(cities_sf)) {
    message("  Computing city distances for grid...")
    nn_g_city <- sf::st_nearest_feature(grid_sf, cities_sf)
    world_grid$dist_city_km <- as.numeric(
      sf::st_distance(grid_sf, cities_sf[nn_g_city, ], by_element = TRUE)) / 1000
  } else {
    world_grid$dist_city_km <- 0  # no cities_sf: constant 0 → MaxEnt drops
  }
  world_grid <- dplyr::filter(world_grid, is.finite(dist_km))
  wg_vect <- terra::vect(as.matrix(world_grid[, c("lon","lat")]), crs = "EPSG:4326")
  wg_depth_raw <- terra::extract(depth_rast, wg_vect)[, 2]
  world_grid$depth <- pmax(0, replace(abs(wg_depth_raw), !is.finite(wg_depth_raw), 0))
  if (!is.null(tidal_rast)) {
    wg_tidal_raw       <- terra::extract(tidal_rast, wg_vect)[, 2]
    world_grid$tidal_m <- pmax(0, replace(wg_tidal_raw, !is.finite(wg_tidal_raw), 0)) / TIDAL_SCALE
  } else {
    world_grid$tidal_m <- 0
  }
  saveRDS(world_grid, GRID_CACHE)
}

world_grid <- world_grid |> dplyr::filter(depth <= MAX_DEPTH_M)
message(sprintf("  world_grid after depth filter (<=%.0f m): %d cells", MAX_DEPTH_M, nrow(world_grid)))
world_grid               <- add_dist_forms(world_grid)
world_grid$log_depth     <- log10(world_grid$depth        + 1)
world_grid$log_upstream  <- log10(world_grid$upstream_km2 + 1)
world_grid$log_discharge <- log10(pmax(world_grid$discharge_cms, 0) + 1)
world_grid$log_city     <- log10(world_grid$dist_city_km + 1)
world_grid$log_tidal    <- log10(world_grid$tidal_m      + 1)
world_grid$dist_km_baseline <- world_grid$dist_km  # for filtering scenario figures
Xg     <- as.matrix(world_grid[, FEAT]);       colnames(Xg)     <- FEAT
Xg_rng <- as.matrix(world_grid[, FEAT_range]); colnames(Xg_rng) <- FEAT_range
for (tp in ECO_LVLS) {
  suit_hab <- predict(all_mx[[tp]], as.data.frame(Xg), type = "cloglog") |> as.numeric()
  suit_rng <- if (!is.null(range_mx[[tp]]))
    predict(range_mx[[tp]], as.data.frame(Xg_rng), type = "cloglog") |> as.numeric()
  else rep(1, nrow(world_grid))
  world_grid[[paste0("suit_range_", tp)]] <- suit_rng          # Stage 1 only
  world_grid[[paste0("suit_",       tp)]] <- suit_hab * suit_rng  # Stage1 × Stage2
}

wg_coast <- world_grid |> dplyr::filter(dist_km_baseline <= coast_thresh)

# Stage 1: biogeographic range suitability (SST + lat)
suit_range_long <- wg_coast |>
  dplyr::select(lon, lat, starts_with("suit_range_")) |>
  tidyr::pivot_longer(starts_with("suit_range_"), names_to = "Type", values_to = "Suitability") |>
  dplyr::mutate(Type = factor(stringr::str_remove(Type, "suit_range_"), levels = ECO_LVLS))

# Stage 1 training AUC, annotated per panel
range_auc_lab <- purrr::map_dfr(ECO_LVLS, \(tp) {
  a <- if (exists("range_auc")) range_auc[[tp]] else NULL
  if (is.null(a) || is.na(a)) return(NULL)
  data.frame(Type = tp, label = sprintf("AUC = %.3f", a))
}) |> dplyr::mutate(Type = factor(Type, levels = ECO_LVLS))

fig_s4_range <- ggplot() +
  geom_sf(data = world, fill = "grey85", colour = "white", linewidth = 0.1) +
  geom_tile(data = suit_range_long, aes(lon, lat, fill = Suitability), alpha = 0.9) +
  geom_sf(data = world, fill = NA, colour = "white", linewidth = 0.1) +
  { if (nrow(range_auc_lab) > 0)
      geom_label(data = range_auc_lab, aes(x = -178, y = -56, label = label),
                 hjust = 0, vjust = 0, size = 2.1, colour = "grey20",
                 fill = alpha("white", 0.78), label.size = 0,
                 label.padding = unit(0.8, "mm"), inherit.aes = FALSE)
    else NULL } +
  ggh4x::facet_wrap2(~Type, ncol = 1, strip = eco_strip()) +
  scale_fill_viridis_c(option = "plasma", direction = -1, name = "Range\nsuitability",
                        limits = c(0, 1), labels = label_percent()) +
  coord_sf(ylim = c(-65, 75), expand = FALSE) +
  theme_pub() + theme(axis.title = element_blank(), axis.text = element_text(size = 7),
                      plot.margin = margin(1, 1, 2, 1, "mm"))
save_fig(fig_s4_range, "FigSEM_Stage1_RangeModel", 89, 145)
message("FigS4a_RangeModel_Stage1 saved")

# Stage 1 × Stage 2: combined habitat suitability
suit_long <- wg_coast |>
  dplyr::select(lon, lat, dplyr::matches("^suit_(Seagrass|Mangrove|Kelp)$")) |>
  tidyr::pivot_longer(-c(lon, lat), names_to = "Type", values_to = "Suitability") |>
  dplyr::mutate(Type = factor(stringr::str_remove(Type, "suit_"), levels = ECO_LVLS))

fig_s4 <- ggplot() +
  geom_sf(data = world, fill = "grey85", colour = "white", linewidth = 0.1) +
  geom_tile(data = suit_long, aes(lon, lat, fill = Suitability), alpha = 0.9) +
  geom_sf(data = world, fill = NA, colour = "white", linewidth = 0.1) +
  ggh4x::facet_wrap2(~Type, ncol = 1, strip = eco_strip()) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Habitat\nsuitability",
                        limits = c(0, 1), labels = label_percent()) +
  coord_sf(ylim = c(-65, 75), expand = FALSE) +
  theme_pub() + theme(axis.title = element_blank(), axis.text = element_text(size = 7),
                      plot.margin = margin(1, 1, 2, 1, "mm"))
save_fig(fig_s4, "FigMain_GlobalSuitability", 89, 145)

# FigSup_SpeciesPred: species-level prediction maps only
# (ecosystem-level maps are in FigS4_GlobalSuitability — not repeated here)
# Layout: column = ecosystem type, row = rank within type (top-10 each)
if (exists("sp_models") && length(sp_models) > 0 &&
    exists("top3_per_type") && !is.null(top3_per_type)) {

  coast_mask <- world_grid$dist_km_baseline <= coast_thresh
  Xg_sp <- as.matrix(world_grid[, FEAT_SP]); colnames(Xg_sp) <- FEAT_SP

  map_theme_base <- list(
    theme_pub(),
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), legend.position = "none")
  )

  # Reorder: column-per-type (rank1_Seagrass | rank1_Mangrove | rank1_Kelp ...)
  top_sp_reordered <- top3_per_type |>
    dplyr::group_by(Type) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::arrange(rank, factor(Type, levels = ECO_LVLS))

  pred_map_panels <- purrr::map(seq_len(nrow(top_sp_reordered)), \(i) {
    tp  <- top_sp_reordered$Type[i]
    spp <- top_sp_reordered$.sp[i]
    key <- paste(tp, spp, sep = "__")
    if (!key %in% names(sp_models)) return(ggplot() + theme_void())

    suit_hab_sp <- predict(sp_models[[key]], as.data.frame(Xg_sp),
                           type = "cloglog") |> as.numeric()
    suit_rng_sp <- if (!is.null(range_mx[[tp]]))
      predict(range_mx[[tp]], as.data.frame(Xg_rng), type = "cloglog") |> as.numeric()
    else rep(1, nrow(world_grid))
    suit <- suit_hab_sp * suit_rng_sp
    suit_df <- data.frame(
      lon         = world_grid$lon[coast_mask],
      lat         = world_grid$lat[coast_mask],
      Suitability = suit[coast_mask]
    )
    auc_val   <- dplyr::filter(sp_vimp, Species == spp)$AUC
    type_col  <- ECO_PAL[tp]
    auc_label <- if (length(auc_val) > 0) sprintf("AUC = %.3f", auc_val[1]) else ""
    ggplot() +
      geom_sf(data = world, fill = "grey85", colour = "white", linewidth = 0.1) +
      geom_tile(data = suit_df, aes(lon, lat, fill = Suitability), alpha = 0.9) +
      geom_sf(data = world, fill = NA, colour = "white", linewidth = 0.08) +
      annotate("text", x = Inf, y = -Inf, label = auc_label,
               hjust = 1.05, vjust = -0.4, size = 2.0, colour = "grey35") +
      scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 1),
                            name = "Suitability", labels = scales::label_percent()) +
      coord_sf(ylim = c(-65, 75), expand = FALSE) +
      labs(title = paste("●", spp)) +
      map_theme_base +
      theme(plot.title = element_text(size = 7, face = "italic",
                                      colour = type_col,
                                      margin = margin(b = 1)))
  })
  sp_grid <- patchwork::wrap_plots(pred_map_panels, ncol = 3)

  n_sp_rows <- ceiling(nrow(top_sp_reordered) / 3)
  fig_ssp_pred <- sp_grid +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom",
          legend.key.width  = unit(1.2, "cm"),
          legend.key.height = unit(0.3, "cm"),
          legend.text = element_text(size = 7))
  save_fig(fig_ssp_pred, "FigSup_SpeciesPred", 170, n_sp_rows * 50 + 25)
  message("FigSup_SpeciesPred saved")
} else {
  message("FigSup_SpeciesPred skipped: sp_models or top3_per_type not available")
}

# =============================================================================
# PART 7 — RIVER FLOW REGULATION SCENARIOS (PATHWAY DECOMPOSITION)
# =============================================================================
# Three pathways tested separately and combined at each intensity:
#   (A) Kd490↓ only : dam traps sediment → clearer coastal water
#                     Kd490_sc = Kd490 × (1 − f × exp(−d/L))
#                     Ecologically BENEFICIAL for Seagrass/Kelp (more light)
#   (B) dist↑  only : dam reduces nutrient/freshwater export → functional isolation
#                     dist_eff = dist × (1 + f × exp(−d/L))
#                     Ecologically HARMFUL for all (dependency on river inputs)
#   (C) Combined A+B: realistic dual-pathway; dominant pathway revealed by decomposition
# f = regulation intensity (0.3/0.5/0.7), L = spatial decay scale (50 km)

SCENARIOS <- list(
  Baseline    = list(f = 0.00, L = 50, kd = FALSE, dk = FALSE, label = "Baseline",
                     grp = "decomp"),
  Kd490_only  = list(f = 0.50, L = 50, kd = TRUE,  dk = FALSE, label = "Kd490↓ only (50%)",
                     grp = "decomp"),
  Dist_only   = list(f = 0.50, L = 50, kd = FALSE, dk = TRUE,  label = "dist↑ only (50%)",
                     grp = "decomp"),
  Combined_50 = list(f = 0.50, L = 50, kd = TRUE,  dk = TRUE,  label = "Combined (50%)",
                     grp = "decomp"),
  Combined_30 = list(f = 0.30, L = 50, kd = TRUE,  dk = TRUE,  label = "30% regulation",
                     grp = "gradient"),
  Combined_70 = list(f = 0.70, L = 50, kd = TRUE,  dk = TRUE,  label = "70% regulation",
                     grp = "gradient")
)

message("\n=== Running scenarios (pathway decomposition + intensity gradient) ===")
sc_results <- list()
for (sc_name in names(SCENARIOS)) {
  sc    <- SCENARIOS[[sc_name]]
  decay <- exp(-world_grid$dist_km / sc$L)
  wg    <- world_grid |>
    dplyr::mutate(
      Kd490   = if (sc$kd) pmax(Kd490   * (1 - sc$f * decay), 1e-4) else Kd490,
      dist_km = if (sc$dk) dist_km * (1 + sc$f * decay)             else dist_km
    )
  Xs <- as.matrix(wg[, FEAT]); colnames(Xs) <- FEAT
  for (tp in ECO_LVLS) {
    wg[[paste0("suit_", tp)]] <-
      predict(all_mx[[tp]], as.data.frame(Xs), type = "cloglog") |> as.numeric()
  }
  sc_results[[sc_name]] <- wg
  message("  Done: ", sc_name)
}

# --- ΔSuitability quantitative summary ---
delta_all <- purrr::map_dfr(setdiff(names(SCENARIOS), "Baseline"), \(sc_name) {
  purrr::map_dfr(ECO_LVLS, \(tp) {
    base_s <- sc_results$Baseline[[paste0("suit_", tp)]]
    scen_s <- sc_results[[sc_name]][[paste0("suit_", tp)]]
    delta  <- scen_s - base_s
    tibble::tibble(
      Scenario    = SCENARIOS[[sc_name]]$label,
      Type        = tp,
      n_cells     = length(delta),
      n_gain_5pct = sum(delta >  0.05, na.rm = TRUE),
      n_loss_5pct = sum(delta < -0.05, na.rm = TRUE),
      pct_gain    = round(100 * mean(delta >  0.05, na.rm = TRUE), 2),
      pct_loss    = round(100 * mean(delta < -0.05, na.rm = TRUE), 2),
      mean_delta  = round(mean(delta, na.rm = TRUE), 4),
      max_gain    = round(max(delta,  na.rm = TRUE), 4),
      max_loss    = round(min(delta,  na.rm = TRUE), 4)
    )
  })
})
write.csv(delta_all, file.path(OUT_DIR, "ScenarioDeltaSummary.csv"), row.names = FALSE)
message("ΔSuitability summary:"); print(delta_all)

# Helper for scenario ΔSuitability data frames (used in FigSc4 below)
make_change_df <- function(base_name, scen_name) {
  sc_results[[base_name]] |>
    dplyr::select(lon, lat, starts_with("suit_")) |>
    dplyr::rename_with(\(x) stringr::str_replace(x, "suit_", "base_"), starts_with("suit_")) |>
    dplyr::left_join(
      sc_results[[scen_name]] |>
        dplyr::select(lon, lat, starts_with("suit_")) |>
        dplyr::rename_with(\(x) stringr::str_replace(x, "suit_", "scen_"), starts_with("suit_")),
      by = c("lon","lat")
    ) |>
    tidyr::pivot_longer(c(starts_with("base_"), starts_with("scen_")),
                         names_to = "key", values_to = "val") |>
    dplyr::mutate(Type   = stringr::str_extract(key, "Seagrass|Mangrove|Kelp"),
                  Period = dplyr::if_else(startsWith(key, "base_"), "base", "scen")) |>
    dplyr::filter(!is.na(Type)) |>
    tidyr::pivot_wider(names_from = Period, values_from = val) |>
    dplyr::mutate(Delta = scen - base, Type = factor(Type, levels = ECO_LVLS))
}

# Fig Sc2: Habitat change vs regulation intensity (% change from baseline, x = intensity %)
grad_names  <- c("Baseline","Combined_30","Combined_50","Combined_70")
reg_pct_map <- c(Baseline = 0, Combined_30 = 30, Combined_50 = 50, Combined_70 = 70)

area_grad <- purrr::map_dfr(grad_names, \(sc_name) {
  wg <- sc_results[[sc_name]] |> dplyr::filter(dist_km_baseline <= coast_thresh)
  purrr::map_dfr(ECO_LVLS, \(tp)
    tibble::tibble(
      sc_name    = sc_name,
      Type       = tp,
      mean_suit  = mean(wg[[paste0("suit_", tp)]], na.rm = TRUE)
    )
  )
}) |>
  dplyr::mutate(
    reg_pct  = reg_pct_map[sc_name],
    Type     = factor(Type, levels = ECO_LVLS),
    pct_chg  = (mean_suit / mean_suit[reg_pct == 0] - 1) * 100, .by = Type
  )

# FigSc2b: same gradient but with ±SD across species (species-level model predictions)
if (exists("sp_models") && length(sp_models) > 0) {
  message("Running species-level scenario gradient for error bands...")
  sp_sc_grad <- purrr::map_dfr(names(sp_models), \(key) {
    parts <- strsplit(key, "__")[[1]]
    tp_s  <- parts[1]; spp_s <- parts[2]
    mx_s  <- sp_models[[key]]
    purrr::map_dfr(grad_names, \(sc_name) {
      sc    <- SCENARIOS[[sc_name]]
      decay <- exp(-world_grid$dist_km / sc$L)
      wg_s  <- world_grid |>
        dplyr::mutate(
          Kd490   = if (sc$kd) pmax(Kd490   * (1 - sc$f * decay), 1e-4) else Kd490,
          dist_km = if (sc$dk) dist_km * (1 + sc$f * decay)             else dist_km
        ) |>
        dplyr::filter(dist_km_baseline <= coast_thresh)
      Xw_s <- as.matrix(wg_s[, FEAT]); colnames(Xw_s) <- FEAT
      suit_s <- predict(mx_s, as.data.frame(Xw_s), type = "cloglog") |> as.numeric()
      tibble::tibble(sc_name = sc_name, Type = tp_s, Species = spp_s,
                     mean_suit = mean(suit_s, na.rm = TRUE))
    })
  }) |>
    dplyr::mutate(reg_pct = reg_pct_map[sc_name], Type = factor(Type, levels = ECO_LVLS)) |>
    dplyr::group_by(Type, Species) |>
    dplyr::mutate(
      base_suit = mean_suit[reg_pct == 0][1],
      pct_chg   = (mean_suit / base_suit - 1) * 100
    ) |>
    dplyr::ungroup()

  sp_grad_summ <- sp_sc_grad |>
    dplyr::group_by(Type, reg_pct) |>
    dplyr::summarize(
      mean_chg = mean(pct_chg, na.rm = TRUE),
      sd_chg   = sd(pct_chg, na.rm = TRUE),
      n_sp     = dplyr::n(),
      .groups  = "drop"
    )

  sp_grad_summ <- sp_grad_summ |>
    dplyr::mutate(
      ci_lo = mean_chg - 1.96 * sd_chg / sqrt(pmax(n_sp, 1)),
      ci_hi = mean_chg + 1.96 * sd_chg / sqrt(pmax(n_sp, 1))
    )

  fig_sc2b <- ggplot(sp_grad_summ, aes(reg_pct, mean_chg, colour = Type, group = Type)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = Type),
                alpha = 0.20, colour = NA) +
    geom_line(linewidth = 1) + geom_point(size = 3) +
    scale_colour_manual(values = ECO_PAL, name = NULL) +
    scale_fill_manual(values = ECO_PAL, name = NULL, guide = "none") +
    # % lives in the axis titles only, so tick labels stay bare numbers
    scale_x_continuous(breaks = c(0, 30, 50, 70)) +
    scale_y_continuous(labels = scales::label_number()) +
    labs(x = "Dam-driven river flow reduction (%)",
         y = "Change in mean suitability (%)") +
    theme_pub() + theme(legend.position = c(0.15, 0.25))
  save_fig(fig_sc2b, "FigSc2b_Gradient_SpeciesCI", 130, 90)
  message("FigSc2b_Gradient_SpeciesCI saved")
}

write.csv(delta_all, file.path(OUT_DIR, "ScenarioDeltaSummary.csv"), row.names = FALSE)

# ── Checkpoint 3: save scenario results ──────────────────────────────────────
if (!exists("RDS_DIR")) { RDS_DIR <- file.path(OUT_DIR, "rds"); dir.create(RDS_DIR, showWarnings = FALSE) }
saveRDS(sc_results, file.path(RDS_DIR, "sc_results.rds"))
if (exists("sp_sc_grad") && !is.null(sp_sc_grad))
  saveRDS(sp_sc_grad, file.path(RDS_DIR, "sp_sc_grad.rds"))
message("Checkpoint 3 saved → ", RDS_DIR)

message("\n========== COMPLETE ==========")
message("Output: ", normalizePath(OUT_DIR))
message("GAM:      Fig1–7")
message("Turb:     FigT1–T3  (+Spearman_CorMatrix.csv)")
message("Dist:     FigD1–D3")
message("SDM:      FigS1–S4  (2-stage: range(SST,lat) × habitat(Sal,Kd490,log_dist) | depth=filter only)")
message("Dist:     FigD4_EmpiricalOccupancy — raw dist-from-river-mouth vs occupancy rate")
message("Scenario: FigSc2–Sc3  (+ScenarioDeltaSummary.csv)")
message("  FigSc2: % change in mean suitability vs regulation intensity")
message("  FigSc3: dist_km pathway only (Kd490 excluded — artifact direction in model)")

# =============================================================================
# PART 8 — SUPPLEMENTARY ANALYSES ADDED IN REVISION
# =============================================================================
# Turbidity decomposition, distance functional-form selection, river-size
# predictor comparison and the one-stage / two-stage contrast. Kept in their own
# file so this one stays readable; it runs in this environment and uses bc,
# bg_df, world_grid, FEAT and FEAT_range as they stand at this point.
if (file.exists("R/02_reviewer_analyses.R")) {
  source("R/02_reviewer_analyses.R", echo = FALSE)
} else {
  message("R/02_reviewer_analyses.R not found — supplementary analyses skipped")
}
