#!/usr/bin/env Rscript
# =============================================================================
# Ecklonia cava species complex — GBIF retrieval under the Akita et al. (2020)
# taxonomy
#
# Akita, S., Hashimoto, K., Hanyuda, T. & Kawai, H. (2020) Molecular phylogeny
# and biogeography of Ecklonia spp. (Laminariales, Phaeophyceae) in Japan
# revealed taxonomic revision of E. kurome and E. stolonifera.
# Phycologia 59(4), 330-339.  doi:10.1080/00318884.2020.1756123
#
# E. cava, E. kurome and E. stolonifera are one species; the latter two become
# subspecies of E. cava. WoRMS and AlgaeBase have adopted this.
#
# WHY THIS SCRIPT EXISTS
# ----------------------
# The GBIF backbone taxonomy carries only TWO species in the genus Ecklonia:
# E. radiata (3196095) and E. maxima (3196111). E. cava, E. kurome and
# E. stolonifera are absent from it, so
#   name_backbone("Ecklonia cava")                 -> matchType = NONE
#   occ_search(scientificName = "Ecklonia cava")   -> 0 records
# and every record bearing one of those names is indexed at GENUS rank with an
# empty `species` field. Any species-level query therefore silently drops the
# entire East Asian Ecklonia material: in v1 the genus contributed a single
# species (E. radiata, 51 thinned records) plus 13 unusable genus-only records,
# which left the kelp group a Southern-Hemisphere and North-Atlantic sample.
#
# The records are reachable through `verbatimScientificName`, which matches the
# name as the publisher supplied it rather than the backbone's reading of it.
#
# SCOPE: only the names the backbone CANNOT resolve. E. radiata and E. maxima
# already enter blue_carbon_full.R correctly through the backbone, so
# retrieving them here would duplicate work and pull ~29,000 redundant records.
# Their counts are reported for the record but not downloaded.
#
# Output: data/ecklonia_records.csv (species_reported, species_accepted, lon, lat, ...)
# Environmental predictors are NOT attached here. blue_carbon_full.R extracts
# them from the same env_rast as every other record, so there is one extraction
# path and no chance of a mismatch.
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite); library(dplyr); library(purrr)
})

OUT_CSV <- "data/ecklonia_records.csv"
dir.create("data", showWarnings = FALSE)

API     <- "https://api.gbif.org/v1/occurrence/search"
PAGE    <- 300L     # GBIF REST maximum per page
RETRIES <- 5L       # per page, with exponential backoff
options(timeout = max(300, getOption("timeout")))

# ── Names the GBIF backbone cannot resolve, and what each maps to ────────────
# The cava complex collapses to one species. Eckloniopsis radicosa was moved
# into Ecklonia by Akita et al. but kept as its OWN species, so it is not merged.
NAMES <- tibble::tribble(
  ~verbatim,                          ~accepted,
  "Ecklonia cava",                    "Ecklonia cava",
  "Ecklonia kurome",                  "Ecklonia cava",
  "Ecklonia stolonifera",             "Ecklonia cava",
  "Ecklonia cava subsp. cava",        "Ecklonia cava",
  "Ecklonia cava subsp. kurome",      "Ecklonia cava",
  "Ecklonia cava subsp. stolonifera", "Ecklonia cava",
  "Ecklonia cava f. stolonifera",     "Ecklonia cava",
  "Ecklonia radicosa",                "Ecklonia radicosa",
  "Eckloniopsis radicosa",            "Ecklonia radicosa",
  "Ecklonia arborea",                 "Ecklonia arborea",
  "Ecklonia bicyclis",                "Ecklonia bicyclis",
  "Ecklonia brevipes",                "Ecklonia brevipes"
)

# Resolved by the backbone already — reported, never downloaded here.
BACKBONE_OK <- c("Ecklonia radiata", "Ecklonia maxima")

# Basis-of-record classes meaning an organism was actually at the coordinate.
# FOSSIL_SPECIMEN and LIVING_SPECIMEN (culture / aquarium material) excluded.
OK_BASIS <- c("PRESERVED_SPECIMEN", "HUMAN_OBSERVATION",
              "MACHINE_OBSERVATION", "MATERIAL_SAMPLE", "OCCURRENCE")

KEEP <- c("verbatimScientificName", "scientificName", "decimalLongitude",
          "decimalLatitude", "basisOfRecord", "year", "countryCode",
          "datasetKey")

build_url <- function(verbatim, offset, limit) {
  sprintf("%s?verbatimScientificName=%s&hasCoordinate=true&hasGeospatialIssue=false&occurrenceStatus=PRESENT&offset=%d&limit=%d",
          API, utils::URLencode(verbatim, reserved = TRUE), offset, limit)
}

#' One page, retried with exponential backoff. GBIF intermittently times out on
#' deep offsets, and a single failure must not discard an otherwise complete
#' retrieval, so every page gets its own retry budget.
get_json <- function(url) {
  for (k in seq_len(RETRIES)) {
    js <- try(jsonlite::fromJSON(url, flatten = TRUE), silent = TRUE)
    if (!inherits(js, "try-error")) return(js)
    wait <- 2^k
    message(sprintf("      attempt %d/%d failed; retrying in %ds", k, RETRIES, wait))
    Sys.sleep(wait)
  }
  stop("GBIF unreachable after ", RETRIES, " attempts: ", url)
}

count_only <- function(verbatim) {
  js <- get_json(build_url(verbatim, 0L, 0L))
  js$count %||% 0L
}

fetch_name <- function(verbatim) {
  total <- count_only(verbatim)
  message(sprintf("  %-34s %s georeferenced", verbatim, format(total, big.mark = ",")))
  if (total == 0) return(NULL)

  out <- list(); offset <- 0L
  while (offset < total && offset < 100000L) {
    js <- get_json(build_url(verbatim, offset, PAGE))
    if (is.null(js$results) || length(js$results) == 0) break
    d <- js$results
    for (k in KEEP) if (is.null(d[[k]])) d[[k]] <- NA   # GBIF omits empty columns
    out[[length(out) + 1L]] <- d[, KEEP, drop = FALSE]
    offset <- offset + PAGE
    if (isTRUE(js$endOfRecords)) break
  }
  if (!length(out)) return(NULL)
  dplyr::bind_rows(out) |> dplyr::mutate(query = verbatim)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

message("Backbone-resolved species (handled by blue_carbon_full.R, not downloaded):")
for (nm in BACKBONE_OK)
  message(sprintf("  %-34s %s georeferenced under the verbatim name",
                  nm, format(count_only(nm), big.mark = ",")))

message("\nQuerying GBIF (verbatimScientificName) for the names the backbone misses ...")
raw <- purrr::map_dfr(NAMES$verbatim, fetch_name)
if (nrow(raw) == 0) stop("GBIF returned nothing — check the network connection.")

# ── Quality control ──────────────────────────────────────────────────────────
# Mirrors the screening in the Methods: coordinates present and free of flagged
# geospatial issues (both applied in the query), a basis of record implying
# presence at the coordinate, not the (0,0) null-island artefact, and one record
# per distinct coordinate per accepted species so that repeat sampling of one
# site does not weight the model. The 1-degree spatial thinning used for the
# models is applied later, in blue_carbon_full.R, with every other taxon.
clean <- raw |>
  dplyr::left_join(NAMES, by = c("query" = "verbatim")) |>
  dplyr::filter(
    basisOfRecord %in% OK_BASIS,
    is.finite(decimalLongitude), is.finite(decimalLatitude),
    !(abs(decimalLongitude) < 0.01 & abs(decimalLatitude) < 0.01),
    abs(decimalLatitude) <= 90, abs(decimalLongitude) <= 180
  ) |>
  dplyr::transmute(
    # the name WE queried, not GBIF's reading of it: for these taxa the backbone
    # collapses everything to the genus, so scientificName says "Ecklonia
    # Hornemann, 1828" for all of them and the provenance would be lost
    species_reported = query,
    species_accepted = accepted,
    lon = decimalLongitude, lat = decimalLatitude,
    year, basisOfRecord, countryCode, datasetKey
  ) |>
  dplyr::distinct(species_accepted, lon, lat, .keep_all = TRUE)

message("\n--- unique coordinates retained per accepted species ---")
print(as.data.frame(dplyr::count(clean, species_accepted, sort = TRUE)))

cava <- dplyr::filter(clean, species_accepted == "Ecklonia cava")
message("\n--- cava complex: contribution of each reported name ---")
print(as.data.frame(dplyr::count(cava, species_reported, sort = TRUE)))
message("\n--- cava complex: top countries ---")
print(as.data.frame(utils::head(dplyr::count(cava, countryCode, sort = TRUE), 10)))
message(sprintf("\nE. cava after merging the complex: %d unique coordinates", nrow(cava)))

write.csv(clean, OUT_CSV, row.names = FALSE)
message("Written: ", OUT_CSV)
