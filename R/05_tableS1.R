# =============================================================================
# R/05_tableS1.R — Supplementary Table S1 (taxa and record counts), revised
#
#   BC_OUT_DIR=output Rscript R/05_tableS1.R
#
# Reads what blue_carbon_full.R saves along the way:
#   OUT_DIR/bc_merged.rds         compiled records incl. the merged E. cava
#   OUT_DIR/bc_depthfiltered.rds  after the ≤100 m depth filter
#   OUT_DIR/rds/bc.rds            after the per-species cap (what the models use)
# and writes data/TableS1_taxa_revised.csv in the layout of the v5 table, with
# two extra record columns. Appended E. cava rows carry no genus field, so the
# genus is taken from the species name wherever it is missing.
# =============================================================================
suppressPackageStartupMessages(library(dplyr))
OUT_DIR <- Sys.getenv("BC_OUT_DIR", "output")

prep <- function(d) {
  d <- as.data.frame(d)
  named <- !is.na(d$species) & nzchar(d$species)
  d$genus <- ifelse(is.na(d$genus) | !nzchar(d$genus),
                    ifelse(named, sub(" .*", "", d$species), NA), d$genus)
  d$named <- named
  d
}
merged <- prep(readRDS(file.path(OUT_DIR, "bc_merged.rds")))
depth  <- prep(readRDS(file.path(OUT_DIR, "bc_depthfiltered.rds")))
used   <- prep(readRDS(file.path(OUT_DIR, "rds", "bc.rds")))

count_by <- function(d, nm) d |> dplyr::count(Type, genus, name = nm)
tab <- merged |>
  dplyr::group_by(Type, genus) |>
  dplyr::summarise(species_n = dplyr::n_distinct(species[named]), records = dplyr::n(),
                   genus_only = sum(!named), .groups = "drop") |>
  dplyr::left_join(count_by(depth, "after_depth"), by = c("Type", "genus")) |>
  dplyr::left_join(count_by(used,  "used"),        by = c("Type", "genus")) |>
  dplyr::mutate(dplyr::across(c(after_depth, used), \(x) dplyr::coalesce(x, 0L)),
                Type = factor(Type, levels = c("Seagrass", "Mangrove", "Kelp"))) |>
  dplyr::arrange(Type, dplyr::desc(records))

rows <- list()
for (tp in levels(tab$Type)) {
  t <- tab[tab$Type == tp, ]
  sp_tot <- dplyr::n_distinct(merged$species[merged$Type == tp & merged$named])
  for (i in seq_len(nrow(t))) rows[[length(rows) + 1]] <- data.frame(
    Ecosystem = if (i == 1) tp else "", Genus = t$genus[i], `Species (n)` = t$species_n[i],
    `Records (n)` = t$records[i], `Genus-only records (n)` = t$genus_only[i],
    `After depth filter (n)` = t$after_depth[i], `Used in models (n)` = t$used[i], check.names = FALSE)
  rows[[length(rows) + 1]] <- data.frame(
    Ecosystem = "", Genus = sprintf("Subtotal (%d genera)", nrow(t)), `Species (n)` = sp_tot,
    `Records (n)` = sum(t$records), `Genus-only records (n)` = sum(t$genus_only),
    `After depth filter (n)` = sum(t$after_depth), `Used in models (n)` = sum(t$used), check.names = FALSE)
  rows[[length(rows) + 1]] <- data.frame(Ecosystem = "", Genus = "", `Species (n)` = NA,
    `Records (n)` = NA, `Genus-only records (n)` = NA, `After depth filter (n)` = NA,
    `Used in models (n)` = NA, check.names = FALSE)
}
rows[[length(rows) + 1]] <- data.frame(
  Ecosystem = "Total", Genus = sprintf("%d genera", dplyr::n_distinct(merged$genus[!is.na(merged$genus)])),
  `Species (n)` = dplyr::n_distinct(merged$species[merged$named]),
  `Records (n)` = nrow(merged), `Genus-only records (n)` = sum(!merged$named),
  `After depth filter (n)` = nrow(depth), `Used in models (n)` = nrow(used), check.names = FALSE)
out <- dplyr::bind_rows(rows)
write.csv(out, "data/TableS1_taxa_revised.csv", row.names = FALSE, na = "")
print(out)
