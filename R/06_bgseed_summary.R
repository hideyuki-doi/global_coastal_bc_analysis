# =============================================================================
# R/06_bgseed_summary.R — headline numbers over background draws
#
#   Rscript R/06_bgseed_summary.R output_bgseeds
#
# Each subdirectory <set>_s<seed> is one run of blue_carbon_full.R with
# BC_FEAT=<set> and BC_BG_SEED=<seed> (see RUN_NOTES §3m). One background
# redraw moves AUC by 0.03-0.06 (§3g), so the numbers the paper quotes are the
# median and range over draws, and the two variable sets are compared within
# each draw.
# =============================================================================
suppressPackageStartupMessages(library(dplyr))
root <- commandArgs(TRUE)[1]; if (is.na(root)) root <- "output_bgseeds"
runs <- list.dirs(root, recursive = FALSE)
runs <- runs[file.exists(file.path(runs, "rds", "sp_sc_grad.rds"))]
message(sprintf("%d complete runs in %s", length(runs), root))

one <- function(d) {
  tag  <- basename(d); set <- sub("_s[0-9]+$", "", tag); seed <- as.integer(sub(".*_s", "", tag))
  rds  <- function(f) readRDS(file.path(d, "rds", f))
  roc  <- rds("all_roc.rds"); rng <- rds("range_auc.rds")
  vimp <- dplyr::bind_rows(rds("all_vimp.rds"))
  st   <- rds("sp_vimp_thinned.rds")
  sc   <- rds("sp_sc_grad.rds")
  nbg  <- nrow(rds("bg_df.rds"))
  dist_var <- grep("dist", unique(vimp$Variable), value = TRUE)[1]
  dplyr::bind_rows(lapply(names(roc), function(tp) {
    a <- roc[[tp]]$mx_cv_auc; a <- a[is.finite(a)]
    data.frame(set = set, seed = seed, Type = tp, n_bg = nbg,
               cv_auc = mean(a), stage1_auc = rng[[tp]],
               imp_dist_eco = vimp$MaxEnt[vimp$Type == tp & vimp$Variable == dist_var],
               imp_dist_sp  = median(st$Importance[st$Type == tp & st$Variable == dist_var]),
               imp_sal_sp   = median(st$Importance[st$Type == tp & st$Variable == "Salinity"]),
               reg70 = mean(sc$pct_chg[sc$Type == tp & sc$reg_pct == 70]))
  }))
}
all <- dplyr::bind_rows(lapply(runs, one))
write.csv(all, file.path(root, "bgseed_runs.csv"), row.names = FALSE)

vars <- c("cv_auc", "stage1_auc", "imp_dist_eco", "imp_dist_sp", "imp_sal_sp", "reg70")
summ <- all |> dplyr::group_by(set, Type) |>
  dplyr::summarise(n_seeds = dplyr::n(),
                   dplyr::across(dplyr::all_of(vars),
                                 list(med = median, min = min, max = max)), .groups = "drop")
write.csv(summ, file.path(root, "bgseed_summary.csv"), row.names = FALSE)

# struct minus current, paired by seed
pair <- all |> dplyr::select(set, seed, Type, cv_auc, reg70) |>
  tidyr::pivot_wider(names_from = set, values_from = c(cv_auc, reg70)) |>
  dplyr::mutate(d_auc = cv_auc_struct - cv_auc_current, d_reg70 = reg70_struct - reg70_current)
pair_s <- pair |> dplyr::group_by(Type) |>
  dplyr::summarise(n = sum(is.finite(d_auc)),
                   d_auc_med = median(d_auc, na.rm = TRUE), d_auc_min = min(d_auc, na.rm = TRUE),
                   d_auc_max = max(d_auc, na.rm = TRUE),
                   d_reg70_med = median(d_reg70, na.rm = TRUE), d_reg70_min = min(d_reg70, na.rm = TRUE),
                   d_reg70_max = max(d_reg70, na.rm = TRUE),
                   share_struct_more_negative = mean(d_reg70 < 0, na.rm = TRUE), .groups = "drop")
write.csv(pair_s, file.path(root, "bgseed_paired.csv"), row.names = FALSE)

fmt <- function(m, lo, hi, k = 3) sprintf(paste0("%.", k, "f [%.", k, "f, %.", k, "f]"), m, lo, hi)
print(as.data.frame(summ |> dplyr::transmute(set, Type, n_seeds,
  cv_auc = fmt(cv_auc_med, cv_auc_min, cv_auc_max),
  stage1 = fmt(stage1_auc_med, stage1_auc_min, stage1_auc_max),
  dist_eco = fmt(100 * imp_dist_eco_med, 100 * imp_dist_eco_min, 100 * imp_dist_eco_max, 0),
  reg70 = fmt(reg70_med, reg70_min, reg70_max, 1))))
message("\nstruct − current, paired by seed:")
print(as.data.frame(pair_s |> dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, 3)))))
message(sprintf("background points after depth filter: %s",
                paste(sort(unique(all$n_bg)), collapse = ", ")))
