# =============================================================================
# R/07_sigk_compare.R — three definitions of σk over ten background draws
#
#   Rscript R/07_sigk_compare.R output_bgseeds
#
# Runs named struct_s<seed> (σk from log Kd490 ~ S + T), structSTx_s<seed>
# (~ S * T) and structT_s<seed> (~ T only), as written by blue_carbon_full.R
# with BC_SIGK. Reports CV AUC, importance, the 70% regulation effect, and the
# shape of the σk response, each paired against ST within a draw.
# =============================================================================
suppressPackageStartupMessages(library(dplyr))
root <- commandArgs(TRUE)[1]; if (is.na(root)) root <- "output_bgseeds"
runs <- list.dirs(root, recursive = FALSE)
runs <- runs[grepl("^struct(STx|T)?_s[0-9]+$", basename(runs)) &
             file.exists(file.path(runs, "rds", "sp_sc_grad.rds"))]
lab  <- function(x) c(struct = "ST", structSTx = "STx", structT = "T")[sub("_s[0-9]+$", "", x)]

one <- function(d) {
  rds  <- function(f) readRDS(file.path(d, "rds", f))
  roc  <- rds("all_roc.rds"); vimp <- dplyr::bind_rows(rds("all_vimp.rds"))
  sc   <- rds("sp_sc_grad.rds"); pd <- rds("sp_pd_summary.rds")
  dplyr::bind_rows(lapply(names(roc), function(tp) {
    a <- roc[[tp]]$mx_cv_auc; a <- a[is.finite(a)]
    p <- pd[pd$Type == tp & pd$Variable == "sigk", ]
    data.frame(sigk = lab(basename(d)), seed = as.integer(sub(".*_s", "", basename(d))), Type = tp,
               cv_auc = mean(a),
               imp_sigk = vimp$MaxEnt[vimp$Type == tp & vimp$Variable == "sigk"],
               imp_dist = vimp$MaxEnt[vimp$Type == tp & vimp$Variable == "log_dist"],
               imp_sal  = vimp$MaxEnt[vimp$Type == tp & vimp$Variable == "Salinity"],
               reg70 = mean(sc$pct_chg[sc$Type == tp & sc$reg_pct == 70]),
               sigk_peak = if (nrow(p)) p$x[which.max(p$suit_mean)] else NA)
  }))
}
all <- dplyr::bind_rows(lapply(runs, one))
write.csv(all, file.path(root, "sigk_runs.csv"), row.names = FALSE)

med <- function(x) sprintf("%.3f [%.3f, %.3f]", median(x), min(x), max(x))
pct <- function(x) sprintf("%.0f [%.0f, %.0f]", 100 * median(x), 100 * min(x), 100 * max(x))
reg <- function(x) sprintf("%.1f [%.1f, %.1f]", median(x), min(x), max(x))
summ <- all |> dplyr::group_by(sigk, Type) |>
  dplyr::summarise(n = dplyr::n(), cv_auc = med(cv_auc), imp_sigk = pct(imp_sigk),
                   imp_dist = pct(imp_dist), imp_sal = pct(imp_sal), reg70 = reg(reg70),
                   sigk_peak = sprintf("%.2f", median(sigk_peak, na.rm = TRUE)), .groups = "drop") |>
  dplyr::arrange(Type, sigk)
write.csv(summ, file.path(root, "sigk_summary.csv"), row.names = FALSE)
print(as.data.frame(summ))

base <- all |> dplyr::filter(sigk == "ST") |> dplyr::select(seed, Type, auc0 = cv_auc, reg0 = reg70)
pair <- all |> dplyr::filter(sigk != "ST") |> dplyr::inner_join(base, by = c("seed", "Type")) |>
  dplyr::group_by(sigk, Type) |>
  dplyr::summarise(n = dplyr::n(),
                   d_auc = sprintf("%+.3f [%+.3f, %+.3f]", median(cv_auc - auc0), min(cv_auc - auc0), max(cv_auc - auc0)),
                   d_reg70 = sprintf("%+.1f [%+.1f, %+.1f]", median(reg70 - reg0), min(reg70 - reg0), max(reg70 - reg0)),
                   .groups = "drop")
write.csv(pair, file.path(root, "sigk_paired.csv"), row.names = FALSE)
message("\npaired against ST, same draw:")
print(as.data.frame(pair))
