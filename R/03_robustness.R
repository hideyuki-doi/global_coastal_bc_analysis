# =============================================================================
# R/03_robustness.R — uncertainty for the revision's two key comparisons
#
# Runs after blue_carbon_full.R, from the checkpoint it writes to OUT_DIR/rds.
#
#   BC_OUT_DIR=output Rscript R/03_robustness.R
#
# PART 8 reports each comparison once: one background draw, one fixed set of
# longitude blocks, importance on the training data. That is enough to see a
# pattern but not to say how firm it is. Here every comparison is repeated
# N_REP times, each time with
#   * presences and background resampled with replacement (bootstrap), and
#   * the four longitude blocks rotated by a random offset, so a fold boundary
#     never sits in the same place twice,
# and scored out of sample on the held-out block. Each replicate fits every
# model on the same resampled data and blocks, so differences between models
# are paired within a replicate.
#
#   R1  one stage vs two stages (8d). Contribution of river distance as
#       (a) drop-variable ΔAUC: CV AUC of the full model minus CV AUC of the
#           same model refitted without distance — the unique, out-of-sample
#           contribution; and
#       (b) permutation ΔAUC on the held-out block.
#   R2  river size (8c) with runoff depth added. Discharge / catchment area is
#       runoff depth (mm per year): low in arid catchments however large they
#       are. Each predictor set is scored against "none" in the same replicate.
#
# Outputs: OUT_DIR/TabS_robust_*.csv
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(mgcv); library(parallel) })

OUT_DIR  <- Sys.getenv("BC_OUT_DIR", "output")
N_REP    <- as.integer(Sys.getenv("BC_N_REP", "50"))
N_CORES  <- as.integer(Sys.getenv("BC_CORES", "16"))
PARTS    <- strsplit(Sys.getenv("BC_ROBUST_PARTS", "R1,R2"), ",")[[1]]
ECO_LVLS <- c("Seagrass", "Mangrove", "Kelp")
DIST_VAR <- "log_dist"
REGMULT  <- 2.0
N_PERM   <- 5

bc    <- readRDS(file.path(OUT_DIR, "rds", "bc.rds"))
bg_df <- readRDS(file.path(OUT_DIR, "rds", "bg_df.rds"))
message(sprintf("03_robustness: %s | bc %d | bg %d | %d replicates on %d cores",
                OUT_DIR, nrow(bc), nrow(bg_df), N_REP, N_CORES))

# ── derived predictors ────────────────────────────────────────────────────────
# Land-attributable turbidity, exactly as in 8a: log Kd490 minus its ocean-
# climate fit (SST + |lat|), fitted over presences and background together.
kd_dat <- dplyr::bind_rows(bc[, c("Kd490", "SST", "abs_lat")],
                           bg_df[, c("Kd490", "SST", "abs_lat")]) |>
  dplyr::filter(Kd490 > 0, is.finite(SST), is.finite(abs_lat)) |>
  dplyr::mutate(log_Kd = log(Kd490))
m_ocean <- mgcv::gam(log_Kd ~ s(SST, k = 8) + s(abs_lat, k = 8), data = kd_dat)

# Runoff depth: 1 m³ s⁻¹ over 1 km² = 31,557.6 mm per year.
add_derived <- function(d) {
  d$kd_resid   <- log(d$Kd490) - as.numeric(predict(m_ocean, newdata = d))
  runoff_mm    <- d$discharge_cms / d$upstream_km2 * 31557.6
  d$log_runoff <- ifelse(d$upstream_km2 > 0 & d$discharge_cms >= 0,
                         log10(runoff_mm + 1), NA_real_)
  d
}
bc    <- add_derived(bc)
bg_df <- add_derived(bg_df)

# ── model sets ────────────────────────────────────────────────────────────────
DESIGNS <- list(   # R1: name = list(features, climate-matched background?)
  "1-stage, ocean-wide bg" = list(c("SST", "lat", "Salinity", "Kd490",    DIST_VAR, "log_discharge"), FALSE),
  "1-stage, Kd residual"   = list(c("SST", "lat", "Salinity", "kd_resid", DIST_VAR, "log_discharge"), FALSE),
  "2-stage (stage 2)"      = list(c("Salinity", "Kd490",    DIST_VAR, "log_discharge"), TRUE),
  "2-stage, Kd residual"   = list(c("Salinity", "kd_resid", DIST_VAR, "log_discharge"), TRUE)
)
SIZE_SETS <- list( # R2: all stage 2, climate-matched
  none               = c("Salinity", "Kd490", DIST_VAR),
  discharge          = c("Salinity", "Kd490", DIST_VAR, "log_discharge"),
  area               = c("Salinity", "Kd490", DIST_VAR, "log_upstream"),
  both               = c("Salinity", "Kd490", DIST_VAR, "log_upstream", "log_discharge"),
  runoff             = c("Salinity", "Kd490", DIST_VAR, "log_runoff"),
  "discharge+runoff" = c("Salinity", "Kd490", DIST_VAR, "log_discharge", "log_runoff")
)
ALL_COLS <- unique(c("lon", "lat", "SST", unlist(lapply(DESIGNS, `[[`, 1)), unlist(SIZE_SETS)))

# ── one replicate ─────────────────────────────────────────────────────────────
fit_mx <- function(y, X) {
  m <- try(maxnet::maxnet(y, X, maxnet::maxnet.formula(y, X, classes = "lqph"),
                          regmult = REGMULT), silent = TRUE)
  if (inherits(m, "try-error")) NULL else m
}
auc_of <- function(y, p) as.numeric(pROC::auc(y, p, quiet = TRUE))

# Resample presences and background, and draw the rotated blocks, for one
# ecosystem and one replicate. Rows are first restricted to those complete in
# every column any model uses, so all models in a replicate see the same rows.
draw_rep <- function(tp, r, match_climate) {
  set.seed(1000L * r + match(tp, ECO_LVLS))
  ok   <- function(d) d[stats::complete.cases(d[, ALL_COLS]) &
                        apply(as.matrix(d[, setdiff(ALL_COLS, "lon")]), 1, \(x) all(is.finite(x))), ]
  pres <- ok(as.data.frame(bc[bc$Type == tp, ]))
  bgp  <- ok(as.data.frame(bg_df))
  if (match_climate) {
    sq <- quantile(pres$SST, c(0.05, 0.95)); lq <- quantile(pres$lat, c(0.05, 0.95))
    m  <- bgp$SST >= sq[1] & bgp$SST <= sq[2] & bgp$lat >= lq[1] & bgp$lat <= lq[2]
    if (sum(m) >= 100) bgp <- bgp[m, ]
  }
  nbg  <- min(nrow(pres) * 10, nrow(bgp), 6000)
  df   <- rbind(cbind(pres[sample(nrow(pres), replace = TRUE), ALL_COLS], presence = 1L),
                cbind(bgp[sample(nrow(bgp), nbg, replace = TRUE),   ALL_COLS], presence = 0L))
  # four equal-count longitude blocks, rotated by a random offset
  lon_r    <- (df$lon + 180 + runif(1, 0, 360)) %% 360
  df$block <- as.integer(cut(lon_r, unique(quantile(lon_r, 0:4 / 4)), include.lowest = TRUE))
  df
}

# CV over the blocks for one feature set; optionally also the model without
# `drop` and the permutation loss of `drop` on the held-out block.
cv_set <- function(df, feats, drop = NULL) {
  res <- lapply(sort(unique(df$block)), function(b) {
    tr <- df[df$block != b, ]; te <- df[df$block == b, ]
    if (length(unique(te$presence)) < 2 || length(unique(tr$presence)) < 2) return(NULL)
    m <- fit_mx(tr$presence, tr[, feats, drop = FALSE]); if (is.null(m)) return(NULL)
    p <- as.numeric(predict(m, te[, feats, drop = FALSE], type = "cloglog"))
    out <- c(auc = auc_of(te$presence, p))
    if (!is.null(drop)) {
      f2 <- setdiff(feats, drop)
      m2 <- fit_mx(tr$presence, tr[, f2, drop = FALSE])
      out["auc_drop"] <- if (is.null(m2)) NA else
        auc_of(te$presence, as.numeric(predict(m2, te[, f2, drop = FALSE], type = "cloglog")))
      out["perm_loss"] <- mean(replicate(N_PERM, {
        tp <- te; tp[[drop]] <- sample(tp[[drop]])
        out["auc"] - auc_of(te$presence, as.numeric(predict(m, tp[, feats, drop = FALSE], type = "cloglog")))
      }))
    }
    out
  })
  res <- do.call(rbind, Filter(Negate(is.null), res))
  if (is.null(res)) return(NULL)
  colMeans(res, na.rm = TRUE)
}

summ <- function(x) c(median = median(x, na.rm = TRUE),
                      lo = unname(quantile(x, 0.025, na.rm = TRUE)),
                      hi = unname(quantile(x, 0.975, na.rm = TRUE)),
                      p_le0 = mean(x <= 0, na.rm = TRUE))
summ_df <- function(d, value, by) d |>
  dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
  dplyr::summarise(n = sum(is.finite(.data[[value]])),
                   median = median(.data[[value]], na.rm = TRUE),
                   lo  = quantile(.data[[value]], 0.025, na.rm = TRUE),
                   hi  = quantile(.data[[value]], 0.975, na.rm = TRUE),
                   p_le0 = mean(.data[[value]] <= 0, na.rm = TRUE), .groups = "drop")

jobs <- expand.grid(r = seq_len(N_REP), tp = ECO_LVLS, stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1")

# ═════════════════════════════════════════════════════════════════════════════
# R1 — one stage vs two stages
# ═════════════════════════════════════════════════════════════════════════════
if ("R1" %in% PARTS) {
  message("\n--- R1: distance contribution by design ---")
  t0 <- Sys.time()
  r1 <- dplyr::bind_rows(parallel::mclapply(seq_len(nrow(jobs)), function(i) {
    tp <- jobs$tp[i]; r <- jobs$r[i]
    dplyr::bind_rows(lapply(names(DESIGNS), function(dn) {
      d  <- DESIGNS[[dn]]
      df <- draw_rep(tp, r, match_climate = d[[2]])
      s  <- cv_set(df, d[[1]], drop = DIST_VAR)
      if (is.null(s)) return(NULL)
      data.frame(Type = tp, rep = r, design = dn, auc = s[["auc"]],
                 auc_drop = s[["auc_drop"]], dAUC_drop = s[["auc"]] - s[["auc_drop"]],
                 dAUC_perm = s[["perm_loss"]])
    }))
  }, mc.cores = N_CORES, mc.preschedule = FALSE))
  message(sprintf("  %d fits in %.1f min", nrow(r1), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  write.csv(r1, file.path(OUT_DIR, "TabS_robust_design_reps.csv"), row.names = FALSE)

  s1 <- dplyr::bind_rows(
    summ_df(r1, "dAUC_drop", c("Type", "design")) |> dplyr::mutate(metric = "dAUC_drop"),
    summ_df(r1, "dAUC_perm", c("Type", "design")) |> dplyr::mutate(metric = "dAUC_perm"),
    summ_df(r1, "auc",       c("Type", "design")) |> dplyr::mutate(metric = "CV AUC"))
  write.csv(s1, file.path(OUT_DIR, "TabS_robust_design_summary.csv"), row.names = FALSE)
  message("\n  drop-variable ΔAUC for river distance (median [95% interval]):")
  print(as.data.frame(s1 |> dplyr::filter(metric == "dAUC_drop") |>
    dplyr::mutate(txt = sprintf("%.3f [%.3f, %.3f]", median, lo, hi)) |>
    dplyr::select(design, Type, txt) |> tidyr::pivot_wider(names_from = Type, values_from = txt)))

  # paired within replicate: two-stage minus one-stage
  pr <- r1 |> dplyr::select(Type, rep, design, dAUC_drop, dAUC_perm) |>
    tidyr::pivot_wider(names_from = design, values_from = c(dAUC_drop, dAUC_perm))
  contr <- function(a, b, lab) dplyr::bind_rows(lapply(ECO_LVLS, function(tp) {
    x <- pr[pr$Type == tp, ]
    dplyr::bind_rows(
      data.frame(Type = tp, contrast = lab, metric = "dAUC_drop",
                 t(summ(x[[paste0("dAUC_drop_", a)]] - x[[paste0("dAUC_drop_", b)]]))),
      data.frame(Type = tp, contrast = lab, metric = "dAUC_perm",
                 t(summ(x[[paste0("dAUC_perm_", a)]] - x[[paste0("dAUC_perm_", b)]]))))
  }))
  c1 <- dplyr::bind_rows(
    contr("2-stage (stage 2)", "1-stage, ocean-wide bg", "2-stage minus 1-stage"),
    contr("2-stage (stage 2)", "1-stage, Kd residual",   "2-stage minus 1-stage Kd residual"),
    contr("1-stage, Kd residual", "1-stage, ocean-wide bg", "1-stage Kd residual minus 1-stage"))
  write.csv(c1, file.path(OUT_DIR, "TabS_robust_design_contrasts.csv"), row.names = FALSE)
  message("\n  paired contrasts (drop-variable ΔAUC; p_le0 = share of replicates <= 0):")
  print(as.data.frame(c1 |> dplyr::filter(metric == "dAUC_drop") |>
                      dplyr::mutate(dplyr::across(c(median, lo, hi), \(x) round(x, 3)))))
}

# ═════════════════════════════════════════════════════════════════════════════
# R2 — river size, with runoff depth
# ═════════════════════════════════════════════════════════════════════════════
if ("R2" %in% PARTS) {
  message("\n--- R2: river size predictors incl. runoff depth ---")
  rr <- cor(bc$log_runoff, bc$log_discharge, use = "complete.obs")
  ra <- cor(bc$log_runoff, bc$log_upstream,  use = "complete.obs")
  message(sprintf("  r(log runoff, log discharge) = %.3f | r(log runoff, log area) = %.3f", rr, ra))
  t0 <- Sys.time()
  r2 <- dplyr::bind_rows(parallel::mclapply(seq_len(nrow(jobs)), function(i) {
    tp <- jobs$tp[i]; r <- jobs$r[i]
    df <- draw_rep(tp, r, match_climate = TRUE)
    dplyr::bind_rows(lapply(names(SIZE_SETS), function(sn) {
      s <- cv_set(df, SIZE_SETS[[sn]])
      if (is.null(s)) return(NULL)
      data.frame(Type = tp, rep = r, set = sn, auc = s[["auc"]])
    }))
  }, mc.cores = N_CORES, mc.preschedule = FALSE))
  message(sprintf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  write.csv(r2, file.path(OUT_DIR, "TabS_robust_riversize_reps.csv"), row.names = FALSE)

  base <- r2 |> dplyr::filter(set == "none") |> dplyr::select(Type, rep, auc_none = auc)
  d2 <- r2 |> dplyr::filter(set != "none") |> dplyr::left_join(base, by = c("Type", "rep")) |>
    dplyr::mutate(gain = auc - auc_none)
  s2 <- dplyr::left_join(
    summ_df(r2, "auc", c("Type", "set")) |> dplyr::rename_with(\(x) paste0("auc_", x), c(median, lo, hi)) |>
      dplyr::select(Type, set, auc_median, auc_lo, auc_hi),
    summ_df(d2, "gain", c("Type", "set")) |> dplyr::select(Type, set, gain_median = median,
                                                             gain_lo = lo, gain_hi = hi, p_le0),
    by = c("Type", "set"))
  s2$r_runoff_discharge <- rr; s2$r_runoff_area <- ra
  write.csv(s2, file.path(OUT_DIR, "TabS_robust_riversize_summary.csv"), row.names = FALSE)
  message("\n  CV AUC gain over 'none' (median [95% interval], p_le0):")
  print(as.data.frame(s2 |> dplyr::filter(set != "none") |>
    dplyr::mutate(txt = sprintf("%+.3f [%+.3f, %+.3f] %.2f", gain_median, gain_lo, gain_hi, p_le0)) |>
    dplyr::select(set, Type, txt) |> tidyr::pivot_wider(names_from = Type, values_from = txt)))
}
message("\n=== 03_robustness complete ===")
