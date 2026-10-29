# =============================================================================
# R/04_structural.R — the co-author's structural model (comment of 2026-10-30)
#
# Runs after blue_carbon_full.R, from the checkpoint in OUT_DIR/rds.
#
#   BC_OUT_DIR=output Rscript R/04_structural.R
#
# The proposal makes the river pathway explicit instead of entering distance
# and river size as SDM predictors:
#
#   (1) salinity carries the land input:   S = S_ref − c · B^β · exp(−D / L)
#       S_ref is the open-ocean salinity nearby, B catchment area, D distance
#       to the mouth. Linear on the log scale:
#           log(S_ref − S) = log c + β·log B − D / L
#   (2) turbidity = salinity + temperature + the rest:
#           log K = a0 + a1·S + a2·T + σk
#       σk, the part of turbidity that neither salinity nor temperature
#       explains, is read as the ocean-derived component.
#   (3) the SDM takes T, S and σk; D and B act only through S.
#
# This script fits (1) and (2), checks how much land signal is left in σk,
# compares the structural SDM with the current one by bootstrap CV (as in
# 03_robustness.R), tests whether D adds anything beyond S, and turns β into a
# contribution of catchment area to suitability via salinity.
#
# Outputs: OUT_DIR/TabS_struct_*.csv
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(terra); library(parallel) })

OUT_DIR  <- Sys.getenv("BC_OUT_DIR", "output")
N_REP    <- as.integer(Sys.getenv("BC_N_REP", "50"))
N_CORES  <- as.integer(Sys.getenv("BC_CORES", "16"))
ECO_LVLS <- c("Seagrass", "Mangrove", "Kelp")
REGMULT  <- 2.0
options(sdmpredictors_datadir = "sdm_data")

bc    <- as.data.frame(readRDS(file.path(OUT_DIR, "rds", "bc.rds")))
bg_df <- as.data.frame(readRDS(file.path(OUT_DIR, "rds", "bg_df.rds")))
grid_cache <- file.path(OUT_DIR, "world_grid_overwater.rds")
wg <- if (file.exists(grid_cache)) as.data.frame(readRDS(grid_cache)) else NULL
message(sprintf("04_structural: %s | bc %d | bg %d | grid %s", OUT_DIR, nrow(bc),
                nrow(bg_df), if (is.null(wg)) "none" else nrow(wg)))

# ═════════════════════════════════════════════════════════════════════════════
# Reference salinity: the open-ocean value nearby
# ═════════════════════════════════════════════════════════════════════════════
# The 90th percentile of salinity in a W-degree block, interpolated back to
# the native grid. High quantile rather than max, so one hypersaline lagoon
# cell does not set the reference; a block wide enough to reach past the
# plumes of all but the largest rivers. W = 2 is the main setting.
sal <- terra::rast(sdmpredictors::load_layers("BO_salinity", equalarea = FALSE))
res_deg <- terra::res(sal)[1]
s_ref_rast <- function(w_deg) {
  f  <- max(1, round(w_deg / res_deg))
  ag <- terra::aggregate(sal, fact = f, na.rm = TRUE,
                         fun = function(x, ...) stats::quantile(x, 0.9, na.rm = TRUE, names = FALSE))
  ag <- terra::focal(ag, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE)  # fill coastal gaps
  terra::resample(ag, sal, method = "bilinear")
}
add_sref <- function(d, r, nm) {
  d[[nm]] <- terra::extract(r, as.matrix(d[, c("lon", "lat")]))[, 1]
  d
}
sref <- list(`1` = s_ref_rast(1), `2` = s_ref_rast(2), `4` = s_ref_rast(4))
for (w in names(sref)) {
  nm <- paste0("S_ref", w)
  bc <- add_sref(bc, sref[[w]], nm); bg_df <- add_sref(bg_df, sref[[w]], nm)
  if (!is.null(wg)) wg <- add_sref(wg, sref[[w]], nm)
}

# ═════════════════════════════════════════════════════════════════════════════
# (1) Salinity equation
# ═════════════════════════════════════════════════════════════════════════════
message("\n--- (1) salinity: log(S_ref − S) = log c + β log B − D / L ---")
# Fitted on coastal grid cells (depth ≤ 100 m, within 500 km of a mouth): the
# physical relation, independent of where species happen to be recorded. Also
# fitted on presences + background as a check. Cells with ΔS ≤ 0.05 PSU carry
# no measurable freshwater signal and cannot be logged; their share is
# reported, and a nonlinear fit on the raw scale keeps them in.
fit_sal <- function(d, w = "2", size = "upstream_km2", label = "") {
  d$dS   <- d[[paste0("S_ref", w)]] - d$Salinity
  d$logB <- log10(d[[size]])
  d <- d[is.finite(d$dS) & is.finite(d$logB) & is.finite(d$dist_km), ]
  pos <- d$dS > 0.05
  m <- lm(log(dS) ~ logB + dist_km, data = d[pos, ])
  co <- coef(m)
  # nonlinear fit on the salinity scale, all cells
  nl <- try(nls(dS ~ exp(lc + beta * logB - dist_km / L), data = d,
                start = list(lc = co[[1]], beta = co[["logB"]] / log(10) * log(10),
                             L = max(5, -1 / co[["dist_km"]])),
                control = nls.control(maxiter = 200, warnOnly = TRUE)), silent = TRUE)
  nlc <- if (inherits(nl, "try-error")) c(lc = NA, beta = NA, L = NA) else coef(nl)
  r2_nl <- if (inherits(nl, "try-error")) NA else
    1 - sum(residuals(nl)^2) / sum((d$dS - mean(d$dS))^2)
  data.frame(fit = label, window_deg = as.numeric(w), size = size, n = nrow(d),
             share_dS_pos = mean(pos),
             beta_log = co[["logB"]], elasticity = co[["logB"]] / log(10),
             L_log_km = -1 / co[["dist_km"]],
             p_beta = summary(m)$coefficients["logB", 4],
             p_D = summary(m)$coefficients["dist_km", 4],
             R2_log = summary(m)$r.squared,
             beta_nls = nlc[["beta"]], L_nls_km = nlc[["L"]], R2_nls_dS = r2_nl)
}
sal_dat_grid <- if (!is.null(wg)) wg[wg$depth <= 100 & wg$dist_km <= 500, ] else NULL
sal_dat_pts  <- dplyr::bind_rows(bc, bg_df)
sal_dat_pts  <- sal_dat_pts[sal_dat_pts$dist_km <= 500, ]
sal_tab <- dplyr::bind_rows(
  if (!is.null(sal_dat_grid)) dplyr::bind_rows(lapply(names(sref), \(w) fit_sal(sal_dat_grid, w, "upstream_km2", "grid"))),
  if (!is.null(sal_dat_grid)) fit_sal(sal_dat_grid, "2", "discharge_cms", "grid"),
  fit_sal(sal_dat_pts, "2", "upstream_km2", "points"),
  fit_sal(sal_dat_pts, "2", "discharge_cms", "points"))
print(sal_tab |> dplyr::mutate(dplyr::across(where(is.numeric), \(x) signif(x, 3))))
write.csv(sal_tab, file.path(OUT_DIR, "TabS_struct_salinity.csv"), row.names = FALSE)

# How much does salinity actually change with distance from a mouth? If the
# freshwater deficit barely varies across the first tens of km, salinity at
# this resolution cannot carry the distance effect, whatever the model.
if (!is.null(sal_dat_grid)) {
  gb <- sal_dat_grid |>
    dplyr::mutate(dS = S_ref2 - Salinity,
                  D_band = cut(dist_km, c(0, 10, 25, 50, 100, 200, 500), include.lowest = TRUE),
                  B_tercile = dplyr::ntile(upstream_km2, 3)) |>
    dplyr::filter(is.finite(dS))
  band_tab <- gb |> dplyr::group_by(D_band) |>
    dplyr::summarise(n = dplyr::n(), dS_median = median(dS), dS_q75 = quantile(dS, .75),
                     dS_q90 = quantile(dS, .9), .groups = "drop")
  band_B <- gb |> dplyr::filter(dist_km <= 50) |> dplyr::group_by(B_tercile) |>
    dplyr::summarise(n = dplyr::n(), B_median_km2 = median(upstream_km2),
                     dS_median = median(dS), dS_q90 = quantile(dS, .9), .groups = "drop")
  message("\n  freshwater deficit ΔS (PSU) by distance band, coastal grid:")
  print(as.data.frame(band_tab |> dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, 3)))))
  message("  ... and by catchment-size tercile within 50 km:")
  print(as.data.frame(band_B |> dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, 3)))))
  write.csv(band_tab, file.path(OUT_DIR, "TabS_struct_dS_by_distance.csv"), row.names = FALSE)
  write.csv(band_B,   file.path(OUT_DIR, "TabS_struct_dS_by_catchment.csv"), row.names = FALSE)
}

main_sal <- sal_tab[sal_tab$fit == (if (!is.null(sal_dat_grid)) "grid" else "points") &
                    sal_tab$window_deg == 2 & sal_tab$size == "upstream_km2", ][1, ]
BETA  <- main_sal$beta_log
L_HAT <- main_sal$L_log_km
message(sprintf("  main: ΔS ∝ B^%.3f | L = %.1f km", BETA / log(10), L_HAT))
if (!is.finite(L_HAT) || L_HAT <= 0) stop("salinity does not decay with distance; L not estimable")

# ═════════════════════════════════════════════════════════════════════════════
# (2) Turbidity equation and σk
# ═════════════════════════════════════════════════════════════════════════════
message("\n--- (2) turbidity: log K = a0 + a1 S + a2 T + σk ---")
kd_dat <- dplyr::bind_rows(bc[, c("Kd490", "Salinity", "SST")], bg_df[, c("Kd490", "Salinity", "SST")]) |>
  dplyr::filter(Kd490 > 0, is.finite(Salinity), is.finite(SST))
m_k <- lm(log(Kd490) ~ Salinity + SST, data = kd_dat)
m_k_lin <- lm(Kd490 ~ Salinity + SST, data = kd_dat)  # the comment's linear form
sig <- function(d) log(d$Kd490) - as.numeric(predict(m_k, newdata = d))
bc$sigk <- sig(bc); bg_df$sigk <- sig(bg_df)
kt <- as.data.frame(summary(m_k)$coefficients)
kt$term <- rownames(kt); rownames(kt) <- NULL
kt$R2 <- summary(m_k)$r.squared; kt$R2_linear_form <- summary(m_k_lin)$r.squared
print(kt)
write.csv(kt, file.path(OUT_DIR, "TabS_struct_turbidity.csv"), row.names = FALSE)

# Is σk free of land signal? Regress it on distance and river size.
chk <- dplyr::bind_rows(bc, bg_df) |>
  dplyr::filter(is.finite(sigk), is.finite(log_dist), is.finite(log_discharge))
m_chk <- lm(sigk ~ log_dist + log_discharge, data = chk)
cc <- summary(m_chk)$coefficients
sig_chk <- data.frame(term = rownames(cc), estimate = cc[, 1], p = cc[, 4],
                      R2 = summary(m_chk)$r.squared,
                      r_sigk_logdist = cor(chk$sigk, chk$log_dist),
                      r_sigk_SST = cor(chk$sigk, chk$SST))
message(sprintf("  σk ~ log_dist + log_discharge: R² = %.3f | r(σk, log dist) = %.3f",
                summary(m_chk)$r.squared, cor(chk$sigk, chk$log_dist)))
write.csv(sig_chk, file.path(OUT_DIR, "TabS_struct_sigk_check.csv"), row.names = FALSE)

# distance on the decay scale estimated from salinity
bc$expD_hat    <- exp(-bc$dist_km / L_HAT)
bg_df$expD_hat <- exp(-bg_df$dist_km / L_HAT)

# Plume index: the comment's salinity equation used directly as a predictor,
# ΔS ∝ size · exp(−D / L), with river size = discharge (or catchment area for
# the comment's B) and L fixed, since salinity cannot estimate it (§1).
PLUME_L <- c(10, 25, 50, 100)
for (L in PLUME_L) {
  bc[[paste0("plumeQ_", L)]]    <- log10(1 + bc$discharge_cms    * exp(-bc$dist_km    / L))
  bg_df[[paste0("plumeQ_", L)]] <- log10(1 + bg_df$discharge_cms * exp(-bg_df$dist_km / L))
}
bc$plumeB_50    <- log10(1 + bc$upstream_km2    * exp(-bc$dist_km    / 50))
bg_df$plumeB_50 <- log10(1 + bg_df$upstream_km2 * exp(-bg_df$dist_km / 50))

# ═════════════════════════════════════════════════════════════════════════════
# (3) SDM comparison, bootstrap CV (as 03_robustness.R)
# ═════════════════════════════════════════════════════════════════════════════
message("\n--- (3) SDM: structural vs current ---")
MODELS <- list(
  # two-stage: T lives in stage 1, so stage 2 takes S and σk
  "2s current (S, K, logD, Q)" = list(c("Salinity", "Kd490", "log_dist", "log_discharge"), TRUE),
  "2s structural (S, σk)"      = list(c("Salinity", "sigk"), TRUE),
  "2s structural + exp(-D/L)"  = list(c("Salinity", "sigk", "expD_hat"), TRUE),
  "2s structural + logD"       = list(c("Salinity", "sigk", "log_dist"), TRUE),
  # one stage, as the comment writes it: T, S, σk against the whole coast
  "1s current (T, lat, S, K, logD, Q)" = list(c("SST", "lat", "Salinity", "Kd490", "log_dist", "log_discharge"), FALSE),
  "1s structural (T, S, σk)"           = list(c("SST", "Salinity", "sigk"), FALSE),
  "1s structural + exp(-D/L)"          = list(c("SST", "Salinity", "sigk", "expD_hat"), FALSE)
)
MODELS5 <- c(
  list("2s current (S, K, logD, Q)" = c("Salinity", "Kd490", "log_dist", "log_discharge"),
       "2s no river size (S, K, logD)" = c("Salinity", "Kd490", "log_dist")),
  setNames(lapply(PLUME_L, \(L) c("Salinity", "Kd490", paste0("plumeQ_", L))),
           paste0("2s plume Q, L=", PLUME_L, " (S, K, P)")),
  list("2s plume B, L=50 (S, K, P)"   = c("Salinity", "Kd490", "plumeB_50"),
       "2s structural + plume Q, L=50 (S, σk, P)" = c("Salinity", "sigk", "plumeQ_50")))
ALL_COLS <- unique(c("lon", "lat", "SST", unlist(lapply(MODELS, `[[`, 1)), unlist(MODELS5)))

fit_mx <- function(y, X) {
  m <- try(maxnet::maxnet(y, X, maxnet::maxnet.formula(y, X, classes = "lqph"),
                          regmult = REGMULT), silent = TRUE)
  if (inherits(m, "try-error")) NULL else m
}
auc_of <- function(y, p) as.numeric(pROC::auc(y, p, quiet = TRUE))
okrows <- function(d) d[stats::complete.cases(d[, ALL_COLS]) &
                        apply(as.matrix(d[, setdiff(ALL_COLS, "lon")]), 1, \(x) all(is.finite(x))), ]
draw_rep <- function(tp, r, match_climate) {
  set.seed(1000L * r + match(tp, ECO_LVLS))
  pres <- okrows(bc[bc$Type == tp, ]); bgp <- okrows(bg_df)
  if (match_climate) {
    sq <- quantile(pres$SST, c(0.05, 0.95)); lq <- quantile(pres$lat, c(0.05, 0.95))
    m  <- bgp$SST >= sq[1] & bgp$SST <= sq[2] & bgp$lat >= lq[1] & bgp$lat <= lq[2]
    if (sum(m) >= 100) bgp <- bgp[m, ]
  }
  nbg <- min(nrow(pres) * 10, nrow(bgp), 6000)
  df  <- rbind(cbind(pres[sample(nrow(pres), replace = TRUE), ALL_COLS], presence = 1L),
               cbind(bgp[sample(nrow(bgp), nbg, replace = TRUE),   ALL_COLS], presence = 0L))
  lon_r    <- (df$lon + 180 + runif(1, 0, 360)) %% 360
  df$block <- as.integer(cut(lon_r, unique(quantile(lon_r, 0:4 / 4)), include.lowest = TRUE))
  df
}
cv_auc <- function(df, feats) {
  a <- sapply(sort(unique(df$block)), function(b) {
    tr <- df[df$block != b, ]; te <- df[df$block == b, ]
    if (length(unique(te$presence)) < 2) return(NA)
    m <- fit_mx(tr$presence, tr[, feats, drop = FALSE]); if (is.null(m)) return(NA)
    auc_of(te$presence, as.numeric(predict(m, te[, feats, drop = FALSE], type = "cloglog")))
  })
  mean(a, na.rm = TRUE)
}
jobs <- expand.grid(r = seq_len(N_REP), tp = ECO_LVLS, stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1")
t0 <- Sys.time()
res <- dplyr::bind_rows(parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  tp <- jobs$tp[i]; r <- jobs$r[i]
  d2 <- draw_rep(tp, r, TRUE); d1 <- draw_rep(tp, r, FALSE)
  dplyr::bind_rows(lapply(names(MODELS), function(mn) {
    m <- MODELS[[mn]]
    data.frame(Type = tp, rep = r, model = mn, auc = cv_auc(if (m[[2]]) d2 else d1, m[[1]]))
  }))
}, mc.cores = N_CORES, mc.preschedule = FALSE))
message(sprintf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
write.csv(res, file.path(OUT_DIR, "TabS_struct_sdm_reps.csv"), row.names = FALSE)

q <- function(x) c(median = median(x, na.rm = TRUE), lo = unname(quantile(x, .025, na.rm = TRUE)),
                   hi = unname(quantile(x, .975, na.rm = TRUE)), p_le0 = mean(x <= 0, na.rm = TRUE))
auc_s <- res |> dplyr::group_by(Type, model) |>
  dplyr::summarise(median = median(auc, na.rm = TRUE), lo = quantile(auc, .025, na.rm = TRUE),
                   hi = quantile(auc, .975, na.rm = TRUE), .groups = "drop")
wide <- tidyr::pivot_wider(res, names_from = model, values_from = auc)
CONTRASTS <- list(
  c("2s structural (S, σk)",      "2s current (S, K, logD, Q)", "structural − current (2-stage)"),
  c("2s structural + exp(-D/L)",  "2s structural (S, σk)",      "D beyond S (2-stage, exp)"),
  c("2s structural + logD",       "2s structural (S, σk)",      "D beyond S (2-stage, log)"),
  c("1s structural (T, S, σk)",   "1s current (T, lat, S, K, logD, Q)", "structural − current (1-stage)"),
  c("1s structural + exp(-D/L)",  "1s structural (T, S, σk)",   "D beyond S (1-stage, exp)"))
con <- dplyr::bind_rows(lapply(CONTRASTS, function(cn) dplyr::bind_rows(lapply(ECO_LVLS, function(tp) {
  x <- wide[wide$Type == tp, ]
  data.frame(Type = tp, contrast = cn[3], t(q(x[[cn[1]]] - x[[cn[2]]])))
}))))
write.csv(auc_s, file.path(OUT_DIR, "TabS_struct_sdm_auc.csv"), row.names = FALSE)
write.csv(con,   file.path(OUT_DIR, "TabS_struct_sdm_contrasts.csv"), row.names = FALSE)
message("\n  CV AUC (median [95%]):")
print(as.data.frame(auc_s |> dplyr::mutate(txt = sprintf("%.3f [%.3f, %.3f]", median, lo, hi)) |>
  dplyr::select(model, Type, txt) |> tidyr::pivot_wider(names_from = Type, values_from = txt)))
message("\n  paired contrasts (median [95%], p≤0):")
print(as.data.frame(con |> dplyr::mutate(txt = sprintf("%+.3f [%+.3f, %+.3f] %.2f", median, lo, hi, p_le0)) |>
  dplyr::select(contrast, Type, txt) |> tidyr::pivot_wider(names_from = Type, values_from = txt)))

# ═════════════════════════════════════════════════════════════════════════════
# (4) Catchment area's contribution through salinity
# ═════════════════════════════════════════════════════════════════════════════
# With ΔS ∝ B^e at fixed D (e = β / ln 10, β per log10 km²), scaling B by k
# scales the freshwater deficit by k^e. Fit the 2-stage structural model on all data, move every presence's
# salinity to S_ref − ΔS·k^β, and report the mean relative change in predicted
# suitability. This is B's effect routed through S, which is what the
# structural model allows it to have.
message("\n--- (4) catchment area via salinity ---")
b_eff <- dplyr::bind_rows(lapply(ECO_LVLS, function(tp) {
  pres <- okrows(bc[bc$Type == tp, ])
  bgp  <- okrows(bg_df)
  sq <- quantile(pres$SST, c(.05, .95)); lq <- quantile(pres$lat, c(.05, .95))
  bgp <- bgp[bgp$SST >= sq[1] & bgp$SST <= sq[2] & bgp$lat >= lq[1] & bgp$lat <= lq[2], ]
  feats <- MODELS[["2s structural (S, σk)"]][[1]]
  y <- c(rep(1L, nrow(pres)), rep(0L, nrow(bgp)))
  X <- rbind(pres[, feats], bgp[, feats])
  m <- fit_mx(y, X)
  base <- as.numeric(predict(m, pres[, feats], type = "cloglog"))
  dplyr::bind_rows(lapply(c(0.5, 2, 10), function(k) {
    # only points fresher than the reference carry a freshwater deficit to scale
    p2 <- pres; dS <- p2$S_ref2 - p2$Salinity; fr <- is.finite(dS) & dS > 0
    p2$Salinity[fr] <- p2$S_ref2[fr] - dS[fr] * exp(BETA * log10(k))   # ΔS ∝ B^(β/ln10)
    new <- as.numeric(predict(m, p2[, feats], type = "cloglog"))
    data.frame(Type = tp, B_factor = k, dS_factor = exp(BETA * log10(k)),
               mean_dSal = mean(p2$Salinity - pres$Salinity, na.rm = TRUE),
               mean_rel_change_pct = 100 * (mean(new, na.rm = TRUE) / mean(base, na.rm = TRUE) - 1))
  }))
}))
print(b_eff |> dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, 3))))
write.csv(b_eff, file.path(OUT_DIR, "TabS_struct_B_via_S.csv"), row.names = FALSE)
# ═════════════════════════════════════════════════════════════════════════════
# (5) Plume index in place of distance + discharge
# ═════════════════════════════════════════════════════════════════════════════
message("\n--- (5) plume index P = log10(1 + size · exp(−D/L)) ---")
t0 <- Sys.time()
res5 <- dplyr::bind_rows(parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  tp <- jobs$tp[i]; r <- jobs$r[i]
  d2 <- draw_rep(tp, r, TRUE)
  dplyr::bind_rows(lapply(names(MODELS5), function(mn)
    data.frame(Type = tp, rep = r, model = mn, auc = cv_auc(d2, MODELS5[[mn]]))))
}, mc.cores = N_CORES, mc.preschedule = FALSE))
message(sprintf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
write.csv(res5, file.path(OUT_DIR, "TabS_struct_plume_reps.csv"), row.names = FALSE)
w5 <- tidyr::pivot_wider(res5, names_from = model, values_from = auc)
con5 <- dplyr::bind_rows(lapply(setdiff(names(MODELS5), "2s current (S, K, logD, Q)"), function(mn)
  dplyr::bind_rows(lapply(ECO_LVLS, function(tp) {
    x <- w5[w5$Type == tp, ]
    data.frame(Type = tp, model = mn, t(q(x[[mn]] - x[["2s current (S, K, logD, Q)"]])))
  }))))
auc5 <- res5 |> dplyr::group_by(Type, model) |>
  dplyr::summarise(median = median(auc, na.rm = TRUE), lo = quantile(auc, .025, na.rm = TRUE),
                   hi = quantile(auc, .975, na.rm = TRUE), .groups = "drop")
write.csv(auc5, file.path(OUT_DIR, "TabS_struct_plume_auc.csv"), row.names = FALSE)
write.csv(con5, file.path(OUT_DIR, "TabS_struct_plume_contrasts.csv"), row.names = FALSE)
message("  CV AUC minus current (median [95%], p≤0):")
print(as.data.frame(con5 |> dplyr::mutate(txt = sprintf("%+.3f [%+.3f, %+.3f] %.2f", median, lo, hi, p_le0)) |>
  dplyr::select(model, Type, txt) |> tidyr::pivot_wider(names_from = Type, values_from = txt)))

message("\n=== 04_structural complete ===")
