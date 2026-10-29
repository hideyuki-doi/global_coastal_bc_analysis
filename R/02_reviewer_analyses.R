#!/usr/bin/env Rscript
# =============================================================================
# PART 8 — Supplementary analyses added in revision
#
# Sourced from the end of blue_carbon_full.R; uses bc, bg_df, ECO_LVLS,
# assign_blk(), FEAT, FEAT_range, OUT_DIR and the plotting helpers as they
# stand at that point.
#
# Four questions, each raised in review:
#
#   8a  Turbidity mixes two processes. Kd490 at a coastal cell is the sum of
#       local ocean optics (phytoplankton, resuspension — both strongly
#       latitudinal) and terrestrial delivery (CDOM, mineral sediment). The
#       first version could not say which half carried the signal, which made
#       the predictor hard to interpret. Here Kd490 is decomposed against
#       climate and against the river terms, and the SDM is refitted on the
#       land-attributable part alone.
#
#   8b  Distance functional form. v1 used log10(d + 1) in the SDM but exp(-d/L)
#       in the flow-regulation scenarios — two different claims about the same
#       process. A plume diluting radially from a point source argues for 1/d.
#       All three are fitted and compared under the same spatial-block CV.
#
#   8c  River size. Catchment area treats an arid basin and a wet one of the
#       same size alike. Discharge integrates area with runoff. Both are
#       available from HydroRIVERS; this compares them head to head.
#
#   8d  Does the two-stage design earn its complexity? If decomposing turbidity
#       removes the climate signal, a single stage might suffice. Fitted both
#       ways, with ocean-wide and with climate-matched background.
#
# Outputs: output/TabS*_*.csv and output/FigR*_*.{pdf,png}
# =============================================================================

message("\n\n=== PART 8: supplementary analyses (revision) ===")

if (!exists("bc") || !exists("bg_df"))
  stop("PART 8 expects bc and bg_df in the environment; source this from blue_carbon_full.R")

suppressPackageStartupMessages({ library(mgcv); library(dplyr) })

set.seed(42)
R8 <- list()   # collected results

# ─────────────────────────────────────────────────────────────────────────────
# Shared helper: fit stage-2 MaxEnt for one ecosystem and return spatial-block
# CV AUC plus permutation importance. Mirrors PART 5 exactly (climate-matched
# background, regmult 2, lqph features, longitudinal-quartile blocks) so the
# numbers here are comparable with the main results.
# ─────────────────────────────────────────────────────────────────────────────
fit_stage2 <- function(tp, feats, data_pres = bc, data_bg = bg_df,
                       match_climate = TRUE, regmult = 2.0) {
  need <- unique(c("lon", "lat", "SST", feats))
  pres <- data_pres |> dplyr::filter(Type == tp) |>
    dplyr::select(dplyr::all_of(need)) |> dplyr::mutate(presence = 1L)

  bgp <- data_bg
  if (match_climate) {
    sq <- quantile(pres$SST, c(0.05, 0.95), na.rm = TRUE)
    lq <- quantile(pres$lat, c(0.05, 0.95), na.rm = TRUE)
    bgp <- dplyr::filter(bgp, SST >= sq[1], SST <= sq[2], lat >= lq[1], lat <= lq[2])
    if (nrow(bgp) < 100) bgp <- data_bg
  }
  nbg <- min(nrow(pres) * 10, nrow(bgp), 6000)
  bgs <- bgp |> dplyr::select(dplyr::all_of(need)) |>
    dplyr::mutate(presence = 0L) |> dplyr::slice_sample(n = nbg)

  df <- dplyr::bind_rows(pres, bgs) |>
    dplyr::filter(dplyr::if_all(dplyr::all_of(feats), is.finite)) |>
    dplyr::mutate(block = assign_blk(lon))
  if (nrow(df) < 50 || length(unique(df$presence)) < 2) return(NULL)

  X <- as.data.frame(df[, feats, drop = FALSE])
  keep <- names(which(apply(X, 2, var) > 1e-10))
  if (length(keep) == 0) return(NULL)
  X <- X[, keep, drop = FALSE]; y <- df$presence

  mod <- maxnet::maxnet(y, X, maxnet::maxnet.formula(y, X, classes = "lqph"),
                        regmult = regmult)

  # spatial-block CV
  auc <- rep(NA_real_, 4)
  for (b in 1:4) {
    tr <- df[df$block != b, ]; te <- df[df$block == b, ]
    if (length(unique(te$presence)) < 2 || nrow(te) < 10) next
    Xtr <- as.data.frame(tr[, keep, drop = FALSE])
    Xte <- as.data.frame(te[, keep, drop = FALSE])
    if (any(apply(Xtr, 2, var) <= 1e-10)) next
    m <- try(maxnet::maxnet(tr$presence, Xtr,
             maxnet::maxnet.formula(tr$presence, Xtr, classes = "lqph"),
             regmult = regmult), silent = TRUE)
    if (inherits(m, "try-error")) next
    p <- as.numeric(predict(m, Xte, type = "cloglog"))
    auc[b] <- as.numeric(pROC::auc(te$presence, p, quiet = TRUE))
  }

  # permutation importance: AUC lost when a variable is shuffled, averaged over
  # N_PERM shuffles. Kept both raw (ΔAUC) and as a share of the total. The share
  # depends on how many other predictors compete for it — a 6-variable
  # one-stage model and a 4-variable stage-2 model split 100% differently even
  # if distance carries the same signal — so design comparisons (8d) should be
  # read on the raw ΔAUC.
  N_PERM <- 10
  base <- as.numeric(pROC::auc(y, as.numeric(predict(mod, X, type = "cloglog")),
                               quiet = TRUE))
  imp_raw <- sapply(keep, function(v) {
    mean(replicate(N_PERM, {
      Xp <- X; Xp[[v]] <- sample(Xp[[v]])
      max(0, base - as.numeric(pROC::auc(y,
        as.numeric(predict(mod, Xp, type = "cloglog")), quiet = TRUE)))
    }))
  })
  imp <- imp_raw / max(sum(imp_raw), 1e-9)
  full     <- setNames(rep(0, length(feats)), feats); full[keep]     <- imp
  full_raw <- setNames(rep(0, length(feats)), feats); full_raw[keep] <- imp_raw

  list(Type = tp, auc_mean = mean(auc, na.rm = TRUE),
       auc_sd = sd(auc, na.rm = TRUE), auc_train = base,
       n_pres = sum(y), importance = full, importance_raw = full_raw)
}

imp_table <- function(fits, label) {
  dplyr::bind_rows(lapply(Filter(Negate(is.null), fits), function(f)
    tibble::tibble(model = label, Type = f$Type,
                   Variable = names(f$importance),
                   Importance = as.numeric(f$importance),
                   dAUC = as.numeric(f$importance_raw),
                   AUC_cv = f$auc_mean, AUC_cv_sd = f$auc_sd,
                   AUC_train = f$auc_train, n_pres = f$n_pres)))
}

# ═════════════════════════════════════════════════════════════════════════════
# 8a — TURBIDITY DECOMPOSITION
# ═════════════════════════════════════════════════════════════════════════════
message("\n--- 8a: decomposing Kd490 into ocean-climate and river components ---")

# Fitted over presences AND background together: that is the domain the SDM
# actually sees, so the decomposition describes the turbidity variation the
# model has to work with, not just the variation at occupied sites.
kd_dat <- dplyr::bind_rows(
    bc    |> dplyr::select(Kd490, SST, abs_lat, log_dist, log_discharge, log_upstream),
    bg_df |> dplyr::select(Kd490, SST, abs_lat, log_dist, log_discharge, log_upstream)
  ) |>
  dplyr::filter(Kd490 > 0, dplyr::if_all(dplyr::everything(), is.finite)) |>
  dplyr::mutate(log_Kd = log(Kd490))

message(sprintf("  n = %s coastal points", format(nrow(kd_dat), big.mark = ",")))

K <- 8
m_full  <- mgcv::gam(log_Kd ~ s(SST, k = K) + s(abs_lat, k = K) +
                              s(log_dist, k = K) + s(log_discharge, k = K), data = kd_dat)
m_ocean <- mgcv::gam(log_Kd ~ s(SST, k = K) + s(abs_lat, k = K), data = kd_dat)
m_river <- mgcv::gam(log_Kd ~ s(log_dist, k = K) + s(log_discharge, k = K), data = kd_dat)

dev <- function(m) summary(m)$dev.expl
part <- tibble::tibble(
  component = c("Total (climate + river)", "Climate alone (SST + |lat|)",
                "River alone (distance + discharge)", "Unique to climate",
                "Unique to river", "Shared climate/river"),
  dev_expl  = c(dev(m_full), dev(m_ocean), dev(m_river),
                dev(m_full) - dev(m_river), dev(m_full) - dev(m_ocean),
                dev(m_ocean) + dev(m_river) - dev(m_full))
)
print(as.data.frame(part |> dplyr::mutate(pct = round(100 * dev_expl, 1))))
write.csv(part, file.path(OUT_DIR, "TabS_Kd490_variance_partition.csv"), row.names = FALSE)
R8$kd_partition <- part

# --- the residual model, as suggested in review -----------------------------
# Kd490 residual after removing the ocean-climate signal, regressed on the two
# river terms. The response is on the LOG scale, so this additive model is
# already the multiplicative one: Kd_land proportional to
# discharge^b1 x (dist + 1)^b2. The coefficients are therefore read as elasticities —
# a tenfold change in discharge multiplies land-attributable turbidity by 10^b1.
kd_dat$kd_resid <- as.numeric(residuals(m_ocean, type = "response"))
m_res <- lm(kd_resid ~ log_discharge + log_dist, data = kd_dat)
co <- summary(m_res)$coefficients
resid_tab <- tibble::tibble(
  term       = rownames(co),
  estimate   = co[, 1], se = co[, 2], t = co[, 3], p = co[, 4],
  multiplicative_per_10x = ifelse(rownames(co) == "(Intercept)", NA, 10^co[, 1])
)
message(sprintf("  residual model R2 = %.3f", summary(m_res)$r.squared))
print(as.data.frame(resid_tab))
write.csv(resid_tab, file.path(OUT_DIR, "TabS_Kd490_residual_model.csv"), row.names = FALSE)
R8$kd_residual <- resid_tab

# Make the land-attributable turbidity available to the SDM comparisons below.
kd_ocean_fit  <- function(d) as.numeric(predict(m_ocean, newdata = d))
bc$kd_resid    <- log(bc$Kd490)    - kd_ocean_fit(bc)
bg_df$kd_resid <- log(bg_df$Kd490) - kd_ocean_fit(bg_df)

# ═════════════════════════════════════════════════════════════════════════════
# 8b — DISTANCE FUNCTIONAL FORM
# ═════════════════════════════════════════════════════════════════════════════
message("\n--- 8b: distance functional form (log vs inverse vs exponential) ---")

dist_sets <- list(
  log_dist = c("Salinity", "Kd490", "log_dist", "log_discharge"),
  inv_dist = c("Salinity", "Kd490", "inv_dist", "log_discharge"),
  exp_dist = c("Salinity", "Kd490", "exp_dist", "log_discharge")
)
dist_res <- dplyr::bind_rows(lapply(names(dist_sets), function(nm) {
  message("  form: ", nm)
  imp_table(lapply(ECO_LVLS, fit_stage2, feats = dist_sets[[nm]]), nm)
}))
write.csv(dist_res, file.path(OUT_DIR, "TabS_distance_form_comparison.csv"), row.names = FALSE)
R8$dist_form <- dist_res

dist_auc <- dist_res |> dplyr::distinct(model, Type, AUC_cv, AUC_cv_sd)
message("\n  spatial-block CV AUC by distance form:")
print(as.data.frame(tidyr::pivot_wider(dist_auc, names_from = model,
                                       values_from = c(AUC_cv, AUC_cv_sd))))

best_form <- dist_auc |> dplyr::group_by(model) |>
  dplyr::summarise(mean_auc = mean(AUC_cv, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(mean_auc))
message("\n  mean CV AUC across ecosystems:")
print(as.data.frame(best_form))
message(sprintf("  => best distance form: %s", best_form$model[1]))
R8$best_dist_form <- best_form$model[1]

# ═════════════════════════════════════════════════════════════════════════════
# 8c — RIVER SIZE: CATCHMENT AREA vs DISCHARGE
# ═════════════════════════════════════════════════════════════════════════════
message("\n--- 8c: river size predictor (catchment area vs discharge) ---")

size_sets <- list(
  area      = c("Salinity", "Kd490", DIST_VAR, "log_upstream"),
  discharge = c("Salinity", "Kd490", DIST_VAR, "log_discharge"),
  both      = c("Salinity", "Kd490", DIST_VAR, "log_upstream", "log_discharge"),
  none      = c("Salinity", "Kd490", DIST_VAR)
)
size_res <- dplyr::bind_rows(lapply(names(size_sets), function(nm) {
  message("  set: ", nm)
  imp_table(lapply(ECO_LVLS, fit_stage2, feats = size_sets[[nm]]), nm)
}))
write.csv(size_res, file.path(OUT_DIR, "TabS_riversize_comparison.csv"), row.names = FALSE)
R8$river_size <- size_res

message("\n  spatial-block CV AUC by river-size predictor:")
print(as.data.frame(size_res |> dplyr::distinct(model, Type, AUC_cv) |>
                    tidyr::pivot_wider(names_from = model, values_from = AUC_cv)))

# correlation between the two measures, to say how much they actually differ
rs_cor <- cor(bc$log_upstream, bc$log_discharge, use = "complete.obs")
message(sprintf("  Pearson r(log area, log discharge) = %.3f", rs_cor))
R8$area_discharge_r <- rs_cor

# ═════════════════════════════════════════════════════════════════════════════
# 8d — ONE STAGE vs TWO STAGES
# ═════════════════════════════════════════════════════════════════════════════
message("\n--- 8d: one-stage vs two-stage ---")

# (i)   one stage, ocean-wide background, climate + habitat together: the
#       conventional global SDM, and the design we argue hides the river signal
# (ii)  one stage, ocean-wide background, turbidity replaced by its
#       land-attributable residual: the version suggested in review
# (iii) stage 2 of the two-stage model, climate-matched background: the paper
one_feats  <- c("SST", "lat", "Salinity", "Kd490", DIST_VAR, "log_discharge")
one_resid  <- c("SST", "lat", "Salinity", "kd_resid", DIST_VAR, "log_discharge")
two_feats  <- c("Salinity", "Kd490", DIST_VAR, "log_discharge")
two_resid  <- c("Salinity", "kd_resid", DIST_VAR, "log_discharge")

stage_res <- dplyr::bind_rows(
  imp_table(lapply(ECO_LVLS, fit_stage2, feats = one_feats, match_climate = FALSE),
            "1-stage, ocean-wide bg"),
  imp_table(lapply(ECO_LVLS, fit_stage2, feats = one_resid, match_climate = FALSE),
            "1-stage, Kd residual"),
  imp_table(lapply(ECO_LVLS, fit_stage2, feats = two_feats, match_climate = TRUE),
            "2-stage (stage 2)"),
  imp_table(lapply(ECO_LVLS, fit_stage2, feats = two_resid, match_climate = TRUE),
            "2-stage, Kd residual")
)
write.csv(stage_res, file.path(OUT_DIR, "TabS_onestage_vs_twostage.csv"), row.names = FALSE)
R8$stages <- stage_res

dist_share <- stage_res |>
  dplyr::filter(Variable %in% c(DIST_VAR, "log_dist", "inv_dist", "exp_dist")) |>
  dplyr::select(model, Type, Importance) |>
  tidyr::pivot_wider(names_from = Type, values_from = Importance)
message("\n  importance share of river distance, by design:")
print(as.data.frame(dist_share))
R8$dist_share <- dist_share

# The same comparison on the raw scale: AUC lost when distance is shuffled.
# Unlike the share, this does not shrink just because SST and lat are added.
dist_dauc <- stage_res |>
  dplyr::filter(Variable %in% c(DIST_VAR, "log_dist", "inv_dist", "exp_dist")) |>
  dplyr::select(model, Type, dAUC) |>
  tidyr::pivot_wider(names_from = Type, values_from = dAUC)
message("\n  AUC lost when river distance is permuted (raw ΔAUC), by design:")
print(as.data.frame(dist_dauc))
R8$dist_dauc <- dist_dauc

# ═════════════════════════════════════════════════════════════════════════════
# FIGURES
# ═════════════════════════════════════════════════════════════════════════════
message("\n--- 8e: figures ---")

# FigR1 — what Kd490 is made of, and how the river part behaves
kd_plot <- kd_dat |> dplyr::slice_sample(n = min(6000, nrow(kd_dat)))
p_a <- ggplot(part |> dplyr::slice(4:6) |>
                dplyr::mutate(component = factor(component, levels = component)),
              aes(component, 100 * dev_expl)) +
  geom_col(fill = "grey45", width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * dev_expl)), vjust = -0.4, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Deviance explained in log Kd490 (%)",
       title = "a  Turbidity: climate vs river") +
  theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_b <- ggplot(kd_plot, aes(10^log_dist - 1, kd_resid)) +
  geom_point(size = 0.25, alpha = 0.12, shape = 16) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), colour = "#0072B2",
              fill = "#0072B2", alpha = 0.2, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50", linewidth = 0.3) +
  scale_x_log10(breaks = c(1, 10, 100, 1000), labels = scales::label_comma()) +
  annotation_logticks(sides = "b", linewidth = 0.2) +
  labs(x = "Distance from river mouth (km)",
       y = "Land-attributable turbidity\n(log Kd490 residual)",
       title = "b  River signal after removing climate") +
  theme_pub()

save_fig(p_a | p_b, "FigR1_Turbidity_Decomposition", 180, 80)

# FigR2 — distance form and design comparison
p_c <- ggplot(dist_auc, aes(model, AUC_cv, colour = Type, group = Type)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  geom_pointrange(aes(ymin = AUC_cv - AUC_cv_sd, ymax = AUC_cv + AUC_cv_sd),
                  size = 0.35) +
  scale_colour_manual(values = ECO_PAL, name = NULL) +
  scale_x_discrete(labels = c(log_dist = "log(d+1)", inv_dist = "1/(d+1)",
                              exp_dist = "exp(-d/50)")) +
  labs(x = "Distance transform", y = "Spatial-block CV AUC",
       title = "a  Distance functional form") +
  theme_pub() + theme(legend.position = "bottom")

ds_long <- dist_share |>
  tidyr::pivot_longer(-model, names_to = "Type", values_to = "Importance") |>
  dplyr::mutate(Type = factor(Type, levels = ECO_LVLS))
p_d <- ggplot(ds_long, aes(model, 100 * Importance, fill = Type)) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  scale_fill_manual(values = ECO_PAL, name = NULL) +
  labs(x = NULL, y = "Importance share of river distance (%)",
       title = "b  Why the range model is needed") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")

save_fig(p_c | p_d, "FigR2_Design_Comparison", 190, 95)

saveRDS(R8, file.path(OUT_DIR, "part8_results.rds"))
message("\n=== PART 8 complete — tables in ", OUT_DIR, "/TabS_*.csv ===")
