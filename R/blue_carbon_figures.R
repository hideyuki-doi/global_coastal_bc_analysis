# =============================================================================
# Blue Carbon × River Proximity — FIGURE GENERATION (standalone)
# Run blue_carbon_full.R at least once first to create output/rds/ checkpoints.
# Then source this file to regenerate / tweak all figures without re-running
# the heavy SDM / scenario computation.
# =============================================================================

# ── 0. Packages ───────────────────────────────────────────────────────────────
library(dplyr);  library(tidyr);  library(stringr); library(tibble); library(purrr)
library(ggplot2); library(patchwork); library(ggridges); library(scales); library(ggh4x)
library(sf);     library(maxnet); library(pROC); library(mgcv)
# shapviz is optional: only FigS2b (a discarded robustness check, not used in the
# manuscript) needs it, and its SHAP input comes from the archived fastshap.
if (requireNamespace("shapviz", quietly = TRUE)) library(shapviz)
library(rnaturalearth); library(rnaturalearthdata)

# ── 1. Paths ──────────────────────────────────────────────────────────────────
OUT_DIR <- "output"          # same as in the analysis script
RDS_DIR <- file.path(OUT_DIR, "rds")
if (!dir.exists(RDS_DIR))
  stop("RDS checkpoints not found. Run blue_carbon_full.R first.\n",
       "Expected: ", RDS_DIR)

# ── 2. Shared constants & theme ───────────────────────────────────────────────
ECO_PAL  <- c(Seagrass = "#009E73", Mangrove = "#E69F00", Kelp = "#0072B2")
ECO_LVLS <- c("Seagrass","Mangrove","Kelp")
eco_strip <- function() {
  ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(fill = alpha(unname(ECO_PAL[ECO_LVLS]), 0.22)),
    text_x       = ggh4x::elem_list_text(colour = unname(ECO_PAL[ECO_LVLS]), face = "bold")
  )
}

# Two-stage SDM: range model (SST+lat) × habitat model (local factors)
FEAT_range <- c("SST","lat")
FEAT       <- c("Salinity","Kd490","log_dist","log_upstream")
FEAT_SP    <- FEAT
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
  log_city     = "log₁₀(Dist. to city + 1) (km)",
  dist_city_km = "Dist. to nearest city (km)",
  log_tidal    = "log₁₀(Tidal range + 1) (m)",
  tidal_m      = "Tidal range (m)"
)

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

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p,
         width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p,
         width = w, height = h, units = "mm", dpi = 300, bg = "white")
  message("Saved: ", name)
}

# ── 3. Load analysis results ──────────────────────────────────────────────────
rds_load <- function(fname) {
  path <- file.path(RDS_DIR, fname)
  if (!file.exists(path)) { warning("Not found: ", fname); return(NULL) }
  readRDS(path)
}

bc           <- rds_load("bc.rds")
bg_df        <- rds_load("bg_df.rds")
bg_coastal   <- rds_load("bg_coastal.rds")
coast_thresh <- rds_load("coast_thresh.rds")
all_mx       <- rds_load("all_mx.rds")
range_auc    <- rds_load("range_auc.rds")
all_auc      <- rds_load("all_auc.rds")
all_vimp     <- rds_load("all_vimp.rds")
all_shap     <- rds_load("all_shap.rds")
all_roc      <- rds_load("all_roc.rds")
all_pdp      <- rds_load("all_pdp.rds")
RF_SHAP      <- rds_load("RF_SHAP.rds")
sp_models    <- rds_load("sp_models.rds")
sp_vimp         <- rds_load("sp_vimp.rds")
sp_vimp_thinned <- rds_load("sp_vimp_thinned.rds")
sp_roc_list     <- rds_load("sp_roc_list.rds")
top3_per_type <- rds_load("top3_per_type.rds")
sp_col        <- rds_load("sp_col.rds")
sp_pd_summary <- rds_load("sp_pd_summary.rds")
sc_results    <- rds_load("sc_results.rds")
sp_sc_grad    <- rds_load("sp_sc_grad.rds")

# world_grid (large; loaded from separate cache)
GRID_CACHE <- file.path(OUT_DIR, "world_grid.rds")
world_grid <- if (file.exists(GRID_CACHE)) readRDS(GRID_CACHE) else NULL
# Restore derived columns if missing from older cache
if (!is.null(world_grid)) {
  if (!"dist_km_baseline" %in% names(world_grid))
    world_grid$dist_km_baseline <- world_grid$dist_km
  if (!"depth" %in% names(world_grid)) {
    BATHY_RAST_F <- file.path(OUT_DIR, "depth_rast.tif")
    if (file.exists(BATHY_RAST_F)) {
      depth_rast_fig <- terra::rast(BATHY_RAST_F)
      wg_v <- terra::vect(as.matrix(world_grid[, c("lon","lat")]), crs = "EPSG:4326")
      wg_depth_raw <- terra::extract(depth_rast_fig, wg_v)[, 2]
      # abs(): MS_bathy_5m uses negative values for ocean (ETOPO convention)
      world_grid$depth <- pmax(0, replace(abs(wg_depth_raw), !is.finite(wg_depth_raw), 0))
      message("  depth added to world_grid from cached raster")
    } else {
      world_grid$depth <- 0  # placeholder; re-run full.R to get real depth
      warning("depth_rast.tif not found — world_grid$depth set to 0. Re-run full script.")
    }
  }
  # Apply same depth domain restriction as full.R (Kelp max ~40 m; ceiling 50 m)
  MAX_DEPTH_M <- 100
  if ("depth" %in% names(world_grid))
    world_grid <- world_grid |> dplyr::filter(depth <= MAX_DEPTH_M)
}

# basemap (tiny package data)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

message("RDS files loaded. Generating figures...")

# =============================================================================
# FIGURE A1 — Predictor correlation matrix (Appendix)
# =============================================================================
if (!is.null(bc)) {
  CORR_LABELS <- c(Salinity     = "Salinity",
                   Kd490        = "Kd490",
                   log_dist     = "log₁₀(Dist+1)",
                   log_upstream = "log₁₀(Wshed+1)")

  cor_mat <- bc |>
    dplyr::select(all_of(FEAT)) |>
    cor(use = "pairwise.complete.obs")

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
}

# =============================================================================
# FIGURE S2 — Variable importance: per-species violin (MaxEnt) +
#              type-level RF-SHAP as diamond markers
# =============================================================================
if (!is.null(sp_vimp)) {
  # Type-level RF-SHAP (reference diamonds)
  rf_shap_type <- dplyr::bind_rows(all_vimp) |>
    tidyr::pivot_longer(dplyr::any_of("RF_SHAP"),
                        names_to = "Model", values_to = "Importance") |>
    dplyr::filter(!is.na(Importance)) |>
    dplyr::mutate(
      Variable = factor(VAR_LABELS[Variable], levels = rev(VAR_LABELS[FEAT])),
      Type     = factor(Type, levels = ECO_LVLS)
    )

  fig_s2 <- ggplot(sp_vimp,
                   aes(Importance, Variable, colour = Type, fill = Type)) +
    geom_jitter(height = 0.18, size = 0.9, alpha = 0.45, shape = 16) +
    geom_boxplot(width = 0.38, fill = NA, outlier.shape = NA, linewidth = 0.5) +
    # RF-SHAP reference diamonds at type level
    geom_point(data = rf_shap_type,
               aes(Importance, Variable),
               shape = 23, size = 3.0, colour = "white",
               fill = ECO_PAL[as.character(rf_shap_type$Type)],
               inherit.aes = FALSE) +
    ggh4x::facet_wrap2(~Type, ncol = 3, strip = eco_strip()) +
    scale_colour_manual(values = ECO_PAL, guide = "none") +
    scale_fill_manual(values   = ECO_PAL, guide = "none") +
    scale_x_continuous(labels = scales::label_percent(),
                       expand = expansion(mult = c(0.02, 0.05))) +
    labs(x = "Importance share (%)", y = NULL,
         caption = "Points + boxes: MaxEnt per species (n ≥ 30) | ◆: RF balanced SHAP (type level)") +
    theme_pub()
  save_fig(fig_s2, "FigS2_VarImp", 170, 100)
}

# =============================================================================
# FIGURE S2b — SHAP beeswarm (type-level RF)
# =============================================================================
if (!is.null(RF_SHAP) && is.list(RF_SHAP) && length(RF_SHAP) > 0) {
  s2b_panels <- purrr::imap(RF_SHAP, \(sv, tp) {
    if (is.null(sv)) return(NULL)
    shapviz::sv_importance(sv, kind = "beeswarm", show_numbers = FALSE,
                           color_bar_title = "Feature value") +
      scale_x_continuous(n.breaks = 4) +
      labs(x = "SHAP value", title = tp) +
      theme_pub() + theme(legend.position = "right")
  })
  s2b_panels <- Filter(Negate(is.null), s2b_panels)
  if (length(s2b_panels) > 0) {
    fig_s2b <- patchwork::wrap_plots(s2b_panels, ncol = 3) +
      plot_annotation(theme = theme_pub())
    save_fig(fig_s2b, "FigS2b_SHAP_Beeswarm", 200, 100)
  }
}

# =============================================================================
# FIGURE S3 — Response curves with CI (±SD across species) for all FEAT_SP vars
# =============================================================================
if (!is.null(sp_pd_summary)) {
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
        labels = scales::label_number(big.mark = " ")(nice_pow10))
    } else if (back_log10) {
      p <- p + scale_x_continuous(
        labels = \(v) scales::label_number()(round(10^v - 1)))
    }
    p
  }

  p_sal  <- make_pd_panel("Salinity",     "Salinity (PSU)",           log_x = FALSE, show_y = TRUE)
  p_kd   <- make_pd_panel("Kd490",        "Kd490 (m⁻¹)",             log_x = TRUE,  show_y = FALSE)
  p_dist <- make_pd_panel("log_dist",     "Distance from river mouth (km)", log_x = FALSE, show_y = TRUE,
                          back_log10 = TRUE)
  p_up    <- make_pd_panel("log_upstream", "Watershed area (km²)",   log_x = FALSE, show_y = FALSE,
                           back_log10 = TRUE, nice_pow10 = c(10, 100, 1000, 10000, 100000))
  fig_s3 <- ((p_sal | p_kd) / (p_dist | p_up)) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  save_fig(fig_s3, "FigS3_ResponseCurves_CI", 130, 110)
}

# =============================================================================
# FIGURE D4 — Occupancy rate vs distance from river mouth
# Mirrors blue_carbon_full.R exactly.
# =============================================================================
if (!is.null(bc) && !is.null(bg_df)) {
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
}

# =============================================================================
# FIGURE ROC — Cross-validation ROC (type level, bold) + per-species training
#              ROC (thin lines) in background
# =============================================================================
if (!is.null(all_roc)) {
  # Per-species MaxEnt training ROC thin background lines
  sp_roc_bg <- NULL
  if (!is.null(sp_roc_list) && length(sp_roc_list) > 0) {
    sp_roc_bg <- purrr::map_dfr(names(sp_roc_list), \(key) {
      d <- sp_roc_list[[key]]
      roc_obj <- tryCatch(pROC::roc(d$y, d$pred, quiet = TRUE), error = function(e) NULL)
      if (is.null(roc_obj)) return(NULL)
      data.frame(Type = d$Type[1], Species = d$Species[1],
                 FPR  = 1 - rev(roc_obj$specificities),
                 TPR  = rev(roc_obj$sensitivities))
    }) |> dplyr::mutate(Type = factor(Type, levels = ECO_LVLS))
  }

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
      pool <- all_roc[[tp]]$mx_pool
      fb   <- if (!is.null(pool)) as.numeric(pROC::auc(pROC::roc(pool$y, pool$pred, quiet = TRUE))) else NA_real_
      return(data.frame(Type = tp, label = sprintf("AUC = %.3f", fb)))
    }
    data.frame(Type = tp, label = sprintf("AUC = %.3f ± %.3f", mean(vals), sd(vals)))
  }) |> dplyr::mutate(Type = factor(Type, levels = ECO_LVLS))

  fig_roc <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey70", linewidth = 0.4) +
    {if (!is.null(sp_roc_bg) && nrow(sp_roc_bg) > 0)
      geom_line(data = sp_roc_bg,
                aes(FPR, TPR, group = Species, colour = Type),
                alpha = 0.35, linewidth = 0.28)
    } +
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
    scale_x_continuous(labels = scales::label_number(scale = 100),
                       breaks = c(0.25, 0.5, 0.75, 1.0),
                       expand = expansion(mult = 0.02), limits = c(0, 1)) +
    scale_y_continuous(labels = scales::label_number(scale = 100),
                       breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                       expand = expansion(mult = 0.02), limits = c(0, 1)) +
    coord_fixed() +
    ggh4x::facet_wrap2(~Type, ncol = 3, strip = eco_strip()) +
    labs(x = "False positive rate (%)",
         y  = "True positive rate (%)",
         caption = "Thin lines: per-species training ROC  |  Ribbon: ecosystem MaxEnt CV (95% CI, 4-fold spatial block)") +
    theme_pub() + theme(plot.caption = element_text(size = 7, colour = "grey50"))
  save_fig(fig_roc, "FigROC_CrossValidation", 170, 80)
}

# =============================================================================
# FIGURE Sp1 — Species-level importance violin (MaxEnt per species)
# =============================================================================
if (!is.null(sp_vimp_thinned)) {
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
    theme_pub() +
    theme(legend.position = "bottom")
  save_fig(fig_sp1, "FigSp1_SpeciesImportance", 170, 90)
}

# =============================================================================
# FIGURE S4 — Global suitability maps
# =============================================================================
if (!is.null(world_grid) && !is.null(all_mx) && !is.null(coast_thresh)) {
  Xg     <- as.matrix(world_grid[, FEAT]);       colnames(Xg)     <- FEAT
  Xg_rng <- as.matrix(world_grid[, FEAT_range]); colnames(Xg_rng) <- FEAT_range
  range_mx_fig <- if (exists("range_mx")) range_mx
                  else if (file.exists(file.path(OUT_DIR, "rds", "range_mx.rds")))
                    readRDS(file.path(OUT_DIR, "rds", "range_mx.rds"))
                  else list()
  for (tp in ECO_LVLS) {
    suit_hab <- predict(all_mx[[tp]], as.data.frame(Xg), type = "cloglog") |> as.numeric()
    suit_rng <- if (!is.null(range_mx_fig[[tp]]))
      predict(range_mx_fig[[tp]], as.data.frame(Xg_rng), type = "cloglog") |> as.numeric()
    else rep(1, nrow(world_grid))
    world_grid[[paste0("suit_range_", tp)]] <- suit_rng
    world_grid[[paste0("suit_",       tp)]] <- suit_hab * suit_rng
  }
  wg_coast <- world_grid |> dplyr::filter(dist_km_baseline <= coast_thresh)

  # Stage 1: range suitability map
  suit_range_long <- wg_coast |>
    dplyr::select(lon, lat, starts_with("suit_range_")) |>
    tidyr::pivot_longer(starts_with("suit_range_"), names_to = "Type", values_to = "Suitability") |>
    dplyr::mutate(Type = factor(stringr::str_remove(Type, "suit_range_"), levels = ECO_LVLS))
  range_auc_lab <- purrr::map_dfr(ECO_LVLS, \(tp) {
    a <- if (!is.null(range_auc)) range_auc[[tp]] else NULL
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

  # Stage 1 × Stage 2: combined suitability map
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
}

# =============================================================================
# FIGURE Sup_SpeciesMap — occurrence maps for ALL qualifying species (Appendix)
# Column order: Seagrass | Mangrove | Kelp
# =============================================================================
if (!is.null(top3_per_type) && !is.null(bc) && !is.null(sp_col)) {
  top_sp_occ <- top3_per_type |>
    dplyr::group_by(Type) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::arrange(rank, factor(Type, levels = ECO_LVLS))

  sp_map_panels <- purrr::map(seq_len(nrow(top_sp_occ)), \(i) {
    tp       <- top_sp_occ$Type[i]
    spp      <- top_sp_occ$.sp[i]
    type_col <- ECO_PAL[tp]
    pres_sf <- bc |>
      dplyr::rename(.sp = dplyr::all_of(sp_col)) |>
      dplyr::filter(Type == tp, .sp == spp) |>
      sf::st_as_sf(coords = c("lon","lat"), crs = 4326)
    # training AUC of this species' own MaxEnt model, as in FigSup_SpeciesPred
    auc_val   <- if (!is.null(sp_vimp)) dplyr::filter(sp_vimp, Species == spp)$AUC else numeric(0)
    auc_label <- if (length(auc_val) > 0) sprintf("AUC = %.3f", auc_val[1]) else ""
    ggplot() +
      geom_sf(data = world, fill = "grey88", colour = "white", linewidth = 0.08) +
      geom_sf(data = pres_sf, colour = type_col,
              size = 0.35, alpha = 0.7, shape = 16) +
      annotate("text", x = Inf, y = -Inf, label = auc_label,
               hjust = 1.05, vjust = -0.4, size = 2.0, colour = "grey35") +
      coord_sf(ylim = c(-65, 75), expand = FALSE) +
      labs(title = paste("●", spp)) +
      theme_pub() +
      theme(axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size = 7, face = "italic",
                                      colour = type_col,
                                      margin = margin(b = 1)))
  })
  n_occ_rows <- ceiling(nrow(top_sp_occ) / 3)
  fig_ssp_map <- patchwork::wrap_plots(sp_map_panels, ncol = 3)
  save_fig(fig_ssp_map, "FigSup_SpeciesMap", 170, n_occ_rows * 50 + 10)

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
}

# =============================================================================
# FIGURE Sup_SpeciesPred — species-level prediction maps only (all qualifying species)
# (ecosystem-level maps in FigS4 — not repeated here)
# Column = ecosystem type, row = rank within type
# =============================================================================
if (!is.null(sp_models) && !is.null(top3_per_type) && !is.null(world_grid) &&
    !is.null(sp_vimp) && !is.null(coast_thresh)) {
  coast_mask <- world_grid$dist_km_baseline <= coast_thresh
  Xg_sp  <- as.matrix(world_grid[, FEAT_SP]); colnames(Xg_sp) <- FEAT_SP

  map_theme_base <- list(
    theme_pub(),
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), legend.position = "none")
  )

  top_sp_pred <- top3_per_type |>
    dplyr::group_by(Type) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::arrange(rank, factor(Type, levels = ECO_LVLS))

  pred_map_panels <- purrr::map(seq_len(nrow(top_sp_pred)), \(i) {
    tp  <- top_sp_pred$Type[i]
    spp <- top_sp_pred$.sp[i]
    key <- paste(tp, spp, sep = "__")
    if (!key %in% names(sp_models)) return(ggplot() + theme_void())
    suit_hab_sp <- predict(sp_models[[key]], as.data.frame(Xg_sp),
                           type = "cloglog") |> as.numeric()
    suit_rng_sp <- if (!is.null(range_mx_fig[[tp]]))
      predict(range_mx_fig[[tp]], as.data.frame(Xg_rng), type = "cloglog") |> as.numeric()
    else rep(1, nrow(world_grid))
    suit <- suit_hab_sp * suit_rng_sp
    auc_val   <- sp_vimp |> dplyr::filter(Species == spp) |>
      dplyr::pull(AUC) |> unique() |> mean(na.rm = TRUE)
    type_col  <- ECO_PAL[tp]
    auc_label <- if (!is.na(auc_val)) sprintf("AUC = %.3f", auc_val) else ""
    suit_df <- data.frame(lon = world_grid$lon[coast_mask],
                          lat = world_grid$lat[coast_mask],
                          Suitability = suit[coast_mask])
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

  n_sp_rows <- ceiling(nrow(top_sp_pred) / 3)
  fig_ssp_pred <- sp_grid +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom",
          legend.key.width  = unit(1.2, "cm"),
          legend.key.height = unit(0.3, "cm"),
          legend.text = element_text(size = 7))
  save_fig(fig_ssp_pred, "FigSup_SpeciesPred", 170, n_sp_rows * 50 + 25)
}

# =============================================================================
# FIGURE Sc2b — Gradient with ±SD across species
# =============================================================================
if (!is.null(sp_sc_grad)) {
  sp_grad_summ <- sp_sc_grad |>
    dplyr::group_by(Type, reg_pct) |>
    dplyr::summarize(mean_chg = mean(pct_chg, na.rm = TRUE),
                     sd_chg   = sd(pct_chg, na.rm = TRUE),
                     n_sp     = dplyr::n(), .groups = "drop")

  # 95% CI across species, matching blue_carbon_full.R
  sp_grad_summ <- sp_grad_summ |>
    dplyr::mutate(ci_lo = mean_chg - 1.96 * sd_chg / sqrt(pmax(n_sp, 1)),
                  ci_hi = mean_chg + 1.96 * sd_chg / sqrt(pmax(n_sp, 1)))

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
}

message("\n========== ALL FIGURES COMPLETE ==========")
message("Output: ", normalizePath(OUT_DIR))
