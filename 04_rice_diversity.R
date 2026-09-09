#!/usr/bin/env Rscript

# 04_rice_diversity.R — Rice Diversity Panel 1 (RiceCGM/RDP1) two-environment mixFunMap
# analysis (reproduction script for the mixFunMap manuscript).
#
# Two parts:
#   1. Analysis (mode = "full"): reruns the three scans from the prepared
#      inputs in inputs/rice_diversity/ (349 accessions, 21 days, 33,697 markers):
#      single-environment logistic mixFunMap scans for Control and Low
#      water (plus ordinary-FunMap and GMMAT-minP benchmarks) and the
#      joint two-environment 3-df parameter G x E scan, following
#      ril/scripts/run_lowwater_logistic_three_methods.R and
#      ril/scripts/run_ricecgm_joint_gxe_pilot.R (full mode).
#   2. Figures (default): redraws every manuscript panel (a-k) from the
#      frozen outputs in results/rice_diversity/ without refitting any model.
#
# Usage:
#   Rscript 04_rice_diversity.R            # figures only (fast; default)
#   Rscript 04_rice_diversity.R full [cores] [output_dir]
#
# Model specification (identical to the manuscript):
#   * Phenotypes: per-environment sample means, analysed at scale x1e-5.
#   * Single environment: mean curve = three-parameter logistic;
#     Q = PC1-PC5 + VanRaden K; mixFunMap 3-df P3D Wald scan (Henderson
#     REML). Benchmarks: ordinary FunMap (additive dosages, SAD(1)) and
#     per-day GMMAT-minP.
#   * Joint G x E: fit_mixfunmap_joint with cov_model = "correlated",
#     kinship_cov = "diagonal", phi fixed to 1; 3-df Wald test of
#     [-I3, I3] beta = 0 on endpoint-scaled (0/1) dosages.
#   * Difference-phenotype (Low water - Control) GMMAT-minP benchmark.
#
# Thresholds: M = 33,697 markers; LD-pruned M_eff = 6,904 (r2 = 0.5);
# Bonferroni 0.05/M_eff and suggestive 1/M lines.
#
# Requirements:
#   figures : R >= 4.4 with ggplot2 and ggrepel
#   full    : additionally the mixFunMap package and GMMAT

options(stringsAsFactors = FALSE)

.args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(.args) >= 1L) tolower(.args[[1L]]) else "figures"

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_argument)) {
  dirname(normalizePath(sub("^--file=", "", file_argument[[1L]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

inputs_dir <- file.path(script_dir, "inputs", "rice_diversity")
benchmark_dir <- file.path(script_dir, "results", "rice_diversity")
figure_dir <- file.path(script_dir, "figures")

control_dir <- file.path(benchmark_dir, "control")
lowwater_dir <- file.path(benchmark_dir, "lowwater")
gxe_dir <- file.path(benchmark_dir, "gxe")
minp_control_dir <- file.path(benchmark_dir, "control_minp")
fm_control_dir <- file.path(benchmark_dir, "control_fm")
minp_lowwater_dir <- file.path(benchmark_dir, "lowwater_minp")
fm_lowwater_dir <- file.path(benchmark_dir, "lowwater_fm")

pheno_control_path <- file.path(benchmark_dir, "phenotype_Control_mean_matrix.tsv.gz")
pheno_lowwater_path <- file.path(benchmark_dir, "phenotype_LowWater_mean_matrix.tsv.gz")
annotation_path <- file.path(benchmark_dir, "three_scans_Meff_locus_annotation_msu7.csv")

# Frozen scan tables are bundled gzip-compressed; fresh rerun tables are
# plain TSV. Resolve whichever exists.
scan_file <- function(dir, name) {
  path <- file.path(dir, name)
  if (file.exists(path)) path else paste0(path, ".gz")
}

# ---------------------------------------------------------------------------
# Part 1: analysis (only runs in mode = "full")
# ---------------------------------------------------------------------------

run_analysis <- function(output_dir = file.path(benchmark_dir, "rerun"),
                         cores = max(1L, parallel::detectCores(logical = FALSE)),
                         n_pcs = 5L, maf_min = 0.05, max_missing = 0.50) {
  if (!requireNamespace("mixFunMap", quietly = TRUE)) {
    stop("The mixFunMap package is required for the analysis mode.", call. = FALSE)
  }
  if (!requireNamespace("GMMAT", quietly = TRUE)) {
    stop("GMMAT is required for the minP scans.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv <- function(x, path) {
    utils::write.table(x, path, sep = "\t", quote = FALSE,
                       row.names = FALSE, na = "NA")
  }
  lambda_of <- function(p, df) {
    p <- p[is.finite(p) & p > 0 & p <= 1]
    if (!length(p)) return(NA_real_)
    stats::median(stats::qchisq(1 - p, df = df)) / stats::qchisq(0.5, df = df)
  }

  ## Prepared inputs (349 samples x 21 days; 33,697 markers) --------------
  x <- readRDS(file.path(inputs_dir, "analysis_inputs_RiceCGM_RDP1_joint.rds"))
  pheno_control <- as.matrix(x$pheno_control) / 1e5     # analysis scale x1e-5
  pheno_low_water <- as.matrix(x$pheno_low_water) / 1e5
  index <- as.numeric(x$index)
  geno <- as.matrix(x$geno)                             # 0/1/2 dosages
  Q <- as.matrix(x$Q)[, seq_len(n_pcs), drop = FALSE]
  K <- as.matrix(x$K)
  marker_map <- data.frame(
    marker = as.character(x$marker_map$marker_id),
    chr = as.character(x$marker_map$chromosome),
    pos = as.numeric(x$marker_map$position_bp),
    maf = as.numeric(x$marker_map$maf),
    stringsAsFactors = FALSE
  )
  marker_map <- marker_map[match(rownames(geno), marker_map$marker), ]
  map <- marker_map[, c("marker", "chr", "pos")]

  ## Single-environment three-method scans ---------------------------------
  scan_environment <- function(pheno, tag) {
    fit <- mixFunMap::fit_mixfunmap(
      pheno = pheno, index = index, Q = Q, K = K, geno = geno,
      mean = mixFunMap::mean_logistic(), engine = "henderson",
      max_outer = 60L, max_reml_optim = 150L, tol_theta = 0.01,
      verbose = FALSE
    )
    if (!isTRUE(fit$usable)) stop(tag, ": null model is not usable.",
                                  call. = FALSE)
    scan_mix <- mixFunMap::scan_mixfunmap(
      fit, geno = geno, test = "wald", cores = cores,
      max_outer_snp = 6L, tol_theta = 1e-1
    )
    beta <- as.data.frame(scan_mix$beta, check.names = FALSE)
    tab_mix <- data.frame(
      map, stat = beta[map$marker, "wald3df"], df = 3L,
      pval = beta[map$marker, "pval"],
      converged = beta[map$marker, "converged"] > 0,
      stringsAsFactors = FALSE
    )
    write_tsv(tab_mix, file.path(output_dir, tag, "mixed_logistic_scan.tsv"))

    null_ofm <- mixFunMap::fit_ordinary_funmap(pheno, index, maxit = 300L,
                                               reltol = 1e-8)
    tab_ofm <- mixFunMap::scan_ordinary_funmap(
      pheno = pheno, geno = geno, index = index, marker_map = map,
      null_fit = null_ofm, genetic_model = "additive",
      min_group_n = 5L, maxit = 300L, reltol = 1e-8
    )
    write_tsv(as.data.frame(tab_ofm),
              file.path(output_dir, paste0(tag, "_fm"),
                        "ordinary_logistic_scan.tsv"))

    tab_minp <- mixFunMap::scan_minp_gmmat(
      pheno = pheno, geno = geno, index = index, Q = Q, K = K,
      marker_map = map, correction = "bonferroni", maxiter = 200L,
      tol = 1e-5, primary_optimizer = "AI", fallback_optimizer = "Brent",
      verbose = FALSE, ncores = min(cores, length(index))
    )
    write_tsv(as.data.frame(tab_minp, check.names = FALSE),
              file.path(output_dir, paste0(tag, "_minp"),
                        "minp_gmmat_scan.tsv"))

    write_tsv(
      data.frame(
        method = "mixed_logistic_qk", markers = nrow(tab_mix),
        valid = sum(is.finite(tab_mix$pval)),
        lambda_gc = lambda_of(tab_mix$pval, 3L)
      ),
      file.path(output_dir, tag, "method_summary.tsv")
    )
    invisible(tab_mix)
  }
  scan_environment(pheno_control, "control")
  scan_environment(pheno_low_water, "lowwater")

  ## Joint two-environment 3-df parameter G x E scan -----------------------
  geno_01 <- geno / 2          # endpoint-scaled dosages in [0, 1]
  fit_joint <- mixFunMap::fit_mixfunmap_joint(
    pheno_env1 = pheno_control, pheno_env2 = pheno_low_water,
    index = index, Q = Q, K = K, geno = NULL,
    mean = mixFunMap::mean_logistic(),
    cov_model = "correlated", kinship_cov = "diagonal",
    phi_mode = "fixed_one", engine = "henderson",
    include_intercept = TRUE, max_outer = 100L, max_reml_optim = 200L,
    tol_theta = 0.01, damp_alpha = 0.2, verbose = FALSE
  )
  scan_gxe <- mixFunMap::scan_mixfunmap(
    fit_joint, geno = geno_01, test = "parameter_gxe",
    genetic_model = "additive", genotype_levels = c(0, 1),
    gxe_statistics = "wald", primary_gxe = "wald", cores = cores,
    max_outer_snp = 60L, tol_parameter = 1e-4, tol_ml = 1e-7,
    tol_curve = 1e-4, retry_failed_alt = TRUE, return_curves = FALSE,
    verbose = FALSE
  )
  tab_gxe <- as.data.frame(scan_gxe$beta, check.names = FALSE)
  tab_gxe$marker <- rownames(scan_gxe$beta)
  tab_gxe <- merge(map, tab_gxe, by = "marker", sort = FALSE)
  write_tsv(tab_gxe, file.path(output_dir, "gxe", "parameter_gxe_scan.tsv"))
  write_tsv(
    data.frame(
      component = "parameter_gxe", markers = nrow(tab_gxe),
      valid = sum(is.finite(tab_gxe$wald_pval)),
      lambda_calibration = lambda_of(tab_gxe$wald_pval, 3L)
    ),
    file.path(output_dir, "gxe", "analysis_gate_summary.tsv")
  )

  ## Difference-phenotype (Low water - Control) GMMAT-minP benchmark -------
  tab_diff <- mixFunMap::scan_minp_gmmat(
    pheno = pheno_low_water - pheno_control, geno = geno, index = index,
    Q = Q, K = K, marker_map = map, correction = "bonferroni",
    maxiter = 200L, tol = 1e-5, primary_optimizer = "AI",
    fallback_optimizer = "Brent", verbose = FALSE,
    ncores = min(cores, length(index))
  )
  write_tsv(as.data.frame(tab_diff, check.names = FALSE),
            file.path(output_dir, "gxe", "difference_minp_scan.tsv"))
  message("Analysis tables written to: ", output_dir)
  invisible(output_dir)
}

if (identical(mode, "full")) {
  cores <- if (length(.args) >= 2L) as.integer(.args[[2L]]) else
    max(1L, parallel::detectCores(logical = FALSE))
  out <- if (length(.args) >= 3L) .args[[3L]] else
    file.path(benchmark_dir, "rerun")
  for (sub in c("control", "lowwater", "gxe", "control_minp", "control_fm",
                "lowwater_minp", "lowwater_fm")) {
    dir.create(file.path(out, sub), recursive = TRUE, showWarnings = FALSE)
  }
  run_analysis(output_dir = out, cores = cores)
  # Redraw the panels below from the fresh run instead of the frozen one.
  control_dir <- file.path(out, "control")
  lowwater_dir <- file.path(out, "lowwater")
  gxe_dir <- file.path(out, "gxe")
  minp_control_dir <- file.path(out, "control_minp")
  fm_control_dir <- file.path(out, "control_fm")
  minp_lowwater_dir <- file.path(out, "lowwater_minp")
  fm_lowwater_dir <- file.path(out, "lowwater_fm")
}

for (path in c(pheno_control_path, pheno_lowwater_path,
               scan_file(control_dir, "mixed_logistic_scan.tsv"),
               file.path(control_dir, "method_summary.tsv"),
               scan_file(lowwater_dir, "mixed_logistic_scan.tsv"),
               file.path(lowwater_dir, "method_summary.tsv"),
               scan_file(gxe_dir, "parameter_gxe_scan.tsv"),
               file.path(gxe_dir, "analysis_gate_summary.tsv"),
               scan_file(minp_control_dir, "minp_gmmat_scan.tsv"),
               scan_file(fm_control_dir, "ordinary_logistic_scan.tsv"),
               scan_file(minp_lowwater_dir, "minp_gmmat_scan.tsv"),
               scan_file(fm_lowwater_dir, "ordinary_logistic_scan.tsv"),
               scan_file(gxe_dir, "difference_minp_scan.tsv"),
               annotation_path)) {
  if (!file.exists(path)) stop("Required input not found: ", path, call. = FALSE)
}

# ---------------------------------------------------------------------------
# Shared style (consistent with the other manuscript figures)
# ---------------------------------------------------------------------------

panel_theme <- function() {
  ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, color = "#E5E7EB"),
      strip.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 11, colour = "black"),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank()
    )
}

save_panel <- function(plot, stem, width, height) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
                  width = width, height = height, units = "in",
                  device = grDevices::cairo_pdf, bg = "white")
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot,
                  width = width, height = height, units = "in", dpi = 600,
                  compression = "lzw", bg = "white")
}

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                    quote = "", comment.char = "")
}

# Task palette: same blue / vermillion / bluish-green family as the other
# manuscript figures.
task_colors <- c(control = "#0072B2", lowwater = "#D55E00", gxe = "#009E73")
task_labels <- c(control = "Control", lowwater = "Low water", gxe = "G \u00d7 E")

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

scan_control <- read_tsv(scan_file(control_dir, "mixed_logistic_scan.tsv"))
scan_lowwater <- read_tsv(scan_file(lowwater_dir, "mixed_logistic_scan.tsv"))
scan_gxe <- read_tsv(scan_file(gxe_dir, "parameter_gxe_scan.tsv"))

summary_control <- read_tsv(file.path(control_dir, "method_summary.tsv"))
summary_lowwater <- read_tsv(file.path(lowwater_dir, "method_summary.tsv"))
gate_gxe <- read_tsv(file.path(gxe_dir, "analysis_gate_summary.tsv"))

lambdas <- c(
  control = summary_control$lambda_gc[summary_control$method == "mixed_logistic_qk"][1],
  lowwater = summary_lowwater$lambda_gc[summary_lowwater$method == "mixed_logistic_qk"][1],
  gxe = gate_gxe$lambda_calibration[gate_gxe$component == "parameter_gxe"][1]
)

# Phenotype matrices: rows = days, columns = samples; the analyses used
# values scaled by 1e-5, and the same scale is used for plotting.
pheno_scale <- 1e-5
read_pheno <- function(path) {
  tab <- read_tsv(path)
  days <- as.numeric(tab$day)
  mat <- as.matrix(tab[, -(1:2), drop = FALSE])
  storage.mode(mat) <- "double"
  list(days = days, mat = mat * pheno_scale)
}
pheno_control <- read_pheno(pheno_control_path)
pheno_lowwater <- read_pheno(pheno_lowwater_path)

# Shared y-axis scale for panels a and b (panel a spans both
# environments; panel b shows their difference on the same scale, with
# the y tick numbers suppressed).
yr <- range(c(pheno_control$mat, pheno_lowwater$mat), na.rm = TRUE)
yr_pad <- 0.05 * diff(yr)
ylim_ab <- c(yr[1] - yr_pad, yr[2] + yr_pad)

# Significance thresholds: Bonferroni on the LD-pruned effective marker
# number (ril/prepared/meff_ldprune_r05.npy, r2 = 0.5) and the suggestive
# 1/M line on the full marker set.
M <- nrow(scan_control)
M_eff <- 6904L
threshold_lines <- data.frame(
  name = factor(
    c(sprintf("Bonferroni (0.05/M_eff = 0.05/%d)", M_eff),
      sprintf("Suggestive (1/M = 1/%d)", M)),
    levels = c(sprintf("Bonferroni (0.05/M_eff = 0.05/%d)", M_eff),
               sprintf("Suggestive (1/M = 1/%d)", M))
  ),
  y = -log10(c(0.05 / M_eff, 1 / M))
)
threshold_colors <- stats::setNames(c("#D62728", "#7F7F7F"),
                                    levels(threshold_lines$name))

# Locus annotation from the benchmark CSV; all G x E loci are labelled on
# the linear Manhattan (panel e): significant (Meff) loci with short gene
# names in black, suggestive (1/M) loci with gene IDs in grey.
annotation <- utils::read.csv(annotation_path, check.names = FALSE,
                              stringsAsFactors = FALSE)
gxe_annotation <- annotation[annotation$scan == "GxE", , drop = FALSE]
gxe_annotation$gene_id <- sub("\\(.*$", "", gxe_annotation$closest_gene)
gxe_annotation$significant <- gxe_annotation$tier == "significant(Meff)"

# Short functional labels; uninformative descriptions (retrotransposon /
# expressed protein) fall back to the gene ID.
gene_short_names <- c(
  LOC_Os01g08700 = "GIGANTEA",
  LOC_Os01g50700 = "Dehydrin",
  LOC_Os05g40810 = "BRCA1 C-term",
  LOC_Os05g48700 = "GA2ox"
)
gxe_annotation$label <- ifelse(gxe_annotation$gene_id %in% names(gene_short_names),
                               unname(gene_short_names[gxe_annotation$gene_id]),
                               gxe_annotation$gene_id)

# ---------------------------------------------------------------------------
# Panel a: per-sample growth dynamics, colored by environment
# ---------------------------------------------------------------------------

make_panel_a <- function() {
  build_curves <- function(pheno, env_key) {
    mat <- pheno$mat
    data.frame(
      day = rep(pheno$days, ncol(mat)),
      value = as.vector(mat),
      sample = rep(colnames(mat), each = nrow(mat)),
      env = env_key
    )
  }
  # Interleave the two environments' sample curves in random order (fixed
  # seed) so neither color clouds over the other; each cloud's density
  # shows its own color where it concentrates.
  curve_df <- rbind(build_curves(pheno_control, "control"),
                    build_curves(pheno_lowwater, "lowwater"))
  set.seed(42)
  curve_df <- curve_df[sample(nrow(curve_df)), , drop = FALSE]
  curve_df$env <- factor(curve_df$env, levels = c("control", "lowwater"),
                         labels = task_labels[c("control", "lowwater")])
  mean_df <- do.call(rbind, lapply(c(control = "control", lowwater = "lowwater"),
    function(env_key) {
      pheno <- if (identical(env_key, "control")) pheno_control else pheno_lowwater
      data.frame(day = pheno$days, value = rowMeans(pheno$mat, na.rm = TRUE),
                 env = task_labels[[env_key]])
    }))
  mean_df$env <- factor(mean_df$env, levels = task_labels[c("control", "lowwater")])

  env_colors <- stats::setNames(unname(task_colors[c("control", "lowwater")]),
                                task_labels[c("control", "lowwater")])
  # Sample curves use medium-strength tints of the environment colors:
  # light enough that the bold mean curves stand out, saturated enough
  # that the two environment clouds stay distinguishable where they
  # overlap. Mean curves are plain bold lines without points.
  env_tints <- stats::setNames(c("#6BAED6", "#F08A3C"),
                               task_labels[c("control", "lowwater")])

  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = curve_df,
      ggplot2::aes(x = .data$day, y = .data$value,
                   group = paste(.data$sample, .data$env),
                   color = .data$env),
      linewidth = 0.2, alpha = 0.35
    ) +
    ggplot2::geom_line(
      data = mean_df[mean_df$env == task_labels[["control"]], , drop = FALSE],
      ggplot2::aes(x = .data$day, y = .data$value),
      color = unname(task_colors[["control"]]), linewidth = 0.9
    ) +
    ggplot2::geom_line(
      data = mean_df[mean_df$env == task_labels[["lowwater"]], , drop = FALSE],
      ggplot2::aes(x = .data$day, y = .data$value),
      color = unname(task_colors[["lowwater"]]), linewidth = 0.9
    ) +
    # Sample curves take the medium tints; the legend keys are overridden
    # to the saturated environment colors of the mean curves.
    ggplot2::scale_color_manual(values = env_tints) +
    ggplot2::coord_cartesian(ylim = ylim_ab) +
    ggplot2::labs(x = NULL, y = NULL) +
    panel_theme() +
    ggplot2::theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(color = unname(env_colors), linewidth = 1.2,
                          alpha = 1)))
}

# ---------------------------------------------------------------------------
# Panel b: difference phenotype (Control - Low water) dynamics
# ---------------------------------------------------------------------------

make_panel_b <- function() {
  common <- intersect(colnames(pheno_control$mat), colnames(pheno_lowwater$mat))
  stopifnot(length(common) > 0, identical(pheno_control$days, pheno_lowwater$days))
  diff_mat <- pheno_control$mat[, common, drop = FALSE] -
    pheno_lowwater$mat[, common, drop = FALSE]
  days <- pheno_control$days

  curve_df <- data.frame(
    day = rep(days, ncol(diff_mat)),
    value = as.vector(diff_mat),
    sample = rep(seq_len(ncol(diff_mat)), each = nrow(diff_mat))
  )
  mean_df <- data.frame(day = days, value = rowMeans(diff_mat, na.rm = TRUE))

  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "#9CA3AF", linewidth = 0.4) +
    ggplot2::geom_line(
      data = curve_df,
      ggplot2::aes(x = .data$day, y = .data$value, group = .data$sample),
      linewidth = 0.25, alpha = 0.25, color = "#9CA3AF"
    ) +
    ggplot2::geom_line(
      data = mean_df, ggplot2::aes(x = .data$day, y = .data$value),
      linewidth = 0.8, color = "#0072B2"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_cartesian(ylim = ylim_ab) +
    panel_theme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
}

# ---------------------------------------------------------------------------
# Panel c: QQ plot of the three tasks (Control and Low water first,
# G x E on top; lambda annotated per task)
# ---------------------------------------------------------------------------

make_panel_c <- function() {
  task_levels <- c("control", "lowwater", "gxe")
  qq_labels <- c(
    control = sprintf("Control (\u03bb = %.2f)", lambdas[["control"]]),
    lowwater = sprintf("Low water (\u03bb = %.2f)", lambdas[["lowwater"]]),
    gxe = sprintf("G \u00d7 E (Plasticity, \u03bb = %.2f)", lambdas[["gxe"]])
  )
  scans <- list(control = scan_control, lowwater = scan_lowwater, gxe = scan_gxe)
  qq_df <- do.call(rbind, lapply(task_levels, function(task) {
    p <- sort(scans[[task]]$pval[is.finite(scans[[task]]$pval) & scans[[task]]$pval > 0])
    n <- length(p)
    data.frame(task = task,
               expected = -log10((seq_len(n) - 0.5) / n),
               observed = -log10(pmax(p, 1e-300)))
  }))
  qq_df$task <- factor(qq_df$task, levels = task_levels, labels = qq_labels[task_levels])
  # Draw order: Control (1), Low water (2), G x E (3) last on top.
  qq_df <- qq_df[order(match(as.integer(qq_df$task), c(1L, 2L, 3L))), ]
  qq_colors <- stats::setNames(unname(task_colors[task_levels]), qq_labels[task_levels])

  ggplot2::ggplot(qq_df, ggplot2::aes(.data$expected, .data$observed,
                                      color = .data$task)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "#333333", linewidth = 0.4) +
    ggplot2::geom_point(size = 0.6, alpha = 0.5) +
    ggplot2::scale_color_manual(values = qq_colors) +
    ggplot2::labs(x = NULL, y = NULL) +
    panel_theme() +
    ggplot2::theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 1.6, alpha = 1)))
}

# ---------------------------------------------------------------------------
# Panel d: circular Manhattan, three rings
#   outer ring = Control, middle ring = Low water, inner ring = G x E
# ---------------------------------------------------------------------------

# Shared chromosome geometry for the Manhattan panels (same marker set in
# all three scans).
chroms <- as.character(1:12)

manhattan_layout <- function(gap_frac = 0.012) {
  chr_max <- tapply(scan_control$pos[scan_control$chr %in% chroms],
                    as.character(scan_control$chr[scan_control$chr %in% chroms]),
                    max)[chroms]
  gap <- gap_frac * sum(as.numeric(chr_max))
  offsets <- c(0, cumsum(as.numeric(chr_max) + gap)[-length(chr_max)])
  names(offsets) <- chroms
  list(chr_max = chr_max, gap = gap, offsets = offsets,
       ticks = offsets + as.numeric(chr_max) / 2)
}

chr_colors <- stats::setNames(rep(c("#1F78B4", "#FF7F00"), length.out = 12),
                              chroms)

# ---------------------------------------------------------------------------
# Panel d: two-ring inward-facing circular Manhattan
#   outer ring = Control, inner ring = Low water; peaks point toward the
#   center so the two environments mirror each other around the ring gap.
# ---------------------------------------------------------------------------

make_panel_d <- function() {
  layout <- manhattan_layout()

  # Distinct color per chromosome (publication-style multi-hue palette).
  chr_colors12 <- stats::setNames(
    c("#4E9BD0", "#F2C12E", "#5D51A5", "#6B7F39", "#D9488B", "#E2670E",
      "#B47BB3", "#66C2A5", "#9E2B25", "#3D8B87", "#E7298A", "#8C6D4F"),
    chroms
  )

  # Widen the gap at the top junction (x = 0) so each ring can carry a
  # small radial -log10(p) axis there, as in the reference style.
  total_len <- sum(as.numeric(layout$chr_max)) + 12 * layout$gap
  start_gap <- 0.030 * total_len
  x_total <- total_len + start_gap
  x_shift <- stats::setNames(layout$offsets + start_gap, chroms)
  chr_ticks <- x_shift + as.numeric(layout$chr_max) / 2

  prep <- function(scan, task) {
    keep <- scan$chr %in% chroms & is.finite(scan$pval) & scan$pval > 0
    df <- scan[keep, c("chr", "pos", "pval")]
    df$chr <- as.character(df$chr)
    df$task <- task
    df$neglog10p <- -log10(pmax(df$pval, 1e-300))
    df$x <- x_shift[df$chr] + df$pos
    df$pt_color <- unname(chr_colors12[df$chr])
    df
  }
  dat <- rbind(prep(scan_control, "control"), prep(scan_lowwater, "lowwater"))

  # Ring geometry: each ring has a baseline circle; points extend INWARD
  # from the baseline. A center hole keeps the inner ring readable.
  hole <- 7
  ring_gap <- 2.2
  ring_max <- c(control = max(dat$neglog10p[dat$task == "control"]),
                lowwater = max(dat$neglog10p[dat$task == "lowwater"]))
  ring_base <- c(
    lowwater = hole + ring_max[["lowwater"]],
    control = hole + ring_max[["lowwater"]] + ring_gap + ring_max[["control"]]
  )
  dat$y <- ring_base[dat$task] - dat$neglog10p

  # Highlight the annotated loci (significant and suggestive alike) of both
  # single-environment scans in red, slightly larger than the cloud points.
  hl <- annotation[annotation$scan %in% c("Control", "LowWater"), , drop = FALSE]
  hl$task <- ifelse(hl$scan == "Control", "control", "lowwater")
  hl$x <- NA_real_
  hl$neglog10p <- NA_real_
  for (i in seq_len(nrow(hl))) {
    sc <- if (identical(hl$task[[i]], "control")) scan_control else scan_lowwater
    hit <- match(hl$lead[[i]], sc$marker)
    pos <- if (is.na(hit)) (hl$start[[i]] + hl$end[[i]]) / 2 else sc$pos[[hit]]
    pv <- if (is.na(hit)) hl$p[[i]] else sc$pval[[hit]]
    hl$x[[i]] <- x_shift[[as.character(hl$chr[[i]])]] + pos
    hl$neglog10p[[i]] <- -log10(max(pv, 1e-300))
  }
  hl$y <- ring_base[hl$task] - hl$neglog10p

  # Threshold circles per ring (Bonferroni red dashed, suggestive grey).
  thr_df <- do.call(rbind, lapply(c("control", "lowwater"), function(task) {
    data.frame(task = task,
               name = levels(threshold_lines$name),
               y = ring_base[[task]] - threshold_lines$y)
  }))
  thr_df$linetype <- ifelse(startsWith(thr_df$name, "Bonferroni"),
                            "dashed", "dotted")

  # Outermost black chromosome ideogram band with white gaps.
  ideo_df <- data.frame(
    xmin = x_shift, xmax = x_shift + as.numeric(layout$chr_max),
    ymin = ring_base[["control"]] + 1.5, ymax = ring_base[["control"]] + 2.9
  )

  # Chromosome labels outside the ideogram band, rotated tangentially.
  chr_label_df <- data.frame(
    x = unname(chr_ticks),
    y = ring_base[["control"]] + 4.6,
    label = paste0("Chr", chroms),
    angle = -360 * (unname(chr_ticks) / x_total)
  )

  # Small radial axis (0, 2, 4, ...) per ring inside the widened top gap.
  axis_df <- do.call(rbind, lapply(c("control", "lowwater"), function(task) {
    v <- seq(0, floor(ring_max[[task]] / 2) * 2, by = 2)
    data.frame(x = start_gap * 0.62, y = ring_base[[task]] - v, label = v)
  }))
  axis_x <- start_gap * 0.62
  axis_line_df <- data.frame(
    task = c("control", "lowwater"),
    y0 = ring_base[c("control", "lowwater")],
    y1 = ring_base[c("control", "lowwater")] - ring_max[c("control", "lowwater")]
  )
  tick_len <- 0.0018 * x_total

  ggplot2::ggplot(dat, ggplot2::aes(.data$x, .data$y)) +
    # ring baseline circles
    ggplot2::geom_hline(
      data = data.frame(y = unname(ring_base)),
      ggplot2::aes(yintercept = .data$y),
      color = "#9CA3AF", linewidth = 0.3
    ) +
    ggplot2::geom_point(color = dat$pt_color, size = 0.4, na.rm = TRUE) +
    ggplot2::geom_point(
      data = hl, ggplot2::aes(.data$x, .data$y),
      color = "#D62728", size = 0.9
    ) +
    ggplot2::geom_hline(
      data = thr_df,
      ggplot2::aes(yintercept = .data$y, linetype = .data$linetype),
      color = ifelse(startsWith(thr_df$name, "Bonferroni"), "#D62728", "#7F7F7F"),
      linewidth = 0.35, show.legend = FALSE
    ) +
    ggplot2::scale_linetype_identity() +
    # radial axis spokes, tick marks and tick numbers in the top gap
    ggplot2::geom_segment(
      data = axis_line_df,
      ggplot2::aes(x = axis_x, xend = axis_x, y = .data$y0, yend = .data$y1),
      color = "#374151", linewidth = 0.55, inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      data = axis_df,
      ggplot2::aes(x = axis_x - tick_len, xend = axis_x,
                   y = .data$y, yend = .data$y),
      color = "#374151", linewidth = 0.45, inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = axis_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 2.6, color = "#111827", hjust = -0.4
    ) +
    # ideogram band and tangential chromosome labels
    ggplot2::geom_rect(
      data = ideo_df,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$ymin, ymax = .data$ymax),
      fill = "black", inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = chr_label_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
                   angle = .data$angle),
      size = 3.4, color = "#111827"
    ) +
    ggplot2::coord_polar(theta = "x") +
    ggplot2::scale_y_continuous(limits = c(0, ring_base[["control"]] + 6.2),
                                expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}

# ---------------------------------------------------------------------------
# Panel e: linear G x E Manhattan with Bonferroni (0.05/M_eff) and
# suggestive (1/M) lines; significant (Meff) loci annotated per the CSV.
# ---------------------------------------------------------------------------

make_panel_e <- function() {
  layout <- manhattan_layout(gap_frac = 0)  # no gaps between chromosomes
  keep <- scan_gxe$chr %in% chroms & is.finite(scan_gxe$pval) & scan_gxe$pval > 0
  mapped <- scan_gxe[keep, c("marker", "chr", "pos", "pval")]
  mapped$chr <- as.character(mapped$chr)
  mapped$neglog10p <- -log10(pmax(mapped$pval, 1e-300))
  mapped$x <- layout$offsets[mapped$chr] + mapped$pos
  mapped$pt_color <- unname(chr_colors[mapped$chr])

  p <- ggplot2::ggplot(mapped, ggplot2::aes(.data$x, .data$neglog10p)) +
    ggplot2::geom_point(color = mapped$pt_color, size = 0.5, na.rm = TRUE) +
    ggplot2::geom_hline(
      data = threshold_lines,
      ggplot2::aes(yintercept = .data$y, color = .data$name),
      linetype = "dashed", linewidth = 0.5
    ) +
    ggplot2::scale_color_manual(values = threshold_colors) +
    ggplot2::scale_x_continuous(breaks = layout$ticks,
                                labels = paste0("Chr", chroms),
                                expand = c(0.005, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    panel_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 7, colour = "black"),
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 8)
    )

  ann <- gxe_annotation
  hit <- match(ann$lead, mapped$marker)
  ann$x <- layout$offsets[as.character(ann$chr)] +
    ifelse(is.na(hit), (ann$start + ann$end) / 2, mapped$pos[hit])
  ann$y <- ifelse(is.na(hit), -log10(ann$p), mapped$neglog10p[hit])
  # All annotated loci (significant and suggestive alike) share the FIP1
  # Manhattan style: red highlight point + black label of the same size.
  p +
    ggplot2::geom_point(
      data = ann, ggplot2::aes(.data$x, .data$y),
      color = "#D62728", size = 1.2
    ) +
    ggrepel::geom_text_repel(
      data = ann, ggplot2::aes(.data$x, .data$y, label = .data$label),
      size = 2.8, color = "#111111", max.overlaps = Inf,
      min.segment.length = 0, segment.size = 0.25, segment.color = "#9CA3AF",
      box.padding = 0.35, point.padding = 0.2, seed = 42
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0, max(mapped$neglog10p, ann$y, na.rm = TRUE) * 1.35),
      clip = "off"
    )
}

# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------

stem <- "Figure4_rice_ril"
panel_a <- make_panel_a()
panel_b <- make_panel_b()
panel_c <- make_panel_c()
panel_d <- make_panel_d()
panel_e <- make_panel_e()

save_panel(panel_a, paste0(stem, "_a_phenotype_dynamics"), 3.6, 2.4)
save_panel(panel_b, paste0(stem, "_b_difference_dynamics"), 3.6, 2.4)
save_panel(panel_c, paste0(stem, "_c_qq"), 3.6, 2.4)
save_panel(panel_d, paste0(stem, "_d_circular_manhattan"), 7.2, 6.4)
save_panel(panel_e, paste0(stem, "_e_manhattan_gxe"), 7.2, 3.4)

message("Rice RIL figures written to: ", figure_dir)

# ---------------------------------------------------------------------------
# Panel f: Venn diagrams of locus-level overlap between the three scans
# ---------------------------------------------------------------------------
# Loci are the LD-clumped intervals listed in
# three_scans_Meff_locus_annotation_msu7.csv. A locus pair from two scans is
# declared the same locus when their intervals overlap on the same
# chromosome; clusters (connected components) of overlapping loci are then
# counted once in the Venn region matching the scans they appear in. This is
# done separately for the Bonferroni (0.05/M_eff, "significant(Meff)") and
# suggestive (1/M, "suggestive(1/M)") tiers.

locus_overlap_clusters <- function(sub) {
  n <- nrow(sub)
  parent <- seq_len(n)
  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  if (n > 1L) {
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        if (sub$chr[i] == sub$chr[j] &&
            sub$start[i] <= sub$end[j] && sub$start[j] <= sub$end[i]) {
          parent[find_root(i)] <- find_root(j)
        }
      }
    }
  }
  roots <- vapply(seq_len(n), find_root, integer(1L))
  lapply(split(seq_len(n), roots), function(idx) unique(sub$scan[idx]))
}

venn_region_counts <- function(tier_name) {
  sub <- annotation[annotation$tier == tier_name, , drop = FALSE]
  patterns <- vapply(
    locus_overlap_clusters(sub),
    function(s) paste(sort(s), collapse = "&"),
    character(1L)
  )
  region_levels <- c(
    "Control", "LowWater", "GxE",
    "Control&LowWater", "Control&GxE", "LowWater&GxE",
    "Control&LowWater&GxE"
  )
  as.integer(table(factor(patterns, levels = region_levels)))
}

venn_scan_totals <- function(tier_name) {
  sub <- annotation[annotation$tier == tier_name, , drop = FALSE]
  stats::setNames(
    as.integer(table(factor(sub$scan, levels = c("Control", "LowWater", "GxE")))),
    c("Control", "LowWater", "GxE")
  )
}

make_panel_f <- function() {
  tiers <- c("significant(Meff)", "suggestive(1/M)")

  venn_sets <- c("Control", "LowWater", "GxE")
  venn_centers <- data.frame(
    set = venn_sets,
    cx = c(0.00, 0.90, 0.45),
    cy = c(0.45, 0.45, -0.33)
  )
  venn_radius <- 0.75
  venn_fill <- stats::setNames(
    unname(task_colors[c("control", "lowwater", "gxe")]), venn_sets
  )
  # Region label anchors inside the three-circle layout, in the order of the
  # region levels used by venn_region_counts(). Every region shows
  # "Bonferroni | suggestive" locus counts at its own anchor.
  region_levels <- c(
    "Control", "LowWater", "GxE",
    "Control&LowWater", "Control&GxE", "LowWater&GxE",
    "Control&LowWater&GxE"
  )
  region_anchor <- data.frame(
    region = region_levels,
    x = c(-0.45, 1.35, 0.45, 0.45, -0.05, 0.95, 0.45),
    y = c(0.58, 0.58, -0.85, 0.66, -0.15, -0.15, 0.20)
  )
  set_anchor <- data.frame(
    set = venn_sets,
    x = c(-0.95, 1.85, 0.45),
    y = c(1.32, 1.32, -1.52)
  )

  bonf_counts <- venn_region_counts(tiers[1L])
  sugg_counts <- venn_region_counts(tiers[2L])
  bonf_totals <- venn_scan_totals(tiers[1L])
  sugg_totals <- venn_scan_totals(tiers[2L])

  circles <- do.call(
    rbind,
    lapply(venn_sets, function(s) {
      theta <- seq(0, 2 * pi, length.out = 240)
      ctr <- venn_centers[venn_centers$set == s, ]
      data.frame(
        x = ctr$cx + venn_radius * cos(theta),
        y = ctr$cy + venn_radius * sin(theta),
        set = s
      )
    })
  )
  labels <- region_anchor
  labels$count <- paste0(bonf_counts, " | ", sugg_counts)
  labels$any_nonzero <- (bonf_counts + sugg_counts) > 0
  names_df <- set_anchor
  names_df$total <- paste0(bonf_totals[names_df$set], " | ",
                           sugg_totals[names_df$set])
  names_df$set_display <- ifelse(
    names_df$set == "GxE", "G\u00d7E (Plasticity)", names_df$set
  )
  title_df <- data.frame(
    x = 0.45, y = 1.95,
    label = 'Bonferroni~(0.05/M[eff])~~"|"~~Suggestive~(1/M)'
  )

  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = circles,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$set,
                   fill = .data$set, color = .data$set),
      alpha = 0.35, linewidth = 0.7
    ) +
    ggplot2::geom_text(
      data = labels[labels$any_nonzero, ],
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$count),
      size = 4.0, fontface = "bold", color = "grey15"
    ) +
    ggplot2::geom_text(
      data = labels[!labels$any_nonzero, ],
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$count),
      size = 3.6, color = "grey55"
    ) +
    ggplot2::geom_text(
      data = names_df,
      ggplot2::aes(x = .data$x, y = .data$y,
                   label = paste0(.data$set_display, " (", .data$total, ")"),
                   color = .data$set),
      size = 4.0, fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = title_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 4.4, parse = TRUE
    ) +
    ggplot2::scale_fill_manual(values = venn_fill) +
    ggplot2::scale_color_manual(values = venn_fill) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(8, 12, 8, 12)
    )
}

panel_f <- make_panel_f()
save_panel(panel_f, paste0(stem, "_f_venn_locus_overlap"), 4.3, 3.4)

message("Rice RIL Venn panel written to: ", figure_dir)

# ---------------------------------------------------------------------------
# Benchmark scan inputs for panels g-j (GMMAT-minP and ordinary Functional
# Mapping; plot-only, frozen outputs)
# ---------------------------------------------------------------------------

scan_minp_control <- read_tsv(scan_file(minp_control_dir, "minp_gmmat_scan.tsv"))
scan_fm_control <- read_tsv(scan_file(fm_control_dir, "ordinary_logistic_scan.tsv"))
scan_minp_lowwater <- read_tsv(scan_file(minp_lowwater_dir, "minp_gmmat_scan.tsv"))
scan_fm_lowwater <- read_tsv(scan_file(fm_lowwater_dir, "ordinary_logistic_scan.tsv"))
scan_minp_gxe <- read_tsv(scan_file(gxe_dir, "difference_minp_scan.tsv"))

# ---------------------------------------------------------------------------
# Panels g-i: inward-facing circular Manhattans of the benchmark scans.
#   scan_list names give the rings, outer ring first. Every ring has its own
#   radial scale (the Functional Mapping scans reach -log10(p) ~ 320 while
#   the minP scans stay below ~6), so each ring carries its own radial axis
#   in the widened top gap; Bonferroni (red dashed) and suggestive (grey
#   dotted) threshold circles are drawn per ring when they fall inside the
#   ring's data range.
# ---------------------------------------------------------------------------

make_circular_benchmark <- function(scan_list) {
  tasks <- names(scan_list)  # outer ring first
  layout <- manhattan_layout()
  chr_max <- layout$chr_max
  total_len <- sum(as.numeric(chr_max)) + length(chroms) * layout$gap
  start_gap <- 0.030 * total_len
  x_total <- total_len + start_gap
  x_shift <- stats::setNames(layout$offsets + start_gap, chroms)
  chr_ticks <- x_shift + as.numeric(chr_max) / 2

  chr_colors12 <- stats::setNames(
    c("#4E9BD0", "#F2C12E", "#5D51A5", "#6B7F39", "#D9488B", "#E2670E",
      "#B47BB3", "#66C2A5", "#9E2B25", "#3D8B87", "#E7298A", "#8C6D4F"),
    chroms
  )

  prep <- function(scan, task) {
    keep <- scan$chr %in% chroms & is.finite(scan$pval) & scan$pval > 0
    df <- scan[keep, c("chr", "pos", "pval")]
    df$chr <- as.character(df$chr)
    df$task <- task
    df$neglog10p <- -log10(pmax(df$pval, 1e-300))
    df$x <- x_shift[df$chr] + df$pos
    df$pt_color <- unname(chr_colors12[df$chr])
    df
  }
  dat <- do.call(rbind, Map(prep, scan_list, tasks))

  # Ring geometry: fixed physical height per ring, points extend INWARD from
  # the baseline circle; radial position is scaled by each ring's own
  # maximum -log10(p).
  hole <- 7
  ring_gap <- 2.2
  ring_h <- 8
  ring_max <- vapply(tasks, function(t) max(dat$neglog10p[dat$task == t]),
                     numeric(1L))
  ring_base <- numeric(length(tasks))
  names(ring_base) <- tasks
  base <- hole
  for (t in rev(tasks)) {
    ring_base[[t]] <- base + ring_h
    base <- base + ring_h + ring_gap
  }
  dat$y <- ring_base[dat$task] -
    dat$neglog10p / ring_max[dat$task] * ring_h
  outer <- tasks[[1L]]

  # Per-ring threshold circles (kept only when inside the ring's range).
  thr_df <- do.call(rbind, lapply(tasks, function(task) {
    keep <- threshold_lines$y <= ring_max[[task]]
    if (!any(keep)) return(NULL)
    data.frame(task = task,
               name = as.character(threshold_lines$name[keep]),
               y = ring_base[[task]] -
                 threshold_lines$y[keep] / ring_max[[task]] * ring_h,
               stringsAsFactors = FALSE)
  }))
  if (!is.null(thr_df)) {
    thr_df$linetype <- ifelse(startsWith(thr_df$name, "Bonferroni"),
                              "dashed", "dotted")
  }

  # Outermost black chromosome ideogram band and tangential labels.
  ideo_df <- data.frame(
    xmin = x_shift, xmax = x_shift + as.numeric(chr_max),
    ymin = ring_base[[outer]] + 1.5, ymax = ring_base[[outer]] + 2.9
  )
  chr_label_df <- data.frame(
    x = unname(chr_ticks),
    y = ring_base[[outer]] + 4.6,
    label = paste0("Chr", chroms),
    angle = -360 * (unname(chr_ticks) / x_total)
  )

  # Per-ring radial axes in the top gap, labelled with actual -log10(p).
  axis_df <- do.call(rbind, lapply(tasks, function(task) {
    v <- pretty(c(0, ring_max[[task]]), n = 4)
    v <- v[v <= ring_max[[task]]]
    data.frame(task = task, x = start_gap * 0.62,
               y = ring_base[[task]] - v / ring_max[[task]] * ring_h,
               label = v)
  }))
  axis_x <- start_gap * 0.62
  axis_line_df <- data.frame(
    task = tasks,
    y0 = ring_base[tasks],
    y1 = ring_base[tasks] - ring_h
  )
  tick_len <- 0.0018 * x_total

  p <- ggplot2::ggplot(dat, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_hline(
      data = data.frame(y = unname(ring_base)),
      ggplot2::aes(yintercept = .data$y),
      color = "#9CA3AF", linewidth = 0.3
    ) +
    ggplot2::geom_point(color = dat$pt_color, size = 0.4, na.rm = TRUE) +
    ggplot2::geom_segment(
      data = axis_line_df,
      ggplot2::aes(x = axis_x, xend = axis_x, y = .data$y0, yend = .data$y1),
      color = "#374151", linewidth = 0.55, inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      data = axis_df,
      ggplot2::aes(x = axis_x - tick_len, xend = axis_x,
                   y = .data$y, yend = .data$y),
      color = "#374151", linewidth = 0.45, inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = axis_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 2.6, color = "#111827", hjust = -0.4
    ) +
    ggplot2::geom_rect(
      data = ideo_df,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$ymin, ymax = .data$ymax),
      fill = "black", inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = chr_label_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
                   angle = .data$angle),
      size = 3.4, color = "#111827"
    ) +
    ggplot2::coord_polar(theta = "x") +
    ggplot2::scale_y_continuous(limits = c(0, ring_base[[outer]] + 6.2),
                                expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(
      legend.position = "none",
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )

  if (!is.null(thr_df) && nrow(thr_df)) {
    p <- p +
      ggplot2::geom_hline(
        data = thr_df,
        ggplot2::aes(yintercept = .data$y, linetype = .data$linetype),
        color = ifelse(startsWith(thr_df$name, "Bonferroni"),
                       "#D62728", "#7F7F7F"),
        linewidth = 0.35, show.legend = FALSE
      ) +
      ggplot2::scale_linetype_identity()
  }
  p
}

# ---------------------------------------------------------------------------
# Panels j-k: QQ plots of the benchmark scans, one per method; color = task
# (Control / Low water / G x E). Panel j = GMMAT-minP (three tasks),
# panel k = Functional Mapping (Control and Low water).
# ---------------------------------------------------------------------------

make_panel_qq_benchmark <- function(method_key) {
  series <- switch(method_key,
    minp = list(
      list(task = "control",  scan = scan_minp_control),
      list(task = "lowwater", scan = scan_minp_lowwater),
      list(task = "gxe",      scan = scan_minp_gxe)
    ),
    fm = list(
      list(task = "control",  scan = scan_fm_control),
      list(task = "lowwater", scan = scan_fm_lowwater)
    ),
    stop("Unknown method key: ", method_key, call. = FALSE)
  )
  qq_df <- do.call(rbind, lapply(series, function(s) {
    p <- sort(s$scan$pval[is.finite(s$scan$pval) & s$scan$pval > 0])
    n <- length(p)
    data.frame(task = s$task,
               expected = -log10((seq_len(n) - 0.5) / n),
               observed = -log10(pmax(p, 1e-300)))
  }))
  qq_labels_j <- c(control = "Control", lowwater = "Low water",
                   gxe = "G \u00d7 E (Plasticity)")
  keep_labels <- qq_labels_j[unique(
    sapply(series, function(s) s$task))]
  qq_df$task <- factor(qq_df$task, levels = names(keep_labels),
                       labels = keep_labels)
  qq_colors_j <- stats::setNames(unname(task_colors[names(keep_labels)]),
                                 keep_labels)
  # minP curves hug the diagonal: legend top-left. The FM curves shoot up
  # the left side and plateau: legend mid-right on a white background.
  legend_at <- if (identical(method_key, "minp")) {
    list(position = c(0.02, 0.98), justification = c(0, 1),
         background = ggplot2::element_blank())
  } else {
    list(position = c(0.97, 0.45), justification = c(1, 0.5),
         background = ggplot2::element_rect(
           fill = scales::alpha("white", 0.9), color = NA))
  }

  ggplot2::ggplot(qq_df, ggplot2::aes(.data$expected, .data$observed,
                                      color = .data$task)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "#333333", linewidth = 0.4) +
    ggplot2::geom_point(size = 0.6, alpha = 0.5) +
    ggplot2::scale_color_manual(values = qq_colors_j) +
    ggplot2::labs(x = NULL, y = NULL) +
    panel_theme() +
    ggplot2::theme(
      legend.position = legend_at$position,
      legend.justification = legend_at$justification,
      legend.background = legend_at$background,
      legend.text = ggplot2::element_text(size = 8),
      legend.spacing.y = ggplot2::unit(0.02, "cm")
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(size = 1.6, alpha = 1)))
}

panel_g <- make_circular_benchmark(list(fm = scan_fm_control,
                                        minp = scan_minp_control))
panel_h <- make_circular_benchmark(list(fm = scan_fm_lowwater,
                                        minp = scan_minp_lowwater))
panel_i <- make_circular_benchmark(list(minp = scan_minp_gxe))
panel_j <- make_panel_qq_benchmark("minp")
panel_k <- make_panel_qq_benchmark("fm")

save_panel(panel_g, paste0(stem, "_g_circular_manhattan_control_fm_minp"),
           7.2, 6.4)
save_panel(panel_h, paste0(stem, "_h_circular_manhattan_lowwater_fm_minp"),
           7.2, 6.4)
save_panel(panel_i, paste0(stem, "_i_circular_manhattan_gxe_minp"), 6.0, 5.4)
save_panel(panel_j, paste0(stem, "_j_qq_minp"), 3.6, 2.4)
save_panel(panel_k, paste0(stem, "_k_qq_fm"), 3.6, 2.4)

# The combined five-series QQ panel is superseded by the per-method panels.
for (ext in c(".pdf", ".tiff")) {
  stale <- file.path(figure_dir, paste0(stem, "_j_qq_fm_minp", ext))
  if (file.exists(stale)) file.remove(stale)
}

message("Rice RIL benchmark FM/minP panels written to: ", figure_dir)
