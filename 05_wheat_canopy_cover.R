#!/usr/bin/env Rscript

# 05_wheat_canopy_cover.R — FIP1 wheat FPWW012 canopy-cover LOP analysis
# (reproduction script for the mixFunMap manuscript).
#
# Two parts:
#   1. Analysis (mode = "full"): reruns the LOP model selection and the
#      three-method genome-wide scan from the analysis-ready inputs in
#      inputs/wheat_canopy_cover/ (plus the 0/1/2 hard-call matrix in
#      inputs/wheat_height/ for the Functional Mapping comparator),
#      following real_data_analysis/run_fip1_cc_lop_gwas.R
#      (criterion = BIC) and run_fip1_cc_bic4_minp_funmap.R.
#   2. Figures (default): redraws every manuscript panel from the frozen
#      outputs in results/wheat_cc/ without refitting any model.
#
# Usage:
#   Rscript 05_wheat_canopy_cover.R            # figures only (default)
#   Rscript 05_wheat_canopy_cover.R full [cores] [output_dir]
#
# Model specification (identical to the manuscript):
#   * Mean curve: Legendre orthogonal polynomial (LOP); degree selected
#     from {2..6} by BIC on the population mean curve -> degree 4.
#   * mixFunMap: Q = PC1-PC5 + VanRaden K, fix_tau3 = FALSE, 5-df P3D Wald
#     joint test on the five LOP coefficients.
#   * Functional Mapping: degree-4 LOP + SAD(1), genotypic 0/1/2 groups,
#     fixed-covariance LRT (no Q/K).
#   * GMMAT-minP: 39 per-date Q + K score tests, min(1, 39 * min_t p_t).
#
# Thresholds: same FIP1 LD-aware thresholds as the height analysis
# (Bonferroni 0.05/M_eff with M_eff = 3,270; suggestive 1/M with
# M = 18,583).
#
# Requirements:
#   figures : R >= 4.4 with ggplot2 and ggrepel
#   full    : additionally the mixFunMap package and GMMAT

options(stringsAsFactors = FALSE)

.args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(.args) >= 1L) tolower(.args[[1L]]) else "figures"

script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
                            grep("^--file=", commandArgs(FALSE),
                                 value = TRUE)[1L]))),
  error = function(...) getwd()
)
inputs_dir <- file.path(script_dir, "inputs", "wheat_canopy_cover")
results_dir <- file.path(script_dir, "results", "wheat_cc")
comparators_dir <- file.path(results_dir, "comparators")
figure_dir <- file.path(script_dir, "figures")
pheno_path <- file.path(
  results_dir, "phenotype_canopy_cover_FPWW012_matrix.tsv.gz"
)

# ---------------------------------------------------------------------------
# Part 1: analysis (only runs in mode = "full")
# ---------------------------------------------------------------------------

run_analysis <- function(output_dir = file.path(results_dir, "rerun"),
                         cores = max(1L, parallel::detectCores(logical = FALSE))) {
  if (!requireNamespace("mixFunMap", quietly = TRUE)) {
    stop("The mixFunMap package is required for the analysis mode.", call. = FALSE)
  }
  if (!requireNamespace("GMMAT", quietly = TRUE)) {
    stop("GMMAT is required for the minP scan.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "comparators"), showWarnings = FALSE)
  write_tsv <- function(x, path) utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )

  inputs <- readRDS(file.path(
    inputs_dir, "analysis_inputs_FPWW012_canopy_cover_lop.rds"
  ))
  pheno <- as.matrix(inputs$pheno)
  index <- as.numeric(inputs$index)
  dates <- inputs$dates
  geno <- as.matrix(inputs$geno)
  Q <- as.matrix(inputs$Q)
  K <- as.matrix(inputs$K)
  marker_map <- as.data.frame(inputs$marker_map, stringsAsFactors = FALSE)
  map <- marker_map[, c("marker", "chr", "pos")]
  geno_hardcall <- as.matrix(readRDS(file.path(
    script_dir, "inputs", "wheat_height",
    "analysis_inputs_FPWW012_height_ordinary3class_v2.rds"
  ))$geno_hardcall)
  genomewide_alpha <- 0.05 / nrow(geno)

  ## 1. mixFunMap LOP: BIC degree selection on the mean curve, then the
  ##    Q + K null fit and the 5-df P3D Wald scan --------------------------
  fit <- mixFunMap::fit_mixfunmap(
    pheno = pheno, index = index, Q = Q, K = K,
    mean = mixFunMap::mean_lop(
      degree = "auto", criterion = "BIC",
      degree_candidates = 2:6, scale_index = TRUE
    ),
    engine = "auto", fix_tau3 = FALSE,
    max_outer = 80L, max_reml_optim = 200L, tol_theta = 0.01,
    reml_every = 1L, verbose = FALSE
  )
  if (!isTRUE(fit$usable)) stop("LOP null model is not usable.", call. = FALSE)
  degree_table <- as.data.frame(fit$mean$model_selection)
  write_tsv(degree_table, file.path(output_dir, "lop_degree_selection.tsv"))
  selected_degree <- as.integer(fit$mean$degree)
  selection_coef <- as.numeric(fit$mean$init(pheno, index))
  population_mean <- rowMeans(pheno)
  selection_fit <- as.numeric(
    mixFunMap::lop_basis(index, selected_degree, scale_index = TRUE) %*%
      selection_coef
  )
  write_tsv(
    data.frame(date = as.character(dates), day = index,
               observed_mean = population_mean, fitted_lop = selection_fit,
               residual = population_mean - selection_fit),
    file.path(output_dir, "lop_population_mean_fit.tsv")
  )

  scan1 <- mixFunMap::scan_mixfunmap(
    fit, geno = geno, test = "wald", cores = cores,
    max_outer_snp = 15L, tol_theta = 0.02
  )
  beta <- as.matrix(scan1$beta)
  retry <- !is.finite(beta[, "pval"]) | beta[, "pval"] <= 0 |
    beta[, "pval"] > 1 | beta[, "converged"] <= 0 | beta[, "outer_converged"] <= 0
  if (any(retry)) {
    scan2 <- mixFunMap::scan_mixfunmap(
      fit, geno = geno[retry, , drop = FALSE], test = "wald", cores = cores,
      max_outer_snp = 45L, tol_theta = 0.01
    )
    beta2 <- as.matrix(scan2$beta)
    beta[match(rownames(beta2), rownames(beta)), colnames(beta2)] <- beta2
  }
  stat_col <- paste0("wald", fit$mean$npar, "df")
  table_mixfunmap <- data.frame(
    marker_map[, c("marker", "chr", "pos", "mapped")],
    stat = as.numeric(beta[marker_map$marker, stat_col]),
    df = fit$mean$npar,
    pval = as.numeric(beta[marker_map$marker, "pval"]),
    converged = as.numeric(beta[marker_map$marker, "converged"]) > 0,
    outer_converged = as.numeric(beta[marker_map$marker, "outer_converged"]) > 0,
    stringsAsFactors = FALSE
  )
  con <- gzfile(file.path(output_dir, "mixfunmap_lop_results.tsv.gz"), "wt")
  utils::write.table(table_mixfunmap, con, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "NA")
  close(con)

  ## 2. Functional Mapping comparator: degree-4 LOP + SAD(1), genotype-
  ##    group fixed-covariance LRT (P3D; = Wald for this linear mean) ------
  B <- mixFunMap::lop_basis(index, 4L, scale_index = TRUE)
  sad1_base <- function(phi, n_time) {
    iv <- numeric(n_time); iv[1L] <- 1
    for (i in 2:n_time) iv[i] <- phi^2 * iv[i - 1L] + 1
    out <- matrix(0, n_time, n_time)
    for (i in seq_len(n_time)) for (j in i:n_time) {
      out[i, j] <- out[j, i] <- phi^(j - i) * iv[i]
    }
    out
  }
  profile_null <- function(phi, retain = FALSE) {
    R <- sad1_base(phi, nrow(pheno))
    chol_R <- tryCatch(chol(R), error = function(e) NULL)
    if (is.null(chol_R)) return(if (retain) NULL else .Machine$double.xmax / 100)
    Bw <- forwardsolve(t(chol_R), B)
    Yw <- forwardsolve(t(chol_R), pheno)
    info_inv <- chol2inv(chol(crossprod(Bw)))
    beta_hat <- info_inv %*% crossprod(Bw, rowMeans(Yw))
    sse <- sum((Yw - as.numeric(Bw %*% beta_hat))^2)
    gamma2 <- sse / (nrow(pheno) * ncol(pheno))
    if (!is.finite(gamma2) || gamma2 <= 0) {
      return(if (retain) NULL else .Machine$double.xmax / 100)
    }
    nll <- 0.5 * (nrow(pheno) * ncol(pheno) * (log(2 * pi) + log(gamma2) + 1) +
                    ncol(pheno) * 2 * sum(log(diag(chol_R))))
    if (!retain) return(nll)
    list(beta = drop(beta_hat), phi = phi, gamma = sqrt(gamma2),
         Bw = Bw, Yw = Yw, info_inv = info_inv, pheno = pheno)
  }
  null_ofm <- profile_null(
    stats::optimize(profile_null, interval = c(1e-4, 0.98),
                    tol = .Machine$double.eps^0.35)$minimum,
    retain = TRUE
  )
  table_funmap <- do.call(rbind, lapply(seq_len(nrow(geno_hardcall)),
    function(j) {
      g <- geno_hardcall[j, ]
      finite <- is.finite(g)
      observed <- (0:2)[vapply(0:2, function(l) any(g[finite] == l), logical(1))]
      groups <- observed[vapply(observed, function(l) sum(g[finite] == l),
                                integer(1)) >= 5L]
      keep <- finite & g %in% groups
      df <- 5L * max(1L, length(groups) - 1L)
      row <- data.frame(map[j, ], mapped = marker_map$mapped[j],
                        stat = NA_real_, df = df, pval = NA_real_,
                        converged = FALSE, stringsAsFactors = FALSE)
      if (length(groups) < 2L) return(row)
      Yw <- null_ofm$Yw[, keep, drop = FALSE]
      g_used <- g[keep]
      beta0 <- null_ofm$info_inv %*% crossprod(null_ofm$Bw, rowMeans(Yw))
      sse0 <- sum((Yw - as.numeric(null_ofm$Bw %*% beta0))^2)
      sse1 <- 0
      for (grp in groups) {
        sel <- g_used == grp
        beta_g <- null_ofm$info_inv %*%
          crossprod(null_ofm$Bw, rowMeans(Yw[, sel, drop = FALSE]))
        sse1 <- sse1 +
          sum((Yw[, sel, drop = FALSE] - as.numeric(null_ofm$Bw %*% beta_g))^2)
      }
      stat <- max(0, sse0 - sse1)
      row$stat <- stat
      row$pval <- stats::pchisq(stat, df = df, lower.tail = FALSE)
      row$converged <- is.finite(row$pval)
      row
    }
  ))
  con <- gzfile(file.path(output_dir, "comparators",
                          "ordinary_funmap_lop4_results.tsv.gz"), "wt")
  utils::write.table(table_funmap, con, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "NA")
  close(con)

  ## 3. GMMAT-minP comparator: 39 per-date Q + K score tests ---------------
  table_minp <- mixFunMap::scan_minp_gmmat(
    pheno = pheno, geno = geno, index = index, Q = Q, K = K,
    marker_map = map, correction = "bonferroni", maxiter = 200L, tol = 1e-5,
    primary_optimizer = "AI", fallback_optimizer = "Brent",
    verbose = FALSE, ncores = cores
  )
  table_minp <- as.data.frame(table_minp, stringsAsFactors = FALSE)
  table_minp$mapped <- marker_map$mapped
  con <- gzfile(file.path(output_dir, "comparators",
                          "minp_gmmat_results.tsv.gz"), "wt")
  utils::write.table(table_minp, con, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "NA")
  close(con)

  ## Shared summary ---------------------------------------------------------
  lambda_of <- function(p, df) {
    ok <- is.finite(p) & p > 0 & p <= 1
    stats::median(stats::qchisq(1 - p[ok], df = df)) /
      stats::qchisq(0.5, df = df)
  }
  write_tsv(
    data.frame(
      dataset = "FIP1_FPWW012_canopy_cover_lop_v1", criterion = "BIC",
      selected_degree = selected_degree, test_df = fit$mean$npar,
      n_samples = ncol(pheno), n_time_points = nrow(pheno),
      n_markers_scanned = nrow(table_mixfunmap),
      lambda_gc = lambda_of(table_mixfunmap$pval, fit$mean$npar),
      top_marker = table_mixfunmap$marker[which.min(table_mixfunmap$pval)],
      top_p = min(table_mixfunmap$pval, na.rm = TRUE),
      genomewide_alpha = genomewide_alpha
    ),
    file.path(output_dir, "analysis_summary.tsv")
  )
  write_tsv(
    data.frame(
      method = c("mixfunmap_lop4_BIC", "ordinary_funmap_lop4",
                 "minp_gmmat"),
      df = c(fit$mean$npar, NA, NA),
      n_valid = c(sum(is.finite(table_mixfunmap$pval)),
                  sum(is.finite(table_funmap$pval)),
                  sum(is.finite(table_minp$pval))),
      lambda_gc = c(lambda_of(table_mixfunmap$pval, fit$mean$npar),
                    lambda_of(table_funmap$pval, table_funmap$df),
                    NA)
    ),
    file.path(output_dir, "comparators", "lambda_gc_audit.tsv")
  )
  message("Analysis tables written to: ", output_dir)
  invisible(output_dir)
}

if (identical(mode, "full")) {
  cores <- if (length(.args) >= 2L) as.integer(.args[[2L]]) else
    max(1L, parallel::detectCores(logical = FALSE))
  out <- if (length(.args) >= 3L) .args[[3L]] else
    file.path(results_dir, "rerun")
  run_analysis(output_dir = out, cores = cores)
  # Redraw the figures below from the fresh run instead of the frozen one.
  results_dir <- out
  comparators_dir <- file.path(out, "comparators")
}

for (path in c(results_dir, pheno_path,
               file.path(results_dir, "analysis_summary.tsv"),
               file.path(comparators_dir, "ordinary_funmap_lop4_results.tsv.gz"),
               file.path(comparators_dir, "minp_gmmat_results.tsv.gz"),
               file.path(comparators_dir, "lambda_gc_audit.tsv"))) {
  if (!file.exists(path)) stop("Required input not found: ", path, call. = FALSE)
}

# ---------------------------------------------------------------------------
# Shared helpers and theme (same style as the other manuscript figure scripts)
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

# Method palette: identical to the other manuscript figure scripts.
method_colors <- c(mixfunmap = "#0072B2", ordinary_funmap = "#D55E00",
                   minp = "#009E73")

save_panel <- function(plot, stem, width, height) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
                  width = width, height = height, units = "in",
                  device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot,
                  width = width, height = height, units = "in", dpi = 600,
                  compression = "lzw")
}

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                    quote = "", comment.char = "")
}

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

pheno_tab <- read_tsv(pheno_path)
mean_fit <- read_tsv(file.path(results_dir, "lop_population_mean_fit.tsv"))
degree_table <- read_tsv(file.path(results_dir, "lop_degree_selection.tsv"))

stopifnot(nrow(degree_table) > 0L, sum(degree_table$selected) == 1L)

# ---------------------------------------------------------------------------
# Panel a: canopy-cover trajectories + fourth-order LOP population-mean fit
# ---------------------------------------------------------------------------

make_panel_a <- function() {
  days <- as.numeric(pheno_tab$day)
  cc <- as.matrix(pheno_tab[, -1, drop = FALSE])
  storage.mode(cc) <- "double"

  curve_df <- data.frame(
    day = rep(days, ncol(cc)),
    cc = as.vector(cc),
    sample = rep(seq_len(ncol(cc)), each = nrow(cc))
  )

  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, color = "#9CA3AF",
                        linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_line(
      data = curve_df,
      ggplot2::aes(x = .data$day, y = .data$cc, group = .data$sample),
      linewidth = 0.25, alpha = 0.25, color = "#9CA3AF"
    ) +
    ggplot2::geom_line(
      data = mean_fit,
      ggplot2::aes(x = .data$day, y = .data$fitted_lop),
      linewidth = 0.9, color = "#0072B2"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    panel_theme()
}

# ---------------------------------------------------------------------------
# BIC inset: LOP degree-selection mini plot embedded in panel a
# (no y-axis numbers, only a "BIC" axis title; optimum at degree 4 in red)
# ---------------------------------------------------------------------------

make_bic_inset <- function() {
  df <- degree_table[order(degree_table$degree), , drop = FALSE]
  best <- df$degree[df$selected][1L]

  ggplot2::ggplot(df, ggplot2::aes(x = .data$degree, y = .data$BIC)) +
    ggplot2::geom_vline(xintercept = best, color = "#D62728",
                        linetype = "dashed", linewidth = 0.35) +
    ggplot2::geom_line(color = "#374151", linewidth = 0.5) +
    ggplot2::geom_point(color = "#374151", size = 1.2) +
    ggplot2::geom_point(
      data = df[df$selected, , drop = FALSE],
      color = "#D62728", size = 1.8
    ) +
    ggplot2::scale_x_continuous(breaks = df$degree) +
    ggplot2::labs(x = NULL, y = "BIC") +
    ggplot2::theme_bw(base_size = 7) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 7, colour = "black"),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 8, colour = "black"),
      plot.background = ggplot2::element_rect(fill = "white", colour = "#9CA3AF",
                                              linewidth = 0.3),
      panel.background = ggplot2::element_rect(fill = "white"),
      plot.margin = ggplot2::margin(2, 4, 2, 2)
    )
}

# ---------------------------------------------------------------------------
# Panel b: three-method QQ plot (Functional Mapping and GMMAT-minP first,
# mixFunMap on top; lambda annotated for Functional Mapping and mixFunMap,
# same convention as Figure 2b)
# ---------------------------------------------------------------------------

scan_mixfunmap <- read_tsv(file.path(results_dir, "mixfunmap_lop_results.tsv.gz"))
scan_fm <- read_tsv(file.path(comparators_dir, "ordinary_funmap_lop4_results.tsv.gz"))
scan_minp <- read_tsv(file.path(comparators_dir, "minp_gmmat_results.tsv.gz"))

# lambda_GC: mixFunMap from the primary run's analysis summary; Functional
# Mapping from the comparator audit (5-df family, the dominant test).
summary_cc <- read_tsv(file.path(results_dir, "analysis_summary.tsv"))
lambda_audit <- read_tsv(file.path(comparators_dir, "lambda_gc_audit.tsv"))
lambdas <- c(
  mixfunmap = summary_cc$lambda_gc[1L],
  ordinary_funmap = lambda_audit$lambda_gc[
    lambda_audit$method == "ordinary_funmap_lop4" & lambda_audit$df == 5][1L]
)

make_panel_b <- function() {
  method_levels <- c("mixfunmap", "ordinary_funmap", "minp")
  method_labels <- c(
    mixfunmap = sprintf("mixFunMap (\u03bb = %.3f)", lambdas[["mixfunmap"]]),
    ordinary_funmap = sprintf("Functional Mapping (\u03bb = %.2f)",
                              lambdas[["ordinary_funmap"]]),
    minp = "GMMAT-minP (index-adjusted p)"
  )
  scans <- list(mixfunmap = scan_mixfunmap, ordinary_funmap = scan_fm,
                minp = scan_minp)
  qq_df <- do.call(rbind, lapply(method_levels, function(meth) {
    p <- sort(scans[[meth]]$pval[is.finite(scans[[meth]]$pval) &
                                 scans[[meth]]$pval > 0])
    n <- length(p)
    data.frame(method = meth,
               expected = -log10((seq_len(n) - 0.5) / n),
               observed = -log10(pmax(p, 1e-300)))
  }))
  qq_df$method <- factor(qq_df$method, levels = method_levels,
                         labels = method_labels[method_levels])
  # Draw order: Functional Mapping (2), GMMAT-minP (3), mixFunMap (1) last.
  qq_df <- qq_df[order(match(as.integer(qq_df$method), c(2L, 3L, 1L))), ]
  qq_colors <- stats::setNames(unname(method_colors[method_levels]),
                               method_labels[method_levels])

  ggplot2::ggplot(qq_df, ggplot2::aes(.data$expected, .data$observed,
                                      color = .data$method)) +
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
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(size = 1.6, alpha = 1)))
}

# ---------------------------------------------------------------------------
# Panel c: mixFunMap Manhattan with Bonferroni (0.05/M_eff, M_eff = 3270
# from LD pruning at r2 = 0.2) and suggestive (1/M, M = 18583) threshold
# lines; the annotated regions of region_annotation.csv are highlighted in
# red and labelled (same style as Figure 2c).
# ---------------------------------------------------------------------------

M <- nrow(scan_mixfunmap)
M_eff <- 3270L
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

chroms21 <- paste0("chr", rep(1:7, each = 3), rep(c("A", "B", "D"), 7))

cc_annotation <- utils::read.csv(
  file.path(results_dir, "region_annotation.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)

# Short functional labels; uninformative descriptions ("uncharacterized",
# "protein enabled", ...) fall back to the gene ID.
gene_short_names_cc <- c(
  TraesCS5B02G430100 = "AO2",
  TraesCS2B02G414300 = "FRO2",
  TraesCS6B02G052900 = "SKIP4",
  TraesCS3A02G453200 = "ATG18a",
  TraesCS1B02G457600 = "SLP2",
  TraesCS1A02G311900 = "bZIP23",
  TraesCS1B02G129200 = "AIG1",
  TraesCS1B02G307500 = "SNX2B",
  TraesCS2A02G098100 = "CaM1",
  TraesCS5B02G257300 = "WRKY2",
  TraesCS2D02G355200 = "DAO5"
)
cc_annotation$label <- ifelse(cc_annotation$gene %in% names(gene_short_names_cc),
                              unname(gene_short_names_cc[cc_annotation$gene]),
                              cc_annotation$gene)

make_manhattan_cc <- function(scan, annotation_df = NULL) {
  mapped <- scan[scan$mapped & scan$chr %in% chroms21, , drop = FALSE]
  mapped$neglog10p <- -log10(pmax(mapped$pval, 1e-300))
  chr_max <- tapply(mapped$pos, mapped$chr, max)[chroms21]
  offsets <- c(0, cumsum(as.numeric(chr_max))[-length(chr_max)])
  names(offsets) <- chroms21
  mapped$x <- offsets[mapped$chr] + mapped$pos
  # Color by subgenome (A / B / D): adjacent chromosomes always differ.
  subgenome_colors <- c(A = "#1F78B4", B = "#FF7F00", D = "#6A3D9A")
  mapped$pt_color <- unname(subgenome_colors[sub("^chr\\d+", "", mapped$chr)])
  ticks <- offsets + as.numeric(chr_max) / 2

  p <- ggplot2::ggplot(mapped, ggplot2::aes(.data$x, .data$neglog10p)) +
    ggplot2::geom_point(color = mapped$pt_color, size = 0.5, na.rm = TRUE) +
    ggplot2::geom_hline(
      data = threshold_lines,
      ggplot2::aes(yintercept = .data$y, color = .data$name),
      linetype = "dashed", linewidth = 0.5
    ) +
    ggplot2::scale_color_manual(values = threshold_colors) +
    ggplot2::scale_x_continuous(breaks = ticks,
                                labels = sub("chr", "", chroms21),
                                expand = c(0.005, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    panel_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 7, colour = "black"),
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 8)
    )

  if (!is.null(annotation_df) && nrow(annotation_df)) {
    ann <- annotation_df
    hit <- match(ann$top_marker, mapped$marker)
    ann$x <- offsets[ann$chr] + ann$top_pos
    ann$y <- ifelse(is.na(hit), -log10(pmax(ann$top_p, 1e-300)),
                    mapped$neglog10p[hit])
    p <- p +
      ggplot2::geom_point(
        data = ann, ggplot2::aes(.data$x, .data$y),
        color = "#D62728", size = 1.2
      ) +
      ggrepel::geom_text_repel(
        data = ann, ggplot2::aes(.data$x, .data$y, label = .data$label),
        size = 2.8, color = "#111111", max.overlaps = Inf,
        min.segment.length = 0, segment.size = 0.25, segment.color = "#6B7280",
        box.padding = 0.35, point.padding = 0.2, seed = 42
      ) +
      ggplot2::coord_cartesian(
        ylim = c(0, max(mapped$neglog10p, ann$y, na.rm = TRUE) * 1.35),
        clip = "off"
      )
  }
  p
}

# ---------------------------------------------------------------------------
# Build and save
# ---------------------------------------------------------------------------

stem <- "Figure5_wheat_cc"

panel_a <- make_panel_a() +
  patchwork::inset_element(make_bic_inset(),
                           left = 0.02, bottom = 0.52,
                           right = 0.40, top = 0.98,
                           align_to = "panel", on_top = TRUE)
panel_b <- make_panel_b()
panel_c <- make_manhattan_cc(scan_mixfunmap, cc_annotation)
panel_d <- make_manhattan_cc(scan_minp)
panel_e <- make_manhattan_cc(scan_fm)

save_panel(panel_a, paste0(stem, "_a_cc_trajectories_lop4_mean"), 3.6, 2.4)
save_panel(panel_b, paste0(stem, "_b_qq"), 3.6, 2.4)
save_panel(panel_c, paste0(stem, "_c_manhattan_mixfunmap"), 7.2, 3.4)
save_panel(panel_d, paste0(stem, "_d_manhattan_minp"), 7.2, 3.4)
save_panel(panel_e, paste0(stem, "_e_manhattan_functional_mapping"), 7.2, 3.4)

# The standalone BIC panel is superseded by the inset; remove stale outputs.
for (ext in c(".pdf", ".tiff")) {
  stale <- file.path(figure_dir, paste0(stem, "_b_lop_degree_bic", ext))
  if (file.exists(stale)) file.remove(stale)
}

cat("Wheat canopy-cover figures written to: ",
    normalizePath(figure_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")

# ---------------------------------------------------------------------------
# Panel fg: genotype-mean curves (top) and conditional genetic SD (bottom)
# for the four manuscript-chosen mixFunMap loci (AO2 5B, FRO2 2B, AIG1 1B,
# WRKY2 5B); same stacked layout as Figure 2fg. Plot-only: reads the frozen
# estimate_genotype_curves() output.
# ---------------------------------------------------------------------------

top4_markers_cc <- c(
  "RFL_Contig2772_1693",      # AO2, chr5B 605.5 Mb; strongest p, headline
  "BS00064448_51",            # FRO2, chr2B 592.0 Mb; two homozygous groups
  "Tdurum_contig42092_348",   # AIG1, chr1B 159.2 Mb
  "wsnp_Ex_c6548_11355524"    # WRKY2, chr5B 439.7 Mb
)
top4_marker_labels_cc <- c(
  RFL_Contig2772_1693 = "AO2 (5B)",
  BS00064448_51 = "FRO2 (2B)",
  Tdurum_contig42092_348 = "AIG1 (1B)",
  wsnp_Ex_c6548_11355524 = "WRKY2 (5B)"
)

# Frozen output of mixFunMap::estimate_genotype_curves(..., background =
# "refit") + estimate_genetic_sd(method = "conditional") for the four loci
# (produced by the original pipeline; refitting takes several minutes).
top4_frozen_cc <- readRDS(file.path(
  script_dir, "results", "wheat_cc", "fip1_cc_top4_interpretability.rds"
))
top4_curves_cc <- top4_frozen_cc$genotype_curves
top4_gsd_cc <- top4_frozen_cc$genetic_sd
stopifnot(all(top4_markers_cc %in% top4_curves_cc$curves$marker))

curve_plot_data <- top4_curves_cc$curves
curve_plot_data$marker_label <- factor(
  top4_marker_labels_cc[curve_plot_data$marker],
  levels = unname(top4_marker_labels_cc[top4_markers_cc])
)
curve_plot_data$genotype_label <- factor(
  paste0("g=", curve_plot_data$genotype),
  levels = paste0("g=", sort(unique(curve_plot_data$genotype)))
)
# Base-genotype labels per marker, decoded from the GABI 90K marker metadata
# (FIP1_308_final_dataset_direct_GABI/FIP1_308_SNP_marker_info_from_GABI90K.csv;
# dosage 0/1/2 = allele_1 homozygote / heterozygote / allele_2 homozygote).
# FRO2 has only the two homozygous groups (274/33).
genotype_base_labels <- list(
  RFL_Contig2772_1693 = c("0" = "CC", "1" = "CA", "2" = "AA"),
  BS00064448_51 = c("0" = "AA", "1" = "AG", "2" = "GG"),
  Tdurum_contig42092_348 = c("0" = "AA", "1" = "AG", "2" = "GG"),
  wsnp_Ex_c6548_11355524 = c("0" = "GG", "1" = "GT", "2" = "TT")
)
curve_plot_data$allele_label <- mapply(
  function(m, g) genotype_base_labels[[m]][[as.character(g)]],
  curve_plot_data$marker, curve_plot_data$genotype,
  USE.NAMES = FALSE, SIMPLIFY = TRUE
)
curve_end_labels <- do.call(
  rbind,
  lapply(
    split(curve_plot_data,
          list(curve_plot_data$marker, curve_plot_data$genotype)),
    function(d) d[which.max(d$index), , drop = FALSE]
  )
)
gsd_plot_data <- top4_gsd_cc$curve
gsd_plot_data$marker_label <- factor(
  top4_marker_labels_cc[gsd_plot_data$marker],
  levels = unname(top4_marker_labels_cc[top4_markers_cc])
)

# Stack the two views in one figure: adjusted genotype-mean curves on top,
# conditional genetic SD below, sharing the marker strips and the x axis.
view_levels <- c("Genotype mean", "Genetic SD")
genotype_levels <- levels(curve_plot_data$genotype_label)
combined_plot_data <- rbind(
  data.frame(
    marker_label = curve_plot_data$marker_label,
    index = curve_plot_data$index,
    y = curve_plot_data$fitted,
    genotype_label = curve_plot_data$genotype_label,
    view = factor(view_levels[1L], levels = view_levels)
  ),
  data.frame(
    marker_label = gsd_plot_data$marker_label,
    index = gsd_plot_data$index,
    y = gsd_plot_data$genetic_sd,
    genotype_label = factor(NA_character_, levels = genotype_levels),
    view = factor(view_levels[2L], levels = view_levels)
  )
)
# Allele labels sit at a fixed spot in the right margin of each facet,
# stacked vertically in genotype order (colors match the curves), so
# converging curve ends can never produce overlapping labels.
y_min <- min(curve_plot_data$fitted)
y_rng <- diff(range(curve_plot_data$fitted))
curve_end_labels$x_label <- max(curve_plot_data$index) * 1.06
curve_end_labels <- curve_end_labels[order(curve_end_labels$marker_label,
                                           curve_end_labels$genotype), ]
# Within-facet stack index via ave(): values stay aligned with the row
# order above (split() would reorder groups alphabetically and misassign).
stack_idx <- ave(seq_len(nrow(curve_end_labels)), curve_end_labels$marker,
                 FUN = seq_along)
curve_end_labels$y <- y_min + y_rng * (0.04 + 0.075 * (stack_idx - 1))
curve_end_labels$view <- factor(view_levels[1L], levels = view_levels)
hline_data <- data.frame(
  yintercept = 0,
  view = factor(view_levels[2L], levels = view_levels)
)
panel_fg <- ggplot2::ggplot(
  combined_plot_data,
  ggplot2::aes(x = .data$index, y = .data$y)
) +
  ggplot2::geom_hline(
    data = hline_data, ggplot2::aes(yintercept = .data$yintercept),
    color = "#9CA3AF", linetype = "dashed", linewidth = 0.35
  ) +
  ggplot2::geom_line(
    data = combined_plot_data[combined_plot_data$view == view_levels[1L], ],
    ggplot2::aes(color = .data$genotype_label, group = .data$genotype_label),
    linewidth = 0.9
  ) +
  ggplot2::geom_line(
    data = combined_plot_data[combined_plot_data$view == view_levels[2L], ],
    color = "#7B2CBF", linewidth = 0.9
  ) +
  # Fixed-position allele stack in the bottom-right margin of each facet.
  ggplot2::geom_text(
    data = curve_end_labels,
    ggplot2::aes(x = .data$x_label, y = .data$y,
                 label = .data$allele_label, color = .data$genotype_label),
    hjust = 0, size = 2.8, fontface = "bold", show.legend = FALSE
  ) +
  ggplot2::facet_grid(view ~ marker_label, scales = "free_y") +
  ggplot2::scale_color_manual(
    values = c("g=0" = "#0072B2", "g=1" = "#E69F00", "g=2" = "#009E73")
  ) +
  ggplot2::scale_x_continuous(
    breaks = function(x) pretty(x, 3),
    expand = ggplot2::expansion(mult = c(0.01, 0.24))
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(x = NULL, y = NULL) +
  panel_theme() +
  ggplot2::theme(
    legend.position = "none",
    strip.text.x = ggplot2::element_text(size = 9, face = "plain"),
    strip.text.y = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "#F3F4F6", color = "#D1D5DB"),
    strip.background.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(5.5, 14, 5.5, 5.5)
  )

save_panel(panel_fg, paste0(stem, "_fg_top4_genotype_curves_and_genetic_sd"),
           7.2, 3.1)

cat("FIP1 canopy-cover top-four interpretability panels written to: ",
    normalizePath(figure_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
