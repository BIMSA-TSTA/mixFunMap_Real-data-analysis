#!/usr/bin/env Rscript

# 03_staph_mic0.R — S. aureus MIC=0 growth-curve GWAS (simplified
# reproduction script for the mixFunMap manuscript).
#
# Two parts:
#   1. Analysis (mode = "full"): reruns the three-method genome-wide scan
#      from the raw strain data in inputs/staph_mic0/ (99 strains, 14 time
#      points, haploid 0/1 SNPs). Compact version of
#      real_data_analysis/run_staph_mic0_logistic_gwas.R with
#      STAPH_MEAN_MODEL="standard_logistic", STAPH_STRUCTURE_MODEL="pca3".
#   2. Figures (default): redraws every manuscript panel from the frozen
#      outputs in results/ without refitting any model.
#
# Usage:
#   Rscript 03_staph_mic0.R            # figures only (fast; default)
#   Rscript 03_staph_mic0.R full [cores] [output_dir]
#
# Model specification (identical to the manuscript):
#   * Mean curve: standard logistic mu(t) = A / (1 + exp(-K (t - T0))).
#   * SNP QC: call rate >= 0.95 and MAF >= 0.05; mean imputation.
#   * Q = PC1-PC3 from LD-pruned SNPs (50-SNP window, r2 <= 0.2);
#     K = VanRaden kinship from all imputed SNPs.
#   * mixFunMap: Q + K, fix_tau3 = FALSE, 3-df P3D Wald joint test.
#   * Functional Mapping: genotypic 0/1 model, SAD(1) residual, no Q/K.
#   * GMMAT-minP: per-time-point Q + K score tests combined as
#     min(1, 14 * min_t p_t).
#
# Thresholds: Bonferroni 0.05/M and suggestive 1/M with M = scanned SNPs.
#
# Requirements:
#   figures : R >= 4.4 with ggplot2, ggrepel, readxl
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

inputs_dir <- file.path(script_dir, "inputs", "staph_mic0")
results_dir <- file.path(script_dir, "results", "staph_mic0")
figure_dir <- file.path(script_dir, "figures")
pheno_path <- file.path(inputs_dir, "growdata_mic0.csv")
annotation_path <- file.path(
  script_dir, "results", "annotation",
  "Table_S_GWAS_significant_loci_annotation_v2.xlsx"
)

# ---------------------------------------------------------------------------
# Part 1: analysis (only runs in mode = "full")
# ---------------------------------------------------------------------------

run_analysis <- function(output_dir = file.path(results_dir, "rerun"),
                         cores = max(1L, parallel::detectCores(logical = FALSE)),
                         maf_min = 0.05) {
  if (!requireNamespace("mixFunMap", quietly = TRUE)) {
    stop("The mixFunMap package is required for the analysis mode.", call. = FALSE)
  }
  if (!requireNamespace("GMMAT", quietly = TRUE)) {
    stop("GMMAT is required for the minP scan.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ## Raw data --------------------------------------------------------------
  phenotype_source <- utils::read.csv(file.path(inputs_dir, "growdata_mic0.csv"),
                                      check.names = FALSE)
  genotype_source <- utils::read.csv(file.path(inputs_dir, "snpdata.csv"),
                                     check.names = FALSE)
  sample_ids <- colnames(genotype_source)[-(1:2)]
  index <- as.numeric(phenotype_source[[2L]])
  pheno <- as.matrix(phenotype_source[, -(1:2), drop = FALSE])
  storage.mode(pheno) <- "double"
  colnames(pheno) <- sample_ids

  snp_order <- as.integer(genotype_source[[1L]])
  snp_position <- as.numeric(genotype_source[[2L]])
  geno_raw <- as.matrix(genotype_source[, -(1:2), drop = FALSE])
  storage.mode(geno_raw) <- "double"
  colnames(geno_raw) <- sample_ids
  stopifnot(all(is.na(geno_raw) | geno_raw %in% c(0, 1)))  # haploid 0/1 calls

  ## SNP QC: call rate >= 0.95, MAF >= 0.05; mean imputation ---------------
  marker_call_rate <- rowMeans(!is.na(geno_raw))
  allele_frequency <- rowMeans(geno_raw, na.rm = TRUE)
  maf <- pmin(allele_frequency, 1 - allele_frequency)
  keep <- marker_call_rate >= 0.95 & is.finite(maf) & maf >= maf_min
  marker_id <- sprintf("SNP_%05d_pos_%d", snp_order, as.integer(snp_position))
  geno_hardcall <- geno_raw[keep, , drop = FALSE]
  rownames(geno_hardcall) <- marker_id[keep]
  mean_impute <- function(x) {
    for (i in which(rowSums(is.na(x)) > 0)) {
      x[i, is.na(x[i, ])] <- mean(x[i, ], na.rm = TRUE)
    }
    x
  }
  geno <- mean_impute(geno_hardcall)                 # scan dosage matrix
  geno_all <- mean_impute(geno_raw)                  # all SNPs, for K
  rownames(geno_all) <- marker_id
  map <- data.frame(marker = marker_id[keep], chr = "chromosome",
                    pos = snp_position[keep], stringsAsFactors = FALSE)

  ## Q = PC1-PC3 from LD-pruned SNPs (window 50, r2 <= 0.2) ----------------
  pca_order <- order(map$pos)
  centered <- geno[pca_order, ] - rowMeans(geno[pca_order, ])
  z <- centered / sqrt(rowSums(centered^2))
  pca_keep <- logical(nrow(z))
  active <- integer()
  for (i in seq_len(nrow(z))) {
    active <- active[active >= i - 50L + 1L]
    correlated <- length(active) > 0L &&
      any((z[i, , drop = FALSE] %*% t(z[active, , drop = FALSE]))^2 > 0.20)
    if (!correlated) {
      pca_keep[i] <- TRUE
      active <- c(active, i)
    }
  }
  pca_fit <- stats::prcomp(t(geno[pca_order, ][pca_keep, ]),
                           center = TRUE, scale. = TRUE)
  Q <- as.matrix(pca_fit$x[, 1:3, drop = FALSE])
  colnames(Q) <- paste0("PC", 1:3)
  rownames(Q) <- sample_ids

  ## K = VanRaden kinship from all imputed SNPs ----------------------------
  K <- mixFunMap::kinship_vanraden(geno_all)

  genomewide_alpha <- 0.05 / nrow(geno)

  ## 1. mixFunMap: standard-logistic Q + K null, 3-df P3D Wald scan --------
  fit <- mixFunMap::fit_mixfunmap(
    pheno = pheno, index = index, Q = Q, K = K,
    mean = mixFunMap::mean_logistic(), engine = "auto",
    include_intercept = TRUE, fix_tau3 = FALSE,
    max_outer = 60L, max_reml_optim = 150L, tol_theta = 0.01,
    reml_every = 1L, verbose = FALSE
  )
  if (!isTRUE(fit$usable)) stop("mixFunMap null model is not usable.", call. = FALSE)
  scan1 <- mixFunMap::scan_mixfunmap(
    fit, geno = geno, test = "wald", cores = cores,
    max_outer_snp = 12L, tol_theta = 0.025
  )
  beta <- as.matrix(scan1$beta)
  retry <- !is.finite(beta[, "pval"]) | beta[, "pval"] <= 0 |
    beta[, "pval"] > 1 | beta[, "converged"] <= 0 | beta[, "outer_converged"] <= 0
  if (any(retry)) {
    scan2 <- mixFunMap::scan_mixfunmap(
      fit, geno = geno[retry, , drop = FALSE], test = "wald", cores = cores,
      max_outer_snp = 40L, tol_theta = 0.01
    )
    beta[match(rownames(as.matrix(scan2$beta)), rownames(beta)),
         colnames(as.matrix(scan2$beta))] <- as.matrix(scan2$beta)
  }
  table_mixfunmap <- data.frame(
    map, stat = as.numeric(beta[map$marker, "wald3df"]), df = 3L,
    pval = as.numeric(beta[map$marker, "pval"]),
    converged = as.numeric(beta[map$marker, "converged"]) > 0,
    outer_converged = as.numeric(beta[map$marker, "outer_converged"]) > 0,
    method = "mixfunmap", stringsAsFactors = FALSE
  )

  ## 2. Functional Mapping: logistic + SAD(1), genotypic 0/1, no Q/K -------
  null_ofm <- mixFunMap::fit_ordinary_funmap(pheno, index, maxit = 300L,
                                             reltol = 1e-8)
  table_funmap <- mixFunMap::scan_ordinary_funmap(
    pheno = pheno, geno = geno_hardcall, index = index, marker_map = map,
    null_fit = null_ofm, genetic_model = "genotypic", genotype_levels = c(0, 1),
    missing_genotype = "omit", sparse_group = "omit", min_group_n = 5L,
    maxit = 300L, reltol = 1e-8
  )
  table_funmap <- as.data.frame(table_funmap, stringsAsFactors = FALSE)
  table_funmap$outer_converged <- NA
  table_funmap$method <- "ordinary_funmap"

  ## 3. GMMAT-minP: 14 per-time-point Q + K score tests --------------------
  table_minp <- mixFunMap::scan_minp_gmmat(
    pheno = pheno, geno = geno, index = index, Q = Q, K = K,
    marker_map = map, correction = "bonferroni", maxiter = 200L, tol = 1e-5,
    primary_optimizer = "AI", fallback_optimizer = "Brent",
    verbose = FALSE, ncores = cores
  )
  table_minp <- as.data.frame(table_minp, stringsAsFactors = FALSE)
  table_minp$converged <- NA
  table_minp$outer_converged <- NA
  table_minp$method <- "minp"

  ## Shared output tables ---------------------------------------------------
  keep_cols <- c("marker", "chr", "pos", "pval", "converged",
                 "outer_converged", "method")
  combined <- rbind(table_mixfunmap[, keep_cols],
                    table_funmap[, keep_cols],
                    table_minp[, keep_cols])
  combined$valid <- is.finite(combined$pval) & combined$pval > 0 &
    combined$pval <= 1 & (is.na(combined$converged) | combined$converged) &
    (combined$method != "mixfunmap" |
       (!is.na(combined$outer_converged) & combined$outer_converged))

  lambda_gc <- function(pval, df) {
    ok <- is.finite(pval) & pval > 0 & pval <= 1
    if (!any(ok) || is.na(df)) return(NA_real_)
    stats::median(stats::qchisq(1 - pval[ok], df = df)) /
      stats::qchisq(0.5, df = df)
  }
  summary_table <- do.call(rbind, lapply(
    list(table_mixfunmap, table_funmap, table_minp),
    function(x) {
      valid <- combined$valid[combined$method == x$method[1L]]
      p <- x$pval[valid]
      top <- which(valid)[which.min(p)]
      df <- if (x$method[1L] == "minp") NA_integer_ else 3L
      data.frame(
        method = x$method[1L], n_requested = nrow(x), n_valid = sum(valid),
        lambda_df = df, lambda_gc = lambda_gc(x$pval[valid], df),
        top_marker = x$marker[top], top_pos = x$pos[top],
        top_p = x$pval[top],
        n_bonferroni_significant = sum(valid & x$pval < genomewide_alpha),
        bonferroni_alpha = genomewide_alpha,
        stringsAsFactors = FALSE
      )
    }
  ))

  con <- gzfile(file.path(output_dir, "three_method_results.tsv.gz"), "wt")
  utils::write.table(combined, con, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "NA")
  close(con)
  utils::write.table(summary_table, file.path(output_dir, "method_summary.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  saveRDS(fit, file.path(output_dir, "mixfunmap_null_fit.rds"),
          compress = "xz", version = 3)
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
}

for (pkg in c("ggplot2", "ggrepel", "readxl")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(pkg, " is required for the Staph MIC=0 figures.", call. = FALSE)
  }
}
for (path in c(pheno_path, annotation_path,
               file.path(results_dir, "three_method_results.tsv.gz"),
               file.path(results_dir, "method_summary.tsv"))) {
  if (!file.exists(path)) stop("Required input not found: ", path, call. = FALSE)
}

# ---------------------------------------------------------------------------
# Shared style (consistent with the wheat and simulation figures)
# ---------------------------------------------------------------------------

method_colors <- c(mixfunmap = "#0072B2", ordinary_funmap = "#D55E00", minp = "#009E73")

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

scan <- read_tsv(file.path(results_dir, "three_method_results.tsv.gz"))
method_summary <- read_tsv(file.path(results_dir, "method_summary.tsv"))
lambdas <- stats::setNames(method_summary$lambda_gc, method_summary$method)

grow <- utils::read.csv(pheno_path, check.names = FALSE)

M <- length(unique(scan$marker))
bonferroni_alpha <- method_summary$bonferroni_alpha[1L]
threshold_lines <- data.frame(
  name = factor(
    c(sprintf("Bonferroni (0.05/M = 0.05/%d)", M),
      sprintf("Suggestive (1/M = 1/%d)", M)),
    levels = c(sprintf("Bonferroni (0.05/M = 0.05/%d)", M),
               sprintf("Suggestive (1/M = 1/%d)", M))
  ),
  y = -log10(c(bonferroni_alpha, 1 / M))
)
threshold_colors <- stats::setNames(c("#D62728", "#7F7F7F"), levels(threshold_lines$name))

read_annotation <- function(sheet, method = NULL) {
  df <- as.data.frame(
    readxl::read_xlsx(annotation_path, sheet = sheet, skip = 4),
    stringsAsFactors = FALSE
  )
  names(df) <- trimws(names(df))
  df <- df[!is.na(df[["Significance tier"]]) & nzchar(df[["Significance tier"]]), ]
  if (!is.null(method) && "Method" %in% names(df)) {
    df <- df[df$Method == method, ]
  }
  df$`Position (bp)` <- as.numeric(df$`Position (bp)`)
  df$`P value` <- as.numeric(df$`P value`)
  df
}

# One label per locus tag: anchor at the marker with the smallest p value.
collapse_by_locus <- function(df) {
  df <- df[order(df$`P value`), , drop = FALSE]
  df[!duplicated(df$`Locus tag`), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Panel a: growth-curve dynamics of all strains
# ---------------------------------------------------------------------------

make_panel_a <- function() {
  hours <- as.numeric(grow$t)
  od <- as.matrix(grow[, grepl("^X", names(grow)), drop = FALSE])
  storage.mode(od) <- "double"

  curve_df <- data.frame(
    time = rep(hours, ncol(od)),
    od = as.vector(od),
    strain = rep(seq_len(ncol(od)), each = nrow(od))
  )
  mean_df <- data.frame(time = hours, od = rowMeans(od, na.rm = TRUE))

  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = curve_df,
      ggplot2::aes(x = .data$time, y = .data$od, group = .data$strain),
      linewidth = 0.25, alpha = 0.3, color = "#9CA3AF"
    ) +
    ggplot2::geom_line(
      data = mean_df, ggplot2::aes(x = .data$time, y = .data$od),
      linewidth = 0.8, color = "#0072B2"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    panel_theme()
}

# ---------------------------------------------------------------------------
# Panel b: three-method QQ plot (Functional Mapping and GMMAT-minP first,
# mixFunMap on top; lambda annotated for Functional Mapping and mixFunMap)
# ---------------------------------------------------------------------------

make_panel_b <- function() {
  method_levels <- c("mixfunmap", "ordinary_funmap", "minp")
  method_labels <- c(
    mixfunmap = sprintf("mixFunMap (\u03bb = %.3f)", lambdas[["mixfunmap"]]),
    ordinary_funmap = sprintf("Functional Mapping (\u03bb = %.2f)",
                              lambdas[["ordinary_funmap"]]),
    minp = "GMMAT-minP (index-adjusted p)"
  )
  qq_df <- do.call(rbind, lapply(method_levels, function(meth) {
    p <- sort(scan$pval[scan$method == meth & scan$valid &
                          is.finite(scan$pval) & scan$pval > 0])
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
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 1.6, alpha = 1)))
}

# ---------------------------------------------------------------------------
# Panels c-e: per-method Manhattan plots (single circular chromosome)
# ---------------------------------------------------------------------------

manhattan_panel <- function(method_key, annotation_df = NULL,
                            window_mb = 0.2) {
  mapped <- scan[scan$method == method_key & scan$valid &
                   is.finite(scan$pval) & scan$pval > 0, , drop = FALSE]
  mapped$neglog10p <- -log10(pmax(mapped$pval, 1e-300))
  # Alternate blocks along the single chromosome in the same palette as
  # the wheat subgenome coloring (A blue / B orange), so the position axis
  # reads in fixed-size blocks.
  block <- floor(mapped$pos / (window_mb * 1e6)) %% 2L
  mapped$pt_color <- c("#1F78B4", "#FF7F00")[block + 1L]

  p <- ggplot2::ggplot(mapped, ggplot2::aes(.data$pos / 1e6, .data$neglog10p)) +
    ggplot2::geom_point(color = mapped$pt_color, size = 0.5, na.rm = TRUE) +
    ggplot2::geom_hline(
      data = threshold_lines,
      ggplot2::aes(yintercept = .data$y, color = .data$name),
      linetype = "dashed", linewidth = 0.5
    ) +
    ggplot2::scale_color_manual(values = threshold_colors) +
    ggplot2::scale_x_continuous(breaks = seq(0, 3, 0.5), expand = c(0.01, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    panel_theme() +
    ggplot2::theme(
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 8)
    )

  if (!is.null(annotation_df) && nrow(annotation_df)) {
    ann <- annotation_df
    ann$x <- ann$`Position (bp)` / 1e6
    hit <- match(ann$Marker, mapped$marker)
    ann$y <- ifelse(is.na(hit), -log10(ann$`P value`), mapped$neglog10p[hit])
    p <- p +
      ggplot2::geom_point(
        data = ann, ggplot2::aes(.data$x, .data$y),
        color = "#D62728", size = 1.2
      ) +
      ggrepel::geom_text_repel(
        data = ann, ggplot2::aes(.data$x, .data$y, label = .data$`Locus tag`),
        size = 2.4, color = "#111111", max.overlaps = Inf,
        min.segment.length = 0, segment.size = 0.25, segment.color = "#6B7280",
        box.padding = 0.45, point.padding = 0.25, seed = 42,
        force_pull = 0.5
      ) +
      ggplot2::coord_cartesian(
        ylim = c(0, max(mapped$neglog10p, ann$y, na.rm = TRUE) * 1.55),
        clip = "off"
      )
  }
  p
}

# ---------------------------------------------------------------------------
# Build and save (widths: a + b == c == d == e)
# ---------------------------------------------------------------------------

stem <- "Figure3_staph_mic0"

panel_a <- make_panel_a()
panel_b <- make_panel_b()
panel_c <- manhattan_panel(
  "mixfunmap",
  collapse_by_locus(read_annotation("Table S2 Staph", "mixFunMap"))
)
panel_d <- manhattan_panel(
  "minp",
  collapse_by_locus(read_annotation("Table S2 Staph", "minP"))
)
panel_e <- manhattan_panel("ordinary_funmap")

save_panel(panel_a, paste0(stem, "_a_growth_dynamics"), 3.6, 2.4)
save_panel(panel_b, paste0(stem, "_b_qq"), 3.6, 2.4)
save_panel(panel_c, paste0(stem, "_c_manhattan_mixfunmap"), 7.2, 3.4)
save_panel(panel_d, paste0(stem, "_d_manhattan_minp"), 7.2, 3.4)
save_panel(panel_e, paste0(stem, "_e_manhattan_functional_mapping"), 7.2, 3.4)

cat("Staph MIC=0 figures written to: ",
    normalizePath(figure_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")

# ---------------------------------------------------------------------------
# Added interpretability panels: four manuscript-chosen mixFunMap loci
# (plot-only; reads the frozen estimate_genotype_curves() output)
# ---------------------------------------------------------------------------

# Manuscript-chosen interpretability loci (not simply the four smallest
# p values): fnbB (balanced Bonferroni locus), acetolactate synthase
# (highest-MAF suggestive locus), gidA (flagship nonsynonymous L->F) and
# PI-PLC (nonsynonymous D->G, same direction as gidA).
top4_markers_interpretability <- c(
  "SNP_22131_pos_2579025",  # fnbB
  "SNP_19912_pos_2291841",  # acetolactate synthase
  "SNP_25160_pos_2818604",  # gidA
  "SNP_00892_pos_54481"     # PI-PLC
)
top4_marker_labels <- c(
  SNP_22131_pos_2579025 = "fnbB",
  SNP_19912_pos_2291841 = "ALS",
  SNP_25160_pos_2818604 = "gidA",
  SNP_00892_pos_54481 = "PI-PLC"
)

# Frozen output of mixFunMap::estimate_genotype_curves(..., background =
# "refit") + estimate_genetic_sd(method = "conditional") for the four loci
# (produced by the original pipeline; refitting takes several minutes).
interpretability_path <- file.path(
  script_dir, "results", "staph_mic0", "staph_mic0_top4_interpretability.rds"
)
top4_frozen <- readRDS(interpretability_path)
top4_curves_interpretability <- top4_frozen$genotype_curves
top4_gsd_interpretability <- top4_frozen$genetic_sd
stopifnot(all(top4_markers_interpretability %in%
                top4_curves_interpretability$curves$marker))

curve_plot_data <- top4_curves_interpretability$curves
curve_plot_data$marker_label <- factor(
  top4_marker_labels[curve_plot_data$marker],
  levels = unname(top4_marker_labels[top4_markers_interpretability])
)
curve_plot_data$genotype_label <- factor(
  paste0("g=", curve_plot_data$genotype),
  levels = paste0("g=", sort(unique(curve_plot_data$genotype)))
)
# Base labels per marker from the reference annotation workbook
# (code/金黄色葡萄球菌（奥维森）.xlsx, column "ref_base<->sample_base");
# genotype 0 = reference base, 1 = alternative base (haploid).
# The ref/alt bases differ between markers, so a single shared legend would
# be ambiguous; the base is printed at the right end of each curve instead.
genotype_base_labels <- list(
  SNP_22131_pos_2579025 = c("0" = "A", "1" = "G"),   # fnbB
  SNP_19912_pos_2291841 = c("0" = "A", "1" = "T"),   # acetolactate synthase
  SNP_25160_pos_2818604 = c("0" = "C", "1" = "A"),   # gidA
  SNP_00892_pos_54481 = c("0" = "A", "1" = "G")      # PI-PLC
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
# Spread labels vertically within a marker when the curve ends nearly
# coincide, otherwise the base texts would overlap.
label_min_gap <- diff(range(curve_plot_data$fitted)) * 0.05
curve_end_labels <- do.call(
  rbind,
  lapply(split(curve_end_labels, curve_end_labels$marker), function(d) {
    d <- d[order(d$fitted), , drop = FALSE]
    if (nrow(d) > 1L) {
      for (i in 2:nrow(d)) {
        if (d$fitted[i] - d$fitted[i - 1L] < label_min_gap) {
          d$fitted[i] <- d$fitted[i - 1L] + label_min_gap
        }
      }
    }
    d
  })
)
gsd_plot_data <- top4_gsd_interpretability$curve
gsd_plot_data$marker_label <- factor(
  top4_marker_labels[gsd_plot_data$marker],
  levels = unname(top4_marker_labels[top4_markers_interpretability])
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
curve_end_labels$y <- curve_end_labels$fitted
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
  ggplot2::geom_text(
    data = curve_end_labels,
    ggplot2::aes(label = .data$allele_label, color = .data$genotype_label),
    hjust = 0, nudge_x = diff(range(curve_plot_data$index)) * 0.02,
    size = 2.8, fontface = "bold", show.legend = FALSE
  ) +
  ggplot2::facet_grid(view ~ marker_label, scales = "free_y") +
  ggplot2::scale_color_manual(
    values = c("g=0" = "#0072B2", "g=1" = "#D55E00")
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

cat("Staph MIC=0 top-four interpretability panels written to: ",
    normalizePath(figure_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")

# ---------------------------------------------------------------------------
# Panel h: two-ring inward-facing circular Manhattan (single circular
#   chromosome) -- outer ring = Functional Mapping, inner ring = GMMAT-minP;
#   peaks point toward the center. The chromosome is colored in alternating
#   0.2 Mb blocks (same blue/orange palette as the wheat subgenome style and
#   the linear Manhattan panels), with a black ideogram band, tangential Mb
#   position labels, per-ring Bonferroni (red dashed) and suggestive (grey
#   dotted) threshold circles, red highlights for loci passing the
#   suggestive threshold, and a small radial -log10(p) axis per ring in the
#   top gap.
# ---------------------------------------------------------------------------

make_panel_h <- function() {
  fm_scan <- scan[scan$method == "ordinary_funmap" & scan$valid &
                    is.finite(scan$pval) & scan$pval > 0, , drop = FALSE]
  minp_scan <- scan[scan$method == "minp" & scan$valid &
                      is.finite(scan$pval) & scan$pval > 0, , drop = FALSE]
  chr_len <- max(scan$pos[scan$valid], na.rm = TRUE)
  start_gap <- 0.030 * chr_len
  x_total <- chr_len + start_gap

  prep <- function(df, task) {
    df <- df[, c("marker", "pos", "pval")]
    df$task <- task
    df$neglog10p <- -log10(pmax(df$pval, 1e-300))
    df$x <- start_gap + df$pos
    block <- floor(df$pos / 2e5) %% 2L
    df$pt_color <- c("#1F78B4", "#FF7F00")[block + 1L]
    df
  }
  dat <- rbind(prep(fm_scan, "fm"), prep(minp_scan, "minp"))

  # Ring geometry: baselines are circles; points extend INWARD.
  hole <- 7
  ring_gap <- 2.2
  ring_max <- c(fm = max(dat$neglog10p[dat$task == "fm"]),
                minp = max(dat$neglog10p[dat$task == "minp"]))
  ring_base <- c(
    minp = hole + ring_max[["minp"]],
    fm = hole + ring_max[["minp"]] + ring_gap + ring_max[["fm"]]
  )
  dat$y <- ring_base[dat$task] - dat$neglog10p

  # Curated minP loci (same annotation source as linear panel d) are
  # highlighted in red in the inner ring; the inflated Functional Mapping
  # scan is left unhighlighted, as in linear panel e.
  ann_minp <- collapse_by_locus(read_annotation("Table S2 Staph", "minP"))
  hl <- data.frame()
  if (nrow(ann_minp)) {
    hit <- match(ann_minp$Marker, minp_scan$marker)
    hl <- data.frame(
      x = start_gap + ann_minp$`Position (bp)`,
      y = ring_base[["minp"]] - ifelse(
        is.na(hit), -log10(pmax(ann_minp$`P value`, 1e-300)),
        -log10(pmax(minp_scan$pval[hit], 1e-300))
      )
    )
  }

  # Threshold circles per ring (Bonferroni red dashed, suggestive grey).
  thr_df <- do.call(rbind, lapply(c("fm", "minp"), function(task) {
    data.frame(task = task,
               name = levels(threshold_lines$name),
               y = ring_base[[task]] - threshold_lines$y)
  }))
  thr_df$linetype <- ifelse(startsWith(thr_df$name, "Bonferroni"),
                            "dashed", "dotted")

  # Single black ideogram band for the whole chromosome.
  ideo_df <- data.frame(
    xmin = start_gap, xmax = start_gap + chr_len,
    ymin = ring_base[["fm"]] + 1.5, ymax = ring_base[["fm"]] + 2.9
  )

  # Tangential Mb position labels every 0.5 Mb outside the ideogram band.
  mb_marks <- seq(0, floor(chr_len / 5e5) * 5e5, by = 5e5)
  pos_label_df <- data.frame(
    x = start_gap + mb_marks,
    y = ring_base[["fm"]] + 4.6,
    label = sprintf("%.1f", mb_marks / 1e6),
    angle = -360 * ((start_gap + mb_marks) / x_total)
  )

  # Small radial axis (0, 2, 4, ...) per ring inside the top gap.
  axis_df <- do.call(rbind, lapply(c("fm", "minp"), function(task) {
    step <- if (ring_max[[task]] > 20) 5 else 2
    v <- seq(0, floor(ring_max[[task]] / step) * step, by = step)
    data.frame(task = task, x = start_gap * 0.62,
               y = ring_base[[task]] - v, label = v)
  }))
  axis_x <- start_gap * 0.62
  axis_line_df <- data.frame(
    task = c("fm", "minp"),
    y0 = ring_base[c("fm", "minp")],
    y1 = ring_base[c("fm", "minp")] - ring_max[c("fm", "minp")]
  )
  tick_len <- 0.0018 * x_total

  ggplot2::ggplot(dat, ggplot2::aes(.data$x, .data$y)) +
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
      data = pos_label_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
                   angle = .data$angle),
      size = 2.8, color = "#111827"
    ) +
    ggplot2::coord_polar(theta = "x") +
    ggplot2::scale_y_continuous(limits = c(0, ring_base[["fm"]] + 6.2),
                                expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(
      legend.position = "none",
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}

panel_h <- make_panel_h()
save_panel(panel_h, paste0(stem, "_h_circular_manhattan"), 7.2, 6.4)

message("Staph MIC=0 circular Manhattan written to: ", figure_dir)
