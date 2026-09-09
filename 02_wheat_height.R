#!/usr/bin/env Rscript

# 02_wheat_height.R — FIP1 wheat FPWW012 plant-height GWAS (reproduction
# script for the mixFunMap manuscript).
#
# Two parts:
#   1. Analysis (mode = "full"): reruns the three-method genome-wide scan
#      from the analysis-ready inputs in inputs/wheat_height/, following
#      real_data_analysis/run_fip1_height_gwas.R.
#   2. Figures (default): redraws every manuscript panel from the frozen
#      outputs in results/ without refitting any model.
#
# Usage:
#   Rscript 02_wheat_height.R            # figures only (fast; default)
#   Rscript 02_wheat_height.R full [cores] [output_dir]
#                                       # rerun the analysis, then figures
#
# Model specification (identical to the manuscript):
#   * Mean curve: three-parameter logistic mu(t) = A / (1 + exp(-K (t - T0))).
#   * mixFunMap: Q = PC1-PC5 (fixed) + VanRaden K (random), random effects on
#     all three curve parameters (fix_tau3 = FALSE), 3-df P3D Wald joint
#     test on (A, K, T0); markers that fail convergence are refitted with
#     stricter iteration controls.
#   * Functional Mapping: ordinary FunMap, genotypic 0/1/2 model with SAD(1)
#     residual covariance, no Q/K (primary test: LRT).
#   * GMMAT-minP: per-time-point Q+K GLMM score tests combined across the
#     22 time points as min(1, 22 * min_t p_t).
#
# Thresholds: LD-aware M_eff = 3,270 (LD pruning r2 = 0.2) ->
# Bonferroni 0.05/M_eff = 1.53e-5; suggestive 1/M = 5.38e-5 with
# M = 18,583 (thresholds.rds, produced by the M_eff pruning step of the
# original pipeline).
#
# Requirements:
#   figures : R >= 4.4 with ggplot2, ggrepel, readxl
#   full    : additionally the mixFunMap package (repository root
#             ./mixFunMap, installed or on .libPaths) and GMMAT

options(stringsAsFactors = FALSE)

.args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(.args) >= 1L) tolower(.args[[1L]]) else "figures"

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_argument)) {
  dirname(normalizePath(sub("^--file=", "", file_argument[[1L]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

inputs_dir <- file.path(script_dir, "inputs", "wheat_height")
results_dir <- file.path(script_dir, "results", "wheat_height")
figure_dir <- file.path(script_dir, "figures")
pheno_path <- file.path(results_dir, "phenotype_height_FPWW012_matrix.tsv.gz")
thresholds_path <- file.path(results_dir, "thresholds.rds")
annotation_path <- file.path(
  script_dir, "results", "annotation",
  "Table_S_GWAS_significant_loci_annotation_v2.xlsx"
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

  # Analysis-ready inputs frozen by the data-curation step (raw GABI 90K
  # delivery -> QC -> dosage/hard-call matrices, PC1-PC5, VanRaden K):
  # pheno 22 x 307, geno 18,583 x 307 (mean-imputed dosages),
  # geno_hardcall 18,583 x 307 (0/1/2 with NA), Q 307 x 5, K 307 x 307.
  inputs <- readRDS(file.path(inputs_dir, "analysis_inputs_FPWW012_height.rds"))
  pheno <- as.matrix(inputs$pheno)
  index <- as.numeric(inputs$index)
  geno <- as.matrix(inputs$geno)
  Q <- as.matrix(inputs$Q)
  K <- as.matrix(inputs$K)
  marker_map <- as.data.frame(inputs$marker_map, stringsAsFactors = FALSE)
  map <- marker_map[, c("marker", "chr", "pos")]
  ordinary_inputs <- readRDS(file.path(
    inputs_dir, "analysis_inputs_FPWW012_height_ordinary3class_v2.rds"
  ))
  geno_hardcall <- as.matrix(ordinary_inputs$geno_hardcall)
  # Use the same LD-aware threshold as the manuscript figures.
  thresholds <- readRDS(thresholds_path)
  primary_threshold <- thresholds$thresholds[[paste0("meff_r2_", thresholds$primary_r2)]]
  genomewide_alpha <- primary_threshold$threshold

  ## 1. mixFunMap: logistic Q + K null, then a 3-df P3D Wald scan ----------
  fit <- mixFunMap::fit_mixfunmap(
    pheno = pheno, index = index, Q = Q, K = K,
    mean = mixFunMap::mean_logistic(), engine = "auto", fix_tau3 = FALSE,
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
    beta[match(rownames(second_beta <- as.matrix(scan2$beta)), rownames(beta)),
         colnames(second_beta)] <- second_beta
  }
  table_mixfunmap <- data.frame(
    marker_map[, c("marker", "chr", "pos", "mapped")],
    stat = as.numeric(beta[marker_map$marker, "wald3df"]),
    df = 3L,
    pval = as.numeric(beta[marker_map$marker, "pval"]),
    converged = as.numeric(beta[marker_map$marker, "converged"]) > 0,
    outer_converged = as.numeric(beta[marker_map$marker, "outer_converged"]) > 0,
    method = "mixfunmap",
    stringsAsFactors = FALSE
  )

  ## 2. Functional Mapping: logistic + SAD(1), genotypic 0/1/2, no Q/K -----
  null_ofm <- mixFunMap::fit_ordinary_funmap(pheno, index, maxit = 300L,
                                             reltol = 1e-8)
  table_funmap <- mixFunMap::scan_ordinary_funmap(
    pheno = pheno, geno = geno_hardcall, index = index, marker_map = map,
    null_fit = null_ofm, genetic_model = "genotypic", genotype_levels = 0:2,
    missing_genotype = "omit", sparse_group = "omit", min_group_n = 5L,
    maxit = 300L, reltol = 1e-8
  )
  table_funmap <- as.data.frame(table_funmap, stringsAsFactors = FALSE)
  table_funmap$mapped <- marker_map$mapped[match(table_funmap$marker,
                                                 marker_map$marker)]
  table_funmap$outer_converged <- NA
  table_funmap$method <- "ordinary_funmap"

  ## 3. GMMAT-minP: 22 per-time-point Q + K score tests --------------------
  table_minp <- mixFunMap::scan_minp_gmmat(
    pheno = pheno, geno = geno, index = index, Q = Q, K = K,
    marker_map = map, correction = "bonferroni", maxiter = 200L, tol = 1e-5,
    primary_optimizer = "AI", fallback_optimizer = "Brent",
    verbose = FALSE, ncores = cores
  )
  table_minp <- as.data.frame(table_minp, stringsAsFactors = FALSE)
  table_minp$mapped <- marker_map$mapped[match(table_minp$marker,
                                               marker_map$marker)]
  table_minp$converged <- NA
  table_minp$outer_converged <- NA
  table_minp$method <- "minp"

  ## Shared output tables ---------------------------------------------------
  keep_cols <- c("marker", "chr", "pos", "mapped", "pval", "converged",
                 "outer_converged", "method")
  combined <- rbind(table_mixfunmap[, keep_cols],
                    table_funmap[, keep_cols],
                    table_minp[, keep_cols])
  combined$valid <- is.finite(combined$pval) & combined$pval > 0 &
    combined$pval <= 1 & (is.na(combined$converged) | combined$converged) &
    (combined$method != "mixfunmap" |
       (!is.na(combined$outer_converged) & combined$outer_converged))

  lambda_gc <- function(pval, df) {
    ok <- is.finite(pval) & pval > 0 & pval <= 1 & is.finite(df)
    if (!any(ok)) return(NA_real_)
    stats::median(stats::qchisq(1 - pval[ok], df = df[ok])) /
      stats::qchisq(0.5, df = df[ok][1L])
  }
  summary_table <- do.call(rbind, lapply(
    list(table_mixfunmap, table_funmap, table_minp),
    function(x) {
      valid <- combined$valid[combined$method == x$method[1L]]
      p <- x$pval[valid]
      top <- if (length(p)) which(valid)[which.min(p)] else NA_integer_
      data.frame(
        method = x$method[1L], n_markers_requested = nrow(x),
        n_valid = sum(valid), valid_fraction = mean(valid),
        lambda_gc = lambda_gc(x$pval[valid], x$df[valid]),
        top_marker = x$marker[top], top_chr = x$chr[top],
        top_pos = x$pos[top], top_p = x$pval[top],
        n_genomewide_significant = sum(valid & x$pval < genomewide_alpha),
        genomewide_alpha = genomewide_alpha,
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
    stop(pkg, " is required for the FIP1 height figures.", call. = FALSE)
  }
}
for (path in c(pheno_path, thresholds_path, annotation_path,
               file.path(results_dir, "three_method_results.tsv.gz"),
               file.path(results_dir, "method_summary.tsv"))) {
  if (!file.exists(path)) stop("Required input not found: ", path, call. = FALSE)
}

# ---------------------------------------------------------------------------
# Shared style (consistent with the simulation figures)
# ---------------------------------------------------------------------------

method_colors <- c(mixfunmap = "#0072B2", ordinary_funmap = "#D55E00", minp = "#009E73")
method_shapes <- c(mixfunmap = 16, ordinary_funmap = 17, minp = 15)
method_linetypes <- c(mixfunmap = "solid", ordinary_funmap = "dashed", minp = "dotdash")

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
thr <- readRDS(thresholds_path)
pheno_tab <- read_tsv(pheno_path)

lambdas <- stats::setNames(method_summary$lambda_gc,
                           method_summary$method)

M <- thr$M
primary <- thr$thresholds[[paste0("meff_r2_", thr$primary_r2)]]
threshold_lines <- data.frame(
  name = factor(
    c(sprintf("Bonferroni (0.05/M_eff = 0.05/%d)", primary$M_eff),
      sprintf("Suggestive (1/M = 1/%d)", M)),
    levels = c(sprintf("Bonferroni (0.05/M_eff = 0.05/%d)", primary$M_eff),
               sprintf("Suggestive (1/M = 1/%d)", M))
  ),
  y = -log10(c(primary$threshold, thr$thresholds$suggestive$threshold))
)
threshold_colors <- stats::setNames(c("#D62728", "#7F7F7F"), levels(threshold_lines$name))

# Significant-locus annotation from the manuscript xlsx (Table S1 Wheat).
annotation <- as.data.frame(
  readxl::read_xlsx(annotation_path, sheet = "Table S1 Wheat", skip = 4),
  stringsAsFactors = FALSE
)
names(annotation) <- trimws(names(annotation))
annotation <- annotation[!is.na(annotation$Method) & nzchar(annotation$Method), ]
annotation <- annotation[annotation$Chr != "unmapped", , drop = FALSE]
annotation$`Position (bp)` <- as.numeric(annotation$`Position (bp)`)
annotation$`P value` <- as.numeric(annotation$`P value`)

# Short functional labels for annotated genes; genes whose annotation is
# "Uncharacterized protein" fall back to the gene ID. When the same short
# name maps to more than one locus (homoeologs such as TMK1 on 4B/4D),
# the chromosome is appended to distinguish the loci.
gene_short_names <- c(
  TraesCS4B02G049800 = "TMK1",
  TraesCS4D02G050200 = "TMK1",
  TraesCS3A02G425500 = "SPS1",
  TraesCS3A02G426200 = "HAC1",
  TraesCS3A02G424900 = "CPR5",
  TraesCS3A02G440200 = "MFAP1",
  TraesCS3A02G440000 = "RHC1A",
  TraesCS5A02G279200 = "bHLH089",
  TraesCS1A02G428100 = "APT1"
)

# One label per gene: anchor at the marker with the smallest p value and
# label with the short functional name (chromosome-suffixed if shared).
collapse_by_gene <- function(df) {
  df <- df[order(df$`P value`), , drop = FALSE]
  df <- df[!duplicated(df$Gene), , drop = FALSE]
  df$label <- ifelse(df$Gene %in% names(gene_short_names),
                     unname(gene_short_names[df$Gene]), df$Gene)
  shared <- df$label %in% df$label[duplicated(df$label)]
  df$label[shared] <- paste0(df$label[shared], " (",
                             sub("chr", "", df$Chr[shared]), ")")
  df
}

# ---------------------------------------------------------------------------
# Panel a: plant-height dynamics of all samples
# ---------------------------------------------------------------------------

make_panel_a <- function() {
  days <- as.numeric(pheno_tab$day)
  height <- as.matrix(pheno_tab[, -(1:2), drop = FALSE])
  storage.mode(height) <- "double"

  curve_df <- data.frame(
    day = rep(days, ncol(height)),
    height = as.vector(height),
    sample = rep(seq_len(ncol(height)), each = nrow(height))
  )
  mean_df <- data.frame(day = days, height = rowMeans(height, na.rm = TRUE))

  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = curve_df,
      ggplot2::aes(x = .data$day, y = .data$height, group = .data$sample),
      linewidth = 0.25, alpha = 0.25, color = "#9CA3AF"
    ) +
    ggplot2::geom_line(
      data = mean_df, ggplot2::aes(x = .data$day, y = .data$height),
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
    p <- sort(scan$pval[scan$method == meth & is.finite(scan$pval) & scan$pval > 0])
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
# Panels c-e: per-method Manhattan plots
# ---------------------------------------------------------------------------

chroms21 <- paste0("chr", rep(1:7, each = 3), rep(c("A", "B", "D"), 7))

manhattan_panel <- function(method_key, annotation_df = NULL) {
  mapped <- scan[scan$method == method_key & scan$mapped & scan$chr %in% chroms21,
                 , drop = FALSE]
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
    ggplot2::scale_x_continuous(breaks = ticks, labels = sub("chr", "", chroms21),
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
    ann$x <- offsets[ann$Chr] + ann$`Position (bp)`
    hit <- match(ann$Marker, mapped$marker)
    ann$y <- ifelse(is.na(hit), -log10(ann$`P value`), mapped$neglog10p[hit])
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
      )
    p <- p + ggplot2::coord_cartesian(
      ylim = c(0, max(mapped$neglog10p, ann$y, na.rm = TRUE) * 1.35),
      clip = "off"
    )
  }
  p
}

# ---------------------------------------------------------------------------
# Build and save (widths: a + b == c == d == e)
# ---------------------------------------------------------------------------

stem <- "Figure2_fip1_height"

panel_a <- make_panel_a()
panel_b <- make_panel_b()
panel_c <- manhattan_panel(
  "mixfunmap",
  collapse_by_gene(annotation[annotation$Method == "mixFunMap", , drop = FALSE])
)
panel_d <- manhattan_panel(
  "minp",
  collapse_by_gene(annotation[annotation$Method == "minP", , drop = FALSE])
)
panel_e <- manhattan_panel("ordinary_funmap")

save_panel(panel_a, paste0(stem, "_a_phenotype_dynamics"), 3.6, 2.4)
save_panel(panel_b, paste0(stem, "_b_qq"), 3.6, 2.4)
save_panel(panel_c, paste0(stem, "_c_manhattan_mixfunmap"), 7.2, 3.4)
save_panel(panel_d, paste0(stem, "_d_manhattan_minp"), 7.2, 3.4)
save_panel(panel_e, paste0(stem, "_e_manhattan_functional_mapping"), 7.2, 3.4)

cat("FIP1 height figures written to: ",
    normalizePath(figure_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")

# ---------------------------------------------------------------------------
# Added interpretability panels: four manuscript-chosen mixFunMap loci
# (plot-only; reads the frozen estimate_genotype_curves() output)
# ---------------------------------------------------------------------------

# Manuscript-chosen interpretability loci (not simply the four smallest
# p values): the TMK1 homoeolog pair on 4B/4D, the negative-T0 SPS1 peak
# on 3A, and the large positive-T0 bHLH089 peak on 5A.
top4_markers_interpretability <- c(
  "RAC875_rep_c105718_304",   # TMK1, chr4B 38.3 Mb
  "RAC875_rep_c105718_585",   # TMK1 homoeolog, chr4D 26.0 Mb
  "Kukri_c31546_66",          # SPS1, chr3A 667.8 Mb
  "wsnp_Ex_c7168_12311649"    # bHLH089, chr5A 488.3 Mb
)
top4_marker_labels <- c(
  RAC875_rep_c105718_304 = "TMK1 (4B)",
  RAC875_rep_c105718_585 = "TMK1 (4D)",
  Kukri_c31546_66 = "SPS1 (3A)",
  wsnp_Ex_c7168_12311649 = "bHLH089 (5A)"
)

interpretability_path <- file.path(
  script_dir, "results", "wheat_height", "fip1_top4_interpretability.rds"
)
# Frozen output of mixFunMap::estimate_genotype_curves(..., background =
# "refit") + estimate_genetic_sd(method = "conditional") for the four loci
# (produced by the original pipeline; refitting takes several minutes).
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
# Base-genotype labels per marker, decoded from the GABI 90K marker metadata
# (FIP1_308_final_dataset_direct_GABI/FIP1_308_SNP_marker_info_from_GABI90K.csv;
# dosage 0/1/2 = allele_1 homozygote / heterozygote / allele_2 homozygote,
# see decode_diplotypes() in real_data_analysis/prepare_fip1_height.py).
# The two alleles differ between markers, so a single shared legend would be
# ambiguous; the base genotype is printed at the right end of each curve.
genotype_base_labels <- list(
  RAC875_rep_c105718_304 = c("0" = "CC", "1" = "CT", "2" = "TT"),
  RAC875_rep_c105718_585 = c("0" = "TT", "1" = "TC", "2" = "CC"),
  Kukri_c31546_66 = c("0" = "CC", "1" = "CT", "2" = "TT"),
  wsnp_Ex_c7168_12311649 = c("0" = "AA", "1" = "AG", "2" = "GG")
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
# coincide (e.g. SPS1), otherwise the base texts would overlap.
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

cat("FIP1 top-four interpretability panels written to: ",
    normalizePath(figure_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")

# ---------------------------------------------------------------------------
# Panel h: two-ring inward-facing circular Manhattan (21 chromosomes)
#   outer ring = Functional Mapping, inner ring = GMMAT-minP; peaks point
#   toward the center. Style follows the rice circular Manhattan (panel d of
#   Figure 4): subgenome colors, black ideogram band, tangential chromosome
#   labels, per-ring Bonferroni (red dashed) and suggestive (grey dotted)
#   threshold circles, red highlights for loci passing the suggestive
#   threshold, and a small radial -log10(p) axis per ring in the top gap.
# ---------------------------------------------------------------------------

make_panel_h <- function() {
  gap_frac <- 0.012
  fm_scan <- scan[scan$method == "ordinary_funmap" & scan$mapped &
                    scan$chr %in% chroms21, , drop = FALSE]
  minp_scan <- scan[scan$method == "minp" & scan$mapped &
                      scan$chr %in% chroms21, , drop = FALSE]
  chr_max <- tapply(fm_scan$pos, as.character(fm_scan$chr), max)[chroms21]
  gap <- gap_frac * sum(as.numeric(chr_max))
  offsets <- c(0, cumsum(as.numeric(chr_max) + gap)[-length(chr_max)])
  names(offsets) <- chroms21
  total_len <- sum(as.numeric(chr_max)) + length(chroms21) * gap
  start_gap <- 0.030 * total_len
  x_total <- total_len + start_gap
  x_shift <- stats::setNames(offsets + start_gap, chroms21)
  chr_ticks <- x_shift + as.numeric(chr_max) / 2

  subgenome_colors <- c(A = "#1F78B4", B = "#FF7F00", D = "#6A3D9A")
  prep <- function(df, task) {
    df <- df[is.finite(df$pval) & df$pval > 0,
             c("marker", "chr", "pos", "pval")]
    df$chr <- as.character(df$chr)
    df$task <- task
    df$neglog10p <- -log10(pmax(df$pval, 1e-300))
    df$x <- x_shift[df$chr] + df$pos
    df$pt_color <- unname(subgenome_colors[sub("^chr\\d+", "", df$chr)])
    df
  }
  dat <- rbind(prep(fm_scan, "fm"), prep(minp_scan, "minp"))

  # Ring geometry: each ring has a baseline circle; points extend INWARD
  # from the baseline. A center hole keeps the inner ring readable.
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
  ann_minp <- collapse_by_gene(
    annotation[annotation$Method == "minP", , drop = FALSE]
  )
  hl <- data.frame()
  if (nrow(ann_minp)) {
    hit <- match(ann_minp$Marker, minp_scan$marker)
    hl <- data.frame(
      x = x_shift[ann_minp$Chr] + ann_minp$`Position (bp)`,
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

  # Outermost black chromosome ideogram band with white gaps.
  ideo_df <- data.frame(
    xmin = x_shift, xmax = x_shift + as.numeric(chr_max),
    ymin = ring_base[["fm"]] + 1.5, ymax = ring_base[["fm"]] + 2.9
  )

  # Chromosome labels outside the ideogram band, rotated tangentially.
  chr_label_df <- data.frame(
    x = unname(chr_ticks),
    y = ring_base[["fm"]] + 4.6,
    label = paste0("Chr", sub("chr", "", chroms21)),
    angle = -360 * (unname(chr_ticks) / x_total)
  )

  # Small radial axis (0, 2, 4, ...) per ring inside the widened top gap.
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
      data = chr_label_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
                   angle = .data$angle),
      size = 3.0, color = "#111827"
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

message("FIP1 circular Manhattan written to: ", figure_dir)
