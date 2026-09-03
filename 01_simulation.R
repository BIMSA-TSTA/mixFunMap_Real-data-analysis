#!/usr/bin/env Rscript

# 01_simulation.R — formal V15 logistic simulation (reproduction script
# for the mixFunMap manuscript).
#
# Two parts:
#   1. Simulation (mode = "full"): reruns the formal V15 design — P3D-Wald,
#      relative-QTL direction (1, 1, -1), K-neutral causal markers, n = 100,
#      target family size = 5, Q/noise ratio = 0.10, K/noise ratio = 0.50,
#      Henderson REML; 28 null cells x 300 replicates + 44 alternative
#      cells x 200 replicates = 17,200 datasets. Implements the design of
#      Fig/run_all_simulations_v15.R (01_simulate_logistic.R +
#      04_mixfunmap_adapter.R + 05_logistic_benchmark.R +
#      07_summarise_results.R) with the same models and estimators.
#   2. Figures (default): redraws Figure 1 (three subplots) and Figures
#      S1-S5 from the frozen results/simulation/simulation_summary.rds
#      without rerunning anything.
#
# Usage:
#   Rscript 01_simulation.R                 # figures only (fast; default)
#   Rscript 01_simulation.R full [cores] [n_rep_null] [n_rep_alt]
#                                           # rerun the simulation, then
#                                           # figures (pass small replicate
#                                           # counts for a pilot run)
#
# Design (identical to the manuscript):
#   * Genotypes: two populations (Fst = 0.01) nested into families of
#     ~5 (family ICC = 0.15); m = 500 scan markers plus m_kinship = 1,000
#     kinship markers; MAF in [0.10, 0.45]. Q = population indicator;
#     K = VanRaden kinship projected off Q.
#   * Phenotypes: three-parameter logistic curves (A, K, T0) = (1, 0.65, 7)
#     at 10 time points in [1, 14]; the causal marker shifts the parameters
#     by theta0 * (1, 1, -1) * rho with rho in the effect grid; nuisance Q
#     and K (family-structured) effects calibrated to 0.10 and 0.50 of the
#     noise variance; SAD(1) residuals with phi = 0.55, SNR = 4.
#   * Methods: mixFunMap (3-df P3D Wald, fixed-tau2 null), ordinary FunMap
#     (genotypic LRT), GMMAT-minP (per-time-point Q + K score tests,
#     Bonferroni-combined).
#   * Metrics per replicate: FWER (any p < 0.05/m under the null), raw and
#     empirical-5%-FWER-matched power at the causal marker, causal-marker
#     rank, dynamic-effect RMSE.
#
# Requirements:
#   figures : R >= 4.4 with ggplot2 and patchwork
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
summary_dir <- file.path(script_dir, "results", "simulation", "summary")
figure_dir <- file.path(script_dir, "figures")

# ---------------------------------------------------------------------------
# Part 1: simulation (only runs in mode = "full")
# ---------------------------------------------------------------------------

## Frozen V15 configuration -------------------------------------------------
sim_cfg <- list(
  seed = 20260812L, alpha = 0.05,
  n = 100L, n_time = 10L, time_range = c(1, 14), target_snr = 4,
  reference_effect = 0.12, m = 500L, m_kinship = 1000L,
  scenarios = c("none", "Q", "K", "QK"),
  grids = list(n = c(50L, 100L, 200L), n_time = c(6L, 10L, 14L),
               snr = c(2, 4, 6), effect = c(0.05, 0.08, 0.10, 0.12, 0.15)),
  n_rep = list(null = 300L, alternative = 200L),
  genotype = list(maf_range = c(0.10, 0.45), causal_maf_range = c(0.25, 0.35),
                  population_fst = 0.01, target_family_size = 5L,
                  family_icc = 0.15, min_group_n = 5L),
  # K-neutral causal-marker balance profile
  balance = list(max_abs_cor_q = 0.05, max_population_maf_gap = 0.02,
                 family_reference_quantiles = c(0.10, 0.90),
                 k_pc_r2_reference_quantiles = c(0.10, 0.30),
                 n_k_pc = 5L, max_attempts = 2000L),
  theta0 = c(A = 1.00, K = 0.65, T0 = 7.00),
  qtl_relative_direction = c(A = 1.00, K = 1.00, T0 = -1.00),
  q_effect_direction = c(A = 0.12, K = 0.03, T0 = 0.45),
  tau2_direction = c(A = 0.010, K = 0.020, T0 = 0.240),
  phi = 0.55, numerical_zero_tau2 = 1e-6,
  calibration = list(target_ratio_q = 0.10, target_ratio_k = 0.50,
                     mc_draws = 500L, seed = 20260806L),
  min_valid_fraction = 0.95, bootstrap_reps = 2000L
)

## Genotype simulation: two populations nested in families -------------------
make_nested_structure <- function(n, target_family_size) {
  sample_ids <- paste0("id", seq_len(n))
  population <- sample(rep(1:2, each = n / 2), n, replace = FALSE)
  family <- character(n)
  for (pop in 1:2) {
    members <- which(population == pop)
    n_family <- max(2L, ceiling(length(members) / target_family_size))
    assignment <- sample(rep(seq_len(n_family), length.out = length(members)),
                         length(members), replace = FALSE)
    family[members] <- paste0("p", pop, "_f", assignment)
  }
  names(population) <- names(family) <- sample_ids
  Q <- matrix(as.numeric(population == 2L), ncol = 1L,
              dimnames = list(sample_ids, "Q1"))
  list(sample_ids = sample_ids, population = population, family = family, Q = Q)
}

center_family_probabilities <- function(raw_probability, family_size,
                                        target_probability, eps = 1e-8) {
  logits <- stats::qlogis(pmin(pmax(raw_probability, eps), 1 - eps))
  offset <- stats::uniroot(
    function(o) stats::weighted.mean(stats::plogis(logits + o), family_size) -
      target_probability,
    interval = c(-50, 50), tol = 1e-12
  )$root
  stats::plogis(logits + offset)
}

draw_nested_marker <- function(population, family, base_maf,
                               population_fst, family_icc) {
  population_levels <- sort(unique(population))
  concentration <- (1 - population_fst) / population_fst
  population_prob <- pmin(pmax(
    stats::rbeta(length(population_levels), base_maf * concentration,
                 (1 - base_maf) * concentration), 0.02), 0.98)
  fam_concentration <- (1 - family_icc) / family_icc
  family_levels <- unique(family)
  family_population <- vapply(family_levels, function(fam) {
    population[match(fam, family)]
  }, numeric(1))
  raw_family_prob <- pmin(pmax(vapply(family_levels, function(fam) {
    p <- population_prob[population_levels == family_population[fam_levels == fam]]
    stats::rbeta(1L, p * fam_concentration, (1 - p) * fam_concentration)
  }, numeric(1)), 0.02), 0.98)
  family_size <- table(family)
  family_prob <- raw_family_prob
  for (pop in population_levels) {
    selected <- which(family_population == pop)
    family_prob[selected] <- center_family_probabilities(
      raw_family_prob[selected],
      as.numeric(family_size[family_levels[selected]]),
      population_prob[population_levels == pop]
    )
  }
  names(family_prob) <- family_levels
  stats::rbinom(length(family), size = 1L, prob = family_prob[family])
}

simulate_nested_markers <- function(population, family, m, maf_range,
                                    population_fst, family_icc, min_group_n,
                                    prefix) {
  geno <- matrix(NA_real_, nrow = m, ncol = length(population),
                 dimnames = list(paste0(prefix, seq_len(m)),
                                 names(population)))
  base_maf <- stats::runif(m, maf_range[1], maf_range[2])
  for (j in seq_len(m)) {
    for (attempt in seq_len(1000L)) {
      marker <- draw_nested_marker(population, family, base_maf[j],
                                   population_fst, family_icc)
      if (sum(marker == 0) >= min_group_n && sum(marker == 1) >= min_group_n) {
        geno[j, ] <- marker
        break
      }
      if (attempt == 1000L) stop("Could not draw a polymorphic marker.")
    }
  }
  geno
}

stabilize_projected_kinship <- function(K, eps = 1e-6) {
  K <- (K + t(K)) / 2
  eig <- eigen(K, symmetric = TRUE)
  K <- eig$vectors %*% (pmax(eig$values, 0) * t(eig$vectors))
  K <- (K + t(K)) / 2
  diag(K) <- diag(K) + eps
  K / mean(diag(K))
}

project_kinship_off_q <- function(K, Q) {
  X <- cbind(1, Q)
  H <- diag(nrow(X)) - X %*% qr.solve(crossprod(X), t(X))
  out <- stabilize_projected_kinship(H %*% K %*% H)
  dimnames(out) <- dimnames(K)
  out
}

simulate_genotypes <- function(n, m, m_kinship, cfg) {
  structure <- make_nested_structure(n, cfg$genotype$target_family_size)
  geno <- simulate_nested_markers(
    structure$population, structure$family, m, cfg$genotype$maf_range,
    cfg$genotype$population_fst, cfg$genotype$family_icc,
    cfg$genotype$min_group_n, "snp"
  )
  kinship_geno <- simulate_nested_markers(
    structure$population, structure$family, m_kinship,
    cfg$genotype$maf_range, cfg$genotype$population_fst,
    cfg$genotype$family_icc, cfg$genotype$min_group_n, "ksnp"
  )
  K_raw <- mixFunMap::kinship_vanraden(kinship_geno)
  dimnames(K_raw) <- list(structure$sample_ids, structure$sample_ids)
  c(structure, list(geno = geno, kinship_geno = kinship_geno,
                    K = project_kinship_off_q(K_raw, structure$Q)))
}

## K-neutral balanced causal marker ------------------------------------------
kinship_pcs <- function(K, n_pc = 5L) {
  eig <- eigen((K + t(K)) / 2, symmetric = TRUE)
  keep <- utils::head(which(eig$values > max(eig$values) * 1e-8), n_pc)
  eig$vectors[, keep, drop = FALSE]
}

make_balanced_causal_marker <- function(g, cfg) {
  maf_range <- cfg$genotype$causal_maf_range
  bal <- cfg$balance
  # Reference ranges from the scan markers: family-MAF range and the R2
  # against the first K PCs (K-neutral: restricted to the lower quantiles).
  reference_geno <- g$geno
  pcs <- kinship_pcs(g$K, bal$n_k_pc)
  marker_maf <- rowMeans(reference_geno)
  eligible <- which(marker_maf >= maf_range[1] & marker_maf <= maf_range[2])
  if (length(eligible) < 10L) eligible <- seq_len(nrow(reference_geno))
  family_range <- vapply(eligible, function(j) {
    diff(range(tapply(reference_geno[j, ], g$family, mean)))
  }, numeric(1))
  k_r2 <- vapply(eligible, function(j) {
    summary(stats::lm(reference_geno[j, ] ~ pcs))$r.squared
  }, numeric(1))
  fam_ref <- stats::quantile(family_range, bal$family_reference_quantiles,
                             names = FALSE, type = 8)
  k_ref <- stats::quantile(k_r2, bal$k_pc_r2_reference_quantiles,
                           names = FALSE, type = 8)
  n_pop <- length(g$population) / 2
  counts <- seq.int(ceiling(maf_range[1] * n_pop), floor(maf_range[2] * n_pop))
  concentration <- (1 - cfg$genotype$family_icc) / cfg$genotype$family_icc
  for (attempt in seq_len(bal$max_attempts)) {
    carrier_count <- sample(counts, 1L)
    marker <- integer(length(g$population))
    for (pop in sort(unique(g$population))) {
      members <- which(g$population == pop)
      fam <- g$family[members]
      target_maf <- carrier_count / length(members)
      weight <- pmin(pmax(stats::rbeta(
        length(unique(fam)), target_maf * concentration,
        (1 - target_maf) * concentration), 1e-6), 1 - 1e-6)
      marker[sample(members, carrier_count, replace = FALSE,
                    prob = weight[match(fam, unique(fam))])] <- 1L
    }
    population_maf <- tapply(marker, g$population, mean)
    family_maf <- tapply(marker, g$family, mean)
    abs_cor_q <- suppressWarnings(abs(stats::cor(marker, g$Q[, 1L])))
    marker_k_r2 <- summary(stats::lm(marker ~ pcs))$r.squared
    if (mean(marker) >= maf_range[1] && mean(marker) <= maf_range[2] &&
        is.finite(abs_cor_q) && abs_cor_q <= bal$max_abs_cor_q &&
        diff(range(population_maf)) <= bal$max_population_maf_gap &&
        diff(range(family_maf)) >= fam_ref[1] &&
        diff(range(family_maf)) <= fam_ref[2] &&
        marker_k_r2 >= k_ref[1] && marker_k_r2 <= k_ref[2]) {
      return(marker)
    }
  }
  stop("Could not construct a balanced K-neutral causal marker.")
}

## Phenotype simulation -------------------------------------------------------
draw_kinship_effect <- function(K, variance) {
  eig <- eigen((K + t(K)) / 2, symmetric = TRUE)
  as.numeric(eig$vectors %*% (sqrt(pmax(eig$values, 0)) *
                                stats::rnorm(nrow(K)))) * sqrt(variance)
}

curve_matrix <- function(theta, index) {
  vapply(seq_len(nrow(theta)), function(i) {
    mixFunMap::logistic_mu(theta[i, ], index)
  }, numeric(length(index)))
}

trajectory_variance <- function(curve_difference) {
  mean(apply(as.matrix(curve_difference), 1L, stats::var))
}

counterfactual_effect <- function(theta_base, parameter_effect, index) {
  mean_fun <- mixFunMap::mean_logistic()
  theta1 <- mean_fun$clip(sweep(theta_base, 2L, parameter_effect, "+"), index)
  rowMeans(curve_matrix(theta1, index) - curve_matrix(theta_base, index))
}

sim_index <- function(n_time) seq(1, 14, length.out = n_time)

## Nuisance calibration: scale the Q direction and the tau2 direction so the
## trajectory-heterogeneity/noise ratios equal 0.10 (Q) and 0.50 (K).
resolve_nuisance_calibration <- function(cfg) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(cfg$calibration$seed)
  index <- sim_index(cfg$n_time)
  theta0 <- cfg$theta0
  reference_curve <- mixFunMap::logistic_mu(theta0, index)
  signal_variance <- stats::var(reference_curve)
  residual_variance <- signal_variance / cfg$target_snr
  reference_matrix <- matrix(reference_curve, nrow = length(index),
                             ncol = cfg$n)
  g <- simulate_genotypes(cfg$n, 2L, cfg$m_kinship, cfg)
  mean_fun <- mixFunMap::mean_logistic()
  theta_reference <- matrix(theta0, cfg$n, 3L, byrow = TRUE)
  solve_scale <- function(ratio_fun, target) {
    upper <- 1
    while (ratio_fun(upper) < target && upper < 1e5) upper <- upper * 2
    stats::uniroot(function(s) ratio_fun(s) - target, c(0, upper),
                   tol = 1e-7)$root
  }
  q_scale <- solve_scale(function(scale) {
    theta_q <- mean_fun$clip(
      theta_reference + g$Q[, 1L] %o% (cfg$q_effect_direction * scale), index
    )
    trajectory_variance(curve_matrix(theta_q, index) - reference_matrix) /
      residual_variance
  }, cfg$calibration$target_ratio_q)

  eig <- eigen((g$K + t(g$K)) / 2, symmetric = TRUE)
  K_root <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)),
                                 nrow = length(eig$values))
  z <- array(stats::rnorm(cfg$n * 3L * cfg$calibration$mc_draws),
             dim = c(cfg$n, 3L, cfg$calibration$mc_draws))
  base_u <- array(0, dim = dim(z))
  for (draw in seq_len(cfg$calibration$mc_draws)) {
    for (par in 1:3) {
      base_u[, par, draw] <- sqrt(cfg$tau2_direction[par]) *
        as.numeric(K_root %*% z[, par, draw])
    }
  }
  tau2_scale <- solve_scale(function(scale) {
    mean(vapply(seq_len(cfg$calibration$mc_draws), function(draw) {
      theta_k <- mean_fun$clip(theta_reference + sqrt(scale) *
                                 base_u[, , draw], index)
      trajectory_variance(curve_matrix(theta_k, index) - reference_matrix)
    }, numeric(1))) / residual_variance
  }, cfg$calibration$target_ratio_k)

  cfg$q_effect <- cfg$q_effect_direction * q_scale
  cfg$tau2 <- cfg$tau2_direction * tau2_scale
  cfg
}

## One dataset ----------------------------------------------------------------
simulate_dataset <- function(job, cfg) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(job$seed)
  index <- sim_index(job$n_time)
  g <- simulate_genotypes(job$n, job$m, job$m_kinship, cfg)
  geno <- g$geno
  qtl_index <- sample(seq_len(job$m), 1L)
  geno[qtl_index, ] <- make_balanced_causal_marker(g, cfg)

  has_Q <- job$scenario %in% c("Q", "QK")
  has_K <- job$scenario %in% c("K", "QK")
  mean_fun <- mixFunMap::mean_logistic()
  theta0 <- cfg$theta0
  delta <- theta0 * cfg$qtl_relative_direction * job$effect_scale
  theta_reference <- matrix(theta0, job$n, 3L, byrow = TRUE)
  q_shift <- if (has_Q) g$Q[, 1L] %o% cfg$q_effect else matrix(0, job$n, 3L)
  k_shift <- matrix(0, job$n, 3L)
  if (has_K) {
    for (j in 1:3) k_shift[, j] <- draw_kinship_effect(g$K, cfg$tau2[j])
  }
  theta_base <- mean_fun$clip(theta_reference + q_shift + k_shift, index)
  theta <- mean_fun$clip(theta_base + geno[qtl_index, ] %o% delta, index)

  signal_variance <- stats::var(mixFunMap::logistic_mu(theta0, index))
  base_sigma <- mixFunMap::sad1_base(cfg$phi, length(index))
  gamma <- sqrt(signal_variance /
                  (job$snr * mean(diag(base_sigma))))
  L <- chol(base_sigma + diag(1e-10, length(index)))
  residual <- gamma * (t(L) %*% matrix(stats::rnorm(length(index) * job$n),
                                       nrow = length(index)))
  pheno <- curve_matrix(theta, index) + residual
  colnames(pheno) <- colnames(geno)
  rownames(pheno) <- paste0("t", seq_along(index))

  chromosome <- rep(seq_len(ceiling(job$m / 50)), each = 50L,
                    length.out = job$m)
  marker_map <- data.frame(
    marker = rownames(geno), chr = chromosome,
    pos = ave(seq_len(job$m), chromosome, FUN = seq_along) * 10000,
    stringsAsFactors = FALSE
  )
  list(
    pheno = pheno, index = index, geno = geno,
    Q = if (has_Q) g$Q else NULL,
    K = if (has_K) g$K else NULL,
    marker_map = marker_map,
    qtl = qtl_index, qtl_id = rownames(geno)[qtl_index],
    true_effect_curve = counterfactual_effect(theta_base, delta, index),
    tau2_fit = if (has_K) cfg$tau2 else
      stats::setNames(rep(cfg$numerical_zero_tau2, 3L), names(cfg$tau2))
  )
}

## mixFunMap benchmark: fixed-tau2 null (DGP variances) + dual P3D tests ------
mixfunmap_internal <- function(name) {
  get(name, envir = asNamespace("mixFunMap"), inherits = FALSE)
}

fit_null_fixed_tau <- function(pheno, index, Q, K, tau2,
                               max_outer = 12L, max_resid_optim = 40L,
                               tol_theta = 0.05, damp_alpha = 0.2) {
  mean_fun <- mixFunMap::mean_logistic()
  K <- mixFunMap::stabilize_kinship(K)
  make_Xfix <- mixfunmap_internal("make_Xfix")
  build_linearized_mats <- mixfunmap_internal("build_linearized_mats")
  solve_lmm <- mixfunmap_internal("solve_lmm")
  update_theta <- mixfunmap_internal("update_theta")
  reml_nll_cpp <- mixfunmap_internal("reml_nll_henderson_cpp")
  henderson_context <- mixfunmap_internal("prepare_henderson_kinship")(K)
  use_henderson <- isTRUE(henderson_context$ok)
  solver <- if (use_henderson) "henderson" else "cpp_dense"

  reml_nll <- function(par, ytil, Xstar, Zstar) {
    phi <- pmin(pmax(stats::plogis(par[1]), 1e-4), 0.98)
    gamma <- exp(par[2])
    Sigma <- gamma^2 * mixFunMap::sad1_base(phi, nrow(pheno))
    val <- tryCatch({
      if (use_henderson) {
        reml_nll_cpp(y = as.numeric(ytil), X = as.matrix(Xstar), Z = Zstar,
                     Kinv = henderson_context$Kinv,
                     logdetK = henderson_context$logdetK,
                     tau2 = as.numeric(tau2), Sigma = as.matrix(Sigma),
                     residual_jitter = 1e-8, fixed_jitter = 1e-10)
      } else {
        build_V <- mixfunmap_internal("build_V")
        V <- build_V(Zstar, K, tau2, Sigma) + diag(1e-8, length(ytil))
        cholV <- chol(V)
        VinvX <- backsolve(cholV, forwardsolve(t(cholV), Xstar))
        Vinvy <- backsolve(cholV, forwardsolve(t(cholV), ytil))
        XtVinvX <- crossprod(Xstar, VinvX) + diag(1e-10, ncol(Xstar))
        cholXt <- chol(XtVinvX)
        beta_hat <- backsolve(cholXt, forwardsolve(t(cholXt),
                                                   crossprod(Xstar, Vinvy)))
        Py <- Vinvy - VinvX %*% beta_hat
        0.5 * (2 * sum(log(diag(cholV))) + 2 * sum(log(diag(cholXt))) +
                 sum(ytil * Py))
      }
    }, error = function(e) 1e30)
    if (!is.finite(val)) 1e30 else val
  }

  Xfix <- make_Xfix(snp_vec = NULL, Q = Q, n = ncol(pheno),
                    include_intercept = TRUE)
  theta <- matrix(mean_fun$init(pheno, index), ncol(pheno), mean_fun$npar,
                  byrow = TRUE, dimnames = list(NULL, mean_fun$par_names))
  theta <- mean_fun$clip(theta, index)
  residual0 <- pheno - matrix(mean_fun$mu(mean_fun$init(pheno, index), index),
                              nrow = nrow(pheno))
  phi <- 0.5
  gamma <- max(stats::sd(as.numeric(residual0), na.rm = TRUE), 1e-4)
  sol <- NULL
  converged <- FALSE
  for (it in seq_len(max_outer)) {
    old_theta <- theta
    lin <- build_linearized_mats(pheno, index, theta, Xfix,
                                 mean_fun = mean_fun)
    opt <- tryCatch(
      stats::optim(c(stats::qlogis(phi), log(gamma)), reml_nll,
                   ytil = lin$ytil, Xstar = lin$Xstar, Zstar = lin$Zstar,
                   method = "L-BFGS-B",
                   lower = c(stats::qlogis(1e-4), log(1e-4)),
                   upper = c(stats::qlogis(0.98), log(10)),
                   control = list(maxit = max_resid_optim)),
      error = function(e) NULL
    )
    if (is.null(opt) || !is.finite(opt$value) || opt$value >= 1e29) break
    phi <- pmin(pmax(stats::plogis(opt$par[1]), 1e-4), 0.98)
    gamma <- exp(opt$par[2])
    sol <- solve_lmm(lin$ytil, lin$Xstar, lin$Zstar, K, phi, gamma, tau2,
                     nrow(pheno), solver = solver, resid_cov = NULL)
    if (!isTRUE(sol$ok)) break
    theta_new <- update_theta(Xfix, sol$b_hat, sol$u_hat, ncol(pheno),
                              ncol(Xfix), mean_fun$npar)
    theta <- mean_fun$clip((1 - damp_alpha) * theta + damp_alpha * theta_new,
                           index)
    colnames(theta) <- mean_fun$par_names
    if (max(abs(theta - old_theta), na.rm = TRUE) < tol_theta) {
      converged <- TRUE
      break
    }
  }
  usable <- !is.null(sol) && isTRUE(sol$ok)
  fit <- list(
    theta = theta, phi = phi, gamma = gamma, tau2 = tau2,
    b_hat = if (usable) sol$b_hat else rep(NA_real_, ncol(Xfix) * 3L),
    u_hat = if (usable) sol$u_hat else rep(NA_real_, ncol(pheno) * 3L),
    Vbeta = if (usable) sol$Vbeta else NULL,
    XtVinvX = if (usable) sol$XtVinvX else NULL,
    ml_nll = if (usable) sol$ml_nll else NA_real_,
    mean = mean_fun, solver = solver, converged = converged, usable = usable,
    Xfix = Xfix, Q = Q, K = K, t = index, pheno = pheno,
    include_intercept = TRUE, damp_alpha = damp_alpha, fix_tau3 = FALSE,
    fast = FALSE, var_method = "fixed_tau2_resid_reml", fix_var = TRUE,
    var_init = list(phi = phi, gamma = gamma, tau2 = tau2),
    max_reml_optim = max_resid_optim, reml_every = 1L,
    reml_stop_after = Inf, resid_cov = NULL
  )
  class(fit) <- c("mixfunmap_fit", "mixfunmap_null", "list")
  fit
}

scan_mixfunmap_sim <- function(sim, cfg) {
  Q <- sim$Q
  K <- if (is.null(sim$K)) {
    diag(ncol(sim$pheno),
         dimnames = list(colnames(sim$pheno), colnames(sim$pheno)))
  } else {
    sim$K
  }
  fit <- fit_null_fixed_tau(sim$pheno, sim$index, Q, K, sim$tau2_fit)
  if (!isTRUE(fit$usable) || !isTRUE(fit$converged)) {
    fit <- fit_null_fixed_tau(sim$pheno, sim$index, Q, K, sim$tau2_fit,
                              max_outer = 20L, max_resid_optim = 100L)
  }
  if (!isTRUE(fit$usable)) stop("Fixed-tau null fit failed.")
  # test = "p3d_lrt" stores both the Wald and the refreshed P3D-LRT p values;
  # the Wald p value is the manuscript primary.
  scan <- mixFunMap::qkfunmap_scan(
    fit_null = fit, geno = sim$geno, cores = 1L, max_outer_snp = 3L,
    fix_tau3 = FALSE, tol_theta = 0.05, test = "p3d_lrt",
    solver = fit$solver
  )
  beta <- as.data.frame(scan$beta)
  retry <- !(is.finite(beta$wald_pval) & beta$wald_pval > 0 &
               beta$wald_pval <= 1 & beta$converged > 0 &
               beta$outer_converged > 0)
  if (any(retry)) {
    retry_scan <- mixFunMap::qkfunmap_scan(
      fit_null = fit, geno = sim$geno[retry, , drop = FALSE], cores = 1L,
      max_outer_snp = 6L, fix_tau3 = FALSE, tol_theta = 0.05,
      test = "p3d_lrt", solver = fit$solver
    )
    retry_beta <- as.data.frame(retry_scan$beta)
    beta[match(rownames(retry_beta), rownames(beta)), names(retry_beta)] <-
      retry_beta
  }
  qtl_row <- match(sim$qtl_id, rownames(beta))
  qtl_beta <- as.numeric(beta[qtl_row, paste0("beta", c("A", "K", "T0"))])
  out <- data.frame(
    marker = rownames(beta), pval = beta$wald_pval,
    wald_pval = beta$wald_pval, lrt_pval = beta$lrt_pval,
    converged = isTRUE(fit$usable) & beta$converged > 0 &
      beta$outer_converged > 0,
    stringsAsFactors = FALSE
  )
  attr(out, "estimated_effect_curve") <-
    counterfactual_effect(fit$theta, qtl_beta, sim$index)
  out
}

## One replicate: simulate, scan with three methods, evaluate -----------------
run_replicate <- function(job, cfg) {
  sim <- simulate_dataset(job, cfg)
  alpha <- cfg$alpha / job$m
  evaluate <- function(result) {
    p <- as.numeric(result$pval)
    valid <- is.finite(p) & p > 0 & p <= 1 & result$converged %in% TRUE
    scan_usable <- mean(valid) >= cfg$min_valid_fraction
    qtl_row <- match(sim$qtl_id, result$marker)
    qtl_valid <- !is.na(qtl_row) && valid[qtl_row]
    min_p <- if (scan_usable && any(valid)) min(p[valid]) else NA_real_
    rank_value <- if (job$experiment == "alternative" && qtl_valid) {
      rank(replace(p, !valid, 1), ties.method = "min")[qtl_row]
    } else if (job$experiment == "alternative") {
      job$m + 1L
    } else {
      NA_real_
    }
    effect_hat <- attr(result, "estimated_effect_curve")
    if (is.null(effect_hat)) {
      retained <- attr(result, "retained_markers")
      effect_hat <- retained[[sim$qtl_id]]
      if (is.list(effect_hat)) effect_hat <- effect_hat$effect_curve
    }
    data.frame(
      n_valid = sum(valid), scan_usable = scan_usable,
      min_scan_p = min_p,
      genomewide_false_positive = if (job$experiment == "null" &&
        scan_usable) min_p < alpha else NA,
      qtl_pval = if (qtl_valid) p[qtl_row] else NA_real_,
      genomewide_detected = job$experiment == "alternative" && qtl_valid &&
        p[qtl_row] < alpha,
      qtl_rank = rank_value,
      effect_rmse = if (job$experiment == "alternative" && qtl_valid &&
        length(effect_hat) == length(sim$index) && all(is.finite(effect_hat))) {
        sqrt(mean((as.numeric(effect_hat) - sim$true_effect_curve)^2))
      } else {
        NA_real_
      },
      converged = scan_usable,
      stringsAsFactors = FALSE
    )
  }

  results <- list(
    mixfunmap = function() scan_mixfunmap_sim(sim, cfg),
    ordinary_funmap = function() mixFunMap::scan_ordinary_funmap(
      pheno = sim$pheno, geno = sim$geno, index = sim$index,
      marker_map = sim$marker_map, min_group_n = cfg$genotype$min_group_n,
      retain_marker = sim$qtl_id, maxit = 300L, reltol = 1e-8
    ),
    minp_gmmat = function() mixFunMap::scan_minp_gmmat(
      pheno = sim$pheno, geno = sim$geno, index = sim$index,
      Q = sim$Q, K = sim$K, marker_map = sim$marker_map,
      correction = "bonferroni", maxiter = 200L, tol = 1e-5,
      primary_optimizer = "AI", fallback_optimizer = "Brent",
      retain_marker = sim$qtl_id
    )
  )
  do.call(rbind, lapply(names(results), function(method) {
    metric <- tryCatch({
      evaluate(as.data.frame(results[[method]](), check.names = FALSE))
    }, error = function(e) {
      data.frame(n_valid = 0L, scan_usable = FALSE, min_scan_p = NA_real_,
                 genomewide_false_positive = NA, qtl_pval = NA_real_,
                 genomewide_detected = FALSE, qtl_rank = NA_real_,
                 effect_rmse = NA_real_, converged = FALSE,
                 stringsAsFactors = FALSE)
    })
    cbind(job[, c("design_id", "job_id", "replicate", "scenario",
                  "experiment", "varied_factor", "factor_level", "n",
                  "n_time", "snr", "effect_scale", "m")],
          method = method, metric, stringsAsFactors = FALSE)
  }))
}

## Design: 28 null cells + 44 alternative cells --------------------------------
build_design <- function(cfg) {
  base <- list(n = cfg$n, n_time = cfg$n_time, snr = cfg$target_snr)
  rows <- list()
  add <- function(scenario, experiment, varied_factor, factor_level,
                  n, n_time, snr, effect_scale, n_rep) {
    rows[[length(rows) + 1L]] <<- data.frame(
      scenario = scenario, experiment = experiment,
      varied_factor = varied_factor, factor_level = factor_level,
      n = n, n_time = n_time, snr = snr, effect_scale = effect_scale,
      n_rep = n_rep, m = cfg$m, stringsAsFactors = FALSE
    )
  }
  for (scenario in cfg$scenarios) {
    add(scenario, "null", "baseline", 0, base$n, base$n_time, base$snr, 0,
        cfg$n_rep$null)
    for (v in setdiff(cfg$grids$n, base$n)) {
      add(scenario, "null", "N", v, v, base$n_time, base$snr, 0, cfg$n_rep$null)
    }
    for (v in setdiff(cfg$grids$n_time, base$n_time)) {
      add(scenario, "null", "T", v, base$n, v, base$snr, 0, cfg$n_rep$null)
    }
    for (v in setdiff(cfg$grids$snr, base$snr)) {
      add(scenario, "null", "SNR", v, base$n, base$n_time, v, 0, cfg$n_rep$null)
    }
    for (v in cfg$grids$effect) {
      add(scenario, "alternative", "effect", v, base$n, base$n_time, base$snr,
          v, cfg$n_rep$alternative)
    }
    for (v in setdiff(cfg$grids$n, base$n)) {
      add(scenario, "alternative", "N", v, v, base$n_time, base$snr,
          cfg$reference_effect, cfg$n_rep$alternative)
    }
    for (v in setdiff(cfg$grids$n_time, base$n_time)) {
      add(scenario, "alternative", "T", v, base$n, v, base$snr,
          cfg$reference_effect, cfg$n_rep$alternative)
    }
    for (v in setdiff(cfg$grids$snr, base$snr)) {
      add(scenario, "alternative", "SNR", v, base$n, base$n_time, v,
          cfg$reference_effect, cfg$n_rep$alternative)
    }
  }
  design <- do.call(rbind, rows)
  design$design_id <- sprintf(
    "%s_%s_%s_n%03d_t%02d_snr%g_eff%03d",
    c(none = "none", Q = "q", K = "k", QK = "qk")[design$scenario],
    design$experiment, tolower(design$varied_factor), design$n,
    design$n_time, design$snr, round(100 * design$effect_scale)
  )
  jobs <- design[rep(seq_len(nrow(design)), design$n_rep), , drop = FALSE]
  jobs$replicate <- unlist(lapply(design$n_rep, seq_len))
  jobs$job_number <- seq_len(nrow(jobs))
  jobs$seed <- as.integer(cfg$seed + jobs$job_number - 1L)
  jobs$job_id <- sprintf("%s_rep%04d", jobs$design_id, jobs$replicate)
  rownames(jobs) <- NULL
  jobs
}

## Summaries -------------------------------------------------------------------
wilson_interval <- function(success, total) {
  if (!is.finite(total) || total <= 0) return(c(lcl = NA_real_, ucl = NA_real_))
  z <- stats::qnorm(0.975)
  p <- success / total
  denominator <- 1 + z^2 / total
  center <- (p + z^2 / (2 * total)) / denominator
  half <- z * sqrt(p * (1 - p) / total + z^2 / (4 * total^2)) / denominator
  c(lcl = pmax(0, center - half), ucl = pmin(1, center + half))
}

bootstrap_interval <- function(x, statistic = "mean", B = 2000L, seed = 1L) {
  x <- as.numeric(x[is.finite(x)])
  if (!length(x)) return(c(lcl = NA_real_, ucl = NA_real_))
  if (length(x) == 1L) return(c(lcl = x, ucl = x))
  fun <- if (statistic == "mean") base::mean else stats::median
  set.seed(seed)
  boot <- replicate(B, fun(sample(x, length(x), replace = TRUE)))
  stats::quantile(boot, c(0.025, 0.975), names = FALSE, na.rm = TRUE) |>
    stats::setNames(c("lcl", "ucl"))
}

summarise_metrics <- function(metrics, cfg) {
  # Empirical 5%-FWER cutoffs: 5th percentile of the null min-scan p per
  # method x scenario x design cell.
  group_fields <- c("scenario", "n", "n_time", "snr", "m", "method")
  null <- metrics[metrics$experiment == "null" &
                    is.finite(metrics$min_scan_p), , drop = FALSE]
  cutoffs <- do.call(rbind, lapply(
    split(null, interaction(null[group_fields], drop = TRUE)),
    function(z) {
      data.frame(z[1L, group_fields, drop = FALSE],
                 empirical_fwer_cutoff = unname(stats::quantile(
                   z$min_scan_p, 0.05, type = 1L, names = FALSE)),
                 stringsAsFactors = FALSE, check.names = FALSE)
    }
  ))
  metrics$empirical_fwer_cutoff <- NULL
  metrics <- merge(metrics, cutoffs, by = group_fields, all.x = TRUE,
                   sort = FALSE)
  metrics$adjusted_genomewide_detected <-
    metrics$experiment == "alternative" &
    is.finite(metrics$qtl_pval) & is.finite(metrics$empirical_fwer_cutoff) &
    metrics$qtl_pval < metrics$empirical_fwer_cutoff

  cell_fields <- c("design_id", "scenario", "experiment", "varied_factor",
                   "factor_level", "n", "n_time", "snr", "effect_scale", "m",
                   "method")
  summary <- do.call(rbind, lapply(
    split(metrics, interaction(metrics[cell_fields], drop = TRUE)),
    function(z) {
      fwer_values <- z$genomewide_false_positive[
        !is.na(z$genomewide_false_positive)]
      fwer_ci <- wilson_interval(sum(fwer_values %in% TRUE),
                                 length(fwer_values))
      power_ci <- wilson_interval(sum(z$genomewide_detected %in% TRUE),
                                  nrow(z))
      adjusted_ci <- wilson_interval(
        sum(z$adjusted_genomewide_detected %in% TRUE), nrow(z))
      rmse_ci <- bootstrap_interval(z$effect_rmse, "mean",
                                    cfg$bootstrap_reps, cfg$seed + 1L)
      rank_ci <- bootstrap_interval(z$qtl_rank, "median",
                                    cfg$bootstrap_reps, cfg$seed + 2L)
      cbind(
        z[1L, cell_fields, drop = FALSE],
        data.frame(
          n_total = nrow(z),
          fwer = if (length(fwer_values)) {
            sum(fwer_values %in% TRUE) / length(fwer_values)
          } else {
            NA_real_
          },
          fwer_lcl = fwer_ci["lcl"], fwer_ucl = fwer_ci["ucl"],
          power = power_ci["lcl"] * 0 + sum(z$genomewide_detected %in% TRUE) /
            nrow(z),
          power_lcl = power_ci["lcl"], power_ucl = power_ci["ucl"],
          adjusted_power = if (z$experiment[1L] == "alternative") {
            sum(z$adjusted_genomewide_detected %in% TRUE) / nrow(z)
          } else {
            NA_real_
          },
          adjusted_power_lcl = adjusted_ci["lcl"],
          adjusted_power_ucl = adjusted_ci["ucl"],
          rmse = mean(z$effect_rmse[is.finite(z$effect_rmse)]),
          rmse_lcl = rmse_ci["lcl"], rmse_ucl = rmse_ci["ucl"],
          rank_median = stats::median(z$qtl_rank[is.finite(z$qtl_rank)]),
          rank_lcl = rank_ci["lcl"], rank_ucl = rank_ci["ucl"],
          convergence_rate = mean(z$converged %in% TRUE),
          stringsAsFactors = FALSE, check.names = FALSE
        )
      )
    }
  ))
  rownames(summary) <- NULL
  names(summary)[names(summary) == "n_time"] <- "T"
  list(summary = summary, metrics = metrics, cutoffs = cutoffs)
}

run_simulation <- function(output_dir, cores = 1L,
                           n_rep_null = NULL, n_rep_alt = NULL) {
  if (!requireNamespace("mixFunMap", quietly = TRUE)) {
    stop("The mixFunMap package is required for the simulation mode.",
         call. = FALSE)
  }
  if (!requireNamespace("GMMAT", quietly = TRUE)) {
    stop("GMMAT is required for the minP scans.", call. = FALSE)
  }
  cfg <- sim_cfg
  if (!is.null(n_rep_null)) cfg$n_rep$null <- as.integer(n_rep_null)
  if (!is.null(n_rep_alt)) cfg$n_rep$alternative <- as.integer(n_rep_alt)
  message("Calibrating nuisance Q / K variance ratios (0.10 / 0.50) ...")
  cfg <- resolve_nuisance_calibration(cfg)
  jobs <- build_design(cfg)
  message(sprintf(
    "Running %d datasets (%d null + %d alternative cells) on %d core(s) ...",
    nrow(jobs), sum(jobs$experiment == "null"),
    sum(jobs$experiment == "alternative"), cores
  ))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  run_one <- function(i) run_replicate(jobs[i, , drop = FALSE], cfg)
  metrics <- if (cores > 1L) {
    cluster <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterEvalQ(cluster, library(mixFunMap))
    parallel::clusterExport(cluster, c(
      "sim_cfg", "make_nested_structure", "center_family_probabilities",
      "draw_nested_marker", "simulate_nested_markers",
      "stabilize_projected_kinship", "project_kinship_off_q",
      "simulate_genotypes", "kinship_pcs", "make_balanced_causal_marker",
      "draw_kinship_effect", "curve_matrix", "trajectory_variance",
      "counterfactual_effect", "sim_index", "simulate_dataset",
      "mixfunmap_internal", "fit_null_fixed_tau", "scan_mixfunmap_sim",
      "run_replicate"
    ), envir = environment())
    do.call(rbind, parallel::parLapply(cluster, seq_len(nrow(jobs)),
                                       run_one))
  } else {
    do.call(rbind, lapply(seq_len(nrow(jobs)), function(i) {
      if (i %% 100L == 0L) message("  dataset ", i, "/", nrow(jobs))
      run_one(i)
    }))
  }

  summarized <- summarise_metrics(metrics, cfg)
  summary_dir_out <- file.path(output_dir, "summary")
  dir.create(summary_dir_out, recursive = TRUE, showWarnings = FALSE)
  saveRDS(summarized$summary,
          file.path(summary_dir_out, "simulation_summary.rds"), version = 3)
  utils::write.csv(summarized$summary,
                   file.path(summary_dir_out, "simulation_summary.csv"),
                   row.names = FALSE)
  utils::write.csv(summarized$metrics,
                   file.path(summary_dir_out,
                             "simulation_replicate_metrics.csv"),
                   row.names = FALSE)
  message("Simulation summary written to: ", summary_dir_out)
  invisible(summary_dir_out)
}

if (identical(mode, "full")) {
  cores <- if (length(.args) >= 2L) as.integer(.args[[2L]]) else 1L
  n_rep_null <- if (length(.args) >= 3L) as.integer(.args[[3L]]) else NULL
  n_rep_alt <- if (length(.args) >= 4L) as.integer(.args[[4L]]) else NULL
  summary_dir <- run_simulation(
    file.path(script_dir, "results", "simulation", "rerun"),
    cores = cores, n_rep_null = n_rep_null, n_rep_alt = n_rep_alt
  )
}

if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("patchwork", quietly = TRUE)) {
  stop("ggplot2 and patchwork are required for the simulation figures.",
       call. = FALSE)
}

summary_path <- file.path(summary_dir, "simulation_summary.rds")
if (!file.exists(summary_path)) {
  stop("simulation_summary.rds was not found at: ", summary_path,
       call. = FALSE)
}

# ---------------------------------------------------------------------------
# Shared style
# ---------------------------------------------------------------------------

simulation_plot_style <- function() {
  list(
    colors = c(mixfunmap = "#0072B2", ordinary_funmap = "#D55E00", minp_gmmat = "#009E73"),
    shapes = c(mixfunmap = 16, ordinary_funmap = 17, minp_gmmat = 15),
    linetypes = c(mixfunmap = "solid", ordinary_funmap = "dashed", minp_gmmat = "dotdash"),
    labels = c(
      mixfunmap = "mixFunMap (Wald)",
      ordinary_funmap = "ordinary FunMap (LRT)",
      minp_gmmat = "GMMAT-minP"
    )
  )
}

prepare_plot_factors <- function(data, style = simulation_plot_style()) {
  data$scenario <- factor(data$scenario, levels = c("none", "Q", "K", "QK"),
                          labels = c("Neither Q nor K", "Q only", "K only", "Q + K"))
  data$method <- factor(data$method, levels = names(style$labels), labels = unname(style$labels))
  data
}

simulation_scales <- function(style = simulation_plot_style()) {
  list(
    ggplot2::scale_color_manual(values = stats::setNames(style$colors, unname(style$labels))),
    ggplot2::scale_shape_manual(values = stats::setNames(style$shapes, unname(style$labels))),
    ggplot2::scale_linetype_manual(values = stats::setNames(style$linetypes, unname(style$labels)))
  )
}

simulation_theme <- function() {
  ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, color = "#E5E7EB"),
      strip.background = ggplot2::element_rect(fill = "#F3F4F6", color = "#9CA3AF"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom", legend.title = ggplot2::element_blank()
    )
}

add_common_aesthetics <- function(plot, style = simulation_plot_style()) {
  scales <- simulation_scales(style)
  plot + scales[[1]] + scales[[2]] + scales[[3]] + simulation_theme()
}

save_simulation_plot <- function(plot, stem, figure_dir, width, height) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
                  width = width, height = height, units = "in", device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot,
                  width = width, height = height, units = "in", dpi = 600,
                  compression = "lzw")
}

# ---------------------------------------------------------------------------
# Sensitivity grids (Figures S1, S2, S4, S5)
# ---------------------------------------------------------------------------

baseline_for_sensitivity <- function(summary, experiment, reference_effect = 0.12) {
  if (experiment == "null") {
    base <- summary[summary$experiment == "null" & summary$varied_factor == "baseline", , drop = FALSE]
  } else {
    base <- summary[summary$experiment == "alternative" & summary$varied_factor == "effect" &
                      abs(summary$effect_scale - reference_effect) < 1e-12, , drop = FALSE]
  }
  copies <- lapply(c("N", "T", "SNR"), function(factor_name) {
    z <- base
    z$varied_factor <- factor_name
    z$factor_level <- switch(factor_name, N = z$n, T = z$T, SNR = z$snr)
    z
  })
  do.call(rbind, copies)
}

plot_sensitivity_grid <- function(summary, experiment, metric, lcl, ucl, y_label,
                                  reference_line = NULL, reference_effect = 0.12,
                                  style = simulation_plot_style(),
                                  add_aes = add_common_aesthetics,
                                  mixfunmap_on_top = FALSE) {
  data <- summary[summary$experiment == experiment & summary$varied_factor %in% c("N", "T", "SNR"), , drop = FALSE]
  data <- rbind(data, baseline_for_sensitivity(summary, experiment, reference_effect))
  data <- prepare_plot_factors(data, style)
  if (mixfunmap_on_top) data <- mixfunmap_drawn_last(data)
  data$varied_factor <- factor(data$varied_factor, levels = c("N", "T", "SNR"),
                               labels = c("Sample size", "Time points", "Signal-to-noise ratio"))
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data$factor_level, y = .data[[metric]], color = .data$method,
                 shape = .data$method, linetype = .data$method, group = .data$method)
  ) +
    ggplot2::geom_line(linewidth = 0.55, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.8, na.rm = TRUE) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data[[lcl]], ymax = .data[[ucl]]),
                           width = 0, linewidth = 0.35, na.rm = TRUE) +
    ggplot2::facet_grid(rows = ggplot2::vars(scenario), cols = ggplot2::vars(varied_factor),
                        scales = "free_x") +
    ggplot2::labs(x = NULL, y = y_label)
  if (!is.null(reference_line)) p <- p + ggplot2::geom_hline(yintercept = reference_line, linetype = 2)
  add_aes(p, style)
}

# ---------------------------------------------------------------------------
# Figure 1: three standalone subplot PDFs
#   a) FWER, b) power at matched 5% FWER, c) median causal-marker rank
# Subplot style: no axis titles, no facet strip labels, larger axis text;
# Functional Mapping and GMMAT-minP are drawn first, mixFunMap on top.
# ---------------------------------------------------------------------------

figure1_plot_style <- function() {
  style <- simulation_plot_style()
  style$labels["ordinary_funmap"] <- "Functional Mapping"
  style
}

# Reorder rows so mixFunMap (factor level 1) is drawn last, i.e. on top of
# Functional Mapping (2) and GMMAT-minP (3); the legend order still follows
# the factor levels.
mixfunmap_drawn_last <- function(data) {
  data[order(match(as.integer(data$method), c(2L, 3L, 1L))), , drop = FALSE]
}

figure1_subplot_theme <- function() {
  simulation_theme() +
    ggplot2::theme(
      strip.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 12, colour = "black")
    )
}

add_figure1_aesthetics <- function(plot, style = figure1_plot_style()) {
  plot +
    ggplot2::scale_color_manual(values = stats::setNames(style$colors, unname(style$labels))) +
    ggplot2::scale_shape_manual(values = stats::setNames(style$shapes, unname(style$labels))) +
    ggplot2::scale_linetype_manual(values = stats::setNames(style$linetypes, unname(style$labels))) +
    figure1_subplot_theme()
}

make_figure1_subplots <- function(
    summary,
    figure_dir,
    stem = "Figure1_v15_p3dwald_relative_qtl_d11m1_neutral_n100_fam5_k050") {
  style <- figure1_plot_style()
  fwer <- mixfunmap_drawn_last(prepare_plot_factors(
    summary[summary$experiment == "null" & summary$varied_factor == "baseline", ],
    style
  ))
  effect <- mixfunmap_drawn_last(prepare_plot_factors(
    summary[summary$experiment == "alternative" & summary$varied_factor == "effect", ],
    style
  ))
  if (!nrow(fwer) || !nrow(effect)) {
    stop("Figure 1 requires baseline null cells and effect-series alternative cells.",
         call. = FALSE)
  }

  # Three evenly spaced x ticks (endpoints + middle grid value) so the
  # enlarged tick labels never collide within or across facets.
  effect_grid <- sort(unique(effect$effect_scale))
  x_breaks <- effect_grid[unique(round(seq(1, length(effect_grid), length.out = 3)))]

  p_fwer <- ggplot2::ggplot(
    fwer, ggplot2::aes(x = .data$method, y = .data$fwer,
                       color = .data$method, shape = .data$method,
                       linetype = .data$method)
  ) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = 2, color = "#4B5563") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$fwer_lcl, ymax = .data$fwer_ucl),
      width = 0.12, linewidth = 0.4, show.legend = FALSE
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_grid(. ~ scenario) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_cartesian(ylim = c(0, max(0.10, fwer$fwer_ucl, na.rm = TRUE)))
  p_fwer <- add_figure1_aesthetics(p_fwer) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())

  p_adjusted_power <- ggplot2::ggplot(
    effect, ggplot2::aes(
      x = .data$effect_scale, y = .data$adjusted_power,
      color = .data$method, shape = .data$method,
      linetype = .data$method, group = .data$method
    )
  ) +
    ggplot2::geom_line(linewidth = 0.55, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.8, na.rm = TRUE) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data$adjusted_power_lcl,
        ymax = .data$adjusted_power_ucl
      ),
      width = 0, linewidth = 0.35, na.rm = TRUE
    ) +
    ggplot2::facet_grid(. ~ scenario) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      expand = ggplot2::expansion(mult = 0.14)
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = NULL, y = NULL)
  p_adjusted_power <- add_figure1_aesthetics(p_adjusted_power)

  p_rank <- ggplot2::ggplot(
    effect, ggplot2::aes(x = .data$effect_scale, y = .data$rank_median,
                         color = .data$method, shape = .data$method,
                         linetype = .data$method, group = .data$method)
  ) +
    ggplot2::geom_hline(yintercept = 10, linetype = 3, color = "#6B7280") +
    ggplot2::geom_line(linewidth = 0.55, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.8, na.rm = TRUE) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(1, .data$rank_lcl), ymax = .data$rank_ucl),
      width = 0, linewidth = 0.35, na.rm = TRUE
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      expand = ggplot2::expansion(mult = 0.14)
    ) +
    ggplot2::facet_grid(. ~ scenario) +
    ggplot2::labs(x = NULL, y = NULL)
  p_rank <- add_figure1_aesthetics(p_rank)

  save_simulation_plot(p_fwer, paste0(stem, "_fwer"), figure_dir, 7.2, 2.8)
  save_simulation_plot(p_adjusted_power, paste0(stem, "_adjusted_power"),
                       figure_dir, 7.2, 3.1)
  save_simulation_plot(p_rank, paste0(stem, "_rank"), figure_dir, 7.2, 3.1)
  utils::write.csv(
    rbind(
      transform(fwer, panel_metric = "fwer"),
      transform(effect, panel_metric = "adjusted_power_and_rank")
    ),
    file.path(figure_dir, paste0(stem, "_plot_data.csv")), row.names = FALSE
  )
  invisible(list(fwer = p_fwer, adjusted_power = p_adjusted_power, rank = p_rank))
}

# ---------------------------------------------------------------------------
# Driver: build every figure from the frozen summary
# ---------------------------------------------------------------------------

make_neutral_n100_simulation_figures <- function(
    summary_dir,
    figure_dir = file.path(dirname(summary_dir), "figures"),
    reference_effect = 0.12) {
  summary <- readRDS(file.path(summary_dir, "simulation_summary.rds"))
  effect <- prepare_plot_factors(
    summary[summary$experiment == "alternative" &
              summary$varied_factor == "effect", ]
  )
  if (!nrow(effect)) {
    stop("Effect-series alternatives are required for the simulation figures.",
         call. = FALSE)
  }

  main <- make_figure1_subplots(
    summary = summary,
    figure_dir = figure_dir,
    stem = "Figure1_v15_p3dwald_relative_qtl_d11m1_neutral_n100_fam5_k050"
  )
  fig1_style <- figure1_plot_style()
  s1 <- plot_sensitivity_grid(
    summary, "null", "fwer", "fwer_lcl", "fwer_ucl", NULL,
    reference_line = 0.05, reference_effect = reference_effect,
    style = fig1_style, add_aes = add_figure1_aesthetics,
    mixfunmap_on_top = TRUE
  )
  s2 <- plot_sensitivity_grid(
    summary, "alternative", "adjusted_power",
    "adjusted_power_lcl", "adjusted_power_ucl", NULL,
    reference_effect = reference_effect,
    style = fig1_style, add_aes = add_figure1_aesthetics,
    mixfunmap_on_top = TRUE
  ) + ggplot2::coord_cartesian(ylim = c(0, 1))

  s3 <- ggplot2::ggplot(
    effect, ggplot2::aes(x = .data$effect_scale, y = .data$rmse,
                         color = .data$method, shape = .data$method,
                         linetype = .data$method, group = .data$method)
  ) +
    ggplot2::geom_line(linewidth = 0.55, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.8, na.rm = TRUE) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$rmse_lcl, ymax = .data$rmse_ucl),
      width = 0, linewidth = 0.35, na.rm = TRUE
    ) +
    ggplot2::facet_grid(. ~ scenario) +
    ggplot2::labs(x = "Relative QTL effect strength (rho)", y = "Dynamic-effect RMSE")
  s3 <- add_common_aesthetics(s3)

  s4 <- plot_sensitivity_grid(
    summary, "alternative", "rmse", "rmse_lcl", "rmse_ucl",
    "Dynamic-effect RMSE", reference_effect = reference_effect
  )
  s5 <- plot_sensitivity_grid(
    summary, "alternative", "rank_median", "rank_lcl", "rank_ucl", NULL,
    reference_line = 10, reference_effect = reference_effect,
    style = fig1_style, add_aes = add_figure1_aesthetics,
    mixfunmap_on_top = TRUE
  ) + ggplot2::scale_y_log10()

  save_simulation_plot(s1, "FigureS1_fwer_sensitivity", figure_dir, 7.2, 8.2)
  save_simulation_plot(s2, "FigureS2_power_sensitivity", figure_dir, 7.2, 8.2)
  save_simulation_plot(s3, "FigureS3_rmse_by_effect", figure_dir, 7.2, 3.1)
  save_simulation_plot(s4, "FigureS4_rmse_sensitivity", figure_dir, 7.2, 8.2)
  save_simulation_plot(s5, "FigureS5_rank_sensitivity", figure_dir, 7.2, 8.2)

  supplementary_data <- rbind(
    transform(summary[summary$experiment == "null", ],
              supplementary_metric = "fwer"),
    transform(summary[summary$experiment == "alternative", ],
              supplementary_metric = "power_rmse_rank")
  )
  utils::write.csv(
    supplementary_data,
    file.path(figure_dir, "Supplementary_figure_plot_data_v15.csv"),
    row.names = FALSE
  )
  invisible(list(main = main, fwer = s1, power = s2,
                 rmse_effect = s3, rmse_sensitivity = s4,
                 rank_sensitivity = s5))
}

figures <- make_neutral_n100_simulation_figures(
  summary_dir = summary_dir,
  figure_dir = figure_dir,
  reference_effect = 0.12
)

cat("Simulation figures written to: ",
    normalizePath(figure_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
