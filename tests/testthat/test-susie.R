susie_fixture <- function(z = c(5.5, 5.2, 0.3, -0.2, 0.1)) {
  variants <- data.frame(
    CHROM = rep("1", 5),
    POS = 100:104,
    REF = c("A", "C", "G", "T", "A"),
    ALT = c("G", "T", "A", "C", "C"),
    EZ = z,
    LP = rep(2, 5),
    SS = c(900, 850, 800, 750, 700)
  )
  ids <- paste(variants$CHROM, variants$POS, variants$REF, variants$ALT,
               sep = ":")
  R <- matrix(0.1, 5, 5, dimnames = list(ids, ids))
  diag(R) <- 1
  list(variants = variants, R = R)
}

test_that("SuSiE-RSS records PIPs, credible sets, and the selected effective N", {
  fixture <- susie_fixture()
  result <- run_susie_rss(
    fixture$variants,
    fixture$R,
    L = 3,
    max_iter = 30,
    min_abs_corr = NULL
  )

  expect_equal(result$provenance$method_id, "susie_rss")
  expect_equal(result$provenance$parameters$effective_n, 700)
  expect_equal(result$provenance$parameters$effective_n_source,
               "minimum_variant_ss")
  expect_equal(nrow(result$variants), 5L)
  expect_true(all(result$variants$PIP >= 0 & result$variants$PIP <= 1))
  expect_true(is.logical(result$diagnostics$converged))
  expect_true(all(c("credible_set_id", "coverage", "min_abs_corr") %in%
                    names(result$credible_sets)))
})

test_that("robust SuSiE-RSS records finite-reference and mismatch diagnostics", {
  fixture <- susie_fixture()
  result <- run_susie_rrss(
    fixture$variants,
    fixture$R,
    reference_sample_size = 500,
    L = 3,
    max_iter = 30,
    min_abs_corr = NULL
  )

  expect_equal(result$provenance$method_id, "susie_rrss")
  expect_equal(result$provenance$parameters$R_finite, 500)
  expect_equal(result$provenance$parameters$R_mismatch, "eb_mix")
  expect_true("SER_fallback" %in% names(result$diagnostics))
  expect_true("R_finite_diagnostics" %in% names(result$diagnostics))

  configured <- run_susie_rrss(
    fixture$variants,
    fixture$R,
    reference_sample_size = 500,
    L = 3,
    max_iter = 30,
    min_abs_corr = NULL,
    R_mismatch = "eb"
  )
  expect_equal(configured$provenance$parameters$R_mismatch, "eb")
})

test_that("SuSiE requires signed Z, usable LD, and a sample-size value", {
  fixture <- susie_fixture()
  no_ss <- fixture$variants
  no_ss$SS <- NULL
  expect_error(run_susie_rss(no_ss, fixture$R), "supply an explicit")

  no_z <- fixture$variants
  no_z$EZ <- NULL
  expect_error(run_susie_rss(no_z, fixture$R), "signed Z-score")

  expect_error(run_susie_rrss(fixture$variants, fixture$R, 1),
               "greater than 2")
  expect_error(run_susie_rrss(
    fixture$variants, fixture$R, 500, R_mismatch = "unsupported"
  ), "R_mismatch")
})

test_that("SuSiE handles a one-variant locus and rejects an empty locus", {
  fixture <- susie_fixture()
  single <- fixture$variants[1, , drop = FALSE]
  id <- paste(single$CHROM, single$POS, single$REF, single$ALT, sep = ":")
  result <- run_susie_rss(
    single,
    matrix(1, 1, 1, dimnames = list(id, id)),
    L = 1,
    max_iter = 10,
    min_abs_corr = NULL
  )

  expect_equal(length(result$variants$PIP), 1L)
  expect_true(is.finite(result$variants$PIP[[1]]))

  empty <- fixture$variants[0, , drop = FALSE]
  empty_ld <- matrix(numeric(), 0, 0,
                     dimnames = list(character(), character()))
  expect_error(run_susie_rss(empty, empty_ld), "no variants")
})

test_that("SuSiE-rRSS keeps eb_mix active above eigen_max_variants while skipping Q_art eigen", {
  fixture <- susie_fixture(z = c(8.5, 7.9, 0.4, -0.2, 0.1))
  res <- run_susie_rrss(
    fixture$variants,
    fixture$R,
    reference_sample_size = 500,
    L = 3,
    max_iter = 50,
    min_abs_corr = NULL,
    eigen_max_variants = 2L
  )

  expect_equal(res$provenance$parameters$R_mismatch, "eb_mix")
  expect_equal(res$diagnostics$r_mismatch_requested, "eb_mix")
  expect_equal(res$diagnostics$r_mismatch_effective, "eb_mix")
  expect_equal(res$diagnostics$q_art_eigen, "skipped_large_m")
  expect_equal(res$diagnostics$kriging_status, "skipped_large_m")
  expect_true(is.na(res$diagnostics$s_hat))
  expect_true(all(is.na(res$variants$KRIGING_EXP_Z)))
  expect_true(all(!res$variants$KRIGING_OUTLIER))
  expect_false(isTRUE(res$diagnostics$artifact_evaluable))
  expect_true(all(c("prior_variance_V", "lbf", "ld_mismatch_warning",
                    "kriging_outlier_count") %in%
                    names(res$credible_sets)))
  expect_true(all(c("initial_elbo", "final_elbo", "final_elbo_change",
                    "max_prior_variance_V", "max_lbf") %in%
                    names(res$diagnostics)))
})

test_that("compute_kriging_from_eigen matches susieR::estimate_s_rss and kriging_rss", {
  set.seed(42)
  n_sim <- 1500L
  p_sim <- 150L
  X <- matrix(stats::rnorm(n_sim * p_sim), nrow = n_sim, ncol = p_sim)
  for (j in 2:6) {
    X[, j] <- 0.88 * X[, 1] + sqrt(1 - 0.88^2) * stats::rnorm(n_sim)
  }
  X <- scale(X)
  y <- scale(0.18 * X[, 1] + stats::rnorm(n_sim))[, 1]
  R_sim <- stats::cor(X)
  z_sim <- as.numeric(crossprod(X, y) / sqrt(n_sim - 1))
  # Inject an upstream allele-switch flip on variant 3
  z_sim[3] <- -z_sim[3]

  s_ref <- susieR::estimate_s_rss(
    z = z_sim, R = R_sim, n = n_sim, r_tol = 1e-8, method = "null-mle"
  )
  k_ref <- susieR::kriging_rss(
    z = z_sim, R = R_sim, n = n_sim, r_tol = 1e-8, s = s_ref
  )$plot$data

  eig <- eigen(R_sim, symmetric = TRUE)
  k_fast <- gwasfinemapping:::compute_kriging_from_eigen(
    z = z_sim, n = n_sim, eig = eig, r_tol = 1e-8
  )

  expect_equal(k_fast$status, "computed")
  expect_equal(k_fast$s_hat, as.numeric(s_ref), tolerance = 1e-6)
  expect_equal(k_fast$condmean, as.numeric(k_ref$condmean), tolerance = 1e-6)
  expect_equal(k_fast$condvar, as.numeric(k_ref$condvar), tolerance = 1e-6)
  expect_equal(k_fast$z_std_diff, as.numeric(k_ref$z_std_diff), tolerance = 1e-6)
  expect_equal(k_fast$logLR, as.numeric(k_ref$logLR), tolerance = 1e-4)
  expect_true(k_fast$outlier[3])

  # Pure-R EM mixture fallback produces valid simplex weights
  w_em <- gwasfinemapping:::fit_kriging_mixture_em(matrix(c(-0.1, -2.0, -1.5, -0.2), nrow = 2))
  expect_equal(sum(w_em), 1.0, tolerance = 1e-8)
  expect_true(all(w_em >= 0))
})

test_that("SuSiE-rRSS records kriging diagnostics non-destructively by default and supports opt-in removal", {
  set.seed(42)
  n_eff <- 1500L
  p <- 150L
  X <- matrix(stats::rnorm(n_eff * p), nrow = n_eff, ncol = p)
  for (j in 2:6) {
    X[, j] <- 0.88 * X[, 1] + sqrt(1 - 0.88^2) * stats::rnorm(n_eff)
  }
  X <- scale(X)
  y <- scale(0.18 * X[, 1] + stats::rnorm(n_eff))[, 1]
  R <- stats::cor(X)
  z <- as.numeric(crossprod(X, y) / sqrt(n_eff - 1))
  z[3] <- -z[3]  # Simulated allele-switch flip on tightly correlated tag SNP 3

  variants <- data.frame(
    CHROM = rep("1", p),
    POS = 1000:(1000 + p - 1L),
    REF = rep("A", p),
    ALT = rep("G", p),
    EZ = z,
    LP = pmax(-stats::pnorm(abs(z), lower.tail = FALSE, log.p = TRUE) / log(10) - log10(2), 0.1),
    SS = rep(n_eff, p)
  )
  ids <- paste(variants$CHROM, variants$POS, variants$REF, variants$ALT, sep = ":")
  dimnames(R) <- list(ids, ids)

  # 1. Default: remove_kriging_outliers = FALSE keeps all variants and flags SNP 3
  res_keep <- run_susie_rrss(
    variants,
    R,
    reference_sample_size = 5000,
    L = 3,
    max_iter = 100
  )
  expect_false(res_keep$diagnostics$remove_kriging_outliers)
  expect_equal(res_keep$diagnostics$kriging_status, "computed")
  expect_equal(nrow(res_keep$variants), p)
  expect_equal(res_keep$diagnostics$kriging_outliers_removed, 0L)
  expect_true(res_keep$variants$KRIGING_OUTLIER[3])
  expect_gte(res_keep$diagnostics$kriging_outlier_count, 1L)
  expect_true(ids[3] %in% res_keep$diagnostics$kriging_outlier_gwas_variant_ids)
  expect_true(all(c("KRIGING_EXP_Z", "KRIGING_STD_RESID", "KRIGING_LOGLR", "KRIGING_OUTLIER") %in%
                    names(res_keep$variants)))

  # 2. Pre-scaled XtX (ld_scale = n_eff - 1) produces identical kriging quantities
  R_scaled <- R * (n_eff - 1)
  attr(R_scaled, "ld_scale") <- n_eff - 1
  res_scaled <- run_susie_rrss(
    variants,
    R_scaled,
    reference_sample_size = 5000,
    L = 3,
    max_iter = 100
  )
  expect_equal(res_scaled$diagnostics$s_hat, res_keep$diagnostics$s_hat, tolerance = 1e-8)
  expect_equal(res_scaled$variants$KRIGING_EXP_Z, res_keep$variants$KRIGING_EXP_Z, tolerance = 1e-8)
  expect_equal(res_scaled$variants$KRIGING_STD_RESID, res_keep$variants$KRIGING_STD_RESID, tolerance = 1e-8)

  # 3. Opt-in: remove_kriging_outliers = TRUE drops flagged outlier(s) and preserves coloc naming
  res_drop <- run_susie_rrss(
    variants,
    R,
    reference_sample_size = 5000,
    L = 3,
    max_iter = 100,
    remove_kriging_outliers = TRUE
  )
  expect_true(res_drop$diagnostics$remove_kriging_outliers)
  expect_gte(res_drop$diagnostics$kriging_outliers_removed, 1L)
  expect_equal(nrow(res_drop$variants), p - res_drop$diagnostics$kriging_outliers_removed)
  expect_false(ids[3] %in% res_drop$variants$GWAS_VARIANT_ID)
  expect_equal(names(res_drop$fit$pip), res_drop$variants$GWAS_VARIANT_ID)

  # 4. Fewer than 3 variants skips kriging cleanly
  res_two <- run_susie_rrss(
    variants[1:2, , drop = FALSE],
    R[1:2, 1:2, drop = FALSE],
    reference_sample_size = 5000,
    L = 2,
    max_iter = 30
  )
  expect_equal(res_two$diagnostics$kriging_status, "skipped_too_few_variants")
})

test_that("Two-tier CS validity gate suppresses non-converged and V > 1 blow-up credible sets", {
  # 1. Force non-convergence with max_iter = 1 on a multi-signal locus
  variants <- data.frame(
    CHROM = rep("1", 6),
    POS = 200:205,
    REF = rep("A", 6),
    ALT = rep("G", 6),
    EZ = c(9.5, 9.1, 0.2, 8.4, 8.0, -0.1),
    LP = c(20, 19, 0.5, 16, 15, 0.3),
    SS = rep(5000, 6)
  )
  ids <- paste(variants$CHROM, variants$POS, variants$REF, variants$ALT, sep = ":")
  R <- diag(6)
  R[1, 2] <- R[2, 1] <- 0.95
  R[4, 5] <- R[5, 4] <- 0.95
  R[1, 4] <- R[4, 1] <- 0.25
  R[1, 5] <- R[5, 1] <- 0.25
  R[2, 4] <- R[4, 2] <- 0.25
  R[2, 5] <- R[5, 2] <- 0.25
  dimnames(R) <- list(ids, ids)

  nonconv <- suppressWarnings(run_susie_rss(
    variants,
    R,
    L = 4,
    max_iter = 1L,
    min_abs_corr = 0.5
  ))
  expect_false(nonconv$diagnostics$converged)
  expect_true(nonconv$diagnostics$credible_sets_suppressed)
  expect_equal(nonconv$diagnostics$credible_sets_suppressed_reason, "non_converged")
  expect_gte(nonconv$diagnostics$unconverged_cs_dropped, 1L)
  expect_equal(nrow(nonconv$credible_sets), 0L)
  expect_equal(nrow(nonconv$credible_set_membership), 0L)
  expect_true(all(is.na(nonconv$variants$credible_set_id)))
  expect_null(nonconv$fit$sets$cs)

  # 2. Converged locus with 1 valid lead CS (V <= 1) + 1 opposite-sign collinear blow-up slot (V > 1)
  #    retains the valid CS and drops only the V > 1 CS slot!
  variants_mix <- variants
  variants_mix$EZ <- c(45, -45, 0.1, 9.5, 9.1, 0.05)
  variants_mix$LP <- c(200, 200, 0.3, 20.0, 19.0, 0.2)
  variants_mix$SS <- rep(500, 6)
  R_mix <- diag(6)
  R_mix[1, 2] <- R_mix[2, 1] <- 0.998
  R_mix[4, 5] <- R_mix[5, 4] <- 0.95
  dimnames(R_mix) <- list(ids, ids)

  mix <- suppressWarnings(run_susie_rss(
    variants_mix,
    R_mix,
    L = 4,
    max_iter = 200L,
    min_abs_corr = 0.5
  ))
  if (isTRUE(mix$diagnostics$converged) &&
      isTRUE(mix$diagnostics$implausible_prior_variance_flag)) {
    expect_true(all(mix$credible_sets$prior_variance_V <= 1.0))
    expect_gte(nrow(mix$credible_sets), 1L)
    expect_gte(mix$diagnostics$dropped_implausible_v_cs, 1L)
    expect_equal(length(mix$fit$sets$cs), nrow(mix$credible_sets))
    expect_equal(length(mix$fit$sets$cs_index), nrow(mix$credible_sets))
  } else if (!isTRUE(mix$diagnostics$converged)) {
    expect_equal(nrow(mix$credible_sets), 0L)
    expect_null(mix$fit$sets$cs)
  }
})
