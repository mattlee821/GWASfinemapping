resolve_effective_n <- function(variants, effective_n, allow_missing = FALSE) {
  sample_sizes <- if ("SS" %in% names(variants) && is.numeric(variants$SS)) {
    variants$SS[!is.na(variants$SS)]
  } else {
    numeric()
  }
  summary <- if (length(sample_sizes) == 0L) {
    NULL
  } else {
    c(
      n_variants_with_ss = length(sample_sizes),
      min_ss = min(sample_sizes),
      median_ss = stats::median(sample_sizes),
      max_ss = max(sample_sizes)
    )
  }

  if (!is.null(effective_n)) {
    if (length(effective_n) != 1L || !is.numeric(effective_n) ||
        is.na(effective_n) || !is.finite(effective_n) || effective_n <= 0) {
      stop("`effective_n` must be NULL or one positive finite number.",
           call. = FALSE)
    }
    return(list(
      value = effective_n,
      source = "explicit_override",
      summary = summary
    ))
  }

  if (length(sample_sizes) == 0L) {
    if (isTRUE(allow_missing)) {
      return(list(
        value = NULL,
        source = "z_only_unreported_n",
        summary = NULL
      ))
    }
    if (!"SS" %in% names(variants)) {
      stop("SS is unavailable; supply an explicit `effective_n` override.",
           call. = FALSE)
    }
    stop("No valid per-variant SS values; supply an explicit `effective_n` override.",
         call. = FALSE)
  }

  list(
    value = min(sample_sizes),
    source = "minimum_variant_ss",
    summary = summary
  )
}

resolve_gwas_variant_ids <- function(variants, fallback_ids) {
  if ("GWAS_VARIANT_ID" %in% names(variants) &&
      all(nzchar(as.character(variants$GWAS_VARIANT_ID))) &&
      !anyNA(variants$GWAS_VARIANT_ID)) {
    return(as.character(variants$GWAS_VARIANT_ID))
  }
  if (all(c("GWAS_CHROM", "GWAS_POS", "GWAS_REF", "GWAS_ALT") %in%
          names(variants))) {
    chrom <- sub("^chr", "", as.character(variants$GWAS_CHROM),
                 ignore.case = TRUE)
    return(paste(
      chrom,
      as.integer(variants$GWAS_POS),
      toupper(as.character(variants$GWAS_REF)),
      toupper(as.character(variants$GWAS_ALT)),
      sep = ":"
    ))
  }
  fallback_ids
}

name_susie_fit_variants <- function(fit, variant_ids) {
  p <- length(variant_ids)
  if (!is.null(fit$pip) && length(fit$pip) == p) {
    names(fit$pip) <- variant_ids
  }
  if (!is.null(fit$XtXr) && length(fit$XtXr) == p) {
    names(fit$XtXr) <- variant_ids
  }
  for (slot in c("alpha", "lbf_variable", "mu", "mu2")) {
    mat <- fit[[slot]]
    if (is.matrix(mat) && ncol(mat) == p) {
      colnames(mat) <- variant_ids
      fit[[slot]] <- mat
    }
  }
  fit$snp <- variant_ids
  if (!is.null(fit$sets) && !is.null(fit$sets$cs) && length(fit$sets$cs) > 0L) {
    for (i in seq_along(fit$sets$cs)) {
      idx <- as.integer(fit$sets$cs[[i]])
      names(idx) <- variant_ids[idx]
      fit$sets$cs[[i]] <- idx
    }
  }
  fit
}

compute_credible_sets_fast <- function(fit, R, coverage = 0.95,
                                       min_abs_corr = 0.5,
                                       max_exact_purity = 500L,
                                       ld_scale = 1) {
  raw_cs <- susieR::susie_get_cs(
    fit,
    Xcorr = NULL,
    coverage = coverage,
    min_abs_corr = NULL,
    check_symmetric = FALSE
  )
  if (is.null(raw_cs$cs) || length(raw_cs$cs) == 0L) {
    return(list(
      cs = NULL,
      coverage = NULL,
      requested_coverage = coverage
    ))
  }
  cs_list <- raw_cs$cs
  effect_indices <- as.integer(sub("^L", "", names(cs_list)))
  claimed_cov <- as.numeric(raw_cs$coverage)

  purity_mat <- matrix(NA_real_, nrow = length(cs_list), ncol = 3L)
  for (i in seq_along(cs_list)) {
    pos <- as.integer(cs_list[[i]])
    if (length(pos) <= 1L) {
      purity_mat[i, ] <- c(1, 1, 1)
    } else {
      pos_eval <- if (length(pos) > max_exact_purity) {
        pos[unique(round(seq.int(1L, length(pos), length.out = 200L)))]
      } else {
        pos
      }
      sub_r <- abs(R[pos_eval, pos_eval, drop = FALSE]) / ld_scale
      vals <- sub_r[upper.tri(sub_r, diag = FALSE)]
      purity_mat[i, ] <- c(
        min(vals),
        sum(vals) / length(vals),
        stats::median(vals)
      )
    }
  }
  purity_df <- as.data.frame(purity_mat)
  colnames(purity_df) <- c("min.abs.corr", "mean.abs.corr", "median.abs.corr")

  keep <- if (!is.null(min_abs_corr)) {
    which(purity_df[, 1L] >= min_abs_corr)
  } else {
    seq_len(nrow(purity_df))
  }
  if (length(keep) == 0L) {
    return(list(
      cs = NULL,
      coverage = NULL,
      requested_coverage = coverage
    ))
  }

  cs_list <- cs_list[keep]
  purity_df <- purity_df[keep, , drop = FALSE]
  claimed_cov <- claimed_cov[keep]
  effect_indices <- effect_indices[keep]

  row_names <- paste0("L", effect_indices)
  names(cs_list) <- row_names
  rownames(purity_df) <- row_names

  ord <- order(purity_df[, 1L], decreasing = TRUE)
  list(
    cs = cs_list[ord],
    purity = purity_df[ord, , drop = FALSE],
    cs_index = effect_indices[ord],
    coverage = claimed_cov[ord],
    requested_coverage = coverage
  )
}

fit_kriging_mixture_em <- function(log_lik_shifted, max_iter = 500L, tol = 1e-8) {
  lik <- exp(log_lik_shifted)
  k <- ncol(lik)
  w <- rep(1 / k, k)
  for (it in seq_len(as.integer(max_iter))) {
    post <- sweep(lik, 2L, w, `*`, check.margin = FALSE)
    rsum <- rowSums(post)
    rsum[rsum <= 0] <- .Machine$double.xmin
    post <- post / rsum
    w_new <- colMeans(post)
    if (max(abs(w_new - w)) < tol) {
      return(w_new)
    }
    w <- w_new
  }
  w
}

# Exact O(M^2)-time, O(M)-extra-memory implementation of susieR::estimate_s_rss(method = "null-mle")
# and susieR::kriging_rss() reusing an already-computed correlation eigendecomposition `eig`.
# Avoids susieR's redundant O(M^3) `eigen(R)` and O(M^3) `V %*% (t(V) * dinv)` (which allocates
# three additional M x M matrices inside `susieR::kriging_rss`), using the exact identities:
#   condvar_i    = 1 / Omega_ii
#   condmean_i   = z_adj_i - (Omega z_adj)_i / Omega_ii
#   z_std_diff_i = (Omega z_adj)_i / sqrt(Omega_ii)
compute_kriging_from_eigen <- function(z, n = NULL, eig, r_tol = 1e-8) {
  p <- length(z)
  empty_res <- function(status_label) {
    list(
      status = status_label,
      s_hat = NA_real_,
      z_adj = rep(NA_real_, p),
      condmean = rep(NA_real_, p),
      condvar = rep(NA_real_, p),
      z_std_diff = rep(NA_real_, p),
      logLR = rep(NA_real_, p),
      outlier = rep(FALSE, p),
      dentist_outlier = rep(FALSE, p)
    )
  }
  if (p < 3L) {
    return(empty_res("skipped_too_few_variants"))
  }
  if (is.null(eig) || is.null(eig$values) || is.null(eig$vectors) ||
      length(eig$values) != p || nrow(eig$vectors) != p || ncol(eig$vectors) != p) {
    return(empty_res("skipped_missing_eigen"))
  }

  tryCatch({
    z_clean <- as.numeric(z)
    z_clean[!is.finite(z_clean)] <- 0
    if (!is.null(n) && is.finite(n) && n > 1) {
      adj <- (n - 1) / (z_clean * z_clean + n - 2)
      z_adj <- sqrt(adj) * z_clean
    } else {
      z_adj <- z_clean
    }

    V <- eig$vectors
    d <- as.numeric(eig$values)
    d[!is.finite(d) | d < r_tol] <- 0

    if (all(is.finite(z_adj))) {
      old_mp <- options(matprod = "blas")
      on.exit(options(old_mp), add = TRUE)
    }

    ztv <- as.numeric(crossprod(V, z_adj))
    ztv2 <- ztv * ztv
    negloglik <- function(s) {
      denom <- (1 - s) * d + s
      0.5 * sum(log(denom)) + 0.5 * sum(ztv2 / denom)
    }
    s_hat <- as.numeric(stats::optim(
      0.5,
      fn = negloglik,
      method = "Brent",
      lower = 0,
      upper = 1
    )$par)

    s_use <- if (s_hat > 1) 0.8 else max(s_hat, 0)
    dinv <- 1 / ((1 - s_use) * d + s_use)
    dinv[!is.finite(dinv)] <- 0

    omega_z <- as.numeric(V %*% (dinv * ztv))
    diag_prec <- if (is.loaded("gwasfm_eigen_diag_quadform", PACKAGE = "gwasfinemapping")) {
      .Call("gwasfm_eigen_diag_quadform", V, as.numeric(dinv), PACKAGE = "gwasfinemapping")
    } else {
      dp <- numeric(p)
      for (b in seq.int(1L, p, by = 512L)) {
        idx <- b:min(p, b + 511L)
        dp <- dp + drop((V[, idx, drop = FALSE]^2) %*% dinv[idx])
      }
      dp
    }

    valid_dp <- is.finite(diag_prec) & (diag_prec > 0)
    condvar <- ifelse(valid_dp, 1 / diag_prec, NA_real_)
    condmean <- ifelse(valid_dp, z_adj - omega_z * condvar, 0)
    z_std_diff <- ifelse(valid_dp, omega_z * sqrt(condvar), 0)

    a_min <- 0.8
    max_z2 <- max(z_std_diff * z_std_diff, na.rm = TRUE)
    a_max <- if (!is.finite(max_z2) || max_z2 < 1) 2 else 2 * sqrt(max_z2)
    npoint <- ceiling(log2(a_max / a_min) / log2(1.05))
    a_grid <- 1.05^(seq(-npoint, 0)) * a_max

    sd_mtx <- outer(sqrt(pmax( ifelse(is.na(condvar), 1, condvar), .Machine$double.eps )), a_grid)
    matrix_llik0 <- stats::dnorm(z_adj - condmean, sd = sd_mtx, log = TRUE)
    lfactors0 <- apply(matrix_llik0, 1L, max)
    matrix_llik0 <- matrix_llik0 - lfactors0

    w <- if (requireNamespace("mixsqp", quietly = TRUE)) {
      mixsqp::mixsqp(matrix_llik0, log = TRUE, control = list(verbose = FALSE))$x
    } else {
      fit_kriging_mixture_em(matrix_llik0)
    }

    logl0mix <- drop(log(exp(matrix_llik0) %*% (w + 1e-15))) + lfactors0
    matrix_llik1 <- stats::dnorm(z_adj + condmean, sd = sd_mtx, log = TRUE)
    lfactors1 <- apply(matrix_llik1, 1L, max)
    matrix_llik1 <- matrix_llik1 - lfactors1
    logl1mix <- drop(log(exp(matrix_llik1) %*% (w + 1e-15))) + lfactors1
    logLR <- as.numeric(logl1mix - logl0mix)
    logLR[!valid_dp] <- NA_real_

    dentist_chisq_thr <- stats::qchisq(5e-8, df = 1, lower.tail = FALSE)
    outlier <- !is.na(logLR) & is.finite(logLR) & (logLR > 2) & (abs(z_adj) > 2)
    dentist_outlier <- valid_dp & is.finite(z_std_diff) & ((z_std_diff * z_std_diff) > dentist_chisq_thr)

    list(
      status = "computed",
      s_hat = s_hat,
      z_adj = as.numeric(z_adj),
      condmean = as.numeric(condmean),
      condvar = as.numeric(condvar),
      z_std_diff = as.numeric(z_std_diff),
      logLR = as.numeric(logLR),
      outlier = as.logical(outlier),
      dentist_outlier = as.logical(dentist_outlier)
    )
  }, error = function(e) {
    res <- empty_res("failed")
    res$error <- conditionMessage(e)
    res
  })
}

fit_susie_from_prescaled_xtx <- function(z, XtX, n, L, coverage, max_iter,
                                         refine, prior_weights = NULL,
                                         R_finite = NULL, R_mismatch = "none",
                                         frob_sq = NULL, ld_scale = 1,
                                         precomputed_eigen = NULL) {
  ns <- asNamespace("susieR")
  apply_pve <- get("apply_pve_adjustment", envir = ns)
  working_fn <- get("summary_stats_working_quantities", envir = ns)
  norm_pw <- get("normalize_prior_weights", envir = ns)
  val_params <- get("validate_and_override_params", envir = ns)
  workhorse <- get("susie_workhorse", envir = ns)
  mix_p_to_z <- get("mixture_reference_p_to_z", envir = ns)

  z_clean <- as.numeric(z)
  z_clean[is.na(z_clean)] <- 0
  p <- length(z_clean)
  if (!is.null(colnames(XtX))) {
    names(z_clean) <- colnames(XtX)
  }
  pve <- apply_pve(z_clean, n, "wald")
  wk <- working_fn(
    z = pve$z,
    n = n,
    pve_adjustment = pve$adjustment,
    scaled_prior_variance = 0.2
  )
  if (!is.null(colnames(XtX))) {
    names(wk$Xty) <- colnames(XtX)
  }

  if (is.null(attr(XtX, "d"))) {
    d_vec <- diag(XtX)
    attr(XtX, "d") <- d_vec
    attr(XtX, "scaled:scale") <- rep(1, p)
  }

  pw <- norm_pw(prior_weights, p)
  conv_method <- if (!is.null(R_finite) || R_mismatch != "none") "pip" else "elbo"

  params_obj <- list(
    L = min(as.integer(L), p),
    scaled_prior_variance = wk$scaled_prior_variance,
    residual_variance = NULL,
    prior_weights = pw,
    null_weight = 0,
    estimate_residual_variance = FALSE,
    estimate_residual_method = "MoM",
    estimate_prior_variance = TRUE,
    estimate_prior_method = "optim",
    prior_variance_grid = NULL,
    mixture_weights = NULL,
    unmappable_effects = "none",
    check_null_threshold = 0,
    prior_tol = 1e-9,
    residual_variance_upperbound = Inf,
    model_init = NULL,
    coverage = coverage,
    min_abs_corr = NULL,
    median_abs_corr = NULL,
    compute_univariate_zscore = FALSE,
    max_iter = as.integer(max_iter),
    tol = NULL,
    convergence_method = conv_method,
    verbose = FALSE,
    track_fit = FALSE,
    residual_variance_lowerbound = 0,
    refine = refine,
    n_purity = "auto",
    cs_extension_corr = NULL,
    alpha0 = NULL,
    beta0 = NULL,
    n = wk$n,
    use_NIG = FALSE,
    intercept = FALSE,
    standardize = FALSE,
    check_prior = TRUE,
    slot_prior = NULL,
    L_greedy = NULL,
    greedy_lbf_cutoff = 0.1
  )
  params_obj <- val_params(params_obj)
  params_obj$coverage <- NULL
  params_obj$min_abs_corr <- NULL

  data_obj <- structure(
    list(
      XtX = XtX,
      X = NULL,
      Xty = wk$Xty,
      yty = wk$yty,
      n = wk$n,
      p = p,
      X_colmeans = rep(NA, p),
      y_mean = NA
    ),
    class = "ss"
  )
  ss <- list(data = data_obj, params = params_obj)

  if (!is.null(R_finite) || R_mismatch != "none") {
    B <- if (is.null(R_finite)) Inf else R_finite
    if (is.null(frob_sq) || !is.finite(frob_sq)) {
      if (is.loaded("gwasfm_frob_and_scale_inplace", PACKAGE = "gwasfinemapping")) {
        frob_sq <- .Call("gwasfm_frob_and_scale_inplace", XtX, 1.0,
                         PACKAGE = "gwasfinemapping") / (ld_scale * ld_scale)
      } else {
        frob_sq <- sum((XtX / ld_scale)^2)
      }
    }
    r_db <- if (is.infinite(B)) frob_sq else (B * frob_sq - p^2) / (B + 1)
    eff_r <- p^2 / max(r_db, 1)
    ss$data$R_finite_B <- B
    ss$data$R_finite_diagnostics <- list(
      B = B,
      p = p,
      R_frob_sq_debiased = r_db,
      effective_rank = eff_r,
      r_over_B = eff_r / B,
      Rhat_diag_deviation = abs(diag(XtX) / ld_scale - 1)
    )
    ss$data$R_mismatch <- R_mismatch
    if (R_mismatch != "none") {
      if (!is.null(precomputed_eigen)) {
        ss$data$eigen_R <- precomputed_eigen
      } else {
        eig <- eigen(XtX, symmetric = TRUE)
        eig$values <- eig$values / ld_scale
        ss$data$eigen_R <- eig
      }
    }
    ss$params$R_mismatch <- R_mismatch
    ss$params$R_mismatch_method <- "mle"
    ss$params$eb_mix_ref <- 5e-8
    ss$data$eb_mix_ref <- 5e-8
    ss$params$eb_mix_z_ref <- mix_p_to_z(5e-8)
    ss$data$eb_mix_z_ref <- ss$params$eb_mix_z_ref
    ss$data$z_marginal <- pve$z
    ss$params$eig_delta_rel <- 1e-3
    ss$params$eig_delta_abs <- 0
    ss$params$artifact_threshold <- 0.1
    ss$params$R_sensitivity_threshold <- 30
  }

  # XtX and z_clean are already verified 100% finite; bypass R's single-threaded O(M^2)
  # ISNAN pre-scan on every `XtX %*% v` so multi-threaded BLAS dgemv scales cleanly on MHC.
  if (all(is.finite(z_clean))) {
    old_matprod <- options(matprod = "blas")
    on.exit(options(old_matprod), add = TRUE)
  }
  model <- workhorse(ss$data, ss$params)
  if (!is.null(model$R_finite_diagnostics)) {
    model$R_finite_diagnostics$SER_fallback <-
      isTRUE(model$R_finite_diagnostics$R_reliability_flag)
  }
  model
}

run_susie_model <- function(variants, R, method_id, effective_n = NULL,
                            reference_sample_size = NULL, L = 10L,
                            coverage = 0.95, min_abs_corr = 0.5,
                            max_iter = 100L, refine = FALSE,
                            R_mismatch = "eb_mix",
                            psd_max_variants = 5000L,
                            eigen_max_variants = 20000L,
                            prior_weights = NULL,
                            allow_missing_n = FALSE,
                            remove_kriging_outliers = FALSE) {
  if (!method_id %in% c("susie_rss", "susie_rrss", "polyfun")) {
    stop("Unsupported SuSiE method identifier.", call. = FALSE)
  }
  if (!is.logical(remove_kriging_outliers) || length(remove_kriging_outliers) != 1L ||
      is.na(remove_kriging_outliers)) {
    stop("`remove_kriging_outliers` must be one non-missing logical value.",
         call. = FALSE)
  }
  if (length(eigen_max_variants) != 1L || !is.numeric(eigen_max_variants) ||
      is.na(eigen_max_variants) || eigen_max_variants < 1L) {
    stop("`eigen_max_variants` must be one positive integer.", call. = FALSE)
  }
  if (length(L) != 1L || !is.numeric(L) || is.na(L) || L < 1 ||
      L != floor(L)) {
    stop("`L` must be one positive integer.", call. = FALSE)
  }
  if (length(coverage) != 1L || !is.numeric(coverage) ||
      is.na(coverage) || !is.finite(coverage) || coverage <= 0 ||
      coverage > 1) {
    stop("`coverage` must be in (0, 1].", call. = FALSE)
  }
  if (!is.null(min_abs_corr) &&
      (length(min_abs_corr) != 1L || !is.numeric(min_abs_corr) ||
       is.na(min_abs_corr) || !is.finite(min_abs_corr) ||
       min_abs_corr < 0 || min_abs_corr > 1)) {
    stop("`min_abs_corr` must be NULL or a number in [0, 1].",
         call. = FALSE)
  }
  if (length(max_iter) != 1L || !is.numeric(max_iter) || is.na(max_iter) ||
      max_iter < 1 || max_iter != floor(max_iter)) {
    stop("`max_iter` must be one positive integer.", call. = FALSE)
  }
  if (!is.logical(refine) || length(refine) != 1L || is.na(refine)) {
    stop("`refine` must be one non-missing logical value.", call. = FALSE)
  }

  variants <- validate_gwas_stats(variants, require_p = TRUE, require_z = TRUE)
  if (nrow(variants) == 0L) {
    stop("SuSiE cannot fit a region with no variants.", call. = FALSE)
  }
  ids <- variant_key(variants)
  gwas_ids <- resolve_gwas_variant_ids(variants, ids)
  n <- resolve_effective_n(variants, effective_n, allow_missing = allow_missing_n)

  r_mismatch_requested <- if (method_id == "susie_rrss") R_mismatch else "none"
  q_art_eigen_status <- "not_applicable"
  if (method_id == "susie_rrss") {
    if (length(reference_sample_size) != 1L ||
        !is.numeric(reference_sample_size) ||
        is.na(reference_sample_size) || !is.finite(reference_sample_size) ||
        reference_sample_size <= 0) {
      stop("SuSiE-rRSS requires a positive LD reference sample size B.",
           call. = FALSE)
    }
    if (reference_sample_size <= 2) {
      stop("LD reference sample size B must be greater than 2.",
           call. = FALSE)
    }
    if (length(R_mismatch) != 1L || !is.character(R_mismatch) ||
        is.na(R_mismatch) ||
        !R_mismatch %in% c("none", "eb", "eb_mix")) {
      stop("`R_mismatch` must be one of 'none', 'eb', or 'eb_mix'.",
           call. = FALSE)
    }
    R_finite <- reference_sample_size
  } else {
    R_finite <- NULL
    R_mismatch <- "none"
  }

  ld_scale <- attr(R, "ld_scale")
  if (is.null(ld_scale) || !is.finite(ld_scale) || ld_scale <= 0) {
    ld_scale <- 1
  }
  frob_sq_attr <- attr(R, "frob_sq")

  will_compute_eigen_r <- (method_id == "susie_rrss") &&
    (nrow(R) <= as.integer(eigen_max_variants))

  ld_check <- validate_ld_matrix(
    R,
    variant_ids = ids,
    check_psd = !will_compute_eigen_r,
    psd_max_variants = psd_max_variants
  )
  psd_status <- attr(ld_check, "psd_status") %||% "passed"

  precomputed_eigen <- NULL
  kriging_res <- NULL
  if (method_id == "susie_rrss") {
    if (will_compute_eigen_r) {
      eig <- eigen(R, symmetric = TRUE)
      eig$values <- eig$values / ld_scale
      psd_tol <- if (nrow(R) <= as.integer(psd_max_variants)) 1e-8 else 1e-5
      if (min(eig$values) < -psd_tol) {
        stop("LD matrix is not positive semidefinite; no repair was applied.",
             call. = FALSE)
      }
      psd_status <- "passed"
      kriging_res <- compute_kriging_from_eigen(
        z = variants$EZ,
        n = n$value,
        eig = eig,
        r_tol = 1e-8
      )
      if (R_mismatch != "none") {
        precomputed_eigen <- eig
        q_art_eigen_status <- "computed"
      } else {
        rm(eig)
      }
    } else {
      kriging_res <- list(
        status = "skipped_large_m",
        s_hat = NA_real_,
        z_adj = rep(NA_real_, nrow(variants)),
        condmean = rep(NA_real_, nrow(variants)),
        condvar = rep(NA_real_, nrow(variants)),
        z_std_diff = rep(NA_real_, nrow(variants)),
        logLR = rep(NA_real_, nrow(variants)),
        outlier = rep(FALSE, nrow(variants)),
        dentist_outlier = rep(FALSE, nrow(variants))
      )
      if (R_mismatch != "none") {
        # Keep R_mismatch ("eb_mix" / "eb") active for SER init + lambda_bias + B_corrected
        # at ALL region sizes while skipping O(M^3) Q_art eigendecomposition above eigen_max_variants.
        precomputed_eigen <- list(
          values = 1.0,
          vectors = matrix(0.0, nrow = nrow(R), ncol = 1L)
        )
        q_art_eigen_status <- "skipped_large_m"
      }
    }
  }

  kriging_initial_outlier_ids <- if (!is.null(kriging_res) && any(kriging_res$outlier)) {
    gwas_ids[kriging_res$outlier]
  } else {
    character()
  }
  kriging_initial_outlier_count <- as.integer(length(kriging_initial_outlier_ids))
  kriging_initial_dentist_count <- if (!is.null(kriging_res) &&
                                       identical(kriging_res$status, "computed")) {
    as.integer(sum(kriging_res$dentist_outlier, na.rm = TRUE))
  } else {
    NA_integer_
  }
  kriging_max_abs_std_resid <- if (!is.null(kriging_res) && any(is.finite(kriging_res$z_std_diff))) {
    max(abs(kriging_res$z_std_diff), na.rm = TRUE)
  } else {
    NA_real_
  }
  kriging_max_loglr <- if (!is.null(kriging_res) && any(is.finite(kriging_res$logLR))) {
    max(kriging_res$logLR, na.rm = TRUE)
  } else {
    NA_real_
  }
  kriging_outliers_removed <- 0L
  kriging_removed_ids <- character()

  if (method_id == "susie_rrss" && isTRUE(remove_kriging_outliers) &&
      !is.null(kriging_res) && any(kriging_res$outlier)) {
    keep_k <- !kriging_res$outlier
    if (sum(keep_k) >= 2L) {
      kriging_outliers_removed <- as.integer(sum(!keep_k))
      kriging_removed_ids <- gwas_ids[!keep_k]
      variants <- variants[keep_k, , drop = FALSE]
      ids <- ids[keep_k]
      gwas_ids <- gwas_ids[keep_k]
      if (!is.null(prior_weights)) {
        prior_weights <- prior_weights[keep_k]
        prior_weights <- prior_weights / sum(prior_weights)
      }
      R <- R[keep_k, keep_k, drop = FALSE]
      attr(R, "ld_scale") <- ld_scale
      frob_sq_attr <- NULL
      if (R_mismatch != "none") {
        if (nrow(R) <= as.integer(eigen_max_variants)) {
          eig_sub <- eigen(R, symmetric = TRUE)
          eig_sub$values <- eig_sub$values / ld_scale
          precomputed_eigen <- eig_sub
          q_art_eigen_status <- "computed"
        } else {
          precomputed_eigen <- list(
            values = 1.0,
            vectors = matrix(0.0, nrow = nrow(R), ncol = 1L)
          )
          q_art_eigen_status <- "skipped_large_m"
        }
      }
      for (k_slot in c("z_adj", "condmean", "condvar", "z_std_diff", "logLR", "outlier", "dentist_outlier")) {
        kriging_res[[k_slot]] <- kriging_res[[k_slot]][keep_k]
      }
    }
  }

  expected_scale <- if (is.null(n$value)) 1 else (n$value - 1)
  if (abs(ld_scale - expected_scale) <= 1e-6 * max(1, expected_scale)) {
    xtx_use <- R
  } else {
    xtx_use <- R * (expected_scale / ld_scale)
    dimnames(xtx_use) <- dimnames(R)
  }

  fit <- fit_susie_from_prescaled_xtx(
    z = variants$EZ,
    XtX = xtx_use,
    n = n$value,
    L = as.integer(L),
    coverage = coverage,
    max_iter = as.integer(max_iter),
    refine = refine,
    prior_weights = prior_weights,
    R_finite = R_finite,
    R_mismatch = R_mismatch,
    frob_sq = frob_sq_attr,
    ld_scale = expected_scale,
    precomputed_eigen = precomputed_eigen
  )
  rm(xtx_use, precomputed_eigen)

  cs <- compute_credible_sets_fast(
    fit,
    R = R,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    ld_scale = ld_scale
  )

  raw_cs_count <- length(cs$cs)
  V_vec <- if (!is.null(fit$V)) as.numeric(fit$V) else numeric()
  lbf_vec <- if (!is.null(fit$lbf)) as.numeric(fit$lbf) else numeric()
  active_idx <- if (raw_cs_count > 0L) as.integer(cs$cs_index) else integer()

  is_converged <- isTRUE(fit$converged)
  max_V_all <- if (length(V_vec) > 0L && any(is.finite(V_vec))) max(V_vec, na.rm = TRUE) else NA_real_
  max_V_cs <- if (length(active_idx) > 0L && length(V_vec) >= max(active_idx) &&
                  any(is.finite(V_vec[active_idx]))) {
    max(V_vec[active_idx], na.rm = TRUE)
  } else {
    NA_real_
  }
  max_lbf_val <- if (length(lbf_vec) > 0L && any(is.finite(lbf_vec))) max(lbf_vec, na.rm = TRUE) else NA_real_

  # Two-tier credible-set validity gate:
  #   Tier 1: If IBSS failed to converge (!is_converged), suppress ALL credible sets at the locus.
  #   Tier 2: If IBSS converged, filter credible sets per-slot to prior variance V[l] <= 1.0
  #           (dropping only collinear opposite-sign blow-up slots while preserving valid CSs
  #           such as the primary cis-pQTL lead in fit$sets$cs and fit$sets$cs_index for coloc.susie).
  implausible_v_flag <- is.finite(max_V_cs) && (max_V_cs > 1.0)
  dropped_implausible_v_cs <- 0L
  if (!is_converged) {
    dropped_cs_count <- as.integer(raw_cs_count)
    suppress_cs <- TRUE
    suppress_reason <- if (implausible_v_flag) {
      "non_converged_and_implausible_prior_variance"
    } else {
      "non_converged"
    }
    cs <- list(
      cs = NULL,
      coverage = NULL,
      requested_coverage = coverage
    )
  } else if (raw_cs_count > 0L && implausible_v_flag) {
    cs_v <- if (length(V_vec) >= max(active_idx)) V_vec[active_idx] else rep(NA_real_, length(active_idx))
    keep_v <- which(is.na(cs_v) | (is.finite(cs_v) & cs_v <= 1.0))
    dropped_implausible_v_cs <- as.integer(raw_cs_count - length(keep_v))
    dropped_cs_count <- dropped_implausible_v_cs
    suppress_cs <- dropped_cs_count > 0L
    if (length(keep_v) == 0L) {
      suppress_reason <- "all_cs_implausible_prior_variance_gt_1"
      cs <- list(
        cs = NULL,
        coverage = NULL,
        requested_coverage = coverage
      )
    } else {
      suppress_reason <- "partial_cs_dropped_prior_variance_gt_1"
      cs$cs <- cs$cs[keep_v]
      cs$purity <- cs$purity[keep_v, , drop = FALSE]
      cs$cs_index <- cs$cs_index[keep_v]
      cs$coverage <- cs$coverage[keep_v]
    }
  } else {
    dropped_cs_count <- 0L
    suppress_cs <- FALSE
    suppress_reason <- NA_character_
  }

  fit$sets <- cs
  fit <- name_susie_fit_variants(fit, gwas_ids)

  r_rel_warn <- isTRUE(fit$R_finite_diagnostics$R_reliability_flag)
  cs_kriging_counts <- if (length(cs$cs) > 0L && !is.null(kriging_res) &&
                           identical(kriging_res$status, "computed")) {
    vapply(cs$cs, function(m_idx) {
      as.integer(sum(kriging_res$outlier[as.integer(m_idx)], na.rm = TRUE))
    }, integer(1))
  } else if (length(cs$cs) > 0L) {
    rep(NA_integer_, length(cs$cs))
  } else {
    integer()
  }
  kriging_cs_outlier_total <- if (!is.null(kriging_res) && identical(kriging_res$status, "computed")) {
    if (length(cs$cs) > 0L) {
      all_cs_members <- unique(unlist(lapply(cs$cs, as.integer), use.names = FALSE))
      as.integer(sum(kriging_res$outlier[all_cs_members], na.rm = TRUE))
    } else {
      0L
    }
  } else {
    NA_integer_
  }

  credible_sets <- if (length(cs$cs) == 0L) {
    base_cs_df <- data.frame(
      credible_set_id = character(),
      effect_index = integer(),
      size = integer(),
      coverage = numeric(),
      min_abs_corr = numeric(),
      mean_abs_corr = numeric(),
      median_abs_corr = numeric(),
      prior_variance_V = numeric(),
      lbf = numeric(),
      ld_mismatch_warning = logical(),
      stringsAsFactors = FALSE
    )
    if (method_id == "susie_rrss") {
      base_cs_df$kriging_outlier_count <- integer()
    }
    base_cs_df
  } else {
    eff_idx <- as.integer(cs$cs_index)
    cs_warn_vec <- if (method_id == "susie_rrss") {
      r_rel_warn | (!is.na(cs_kriging_counts) & cs_kriging_counts > 0L)
    } else {
      rep(r_rel_warn, length(eff_idx))
    }
    base_cs_df <- data.frame(
      credible_set_id = names(cs$cs),
      effect_index = eff_idx,
      size = as.integer(lengths(cs$cs)),
      coverage = as.numeric(cs$coverage),
      min_abs_corr = cs$purity$min.abs.corr,
      mean_abs_corr = cs$purity$mean.abs.corr,
      median_abs_corr = cs$purity$median.abs.corr,
      prior_variance_V = if (length(V_vec) >= max(eff_idx)) as.numeric(V_vec[eff_idx]) else rep(NA_real_, length(eff_idx)),
      lbf = if (length(lbf_vec) >= max(eff_idx)) as.numeric(lbf_vec[eff_idx]) else rep(NA_real_, length(eff_idx)),
      ld_mismatch_warning = cs_warn_vec,
      stringsAsFactors = FALSE
    )
    if (method_id == "susie_rrss") {
      base_cs_df$kriging_outlier_count <- cs_kriging_counts
    }
    base_cs_df
  }
  membership <- if (length(cs$cs) == 0L) {
    data.frame(
      credible_set_id = character(),
      effect_index = integer(),
      variant_id = character(),
      gwas_variant_id = character(),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, lapply(seq_along(cs$cs), function(i) {
      members <- as.integer(cs$cs[[i]])
      data.frame(
        credible_set_id = names(cs$cs)[[i]],
        effect_index = as.integer(cs$cs_index[[i]]),
        variant_id = ids[members],
        gwas_variant_id = gwas_ids[members],
        stringsAsFactors = FALSE
      )
    }))
  }

  cs_lookup <- rep(NA_character_, nrow(variants))
  if (length(cs$cs) > 0L) {
    for (i in seq_along(cs$cs)) {
      members <- as.integer(cs$cs[[i]])
      cs_lookup[members] <- names(cs$cs)[[i]]
    }
  }

  elbo_len <- length(fit$elbo)
  initial_elbo_val <- if (elbo_len >= 1L) as.numeric(fit$elbo[[1L]]) else NA_real_
  final_elbo_val <- if (elbo_len >= 1L) as.numeric(utils::tail(fit$elbo, 1L)) else NA_real_
  final_elbo_chg <- if (elbo_len >= 2L) {
    as.numeric(fit$elbo[[elbo_len]] - fit$elbo[[elbo_len - 1L]])
  } else {
    NA_real_
  }

  diagnostics <- list(
    converged = is_converged,
    convergence_status = if (is.null(fit$converged)) NA else fit$converged,
    convergence_reason = fit$convergence_reason %||% NA_character_,
    iterations = fit$niter %||% NA_integer_,
    initial_elbo = initial_elbo_val,
    final_elbo = final_elbo_val,
    final_elbo_change = final_elbo_chg,
    max_prior_variance_V = max_V_all,
    max_cs_prior_variance_V = max_V_cs,
    max_lbf = max_lbf_val,
    implausible_prior_variance_flag = implausible_v_flag,
    credible_sets_suppressed = suppress_cs,
    credible_sets_suppressed_reason = suppress_reason,
    unconverged_cs_dropped = dropped_cs_count,
    dropped_implausible_v_cs = dropped_implausible_v_cs,
    psd_check = psd_status,
    q_art_eigen = q_art_eigen_status,
    r_mismatch_requested = r_mismatch_requested,
    r_mismatch_effective = R_mismatch,
    lambda_bias = fit$R_finite_diagnostics$lambda_bias %||% NA_real_,
    B_corrected = fit$R_finite_diagnostics$B_corrected %||% NA_real_,
    Q_art = fit$R_finite_diagnostics$Q_art %||% NA_real_,
    artifact_evaluable = fit$R_finite_diagnostics$artifact_evaluable %||% NA,
    R_reliability_flag = fit$R_finite_diagnostics$R_reliability_flag %||% NA,
    SER_fallback = fit$R_finite_diagnostics$SER_fallback %||% NA,
    artifact_flag = fit$R_finite_diagnostics$artifact_flag %||% NA,
    R_sensitivity_flag = fit$R_sensitivity_flag %||% NA,
    R_finite_diagnostics = fit$R_finite_diagnostics
  )
  if (method_id == "susie_rrss" && !is.null(kriging_res)) {
    diagnostics$kriging_status <- kriging_res$status
    diagnostics$s_hat <- kriging_res$s_hat
    diagnostics$kriging_outlier_count <- if (identical(kriging_res$status, "computed")) {
      kriging_initial_outlier_count
    } else {
      NA_integer_
    }
    diagnostics$kriging_cs_outlier_count <- kriging_cs_outlier_total
    diagnostics$dentist_outlier_count <- kriging_initial_dentist_count
    diagnostics$max_abs_kriging_std_resid <- kriging_max_abs_std_resid
    diagnostics$max_kriging_loglr <- kriging_max_loglr
    diagnostics$remove_kriging_outliers <- isTRUE(remove_kriging_outliers)
    diagnostics$kriging_outliers_removed <- kriging_outliers_removed
    diagnostics$kriging_outlier_gwas_variant_ids <- utils::head(kriging_initial_outlier_ids, 50L)
    if (kriging_outliers_removed > 0L) {
      diagnostics$kriging_removed_gwas_variant_ids <- utils::head(kriging_removed_ids, 50L)
    }
  }

  lead_idx <- which.max(
    ifelse(is.na(variants$LP), -Inf, as.numeric(variants$LP))
  )
  ld_r2_to_lead <- if ("LD_R2_TO_LEAD" %in% names(variants) &&
                       !anyNA(variants$LD_R2_TO_LEAD)) {
    as.numeric(variants$LD_R2_TO_LEAD)
  } else if (length(lead_idx) == 1L && nrow(R) == nrow(variants)) {
    (as.numeric(R[, lead_idx]) / ld_scale)^2
  } else {
    rep(NA_real_, nrow(variants))
  }
  expected_z <- if (!is.null(fit$XtXr) && length(fit$XtXr) == nrow(variants) &&
                    !is.null(n$value) && n$value > 1) {
    as.numeric(fit$XtXr / sqrt(n$value - 1))
  } else if (!is.null(fit$XtXr) && length(fit$XtXr) == nrow(variants) && is.null(n$value)) {
    as.numeric(fit$XtXr)
  } else {
    rep(NA_real_, nrow(variants))
  }

  variant_results <- variants
  if (!"GWAS_VARIANT_ID" %in% names(variant_results)) {
    variant_results$GWAS_VARIANT_ID <- gwas_ids
  }
  if (!"LD_R2_TO_LEAD" %in% names(variant_results)) {
    variant_results$LD_R2_TO_LEAD <- ld_r2_to_lead
  }
  variant_results$method_id <- method_id
  variant_results$PIP <- as.numeric(fit$pip)
  variant_results$EXPECTED_Z <- expected_z
  if (method_id == "susie_rrss" && !is.null(kriging_res)) {
    variant_results$KRIGING_EXP_Z <- kriging_res$condmean
    variant_results$KRIGING_STD_RESID <- kriging_res$z_std_diff
    variant_results$KRIGING_LOGLR <- kriging_res$logLR
    variant_results$KRIGING_OUTLIER <- kriging_res$outlier
  }
  variant_results$credible_set_id <- cs_lookup

  prov_params <- list(
    L = as.integer(L),
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    max_iter = as.integer(max_iter),
    refine = refine,
    effective_n = n$value,
    effective_n_source = n$source,
    effective_n_summary = n$summary,
    R_finite = R_finite,
    R_mismatch = R_mismatch,
    R_mismatch_requested = r_mismatch_requested,
    psd_max_variants = as.integer(psd_max_variants),
    eigen_max_variants = as.integer(eigen_max_variants)
  )
  if (method_id == "susie_rrss") {
    prov_params$remove_kriging_outliers <- isTRUE(remove_kriging_outliers)
  }

  list(
    variants = variant_results,
    credible_sets = credible_sets,
    credible_set_membership = membership,
    diagnostics = diagnostics,
    provenance = list(
      method_id = method_id,
      method_version = as.character(utils::packageVersion("susieR")),
      parameters = prov_params
    ),
    fit = fit
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

run_susie_rss <- function(variants, R, effective_n = NULL, L = 10L,
                          coverage = 0.95, min_abs_corr = 0.5,
                          max_iter = 100L, refine = FALSE,
                          psd_max_variants = 5000L,
                          allow_missing_n = FALSE) {
  run_susie_model(
    variants = variants,
    R = R,
    method_id = "susie_rss",
    effective_n = effective_n,
    L = L,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    max_iter = max_iter,
    refine = refine,
    psd_max_variants = psd_max_variants,
    allow_missing_n = allow_missing_n
  )
}

run_susie_rrss <- function(variants, R, reference_sample_size,
                           effective_n = NULL, L = 10L, coverage = 0.95,
                           min_abs_corr = 0.5, max_iter = 100L,
                           refine = FALSE, R_mismatch = "eb_mix",
                           psd_max_variants = 5000L,
                           eigen_max_variants = 20000L,
                           allow_missing_n = FALSE,
                           remove_kriging_outliers = FALSE) {
  run_susie_model(
    variants = variants,
    R = R,
    method_id = "susie_rrss",
    effective_n = effective_n,
    reference_sample_size = reference_sample_size,
    L = L,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    max_iter = max_iter,
    refine = refine,
    R_mismatch = R_mismatch,
    psd_max_variants = psd_max_variants,
    eigen_max_variants = eigen_max_variants,
    allow_missing_n = allow_missing_n,
    remove_kriging_outliers = remove_kriging_outliers
  )
}
