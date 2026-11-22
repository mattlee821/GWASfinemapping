#!/usr/bin/env Rscript
# Serialize a MultiSuSiE cis_gene fit into GWASfinemapping's coloc-ready .rda contract.

parse_args <- function(args) {
  if (length(args) %% 2L != 0L) {
    stop("Options must be supplied as --key value pairs.", call. = FALSE)
  }
  opts <- list(compress = "xz")
  for (i in seq.int(1L, length(args), by = 2L)) {
    k <- gsub("-", "_", sub("^--", "", args[[i]]), fixed = TRUE)
    opts[[k]] <- args[[i + 1L]]
  }
  opts
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

opts <- parse_args(commandArgs(trailingOnly = TRUE))
work_dir <- normalizePath(opts$work_dir, mustWork = TRUE)
out_rda <- opts$output
out_format <- tolower(opts$format %||% "rda")

meta <- jsonlite::fromJSON(file.path(work_dir, "fit_meta.json"), simplifyVector = FALSE)
variants <- utils::read.delim(
  file.path(work_dir, "variants.tsv"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
)

P <- nrow(variants)
L <- as.integer(meta$L)
K <- as.integer(meta$K)
gwas_ids <- as.character(variants$GWAS_VARIANT_ID)
panel_ids <- as.character(variants$VARIANT_ID)

read_f64 <- function(fname, len) {
  con <- file(file.path(work_dir, fname), "rb")
  on.exit(close(con), add = TRUE)
  readBin(con, what = "numeric", n = len, size = 8L, endian = "little")
}

pip <- read_f64("pip.bin", P)
names(pip) <- gwas_ids

# Stored row-major (L x P) in Python -> read into P x L and transpose to L x P
alpha <- t(matrix(read_f64("alpha.bin", L * P), nrow = P, ncol = L))
colnames(alpha) <- gwas_ids

lbf_variable <- t(matrix(read_f64("lbf_variable.bin", L * P), nrow = P, ncol = L))
colnames(lbf_variable) <- gwas_ids

V_mat <- t(matrix(read_f64("V.bin", K * L), nrow = L, ncol = K))
lbf_vec <- as.numeric(unlist(meta$lbf))

cs_raw <- meta$credible_sets %||% list()
if (length(cs_raw) == 0L) {
  cs_list <- NULL
  cs_index <- integer()
  cov_vec <- numeric()
  purity_df <- data.frame(
    min.abs.corr = numeric(),
    mean.abs.corr = numeric(),
    median.abs.corr = numeric()
  )
  credible_sets_df <- data.frame(
    credible_set_id = character(),
    effect_index = integer(),
    size = integer(),
    coverage = numeric(),
    min_abs_corr = numeric(),
    mean_abs_corr = numeric(),
    median_abs_corr = numeric(),
    prior_variance_V = numeric(),
    lbf = numeric(),
    stringsAsFactors = FALSE
  )
  membership_df <- data.frame(
    credible_set_id = character(),
    effect_index = integer(),
    variant_id = character(),
    gwas_variant_id = character(),
    stringsAsFactors = FALSE
  )
} else {
  cs_names <- vapply(cs_raw, function(x) as.character(x$credible_set_id), character(1))
  cs_index <- vapply(cs_raw, function(x) as.integer(x$effect_index), integer(1))
  cov_vec <- vapply(cs_raw, function(x) as.numeric(x$coverage), numeric(1))
  min_c <- vapply(cs_raw, function(x) as.numeric(x$min_abs_corr), numeric(1))
  mean_c <- vapply(cs_raw, function(x) as.numeric(x$mean_abs_corr), numeric(1))
  med_c <- vapply(cs_raw, function(x) as.numeric(x$median_abs_corr), numeric(1))
  v_max <- vapply(cs_raw, function(x) as.numeric(x$max_std_prior_variance_V), numeric(1))
  lbf_cs <- vapply(cs_raw, function(x) as.numeric(x$lbf), numeric(1))

  cs_list <- stats::setNames(vector("list", length(cs_raw)), cs_names)
  mem_rows <- vector("list", length(cs_raw))
  for (i in seq_along(cs_raw)) {
    idx_1based <- as.integer(unlist(cs_raw[[i]]$indices_1based))
    names(idx_1based) <- gwas_ids[idx_1based]
    cs_list[[i]] <- idx_1based
    mem_rows[[i]] <- data.frame(
      credible_set_id = cs_names[[i]],
      effect_index = cs_index[[i]],
      variant_id = panel_ids[idx_1based],
      gwas_variant_id = gwas_ids[idx_1based],
      stringsAsFactors = FALSE
    )
  }
  purity_df <- data.frame(
    min.abs.corr = min_c,
    mean.abs.corr = mean_c,
    median.abs.corr = med_c,
    row.names = cs_names
  )
  credible_sets_df <- data.frame(
    credible_set_id = cs_names,
    effect_index = cs_index,
    size = lengths(cs_list),
    coverage = cov_vec,
    min_abs_corr = min_c,
    mean_abs_corr = mean_c,
    median_abs_corr = med_c,
    prior_variance_V = v_max,
    lbf = lbf_cs,
    stringsAsFactors = FALSE
  )
  membership_df <- do.call(rbind, mem_rows)
}

fit_obj <- list(
  alpha = alpha,
  lbf_variable = lbf_variable,
  lbf = lbf_vec,
  V = V_mat,
  pip = pip,
  snp = gwas_ids,
  converged = isTRUE(meta$converged),
  niter = as.integer(meta$iterations %||% NA_integer_),
  sets = list(
    cs = cs_list,
    purity = purity_df,
    cs_index = cs_index,
    coverage = cov_vec,
    requested_coverage = as.numeric(meta$requested_coverage %||% 0.95)
  )
)
class(fit_obj) <- c("multisusie_fit", "susie")

region_entry <- list(
  region_id = as.character(meta$region_id),
  region_type = "cis_gene",
  result = list(
    variants = variants,
    credible_sets = credible_sets_df,
    credible_set_membership = membership_df,
    diagnostics = meta$diagnostics,
    provenance = meta$provenance,
    fit = fit_obj
  )
)

regions_named <- stats::setNames(list(region_entry), as.character(meta$region_id))
method_object <- list(
  schema_version = 1L,
  method_id = "multisusie",
  regions = regions_named
)

dir.create(dirname(out_rda), recursive = TRUE, showWarnings = FALSE)
stem_path <- sub("\\.(rda|rds)$", "", out_rda, ignore.case = TRUE)
if (out_format %in% c("rda", "both")) {
  env_out <- new.env(parent = emptyenv())
  assign("multisusie", method_object, envir = env_out)
  save(list = "multisusie", envir = env_out, file = paste0(stem_path, ".rda"), compress = opts$compress)
}
if (out_format %in% c("rds", "both")) {
  saveRDS(method_object, file = paste0(stem_path, ".rds"), compress = opts$compress)
}
