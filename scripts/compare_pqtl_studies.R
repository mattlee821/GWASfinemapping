#!/usr/bin/env Rscript
# ==============================================================================
# Cross-Study, Cross-Platform (Olink vs SomaScan), and Cross-Ancestry (EUR/AFR/ALL)
# Fine-Mapping & Colocalisation Comparison Report Generator
# ==============================================================================
# Usage:
#   Rscript scripts/compare_pqtl_studies.R \
#     --finemap-dir /path/to/finemapping \
#     --outdir /path/to/comparison_report
# ==============================================================================

parse_cli <- function(args) {
  opts <- list(
    finemap_dir = "finemapping",
    outdir = "comparison_report",
    p1 = 1e-4,
    p2 = 1e-4,
    p12 = 5e-6
  )
  if (length(args) %% 2L != 0L) {
    stop("Supply CLI options as --key value pairs.", call. = FALSE)
  }
  for (i in seq.int(1L, length(args), by = 2L)) {
    k <- gsub("-", "_", sub("^--", "", args[[i]]), fixed = TRUE)
    opts[[k]] <- args[[i + 1L]]
  }
  opts$p1 <- as.numeric(opts$p1)
  opts$p2 <- as.numeric(opts$p2)
  opts$p12 <- as.numeric(opts$p12)
  opts
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

html_esc <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

log_sum_exp <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0L) return(-Inf)
  m <- max(v)
  m + log(sum(exp(v - m)))
}

# Exact analytical coloc.bf_bf over validated SuSiE credible-set rows (fit$sets$cs_index)
# Matches Wallace (2021, PLOS Genetics) without requiring the external `coloc` R package.
compute_susie_coloc_summary <- function(fit1, fit2, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6) {
  if (is.null(fit1) || is.null(fit2) ||
      is.null(fit1$sets$cs) || length(fit1$sets$cs) == 0L ||
      is.null(fit2$sets$cs) || length(fit2$sets$cs) == 0L ||
      !is.matrix(fit1$lbf_variable) || !is.matrix(fit2$lbf_variable)) {
    return(list(n_pairs = 0L, max_pph4 = NA_real_, max_pph3 = NA_real_, top_hit_pair = NA_character_))
  }
  common_snps <- intersect(colnames(fit1$lbf_variable), colnames(fit2$lbf_variable))
  if (length(common_snps) < 5L) {
    return(list(n_pairs = 0L, max_pph4 = NA_real_, max_pph3 = NA_real_, top_hit_pair = NA_character_))
  }
  idx1 <- as.integer(fit1$sets$cs_index)
  idx2 <- as.integer(fit2$sets$cs_index)
  m1 <- fit1$lbf_variable[idx1, common_snps, drop = FALSE]
  m2 <- fit2$lbf_variable[idx2, common_snps, drop = FALSE]
  cs_names1 <- names(fit1$sets$cs)
  cs_names2 <- names(fit2$sets$cs)

  best_h4 <- -1
  best_h3 <- -1
  best_lbl <- NA_character_
  n_pairs <- 0L

  for (i in seq_len(nrow(m1))) {
    b1 <- as.numeric(m1[i, ])
    l1_sum <- log_sum_exp(b1)
    for (j in seq_len(nrow(m2))) {
      b2 <- as.numeric(m2[j, ])
      l2_sum <- log_sum_exp(b2)
      l12_sum <- log_sum_exp(b1 + b2)

      lh0 <- 0
      lh1 <- log(p1) + l1_sum
      lh2 <- log(p2) + l2_sum
      # log(exp(l1_sum + l2_sum) - exp(l12_sum)) guarded against both cancellation and >709 overflow
      diff_arg <- (l1_sum + l2_sum) - l12_sum
      lh3 <- if (diff_arg > 35) {
        log(p1) + log(p2) + l12_sum + diff_arg
      } else if (diff_arg > 1e-12) {
        log(p1) + log(p2) + l12_sum + log(expm1(diff_arg))
      } else {
        -Inf
      }
      lh4 <- log(p12) + l12_sum
      all_lh <- c(lh0, lh1, lh2, lh3, lh4)
      denom <- log_sum_exp(all_lh)
      pp <- exp(all_lh - denom)
      n_pairs <- n_pairs + 1L
      if (pp[[5L]] > best_h4) {
        best_h4 <- pp[[5L]]
        best_h3 <- pp[[4L]]
        top_snp <- common_snps[[which.max(b1 + b2)]]
        best_lbl <- sprintf("%s×%s (%s)", cs_names1[[i]], cs_names2[[j]], top_snp)
      }
    }
  }
  list(
    n_pairs = n_pairs,
    max_pph4 = if (best_h4 >= 0) round(best_h4, 4) else NA_real_,
    max_pph3 = if (best_h3 >= 0) round(best_h3, 4) else NA_real_,
    top_hit_pair = best_lbl
  )
}

load_method_container <- function(methods_dir, stem, obj_name) {
  rds_p <- file.path(methods_dir, paste0(stem, ".rds"))
  if (file.exists(rds_p)) return(readRDS(rds_p))
  rda_p <- file.path(methods_dir, paste0(stem, ".rda"))
  if (file.exists(rda_p)) {
    e <- new.env(parent = emptyenv())
    load(rda_p, envir = e)
    if (exists(obj_name, envir = e, inherits = FALSE)) return(get(obj_name, envir = e, inherits = FALSE))
    nms <- ls(e)
    if (length(nms) == 1L) return(get(nms[[1L]], envir = e))
  }
  NULL
}

extract_cis_method <- function(obj) {
  if (is.null(obj) || is.null(obj$regions) || length(obj$regions) == 0L) return(NULL)
  cis <- Filter(function(x) identical(x$region_type, "cis_gene"), obj$regions)
  if (length(cis) == 0L) return(NULL)
  cis[[1L]]$result
}

opts <- parse_cli(commandArgs(trailingOnly = TRUE))
finemap_dir <- normalizePath(opts$finemap_dir, mustWork = TRUE)
outdir <- opts$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

target_genes <- c(
  Q9BYZ8 = "REG4", Q03403 = "TFF2", O95994 = "AGR2", P08727 = "KRT19",
  Q07654 = "TFF3", Q99988 = "GDF15", P39900 = "MMP12", O95407 = "TNFRSF6B",
  P04155 = "TFF1", Q16790 = "CA9", P50120 = "RBP2", Q8IZP9 = "ADGRG2",
  P10145 = "CXCL8", Q9P0M4 = "IL17C", O15444 = "CCL25", O14836 = "TNFRSF13B",
  O00220 = "TNFRSF10A", P80098 = "CCL7", P01589 = "IL2RA", P02745 = "C1QA",
  Q9ULX7 = "CA14", O43508 = "TNFSF12", P01579 = "IFNG", P78556 = "CCL20",
  P10147 = "CCL3", P13500 = "CCL2", P09237 = "MMP7", P09603 = "CSF1"
)

pqtl_studies <- c(
  "sun_2023_37794186_EUR",
  "sun_2023_37794186_ALL",
  "ferkingstad_2021_34857953",
  "pietzner_2021_34648354",
  "zhang_2022_35501419_EUR",
  "zhang_2022_35501419_AFR",
  "somascan_multisusie_3cohort",
  "somascan_multisusie_4cohort"
)

records <- list()
cis_cache <- list()
study_cs_sizes <- list()

for (s_id in pqtl_studies) {
  s_path <- file.path(finemap_dir, s_id)
  if (!dir.exists(s_path)) next
  g_dirs <- sort(list.dirs(s_path, full.names = TRUE, recursive = FALSE))
  for (gd in g_dirs) {
    qcp <- file.path(gd, "qc.json")
    prp <- file.path(gd, "provenance.json")
    if (!file.exists(qcp) || !file.exists(prp)) next
    gid <- basename(gd)
    u_prefix <- sub("_.*$", "", gid)
    qc <- jsonlite::fromJSON(qcp, simplifyVector = FALSE)
    prov <- jsonlite::fromJSON(prp, simplifyVector = FALSE)
    cis_meta <- qc$cis_gene %||% list()
    uniprot <- as.character(cis_meta$uniprot_id %||% (prov$gene_annotation$uniprot_id %||% u_prefix))
    if (!uniprot %in% names(target_genes)) next
    gene_sym <- target_genes[[uniprot]]

    cis_loci <- Filter(function(l) identical(l$region_type, "cis_gene"), qc$loci %||% list())
    has_cis <- length(cis_loci) > 0L && identical(cis_loci[[1L]]$status, "analyzed")
    loc <- if (has_cis) cis_loci[[1L]] else list()
    m <- loc$methods %||% list()
    is_msusie <- !is.null(m$multisusie)
    rrss <- (if (is_msusie) m$multisusie else m$susie_rrss) %||% list()
    rss <- m$susie_rss %||% list()
    crm <- m$carma %||% list()
    pf <- m$polyfun %||% list()
    pld_r2 <- m$p_ld$selected_by_r2 %||% list()

    rec_key <- paste0(s_id, "/", gid)
    platform <- if (grepl("^sun_2023", s_id)) "Olink_Explore" else "SomaScan_v4"

    rrss_res <- NULL
    crm_res <- NULL
    cs_sizes_here <- integer()
    max_pip_here <- NA_real_
    if (has_cis) {
      mdir <- file.path(gd, "methods")
      rrss_res <- if (is_msusie) {
        extract_cis_method(load_method_container(mdir, "multisusie", "multisusie"))
      } else {
        extract_cis_method(load_method_container(mdir, "susie-rrss", "susie_rrss"))
      }
      crm_res <- extract_cis_method(load_method_container(mdir, "carma", "carma"))
      if (!is.null(rrss_res$credible_sets) && nrow(rrss_res$credible_sets) > 0L) {
        cs_sizes_here <- as.integer(rrss_res$credible_sets$size)
        study_cs_sizes[[s_id]] <- c(study_cs_sizes[[s_id]] %||% integer(), cs_sizes_here)
      }
      if (!is.null(rrss_res$variants) && "PIP" %in% names(rrss_res$variants) && nrow(rrss_res$variants) > 0L) {
        max_pip_here <- round(max(as.numeric(rrss_res$variants$PIP), na.rm = TRUE), 4)
      }
      cis_cache[[rec_key]] <- list(
        study_id = s_id,
        platform = platform,
        ld_pop = as.character(prov$ld_population %||% ""),
        gwas_id = gid,
        uniprot_id = uniprot,
        gene_symbol = gene_sym,
        lead_id = as.character(loc$lead_gwas_variant_id %||% ""),
        rrss = rrss_res,
        carma = crm_res
      )
    }

    records[[length(records) + 1L]] <- data.frame(
      study_id = s_id,
      platform = platform,
      gwas_ancestry = as.character(prov$gwas_ancestry %||% ""),
      ld_population = as.character(prov$ld_population %||% ""),
      gwas_id = gid,
      uniprot_id = uniprot,
      gene_symbol = gene_sym,
      overall_status = as.character(qc$overall_status %||% ""),
      cis_status = as.character(cis_meta$status %||% ""),
      total_loci = as.integer(qc$analyzed_regions_count %||% 0L),
      cis_variants_matched = as.integer(loc$variants_ld_matched %||% 0L),
      cis_removed_min_maf = as.integer(loc$variants_removed_min_maf %||% 0L),
      lead_variant_grch38 = as.character(loc$lead_gwas_variant_id %||% NA_character_),
      lead_lp = as.numeric(loc$lead_lp %||% cis_meta$max_lp %||% NA_real_),
      pld_r2_0.001 = as.integer(pld_r2[["0.001"]] %||% 0L),
      rss_cs = if (has_cis && !is_msusie) as.integer(rss$credible_sets %||% 0L) else NA_integer_,
      rrss_cs = if (has_cis) as.integer(rrss$credible_sets %||% 0L) else NA_integer_,
      rrss_single_snp_cs = if (has_cis) sum(cs_sizes_here == 1L) else NA_integer_,
      rrss_cs_le_5 = if (has_cis) sum(cs_sizes_here <= 5L) else NA_integer_,
      rrss_median_cs_size = if (length(cs_sizes_here) > 0L) stats::median(cs_sizes_here) else NA_real_,
      rrss_max_pip = max_pip_here,
      rrss_dropped_cs = if (has_cis) as.integer(rrss$unconverged_cs_dropped %||% 0L) else NA_integer_,
      rrss_converged = if (has_cis) as.logical(rrss$converged %||% NA) else NA,
      rrss_iterations = if (has_cis) as.integer(rrss$iterations %||% NA_integer_) else NA_integer_,
      rrss_max_V = if (has_cis) as.numeric(rrss$max_prior_variance_V %||% NA_real_) else NA_real_,
      rrss_lambda_bias = if (has_cis && !is_msusie) as.numeric(rrss$lambda_bias %||% NA_real_) else NA_real_,
      rrss_Q_art = if (has_cis && !is_msusie) as.numeric(rrss$Q_art %||% NA_real_) else NA_real_,
      carma_cs = if (has_cis && !is_msusie) as.integer(crm$credible_sets %||% 0L) else NA_integer_,
      carma_outliers = if (has_cis && !is_msusie) as.integer(crm$outlier_count %||% 0L) else NA_integer_,
      polyfun_cs = if (has_cis && !is_msusie) as.integer(pf$credible_sets %||% 0L) else NA_integer_,
      runtime_seconds = round(as.numeric(prov$runtime_seconds %||% NA_real_), 1),
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- if (length(records) > 0L) do.call(rbind, records) else data.frame()
summary_df <- summary_df[order(summary_df$gene_symbol, summary_df$study_id, summary_df$gwas_id), , drop = FALSE]
utils::write.table(
  summary_df,
  file = file.path(outdir, "01_cis_method_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Pairwise concordance & coloc.susie evaluation within each protein (including CXCL8 replicates & cross-study)
pair_rows <- list()
keys <- names(cis_cache)
if (length(keys) >= 2L) {
  for (a in seq_len(length(keys) - 1L)) {
    ca <- cis_cache[[keys[[a]]]]
    for (b in (a + 1L):length(keys)) {
      cb <- cis_cache[[keys[[b]]]]
      if (!identical(ca$uniprot_id, cb$uniprot_id)) next

      va <- ca$rrss$variants
      vb <- cb$rrss$variants
      shared_ids <- if (!is.null(va) && !is.null(vb)) {
        intersect(va$GWAS_VARIANT_ID, vb$GWAS_VARIANT_ID)
      } else {
        character()
      }
      pip_cor_rrss <- NA_real_
      if (length(shared_ids) >= 10L) {
        pa <- va$PIP[match(shared_ids, va$GWAS_VARIANT_ID)]
        pb <- vb$PIP[match(shared_ids, vb$GWAS_VARIANT_ID)]
        if (stats::sd(pa, na.rm = TRUE) > 0 && stats::sd(pb, na.rm = TRUE) > 0) {
          pip_cor_rrss <- round(stats::cor(pa, pb, use = "complete.obs"), 4)
        }
      }

      cva <- ca$carma$variants
      cvb <- cb$carma$variants
      pip_cor_carma <- NA_real_
      if (!is.null(cva) && !is.null(cvb)) {
        sh_c <- intersect(cva$GWAS_VARIANT_ID, cvb$GWAS_VARIANT_ID)
        if (length(sh_c) >= 10L) {
          pca <- cva$PIP[match(sh_c, cva$GWAS_VARIANT_ID)]
          pcb <- cvb$PIP[match(sh_c, cvb$GWAS_VARIANT_ID)]
          if (stats::sd(pca, na.rm = TRUE) > 0 && stats::sd(pcb, na.rm = TRUE) > 0) {
            pip_cor_carma <- round(stats::cor(pca, pcb, use = "complete.obs"), 4)
          }
        }
      }

      mem_a <- ca$rrss$credible_set_membership$gwas_variant_id %||% character()
      mem_b <- cb$rrss$credible_set_membership$gwas_variant_id %||% character()
      shared_cs_snps <- length(intersect(unique(mem_a), unique(mem_b)))

      cmem_a <- ca$carma$credible_set_membership$gwas_variant_id %||% character()
      cmem_b <- cb$carma$credible_set_membership$gwas_variant_id %||% character()
      shared_carma_cs_snps <- length(intersect(unique(cmem_a), unique(cmem_b)))

      coloc_out <- compute_susie_coloc_summary(
        if (!is.null(ca$rrss)) ca$rrss$fit else NULL,
        if (!is.null(cb$rrss)) cb$rrss$fit else NULL,
        p1 = opts$p1, p2 = opts$p2, p12 = opts$p12
      )

      comparison_type <- if (identical(ca$study_id, cb$study_id)) {
        "within_study_replicate"
      } else if (grepl("^sun_2023", ca$study_id) && grepl("^sun_2023", cb$study_id)) {
        if (identical(ca$gwas_id, cb$gwas_id)) "sun_eur_vs_all_same_assay" else "sun_eur_vs_all_cross_panel"
      } else if (!identical(ca$platform, cb$platform)) {
        "cross_platform_olink_vs_somascan"
      } else if (!identical(ca$ld_pop, cb$ld_pop)) {
        "cross_ancestry_somascan"
      } else {
        "within_platform_somascan"
      }

      pair_rows[[length(pair_rows) + 1L]] <- data.frame(
        gene_symbol = ca$gene_symbol,
        uniprot_id = ca$uniprot_id,
        comparison_type = comparison_type,
        study_a = ca$study_id,
        gwas_a = ca$gwas_id,
        study_b = cb$study_id,
        gwas_b = cb$gwas_id,
        same_lead_variant = identical(ca$lead_id, cb$lead_id) && nzchar(ca$lead_id),
        shared_cis_variants = length(shared_ids),
        shared_rrss_cs_variants = shared_cs_snps,
        shared_carma_cs_variants = shared_carma_cs_snps,
        pip_cor_rrss = pip_cor_rrss,
        pip_cor_carma = pip_cor_carma,
        coloc_cs_pairs_tested = coloc_out$n_pairs,
        coloc_max_PPH4 = coloc_out$max_pph4,
        coloc_max_PPH3 = coloc_out$max_pph3,
        coloc_top_pair = coloc_out$top_hit_pair,
        stringsAsFactors = FALSE
      )
    }
  }
}

pairs_df <- if (length(pair_rows) > 0L) do.call(rbind, pair_rows) else data.frame()
if (nrow(pairs_df) > 0L) {
  pairs_df <- pairs_df[order(pairs_df$gene_symbol, pairs_df$comparison_type, pairs_df$study_a), , drop = FALSE]
}
utils::write.table(
  pairs_df,
  file = file.path(outdir, "02_pairwise_coloc_and_concordance.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cxcl8_df <- if (nrow(pairs_df) > 0L) {
  pairs_df[pairs_df$gene_symbol == "CXCL8" & pairs_df$comparison_type == "within_study_replicate", , drop = FALSE]
} else {
  data.frame()
}
utils::write.table(
  cxcl8_df,
  file = file.path(outdir, "03_cxcl8_technical_replicates.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Study-level aggregate summary
studies_present <- unique(summary_df$study_id)
agg_rows <- lapply(studies_present, function(sid) {
  sub <- summary_df[summary_df$study_id == sid, , drop = FALSE]
  act <- sub[sub$cis_variants_matched > 0L, , drop = FALSE]
  act_cs <- act[!is.na(act$rrss_cs) & act$rrss_cs > 0L, , drop = FALSE]
  all_sz <- study_cs_sizes[[sid]] %||% integer()
  data.frame(
    study_id = sid,
    platform = sub$platform[[1L]],
    gwas_ancestry = sub$gwas_ancestry[[1L]],
    ld_population = sub$ld_population[[1L]],
    assays_total = nrow(sub),
    cis_activated = nrow(act),
    rss_with_cs = sum(!is.na(act$rss_cs) & act$rss_cs > 0L),
    rrss_converged = sum(!is.na(act$rrss_converged) & act$rrss_converged),
    rrss_with_cs = nrow(act_cs),
    total_rrss_cs = sum(act$rrss_cs, na.rm = TRUE),
    single_snp_cs = sum(all_sz == 1L),
    cs_le_5 = sum(all_sz <= 5L),
    median_cs_size = if (length(all_sz) > 0L) stats::median(all_sz) else NA_real_,
    median_max_pip = if (nrow(act_cs) > 0L) round(stats::median(act_cs$rrss_max_pip, na.rm = TRUE), 4) else NA_real_,
    carma_with_cs = sum(!is.na(act$carma_cs) & act$carma_cs > 0L),
    polyfun_with_cs = sum(!is.na(act$polyfun_cs) & act$polyfun_cs > 0L),
    loci_max_V_gt_1 = sum(!is.na(act$rrss_max_V) & act$rrss_max_V > 1.0),
    median_rrss_iter = if (nrow(act) > 0L) stats::median(act$rrss_iterations, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
})
agg_df <- if (length(agg_rows) > 0L) do.call(rbind, agg_rows) else data.frame()
utils::write.table(
  agg_df,
  file = file.path(outdir, "00_study_level_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Render self-contained HTML dashboard
make_html_table <- function(df, max_rows = 500L) {
  if (nrow(df) == 0L) return("<p class='muted'>No rows available yet.</p>")
  df_show <- utils::head(df, max_rows)
  ths <- paste0("<th>", vapply(names(df_show), html_esc, character(1)), "</th>", collapse = "")
  trs <- vapply(seq_len(nrow(df_show)), function(i) {
    tds <- vapply(df_show[i, , drop = FALSE], function(col) {
      val <- col[[1L]]
      if (is.na(val)) return("<td class='muted'>&mdash;</td>")
      if (is.numeric(val)) return(sprintf("<td class='num'>%s</td>", html_esc(format(val, digits = 4, trim = TRUE))))
      sprintf("<td>%s</td>", html_esc(val))
    }, character(1))
    paste0("<tr>", paste(tds, collapse = ""), "</tr>")
  }, character(1))
  paste0("<div class='table-wrap'><table><thead><tr>", ths, "</tr></thead><tbody>",
         paste(trs, collapse = "\n"), "</tbody></table></div>")
}

html_out <- sprintf(
'<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>Cross-Study pQTL Fine-Mapping &amp; Colocalisation Comparison</title>
<style>
:root { --bg:#f8fafc; --card:#fff; --fg:#0f172a; --muted:#64748b; --border:#e2e8f0; }
@media (prefers-color-scheme: dark) {
  :root { --bg:#0f172a; --card:#1e293b; --fg:#f1f5f9; --muted:#94a3b8; --border:#334155; }
}
body { margin:0; padding:24px; background:var(--bg); color:var(--fg); font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; font-size:13.5px; line-height:1.45; }
.container { max-width:1480px; margin:0 auto; }
h1 { margin:0 0 6px 0; font-size:22px; }
h2 { margin:24px 0 8px 0; font-size:16px; }
.muted { color:var(--muted); }
.panel { background:var(--card); border:1px solid var(--border); border-radius:8px; padding:14px 16px; margin-bottom:16px; }
.table-wrap { overflow-x:auto; max-height:620px; overflow-y:auto; }
table { width:100%%; border-collapse:collapse; font-size:12.5px; }
th, td { padding:6px 9px; border-bottom:1px solid var(--border); text-align:left; white-space:nowrap; }
th { position:sticky; top:0; background:var(--card); font-weight:600; color:var(--muted); }
.num { text-align:right; font-variant-numeric:tabular-nums; }
</style>
</head>
<body>
<div class="container">
  <h1>Cross-Study pQTL Fine-Mapping &amp; Colocalisation Report (28 CRC Target Proteins)</h1>
  <p class="muted">Studies analyzed: %s &middot; Generated from <code>%s</code></p>

  <h2>1. Cohort-Level Fine-Mapping Summary (cis_gene Windows)</h2>
  <div class="panel">%s</div>

  <h2>2. CXCL8 Technical Replicates Across 4 Olink Explore Panels (EUR vs ALL)</h2>
  <div class="panel">%s</div>

  <h2>3. Pairwise Cross-Study / Cross-Platform / Cross-Ancestry Colocalisation &amp; PIP Concordance</h2>
  <div class="panel">%s</div>

  <h2>4. Per-Assay cis_gene Fine-Mapping Metrics Across All Studies</h2>
  <div class="panel">%s</div>
</div>
</body>
</html>',
  html_esc(paste(studies_present, collapse = ", ")),
  html_esc(finemap_dir),
  make_html_table(agg_df),
  make_html_table(cxcl8_df),
  make_html_table(pairs_df, max_rows = 1000L),
  make_html_table(summary_df, max_rows = 1000L)
)

html_path <- file.path(outdir, "comparison_report.html")
writeLines(html_out, con = html_path, useBytes = TRUE)
cat(sprintf(
  "[compare_pqtl_studies] Wrote %d assay rows (%d studies) and %d pairwise comparisons -> %s\n",
  nrow(summary_df), length(studies_present), nrow(pairs_df), outdir
))
