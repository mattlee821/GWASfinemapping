# Outputs, Colocalisation & Interactive R

**GWASfinemapping** is engineered around a **storage-light, `coloc.susie`-ready output contract** (~25–40 MB per active GWAS with `xz` compression, and <50 KB for traits with zero genome-wide significant hits) that requires **no stored LD matrices** for downstream multi-signal colocalisation or Mendelian randomisation.

---

## 1. Output Directory Hierarchy

```text
<outdir>/
├── manifest.html                       # Batch dashboard across all single-study & MultiSuSiE runs
├── .finemap_state.tsv                  # Idempotent completion ledger keyed by settings_hash
├── <STUDY_ID>/                         # Standard single-study output (e.g. sun_2023_37794186_EUR)
│   └── <GWAS_ID>/
│       ├── metadata.json               # Upstream study metadata + per-method 95% CS totals
│       ├── provenance.json             # Git SHA, CLI args, chain SHA-256, Pan-UKB N, runtime
│       ├── qc.json                     # Per-locus liftover/LD attrition, flags & step_telemetry
│       ├── report.html                 # Self-contained HTML report (tables, QC badges, base64 plots)
│       ├── logs/                       # Per-stage logs, Nextflow trace.tsv, slurm-*.out
│       └── methods/
│           ├── p-ld.rda                # Object `p_ld`: r² < 0.1, 0.01, 0.001 heuristic sets
│           ├── susie-rss.rda           # Object `susie_rss`: PIPs, gated 95% CSs, coloc `fit`
│           ├── susie-rrss.rda          # Object `susie_rrss`: eb_mix robust PIPs, CSs, coloc `fit`
│           ├── carma.rda               # Object `carma`: shotgun PIPs, CSs, OUTLIER flags
│           └── polyfun.rda             # Object `polyfun`: SNPVAR-weighted PIPs, CSs, coloc `fit`
└── <GROUP_STUDY_ID>/                   # Opt-in MultiSuSiE output (e.g. somascan_multisusie_3cohort)
    └── <GROUP_ID>/
        ├── metadata.json               # Group metadata, K cohorts, combined N, CS count
        ├── provenance.json             # Per-cohort VCF paths, LD panels, solver settings
        ├── qc.json                     # Per-cohort variant overlap, rho matrix, max_V, gate status
        ├── report.html                 # Multi-cohort convergence, CS lead table & cohort Z summary
        └── methods/
            └── multisusie.rda          # Object `multisusie` (plus `susie_rrss` alias): coloc `fit`
```

---

## 2. R Container Schema (`methods/*.rda` / `*.rds`)

Each Bayesian container (`susie_rss`, `susie_rrss`, `carma`, `polyfun`, `multisusie`) stores:

- **`$method_id`:** `"susie_rss"`, `"susie_rrss"`, `"carma"`, `"polyfun"`, or `"multisusie"`.
- **`$regions[[region_id]]`:** Named list keyed by region identifier (e.g. `"cis_MMP12_chr11_101862586_103874982"`), containing `region_id`, `region_type` (`"cis_gene"`, `"lead_locus"`, or `"mhc"`), and `result`:
  - **`result$variants`:** Data frame with canonical GRCh38 `GWAS_VARIANT_ID` (`chr:pos:ref:alt`), `GWAS_CHROM`, `GWAS_POS`, `GWAS_REF`, `GWAS_ALT`, `ID` (rsID), `LP`, `ES`, `SE`, `SS`, `AF`, `EZ`, lifted GRCh37 `LD_CHROM`, `LD_POS`, `LD_REF`, `LD_ALT`, `PANEL_IDX`, `PANEL_ALLELE_SWAPPED`, `LD_R2_TO_LEAD`, `cis_trans`, `mhc_region`, `mhc_core_grc`, **`PIP`**, and **`credible_set_id`** (plus `OUTLIER` in `carma`; `SNPVAR`, `PRIOR_WEIGHT`, `PRIOR_MATCHED` in `polyfun`; per-cohort `Z_<cohort>` and `AF_<cohort>` in `multisusie`).
  - **`result$credible_sets`:** One row per purity-passing ($\min |r| \ge 0.50$), validity-gated ($V_l \le 1.0$) 95% credible set (`credible_set_id`, `effect_index`, `size`, `coverage`, `min_abs_corr`, `mean_abs_corr`, `median_abs_corr`).
  - **`result$credible_set_membership`:** Long-format mapping of `credible_set_id` to `GWAS_VARIANT_ID` and variant `PIP`.
  - **`result$diagnostics`:** Convergence flag, IBSS iteration count, ELBO, `max_V`, `n_cs_dropped_high_V`, `eb_mix` (`lambda_bias`, `R_reliability_flag`, `artifact_flag`), or `MultiSuSiE` cross-cohort effect-correlation matrix `rho`.
  - **`result$fit`** (`susie_rss`, `susie_rrss`, `polyfun`, `multisusie`): Compact S3 `"susie"` object ($O(L \times M)$, no $M \times M$ LD matrix stored) with `colnames(fit$alpha)`, `colnames(fit$lbf_variable)`, `colnames(fit$mu)`, `colnames(fit$mu2)`, and `names(fit$pip)` all keyed to GRCh38 `GWAS_VARIANT_ID`.

---

## 3. Multi-Signal Colocalisation in R (`coloc.susie` & `coloc.bf_bf`)

### Direct `coloc::coloc.susie()` Between Single-Study or `MultiSuSiE` Fits

Because both single-study (`susie-rrss.rda`) and multi-cohort (`multisusie.rda`) outputs share identical GRCh38 `chr:pos:ref:alt` variant naming and validity-gated `fit$sets$cs_index`, `coloc::coloc.susie()` works out of the box across platforms, ancestries, and diseases without loading any LD matrix:

```r
library(coloc)

e_olink <- new.env(parent = emptyenv())
e_soma  <- new.env(parent = emptyenv())

load("finemapping/sun_2023_37794186_EUR/Q03403_OID20787/methods/susie-rrss.rda", envir = e_olink)
load("finemapping/somascan_multisusie_3cohort/Q03403_9191-8/methods/multisusie.rda", envir = e_soma)

fit_olink <- e_olink$susie_rrss$regions[[1]]$result$fit
fit_soma  <- e_soma$multisusie$regions[[1]]$result$fit

# Pairwise multi-signal colocalisation across all active credible-set slots
res <- coloc::coloc.susie(fit_olink, fit_soma)
print(res$summary[, c("idx1", "idx2", "nsnps", "PP.H3.abf", "PP.H4.abf", "hit1", "hit2")])
```

### Hybrid SuSiE × LD-Free Wakefield ABF (`coloc.bf_bf`) for Multi-Ancestry Disease Outcomes

When colocalising a well-powered, fine-mapped molecular trait (`susie-rrss` or `multisusie` log-BF rows `fit$lbf_variable[fit$sets$cs_index, ]`) against a modest-$|Z|$ multi-ancestry case-control outcome GWAS where an exact matched LD panel is unavailable (Wallace 2021, *PLOS Genetics*), compute single-variant Wakefield approximate Bayes factors directly from `(ES, SE)` on the outcome and pass both log-BF matrices to `coloc::coloc.bf_bf()`:

```r
# Extract fine-mapped pQTL per-CS log Bayes factors (L_active x M)
cs_idx <- fit_olink$sets$cs_index
lbf_pqtl <- fit_olink$lbf_variable[cs_idx, , drop = FALSE]

# Compute LD-free Wakefield log-ABF on outcome summary stats (prior sd W = 0.20 for log-OR)
W2 <- 0.20^2
v  <- crc_df$SE^2
r  <- W2 / (v + W2)
z2 <- (crc_df$ES / crc_df$SE)^2
lbf_crc <- matrix(0.5 * (log(1 - r) + r * z2), nrow = 1, dimnames = list("crc", crc_df$GWAS_VARIANT_ID))

coloc_hybrid <- coloc::coloc.bf_bf(lbf_pqtl, lbf_crc)
```

### Batch Cross-Study & Cross-Platform Evaluation (`compare_pqtl_studies.R`)

To compare multiple completed study/model directories in `<outdir>` (producing `00_study_level_summary.tsv`, `01_cis_method_summary.tsv`, `02_pairwise_coloc_and_concordance.tsv`, and `03_gene_platform_matrix.tsv` with exact IEEE-754 overflow-safe `coloc` evaluation even when $\text{LBF}_1 + \text{LBF}_2 > 709.8$):

```bash
source scripts/source.sh
Rscript scripts/compare_pqtl_studies.R \
  --finemap-dir /scratch/user/finemapping \
  --outdir /scratch/user/comparison_report \
  --studies sun_2023_37794186_EUR,pietzner_2021_34648354,somascan_multisusie_3cohort
```

---

## 4. Interactive Line-by-Line Debugging in RStudio or Positron

Every method script in `GWASfinemapping` is structured so you can step through Stage 1 screening, Stage 2 liftover & Pan-UKB LD extraction, and Stage 3 model fitting interactively in RStudio or Positron without Nextflow:

1. Open `scripts/interactive_prep_locus.R` to prepare a single locus (with `--keep-intermediates` layout) and inspect `matched_variants.tsv`, `prep_qc.json`, and `matrix.npy`.
2. Open any standalone method runner and execute line-by-line against that prepared locus directory:
   - `scripts/methods/run_p_ld.R`
   - `scripts/methods/run_susie_rss.R`
   - `scripts/methods/run_susie_rrss.R`
   - `scripts/methods/run_carma.R`
   - `scripts/methods/run_polyfun.R`
