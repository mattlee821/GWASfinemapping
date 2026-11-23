# Protein `cis-pQTL` Guide: Why Active Kriging Filtering Is Off by Default

Proteogenomic fine-mapping (`cis-pQTLs` on Olink Explore and SomaScan v4/v4.1 across cohorts such as UKB-PPP, Fenland, ARIC, and deCODE) operates in a **qualitatively different statistical regime** from polygenic complex-trait GWAS. Heuristic pre-filtering rules that are popular in complex-trait pipelines — specifically **active leave-one-out kriging outlier deletion** (`kriging_rss` `logLR > 2`) and **null-MLE LD shrinkage** (`lambda = s_hat`) — systematically damage well-powered `cis-pQTL` loci.

This guide explains **why active kriging removal is disabled by default (`--remove-kriging-outliers false`)**, walks through real protein failure modes from our [28-Protein Cross-Platform & MultiSuSiE Benchmark](sun_eur_vs_all_comparison_report.md), and shows concrete CLI and R workflows for single-study and multi-cohort protein fine-mapping.

---

## 1. Complex-Trait GWAS vs Protein `cis-pQTL` Regimes

| Property | Typical Complex-Trait GWAS Locus | Typical `cis-pQTL` Window (`gene ± 1 Mb`) |
|---|---|---|
| **Peak $|Z|$ / $-\log_{10} P$** | $|Z| \approx 5.5\text{–}15$ (`LP = 7.3–50`) | $|Z| = 20\text{–}109$ (**`LP = 100–2,595`**, e.g. `CCL25`, `MMP12`, `TNFRSF10A`, `GDF15`) |
| **Causal Architecture per Locus** | Usually 1–2 causal variants | **Extreme allelic heterogeneity:** `64–95%` of active `cis` windows harbour **2–10 independent signals** (regulatory eQTLs + coding PAVs + sQTLs) |
| **Coding & Structural Variants** | Mostly common non-coding enhancers well tagged by reference LD | Frequent **missense / PTV protein-altering variants (PAVs)**, low-frequency coding alleles, and tandem paralogue / CNV haplotypes (`FCGR2A`, `CFHR5`, `CCL3`) |
| **Allele Harmonisation Quality** | Often heterogeneous fixed-effect meta-analyses or legacy rsID tables (`Case X`) | Standardised per-cohort GRCh38 GWAS-VCFs left-aligned against reference FASTA (`Case Y`) |

---

## 2. Three Ways Active Kriging / Null-MLE Shrinkage Breaks `cis-pQTLs`

### Failure Mode 1: Opposite-Sign Secondary `cis` Signals Trigger False `logLR > 2` Flips (`MMP12` & `GDF15`)

Leave-one-out kriging (`susieR::kriging_rss`) is a **single-variant conditional null test**: it predicts variant $j$'s marginal score $\hat{\mu}_j = \mathbb{E}[z_{\text{adj},j} \mid z_{\text{adj},-j}]$ assuming no causal effect sits on variant $j$ itself.

Consider a multi-causal `cis-pQTL` window such as **`MMP12` (`P39900`, 6 resolved 95% CSs)** or **`GDF15` (`Q99988`, 6 resolved 95% CSs)**:

- Suppose primary causal variant $1$ (e.g. a high-impact coding PAV or core promoter SNP) has $z_{\text{adj},1} = +55$.
- Suppose a genuine, independent secondary causal variant $2$ (e.g. a distal enhancer eQTL that lowers transcription on the haplotype carrying the minor allele) has true marginal association $z_{\text{adj},2} = -6.5$, and happens to sit in modest positive LD ($r_{12} = +0.25$) with variant $1$ and its tags.
- Before conditioning on causal effects, leave-one-out kriging predicts $\hat{\mu}_2 \approx +13.8$ from the massive $+55$ lead haplotype, whereas the observed score is $z_{\text{adj},2} = -6.5$.
- Because $\hat{\mu}_2$ and $z_{\text{adj},2}$ have **opposite signs** and large magnitudes, `KRIGING_LOGLR` easily exceeds `+2` with $|z_{\text{adj},2}| > 2$ — **even when allele orientation is 100% exact and LD is in-sample!**

If `--remove-kriging-outliers true` were enabled by default, variant $2$ would be deleted **before** `SuSiE-rRSS` ever ran — destroying the clean non-PAV secondary `cis` instrument needed for downstream Mendelian randomisation (`A1` vs coding `A2`/`C` separation). When kept in the model (`--remove-kriging-outliers false`), `SuSiE-rRSS` jointly fits variant $1$ in slot `L1`, subtracts $+55 \cdot r_{12}$ from the residual, and cleanly recovers variant $2$ as a single-SNP 95% credible set!

### Failure Mode 2: Deleting Poorly Tagged Coding PAVs & Founder-Drifted Leads (`DENTIST` Trap)

Coding missense variants, indels, and variants in segmental duplications or founder-drifted haplotypes (such as Icelandic `deCODE` `Ferkingstad` mapped against UK `Pan-UKB EUR`) often have weaker $r^2$ to flanking common array/reference tags than expected for their true causal effect:

- Their conditional expectation from neighbours is near zero ($|\hat{\mu}_j| \ll |z_{\text{adj},j}|$), driving the standardized residual $|t_j| = |z_{\text{adj},j} - \hat{\mu}_j| / \sqrt{\text{condvar}_j} \gg 5.45$ ($P_{\text{DENTIST}} \ll 5 \times 10^{-8}$).
- Dropping `DENTIST` / kriging outliers on `cis-pQTLs` selectively strips the **exact causal missense PAV** (creating a false "non-PAV" credible set on a flanking partial tag and risking misclassification of an epitope-binding artifact as a true abundance QTL!).

### Failure Mode 3: Null-MLE `s_hat` Inflation Spawns False Credible Sets Under Legacy `lambda = s_hat`

`estimate_s_rss(method = "null-mle")` estimates $\hat{s}$ under $z_{\text{adj}} \sim \mathcal{N}(0, (1 - s)R + sI)$ with **zero causal effects**. On `Sun EUR` UK Biobank Olink assays where LD is **in-sample** (`Pan-UKB EUR`, $B = 420{,}542$), extreme true `cis` signals (`LP = 300–2,595`) routinely drive pre-fit null-MLE $\hat{s}$ to `0.15–0.65`.

- **Never filter loci on `s_hat > 0.1` alone**, and **never pass `lambda = s_hat` into `susie_rss`**: uniform shrinkage $(1 - \hat{s})R + \hat{s}I$ weakens off-diagonal LD between the true `LP = 1,000` lead and its $r = 0.95$ tags, causing IBSS to allocate spurious extra credible sets to residual tag leakage.
- Instead, **GWASfinemapping** uses `susieR` 2.0 `SuSiE-rRSS` (`R_finite` + always-on `eb_mix`), which estimates `lambda_bias` **after** removing active sparse effects and gates unresolvable slots via Two-Tier `converged + V_l <= 1.0` pruning.

---

## 3. Recommended Protein `cis-pQTL` Workflows

### 3.1 Single-Cohort Olink / SomaScan Fine-Mapping (Default Mode)

Run `scripts/finemap.sh` with **default** `--remove-kriging-outliers false` and ancestry-matched Pan-UKB LD (`EUR` for `ukbppp Sun EUR`, `Fenland`, and `ARIC EUR`; `AFR` for `ARIC AFR`):

```bash
# 1. UKB-PPP Olink Explore (always use Sun EUR, never multi-ancestry Sun ALL against EUR LD!)
sbatch scripts/finemap.sh \
  --input-dir /path/to/GWASprocessing/GWAS/sun_2023_37794186_EUR \
  --ld-pop EUR \
  --max-forks 16 \
  --slurm

# 2. ARIC African-ancestry SomaScan v4 matched to Pan-UKB AFR (B = 6,637)
sbatch scripts/finemap.sh \
  --input-dir /path/to/GWASprocessing/GWAS/zhang_2022_35501419_AFR \
  --ld-pop AFR \
  --max-forks 8 \
  --slurm
```

> [!WARNING]
> **Avoid Multi-Ancestry `Sun ALL` Against Single-Population LD:** In our [28-Protein Benchmark](sun_eur_vs_all_comparison_report.md), running multi-ancestry `Sun ALL` ($N = 52.4\text{k}$) against `Pan-UKB EUR` collapsed standard `SuSiE-RSS` from `87.1%` to `43.3%` and inflated prior variances up to $V = 398.6$ due to $Z_{\text{phantom}} = \Delta r \cdot Z_{\text{causal}}$, whereas ancestry-matched `Sun EUR` ($N = 34.6\text{k}$) achieved **`96.8%` `SuSiE-rRSS` convergence (`130` CSs, `53` 1-SNP CSs)**.

### 3.2 Within-Platform Multi-Cohort / Cross-Ancestry `MultiSuSiE` (Opt-In)

For multi-cohort SomaScan (or multi-cohort Olink sharing a target gene), combine **non-isolate cohorts within the same platform** via opt-in `MultiSuSiE`. Every cohort is automatically profiled with per-cohort kriging diagnostics (`kriging_by_cohort`) without dropping variants:

```bash
# Jointly fine-map SomaScan v4 aptamers across Fenland EUR + ARIC EUR + ARIC AFR
sbatch scripts/finemap.sh \
  --multisusie-studies /path/to/pietzner_2021_34648354,/path/to/zhang_2022_35501419_EUR,/path/to/zhang_2022_35501419_AFR \
  --group-study-id somascan_multisusie_3cohort \
  --cis-significance-lp 0 \
  --max-forks 8 \
  --slurm
```

- **Why keep `deCODE` (`Ferkingstad`) as a separate single-study run:** Icelandic founder LD drift against UK `Pan-UKB EUR` is safely handled on single-study `deCODE` via `SuSiE-rRSS` (`78.6%` converged) + `CARMA` (`100%` converged) + explicit `!converged` flags, whereas pooling `deCODE` ($N = 35.6\text{k}$) inside 4-cohort `MultiSuSiE` imports founder-LD mismatch into otherwise clean `Fenland + ARIC` high-$|Z|$ loci (`MMP12`, `GDF15`, `REG4`, `TFF1`).
- **Using `MultiSuSiE` per-cohort `s_hat` & `kriging_cs_outlier_count` to screen multi-country cohorts:** Before pooling a multi-country or heterogeneous cohort into `MultiSuSiE`, inspect `kriging_by_cohort` in `qc.json` or `report.html` across your pilot proteins to verify its `s_hat` and `kriging_cs_outlier_count` track UK/US `EUR` cohorts rather than founder-drifted `deCODE`.

---

## 4. Inspecting Non-Destructive Diagnostics & Gating `ld_qc_fail` in R

For downstream colocalisation, PAV/epitope attribution, and cis-MR, load **uniform-prior `susie-rrss.rda`** (never `polyfun.rda`, whose functional coding priors up-weight missense PAVs *a priori*) and apply a calibrated post-fit QC check:

```r
e <- new.env(parent = emptyenv())
load("finemapping/ferkingstad_2021_34857953/P39900_4496-60/methods/susie-rrss.rda", envir = e)

reg  <- e$susie_rrss$regions[[1]]$result
diag <- reg$diagnostics
cs   <- reg$credible_sets
vars <- reg$variants

# 1. Calibrated locus-level LD / solver QC gate (do NOT gate on s_hat > 0.1 alone!)
ld_qc_fail <- !isTRUE(diag$converged) ||
  isTRUE(diag$credible_sets_suppressed) ||
  isTRUE(diag$R_reliability_flag) ||
  isTRUE(diag$artifact_flag) ||
  (!is.na(diag$kriging_cs_outlier_count) && diag$kriging_cs_outlier_count > 0L)

cat(sprintf(
  "converged=%s | active_CS=%d | s_hat=%.3f | lambda_bias=%.4f | krig_CS_out=%s | ld_qc_fail=%s\n",
  diag$converged, nrow(cs), diag$s_hat, diag$lambda_bias,
  as.character(diag$kriging_cs_outlier_count), ld_qc_fail
))

# 2. Inspect per-credible-set purity, prior variance V, and in-CS kriging outliers
print(cs[, c("credible_set_id", "size", "min_abs_corr", "prior_variance_V",
             "lbf", "kriging_outlier_count", "ld_mismatch_warning")])

# 3. If kriging_outlier_count > 0 at a multi-causal cis-pQTL, check whether flagged SNPs
#    are opposite-sign secondary signals vs true isolated allele switches:
subset(vars, KRIGING_OUTLIER,
       select = c("GWAS_VARIANT_ID", "EZ", "KRIGING_EXP_Z", "KRIGING_LOGLR",
                  "PIP", "credible_set_id"))
```

### Summary of Downstream Actions for Protein Loci

1. **`converged == TRUE`, `kriging_cs_outlier_count == 0`, `!R_reliability_flag` (even if `s_hat` is `0.2–0.6`):** High-$|Z|$ `cis-pQTL` converged cleanly under `eb_mix`. Proceed directly with `reg$fit` in `coloc::coloc.susie()` and per-CS PAV annotation.
2. **`kriging_cs_outlier_count > 0` inside an active 95% CS:** Open `report.html` Panel 3 (`Kriging Expected Z vs Observed Z`) and cross-check `methods/carma.rda` (`OUTLIER` column). If the flagged variant is a secondary opposite-sign `cis` signal (`PIP > 0.5` in its own CS, concordant across cohorts), retain it; if it is an isolated sign-flipped SNP in a single external summary-stat file, flag that cohort's CS as `ld_mismatch_warning` (`ld_qc_fail = TRUE`).
3. **`!converged` or `SER_fallback == TRUE` (e.g. extreme founder drift on `deCODE` `CCL25` / `TNFSF12`):** Two-Tier gating has already suppressed unreliable `SuSiE` CSs (`nrow(cs) == 0`). Route the single-cohort `SuSiE` comparison to `ld_qc_fail = TRUE` (unclassified / external-LD limited) or inspect `methods/carma.rda` for spike-and-slab lead recovery.
