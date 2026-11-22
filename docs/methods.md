# Statistical Methods & Edge-Case QC

This page documents the statistical models, coordinate/LD hygiene checks, and numerical safeguards implemented across all stages of **GWASfinemapping**.

---

## 1. Stage 1: Region Screening & Window Construction (`screen_gwas.py`)

`scripts/screen_gwas.py` streams `bcftools view` over the indexed GRCh38 GWAS-VCF, filtering to variants with $\min(\text{AF}, 1 - \text{AF}) \ge \text{min\_maf}$ (default `0.01`) and $-\log_{10} P \ge \text{significance\_lp}$ (default `7.30103`, $P \le 5 \times 10^{-8}$). Traits with zero passing variants finish in ~30 seconds with `overall_status = "NO_SIGNIFICANT_LOCI"` without invoking liftover, S3, or R.

Passing associations are partitioned into three region types:

1. **`cis_gene` (`gene_start - 1 Mb` to `gene_end + 1 Mb`):**
   - For molecular traits (`UNIPROT_OID`, `UNIPROT_SEQID`, with UniProt isoform suffixes such as `O43521-2` stripped automatically, or explicit `--gene-id`), coordinates are resolved against `gene_regions_GRCh38.tsv` (`11,105` unambiguous human genes from `usefuldata::mapping_GRCh38_p14_uniprot`; multi-gene complexes in `gene_regions_GRCh38.ambiguous.tsv` fall back safely to `lead_locus`).
   - Activated when $\max(\text{LP})$ inside the cis window meets `--cis-significance-lp` (or unconditionally when `--cis-significance-lp 0`, e.g. for threshold-free `MultiSuSiE` evaluation).
2. **`mhc` (`chr6:25,500,000–34,000,000` GRCh38):**
   - Following Sun et al. (2023), any significant association inside the extended MHC activates **one merged 8.50 Mb locus** (~8.47 Mb on GRCh37, strictly inside Pan-UKB's 10 Mb stored BlockMatrix band) rather than splitting long-range HLA haplotypes across adjacent windows.
   - Variants carry both `mhc_region` (`25.5–34.0 Mb`) and `mhc_core_grc` (NCBI GRC core MHC `chr6:28,510,120–33,480,577` GRCh38) flags.
3. **`lead_locus` (Kanai greedy `±1.5 Mb` = 3 Mb merged windows):**
   - Remaining genome-wide significant variants outside `mhc` are ranked by descending `LP`, buffered by `±1.5 Mb`, and merged on overlap.
   - Any merged chain wider than `--max-region-width-bp` (`9,500,000` bp) is recursively split at its deepest $-\log_{10} P$ valley between leads so every variant pair stays within Pan-UKB's 10 Mb LD band.

---

## 2. Stage 2: Liftover, Allele Orientation & Pan-UKB LD Extraction (`prep_locus.py`)

Summary-statistic fine-mapping breaks down silently if even a handful of variants suffer reference mismatches, cross-contig liftover artifacts, or allele-flip sign errors. Every locus passes through five sequential guards:

1. **Per-base `pyliftover` + Target-Contig Guard (`lifted_to_different_chromosome`):**
   - Lifts every nucleotide of `REF` from GRCh38 to GRCh37 (complementing reverse-strand chain blocks).
   - Whole-genome sequencing studies (e.g. deCODE 49x WGS) call variants inside pericentromeric segmental duplications whose UCSC `hg38ToHg19` chain target is an unlocalized GRCh37 contig (such as `chr1_gl000192_random` near `REG4`). Any variant lifting to a contig other than the target chromosome is filtered **prior** to `bcftools norm` and tallied in `qc.json`.
2. **GRCh37 FASTA Verification & Indel Left-Alignment:**
   - `bcftools norm --check-ref s -f GRCh37.fa` left-aligns indels and drops any record whose lifted `REF` disagrees with the Ensembl GRCh37 primary assembly.
3. **Pericentromeric `<=9.6 Mb` GRCh37 Band Clamping:**
   - Pericentromeric GRCh38 loci (`1p12` `REG4`, `9p13` `CA9`) contain `NOTCH2NL`/`NBPF` paralogous blocks that lift >29 Mb across the centromere onto the opposite GRCh37 arm, exceeding Pan-UKB's 10 Mb stored BlockMatrix radius. `_clamp_to_grch37_ld_band()` anchors on the lifted target gene / lead variant and clamps any outlier coordinates beyond a `9.6 Mb` GRCh37 span.
4. **Offline PyArrow Index Cut & Sign-Invariant Allele Matching:**
   - Slices the locus row groups from `UKBB.<POP>.variant_index.parquet` in milliseconds offline.
   - Resolves both exact `(REF, ALT)` and swapped `(ALT, REF)` matches (`PANEL_ALLELE_SWAPPED = TRUE`). When swapped, row and column signs of $R$ are flipped in place while `ES` and `EZ` remain oriented to the GWAS `ALT` allele — preserving the likelihood $z^\top R^{-1} z$ identically.
5. **Bounded BlockMatrix Read & Degenerate Row Pruning:**
   - Streams only the required blocks from `s3://pan-ukb-us-east-1/ld_release/UKBB.<POP>.ldadj.bm` with bounded in-flight memory, prunes rare Pan-UKB rows with degenerate zero diagonals or unstored band `NaN`s, symmetrises and unit-normalises $R$ in place, and writes a Fortran-order float64 `matrix.npy`.

---

## 3. Single-Study Prioritisation & Bayesian Fine-Mapping Methods

### 3.1 `p-ld` (`R/prioritize_pld.R`)
Deterministic Selected-Push Sieve ranking variants by descending `LP` (tie-broken by coordinate and `GWAS_VARIANT_ID`) and greedily retaining variants with $r^2 < t$ against all previously selected leads at $t \in \{0.1, 0.01, 0.001\}$. Paired with `--min-maf 0.01`, this eliminates false pseudo-independent lead inflation caused by ultra-rare variants with near-zero sample $r^2$.

### 3.2 `susie-rss` (`R/finemap_susie.R`)
Fits `susieR` 2.0.0 (`L = 10`, `max_iter = 100` genome-wide / `200` `cis_gene`, $n_{\text{eff}} = \min(\text{SS})$, 95% credible sets filtered at exact submatrix purity $\min |r| \ge 0.50$).
- **Zero-copy C scaling:** `gwasfm_frob_and_scale_inplace()` in `src/fast_matrix.c` converts $R \to X^\top X = (n_{\text{eff}} - 1)R$ in place without allocating extra $M \times M$ R matrices.
- **MHC BLAS acceleration:** Once finiteness is verified in C, `options(matprod = "blas")` bypasses R's single-threaded `O(M^2)` `ISNAN` pre-scan on every matrix-vector multiply, accelerating 35,000-variant MHC loci by 3x–9x across 16 OpenBLAS threads.

### 3.3 `susie-rrss` (`R/finemap_susie.R` — Recommended Primary Single-Study Method)
Extends `susie-rss` in `susieR` 2.0.0 with two external-LD protections:
- **Exact finite-sample shrinkage (`R_finite`):** Uses the exact Pan-UKB sample size ($B = 420{,}542$ for `EUR`, $6{,}637$ for `AFR`, $8{,}876$ for `CSA`, $2{,}709$ for `EAS`, $1{,}599$ for `MID`, $982$ for `AMR`).
- **Always-on empirical-Bayes residual mixture (`R_mismatch = "eb_mix"`):** Kept active across **all** locus sizes $M$ so $O(M)$ Maller single-effect initialization and $\lambda_{\text{bias}}$ / `B_corrected` variance inflation protect wide 2 Mb `cis_gene` windows ($M \approx 5{,}000\text{–}20{,}000$), while optional $O(M^3)$ `Q_art` subspace eigendecomposition runs up to `eigen_max_variants = 20000L`.

### 3.4 `carma` (`R/finemap_carma.R`)
Fits shotgun stochastic search Spike-and-Slab Bayesian fine-mapping (`ZikunY/CARMA`, Poisson model-size prior, `L = 10`, deterministic seed) with Bayes-factor LD-discrepancy outlier detection (`outlier.switch = TRUE`, `outlier.BF.index = 1/3.2`), accelerated by a vectorised C/LAPACK `dsyev` batch evaluator (`gwasfm_carma_spike_slab_batch`). Especially resilient as an orthogonal rescue method when extreme $|Z|$ values trigger IBSS non-convergence.

### 3.5 `polyfun` (`R/finemap_polyfun.R`)
Matches autosomal GRCh37 variants against `19,476,620` UK Biobank `baselineLF2.2` meta-analyzed per-SNP heritabilities (`snpvar_meta`, Weissbrod et al. 2020 Approach 1, $q = 100$ max/min locus floor) and fits prior-weighted `SuSiE-RSS` (`prior_weights` $\propto \sigma^2_j$).

> [!NOTE]
> **Choosing `susie-rrss` vs `polyfun` Downstream:** Use `polyfun` when maximising complex-trait variant discovery under functional enrichment; use **uniform-prior `susie-rrss`** whenever testing whether posterior mass sits on a coding protein-altering variant (PAV) vs a non-coding cis-eQTL (e.g. pQTL binding-artifact attribution), because `polyfun` up-weights coding annotations *a priori* by ~10x–50x.

---

## 4. Two-Tier Credible-Set Validity Gate (`converged` + Per-Slot $V_l \le 1.0$)

Under external or mismatched LD ($\Delta r = R^{\text{GWAS}} - R^{\text{LD}}$), a causal variant with large $|Z_{\text{causal}}|$ induces phantom residual scores $Z_{\text{phantom}} = \Delta r \cdot Z_{\text{causal}}$. Standard scalar prior-variance optimisation can drive opposite-sign collinear variants to unphysical prior variances $V_l \gg 1$ (observed up to $V = 398.6$ in multi-ancestry `Sun ALL` against `EUR` LD).

Both `finemap_susie.R` and `run_multisusie.py` enforce a strict **two-tier validity gate**:

1. **Tier 1 (Non-convergence suppression):** If `!isTRUE(fit$converged)`, **all** credible sets at the locus are suppressed (`credible_sets` emptied and `fit$sets$cs <- NULL`), preventing downstream `coloc.susie` from evaluating unconverged IBSS components.
2. **Tier 2 (Surgical per-slot $V_l \le 1.0$ pruning):** On converged fits over unit-variance standardized summary statistics (where a single variant explaining 100% of phenotypic variance has $V \le 1$), any individual component $l$ with $V_l > 1.0$ is surgically removed from `credible_sets`, `fit$sets$cs`, and `fit$sets$cs_index` while preserving valid primary credible sets ($V_l \le 1.0$). All suppressed counts are recorded in `diagnostics$n_cs_dropped_high_V` and `qc.json`.

---

## 5. Opt-In `MultiSuSiE`: Unconditionally SPD 64-Bit Woodbury Kernel (`run_multisusie.py`)

Upstream `MultiSuSiE` (`jordanero/MultiSuSiE`) fits a joint $K$-population Sum of Single Effects model where causal indicator weights $\alpha_l$ are shared across cohorts while effect sizes follow a multivariate normal prior across populations with cross-ancestry correlation matrix $\rho$ and per-population prior variances $V = (V_1, \dots, V_K)$.

During real-world proteogenomic testing on high-$|Z|$ cis-pQTLs (`LP > 500`), upstream `MultiSuSiE` crashed with `LinAlgError: Posterior precision M is not positive definite` because:
1. Near-singular prior covariance $A = \rho \odot (\sqrt{V}\sqrt{V}^\top)$ ($\rho_{kk'} = 0.75$, small $V_k$) cannot be safely inverted directly via $A^{-1} + Q_j$, and
2. Residual variance updates $\hat{\sigma}^2_k = \text{ER2}_k / n_k$ lacked a positive floor when minor LD mismatch drove expected RSS slightly negative ($\text{ER2}_k < 0 \implies Q_{jk} < 0$).

`scripts/run_multisusie.py` replaces both `compute_lbf` kernels in 64-bit floating point with:

1. **Residual variance floor:** $\sigma^2_k \leftarrow \max(\hat{\sigma}^2_k, 10^{-2})$, guaranteeing $q_{jk} = (X^\top X)_{jj,k} / \sigma^2_k \ge 0$.
2. **Unconditionally symmetric positive-definite Woodbury / Cholesky reformulation:** Writing $Q_j = \operatorname{diag}(q_{j1}, \dots, q_{jK})$ and $s_j = \sqrt{q_j}$, the symmetrised matrix
   $$B_j = I_K + Q_j^{1/2} A Q_j^{1/2}$$
   satisfies $\lambda_{\min}(B_j) \ge 1$ for **every** PSD prior $A \succeq 0$ (including rank-deficient $\rho_{kk'} = 1$ or boundary $V_k = 0$), yielding exact log-Bayes factors and posterior moments without ever inverting $A$:
   $$\operatorname{LBF}_j = \frac{1}{2} z_j^\top \left(A z_j - A Q_j^{1/2} B_j^{-1} Q_j^{1/2} A z_j\right) - \frac{1}{2} \log \det B_j, \qquad \Sigma_{\text{post},j} = A - A Q_j^{1/2} B_j^{-1} Q_j^{1/2} A.$$
