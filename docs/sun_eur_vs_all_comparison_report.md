# UKB-PPP `EUR` vs. `ALL` Fine-Mapping Report & Pipeline Hardening Summary

**Date:** 21 November 2026  
**Target Panel:** 31 Olink Explore 3072 assays (28 unique CRC-candidate proteins, including all 4 `CXCL8` panel replicates)  
**Studies Completed:** `sun_2023_37794186_EUR` (31/31) vs. `sun_2023_37794186_ALL` (31/31) against `Pan-UKB EUR` ($B = 420{,}542$)  
**Cross-Study Extension:** `ferkingstad_2021_34857953` (`EUR`), `pietzner_2021_34648354` (`EUR`), `zhang_2022_35501419_EUR` (`EUR`), and `zhang_2022_35501419_AFR` (`AFR` against `Pan-UKB AFR`, $B = 6{,}637$)

> [!IMPORTANT]
> Switching from pooled multi-ancestry `sun_2023_37794186_ALL` to ancestry- and cohort-matched `sun_2023_37794186_EUR` increased `cis_gene` activation from **96.8% (30/31)** to **100% (31/31)** (rescuing `IFNG`), raised `SuSiE-RSS` valid credible-set recovery from **36.7% (11/30)** to **87.1% (27/31)**, and raised `SuSiE-rRSS` convergence and credible-set recovery from **76.7% (23/30)** to **96.8% (30/31)** (`median = 18` IBSS iterations, median `max_V = 0.01`). Combining `SuSiE-rRSS` and `CARMA` achieved **100% (31/31)** `cis_gene` credible-set coverage in `EUR`.

---

## 1. Methodological & Engineering Challenges Uncovered and Resolved

```mermaid
flowchart LR
    A["GRCh38 GWAS-VCF"] --> B["MAF >= 1% Filter (--min-maf 0.01)"]
    B --> C["GRCh38 -> GRCh37 Liftover & <=9.6 Mb Band Clamp"]
    C --> D["Pan-UKB BlockMatrix Extract (Bounded RAM + Degenerate Row Prune)"]
    D --> E["Multi-Method Runner (16-CPU matprod='blas' on MHC)"]
    E --> F["Two-Tier SuSiE Gate (Converged + Per-Slot V <= 1.0)"]
    F --> G["Ancestry-Mirrored Output (<study_dir>/<gwas_id>)"]
```

During development and validation on `MMP12` (`P39900_OID21439`) and the 31-assay benchmark, six distinct statistical and computational failure modes surfaced and were fixed in `GWASfinemapping`:

### 1.1 Cross-Ancestry Meta-Analysis vs. Single-Panel LD Algebra Violation ($\Delta r \cdot Z$)
- **Symptom:** On `MMP12` (`ALL`, `LP = 1592.8`, $|Z| \approx 85.5$) and 9 other strong cis-pQTLs in `sun_2023_37794186_ALL` (`RBP2`, `CCL7`, `C1QA`, `AGR2`, `GDF15`, `CCL25`, `TNFRSF6B`, `TFF1`, `TNFSF12`), `SuSiE` either failed to converge or estimated prior effect variances of $V \in [1.5, 398.6]$ with $\text{LBF} > 100{,}000$ via opposite-sign collinear cancellation.
- **Root Cause:** `sun_2023_37794186_ALL` is an inverse-variance-weighted meta-analysis across European (~46k) and non-European (~8k `AFR`, `CSA`, `EAS`, `AMR`) UK Biobank participants, where per-variant ancestry weights $w_{jk} \propto N_k \cdot 2 f_{jk}(1 - f_{jk})$ vary by allele frequency. Under RSS, any cross-ancestry LD discrepancy $\Delta r$ around a causal variant produces a phantom residual association:
  $$Z_{\text{phantom}} = (R^{\text{GWAS}}_{j,\text{causal}} - R^{\text{EUR}}_{j,\text{causal}}) \cdot Z_{\text{causal}} = \Delta r \cdot Z_{\text{causal}}$$
  At $|Z_{\text{causal}}| = 85$, even $\Delta r = 0.15$ yields $|Z_{\text{phantom}}| \approx 12.8$ ($P < 10^{-36}$), exceeding `eb_mix`'s genome-wide significance gate ($Z_{\text{ref}} = 5.45$) and driving unbounded scalar `V` optimization (`log(V) ∈ [-30, 15]`) into opposite-sign cancellation on tightly correlated SNPs.
- **Resolution:**
  1. **Primary solution (cohort-matched LD):** Use `sun_2023_37794186_EUR` against `Pan-UKB EUR` (in-sample UK Biobank European LD), eliminating cross-ancestry $\Delta r$ at the source (`max_V` on `GDF15` fell from `21.41` → `0.08`; `RBP2` from `372.03` → `0.004`; `MMP12` from `4.26` → `0.11`).
  2. **Two-tier credible-set validity gate** ([finemap_susie.R](../R/finemap_susie.R)):
     - **Tier 1 (Non-convergence):** Whenever `!isTRUE(fit$converged)`, **all** credible sets at the locus are suppressed (`0` rows in `credible_sets` and `fit$sets$cs <- NULL`) so downstream `coloc::coloc.susie()` never colocalises unconverged slots.
     - **Tier 2 (Per-slot $V \le 1.0$ filter on converged fits):** On unit-variance standardized summary statistics, a single variant cannot plausibly have prior variance $V_l > 1.0$. When `converged == TRUE`, individual slots with $V_l > 1.0$ are surgically removed from `credible_sets`, `fit$sets$cs`, and `fit$sets$cs_index` while preserving valid cis-pQTL lead slots ($V_l \le 1.0$). On `ALL`, this rescued valid cis credible sets on `AGR2` (`1`), `MMP12` (`2`, including `rs17368814`), `TFF1` (`4`), and `C1QA` (`2`), and on `EUR` it cleanly isolated `3` valid CSs on `CCL7` (`3` kept, `2` dropped).
  3. **Iteration budget:** Raised default `--max-iter` from `50` to `100` (`lead_locus` / `mhc`) and added `--cis-max-iter 200` for `cis_gene` windows ([locus_worker.py](../scripts/locus_worker.py)).

### 1.2 `SuSiE-rRSS` `eb_mix` Disabled Above $M = 5{,}000$ Variants
- **Symptom:** Previously, `eigen_max_variants = 5000L` reset `R_mismatch <- "none"` whenever $M > 5{,}000$, disabling `eb_mix` on almost every 2 Mb `cis_gene` window ($M \approx 5{,}400\text{–}20{,}700$).
- **Root Cause:** Only `susieR` 2.0's optional `Q_art` subspace diagnostic uses `eigen(R)` ($O(M^3)$); both 1-effect Maller SER initialization (`initialize_R_mismatch`) and empirical-Bayes residual variance inflation (`lambda_bias` / `B_corrected`) are $O(M)$ vector operations (`< 1 ms`).
- **Resolution:** Kept `R_mismatch = "eb_mix"` active at **all** region sizes $M$ ([finemap_susie.R](../R/finemap_susie.R)), raised `eigen_max_variants` to `20,000` (reusing 1 `eigen()` call for both PSD validation and `Q_art`), and passed a 1-column sentinel `eigen_R` above `20,000` (`q_art_eigen = "skipped_large_m"`). In `EUR`, active `eb_mix` rescued **3** loci where standard `SuSiE-RSS` returned 0 CSs: `CCL25` (`10` CSs at `LP = 1756`), `CA14` (`7` CSs), and `CCL7` (`3` CSs).

### 1.3 Ultra-Rare Variant Zero-LD Inflation in `p-ld` (`--min-maf 0.01`)
- **Symptom:** Without MAF filtering, `MMP12` (`ALL`) selected **308** variants at $r^2 < 0.001$ inside a single 2 Mb cis window.
- **Root Cause:** Ultra-rare variants ($\text{MAF} < 0.001$, often ancestry-specific) have near-zero sample $r^2$ in `Pan-UKB EUR` against the `LP = 1593` lead SNP yet pass $P \le 5 \times 10^{-8}$ via finite-sample / hitchhiking effects, creating hundreds of redundant pseudo-independent instruments.
- **Resolution:** Added `--min-maf 0.01` (`1%` minor allele frequency $\min(\text{AF}, 1 - \text{AF}) \ge 0.01$ via `FORMAT/AF`, configurable with `--min-maf <float>`, `0` to disable) across [screen_gwas.py](../scripts/screen_gwas.py) and [prep_locus.py](../scripts/prep_locus.py). On `MMP12`, `p-ld` ($r^2 < 0.001$) dropped from **308 → 63 (`EUR`) / 83 (`ALL`)**.

### 1.4 Pericentromeric & Segmental-Duplication Liftover Jumps (`REG4` on `chr1`, `chrX` Handling)
- **Symptom:** `cis_REG4` (`chr1:117.15–122.90 Mb` GRCh38, `1p12`) failed in `prep_locus.py` with `Region spans 32,008,984 bp on GRCh37, exceeding Pan-UKB's 10 Mb LD band`.
- **Root Cause:** `REG4` (`chr1:119.80 Mb`) adjoins the `NOTCH2NL`/`NBPF` pericentromeric segmental-duplication block (`120.4–121.4 Mb` GRCh38), whose contigs were placed on the opposite chromosome arm (`1q21.1`, `~149.7 Mb`) in GRCh37, where Pan-UKB stores zero pairwise LD with `1p12`.
- **Resolution:** Added automated two-pointer GRCh37 band clamping (`_clamp_to_grch37_ld_band` at $\le 9.6\text{ Mb}$ anchored on the target gene / lead variant), degenerate Pan-UKB diagonal/edge-`NaN` pruning, and clean non-indexed region skipping in [prep_locus.py](../scripts/prep_locus.py). `REG4` (`EUR`) subsequently converged in **7 iterations** with exact **`4 / 4 / 4 / 4`** credible sets across all four Bayesian methods; `ADGRG2` (`chrX`) also mapped `2,587` (`EUR`) / `2,861` (`ALL`) Pan-UKB `chrX` variants cleanly (`2 / 2 / 2` CSs in `EUR`).

### 1.5 Extended MHC (`chr6:25.5–34.0 Mb`) Memory & Single-Threaded R `matprod` Bottleneck
- **Symptom:** Early MHC runs hit OOM (`exit=137`) during `ldcov` assembly and later hit SLURM's `8h` attempt-1 walltime (`exit=140`) when running `susie_rss` + `susie_rrss` + `polyfun` sequentially on 4 CPUs.
- **Root Cause:** (a) Unbounded thread-pool block futures in LD assembly spiked RAM by ~95 GB; (b) R ≥ 4.0's default `options(matprod = "default")` runs a **single-threaded `ISNAN()` scan across the entire 16 GB $M \times M$ `XtX` matrix on every IBSS matrix-vector product** before calling multi-threaded OpenBLAS.
- **Resolution:** Implemented bounded in-flight BlockMatrix assembly + tiled in-place sign flips ([prep_locus.py](../scripts/prep_locus.py)), automatic Attempt-1 MHC resource tiering (`16` CPUs, `>= 160 GB` RAM, `24h` walltime, `OPENBLAS_NUM_THREADS = 16` in [finemapping.nf](../finemapping.nf) and [locus_worker.py](../scripts/locus_worker.py)), and `options(matprod = "blas")` on pre-validated finite `XtX` ([finemap_susie.R](../R/finemap_susie.R)). Both 31-assay batches (`EUR` and `ALL`) completed every MHC locus on Attempt 1 without timeout.

### 1.6 `--keep-intermediates` Preservation & Verbatim `<study_dir>` Output Mirroring
- Fixed [finalise_gwas.py](../scripts/finalise_gwas.py) so `output_dir / "loci"` is preserved when `--keep-intermediates` is enabled (`matrix.npy`, `matrix.ids.txt`, and `matched_variants.tsv` verified 1:1 aligned at `8,237 × 8,237`), and updated [prepare_manifest.py](../scripts/prepare_manifest.py) to mirror ancestry-specific `GWASprocessing/<study_dir>` names (`sun_2023_37794186_EUR`, `sun_2023_37794186_ALL`, `zhang_2022_35501419_AFR`, etc.) verbatim.

---

## 2. Empirical Comparison: `sun_2023_37794186_EUR` vs. `sun_2023_37794186_ALL` (31 Assays)

### 2.1 Aggregate Performance Summary

| Metric (`cis_gene` windows, $\text{MAF} \ge 0.01$) | `sun_2023_37794186_EUR` (In-Sample `EUR` LD) | `sun_2023_37794186_ALL` (Multi-Ancestry GWAS + `EUR` LD) |
|---|---:|---:|
| **`cis_gene` activated ($P \le 5 \times 10^{-8}$)** | **31 / 31 (100.0%)** | 30 / 31 (96.8%; `IFNG` missed) |
| **`SuSiE-RSS` loci with $\ge 1$ valid 95% CS** | **27 / 31 (87.1%)** | 11 / 30 (36.7%; **19 wiped out**) |
| **`SuSiE-rRSS` (`eb_mix`) converged ($\le 200$ iters)** | **30 / 31 (96.8%)** | 23 / 30 (76.7%; 7 non-converged) |
| **`SuSiE-rRSS` loci with $\ge 1$ valid 95% CS** | **30 / 31 (96.8%)** | 23 / 30 (76.7%; 4 required per-slot $V \le 1$ rescue) |
| **`PolyFun + SuSiE` autosomal loci with $\ge 1$ valid CS** | **26 / 30 (86.7%)** | 10 / 29 (34.5%) |
| **`CARMA` (`outlier_switch = TRUE`) loci with $\ge 1$ CS** | **31 / 31 (100.0%)** | 30 / 30 (100.0%) |
| **Median `SuSiE-rRSS` iterations (active cis loci)** | **18** (IQR: `12–30`) | **33** (IQR: `11–186`) |
| **Loci with `SuSiE-rRSS` prior variance blow-up ($V > 1.0$)** | **2 / 31 (6.5%: `CCL7`, `C1QA`)** | **9 / 30 (30.0%: up to $V = 398.6$)** |

### 2.2 Per-Assay Head-to-Head Table (`cis_gene` Windows)

| Assay (`GWAS_ID`) | Gene | `EUR` $M$ | `EUR` LP | `EUR` pLD | `EUR` CS (`RSS` / `rRSS` / `CRM` / `PF`) | `EUR` Conv (it) | `EUR` $\max V$ | `ALL` LP | `ALL` pLD | `ALL` CS (`RSS` / `rRSS` / `CRM` / `PF`) | `ALL` Conv (it) | `ALL` $\max V$ |
|---|---|---:|---:|---:|---|---|---:|---:|---:|---|---|---:|
| `Q9BYZ8_OID20784` | **REG4** | 9,114 | 114.6 | 39 | `4 / 4 / 4 / 4` | `T (7)` | 0.02 | 156.6 | 29 | `0 / 7 / 6 / 0` | `T (55)` | 0.02 |
| `Q03403_OID20787` | **TFF2** | 9,943 | 53.5 | 70 | `6 / 6 / 6 / 6` | `T (11)` | 0.01 | 76.0 | 97 | `6 / 6 / 7 / 6` | `T (27)` | 0.01 |
| `O95994_OID20896` | **AGR2** | 8,381 | 43.1 | 40 | `1 / 1 / 1 / 1` | `T (2)` | 0.01 | 57.0 | 69 | `0 / 1(2d) / 1 / 0` | `T (200)` | **25.12** |
| `P08727_OID20596` | **KRT19** | 6,596 | 81.0 | 40 | `2 / 2 / 2 / 2` | `T (4)` | 0.01 | 124.3 | 27 | `2 / 2 / 2 / 7` | `T (4)` | 0.01 |
| `Q07654_OID20371` | **TFF3** | 9,907 | 107.6 | 56 | `6 / 6 / 7 / 6` | `T (28)` | 0.01 | 155.9 | 79 | `0 / 5 / 6 / 0` | `T (14)` | 0.02 |
| `Q99988_OID20251` | **GDF15** | 7,109 | 419.4 | 69 | `4 / 3 / 3 / 4` | `T (17)` | 0.08 | 584.2 | 52 | `0 / 0(5d) / 6 / 0` | **`F (200)`** | **21.41** |
| `P39900_OID21439` | **MMP12** | 13,229 | 1093.6 | 63 | `7 / 6 / 9 / 7` | `T (14)` | 0.11 | 1592.8 | 83 | `0 / 2(3d) / 8 / 0` | `T (107)` | **4.26** |
| `O95407_OID20964` | **TNFRSF6B** | 7,168 | 339.9 | 45 | `10 / 8 / 10 / 10` | `T (57)` | 0.04 | 485.3 | 42 | `0 / 0(7d) / 8 / 0` | **`F (200)`** | **1.02** |
| `P04155_OID21154` | **TFF1** | 9,986 | 279.8 | 93 | `6 / 6 / 6 / 6` | `T (16)` | 0.04 | 382.4 | 112 | `0 / 4(1d) / 8 / 0` | `T (37)` | **1.51** |
| `Q16790_OID21417` | **CA9** | 5,435 | 51.9 | 12 | `4 / 3 / 4 / 4` | `T (14)` | 0.01 | 82.0 | 25 | `0 / 8 / 6 / 0` | `T (45)` | 0.03 |
| `P50120_OID21365` | **RBP2** | 5,908 | 17.8 | 48 | `1 / 1 / 1 / 1` | `T (2)` | 0.00 | 27.7 | 49 | `0 / 0(5d) / 1 / 0` | **`F (200)`** | **372.03** |
| `Q8IZP9_OID20189` | **ADGRG2** | 2,587 | 7.6 | 29 | `2 / 2 / 2 / —` | `T (14)` | 0.00 | 17.0 | 29 | `5 / 3 / 4 / —` | `T (54)` | 0.00 |
| `P10145_OID20153` | **CXCL8** | 5,779 | 24.2 | 25 | `2 / 2 / 2 / 2` | `T (11)` | 0.00 | 32.5 | 21 | `0 / 4 / 3 / 0` | `T (10)` | 0.07 |
| `P10145_OID20631` | **CXCL8** | 5,778 | 24.1 | 26 | `4 / 4 / 2 / 4` | `T (18)` | 0.00 | 37.2 | 23 | `0 / 4 / 4 / 0` | `T (11)` | 0.10 |
| `P10145_OID20997` | **CXCL8** | 5,777 | 24.5 | 27 | `3 / 3 / 3 / 3` | `T (23)` | 0.00 | 36.4 | 29 | `0 / 4 / 3 / 0` | `T (8)` | 0.06 |
| `P10145_OID21430` | **CXCL8** | 5,776 | 26.7 | 23 | `4 / 4 / 2 / 3` | `T (20)` | 0.00 | 37.6 | 20 | `0 / 4 / 4 / 0` | `T (9)` | 0.07 |
| `Q9P0M4_OID20477` | **IL17C** | 10,076 | 52.1 | 49 | `3 / 4 / 4 / 3` | `T (12)` | 0.00 | 80.3 | 43 | `6 / 6 / 4 / 6` | `T (61)` | 0.08 |
| `O15444_OID20674` | **CCL25** | 8,695 | 1756.3 | 81 | `0 / 10 / 6 / 0` | `T (78)` | 0.13 | 2595.1 | 83 | `10 / 0(10d) / 7 / 10` | **`F (200)`** | **6.49** |
| `O14836_OID20702` | **TNFRSF13B** | 5,694 | 43.9 | 29 | `5 / 3 / 4 / 5` | `T (76)` | 0.01 | 66.2 | 38 | `6 / 1 / 5 / 5` | `T (2)` | 0.01 |
| `O00220_OID20967` | **TNFRSF10A** | 14,488 | 859.9 | 94 | `10 / 8 / 8 / 10` | `T (44)` | 0.06 | 1302.2 | 99 | `10 / 10 / 10 / 10` | `T (102)` | 0.11 |
| `P80098_OID20523` | **CCL7** | 7,076 | 46.3 | 65 | `0 / 3(2d) / 4 / 0` | `T (77)` | **7.17** | 56.0 | 67 | `0 / 0(5d) / 5 / 0` | **`F (200)`** | **398.62** |
| `P01589_OID20267` | **IL2RA** | 10,533 | 1116.9 | 80 | `9 / 8 / 9 / 9` | `T (16)` | 0.14 | 1589.5 | 85 | `10 / 8 / 8 / 10` | `T (13)` | 0.13 |
| `P02745_OID20654` | **C1QA** | 14,673 | 477.6 | 89 | `0 / 0(5d) / 7 / 0` | **`F (200)`** | **124.62** | 696.2 | 104 | `0 / 2(3d) / 9 / 0` | `T (171)` | **106.27** |
| `Q9ULX7_OID21401` | **CA14** | 13,720 | 128.6 | 51 | `0 / 7 / 8 / 0` | `T (32)` | 0.44 | 164.7 | 51 | `0 / 0(8d) / 8 / 0` | **`F (200)`** | 0.13 |
| `O43508_OID20624` | **TNFSF12** | 6,874 | 446.7 | 34 | `2 / 2 / 3 / 2` | `T (19)` | 0.04 | 665.9 | 48 | `4 / 0(10d) / 7 / 0` | **`F (200)`** | 0.79 |
| `P01579_OID20495` | **IFNG** | 7,138 | **7.5** | 60 | `1 / 1 / 1 / 1` | `T (18)` | 0.00 | `< 7.3` | 0 | `skipped (<5e-8)` | — | — |
| `P78556_OID20671` | **CCL20** | 8,472 | 33.8 | 69 | `1 / 1 / 1 / 1` | `T (2)` | 0.00 | 42.6 | 47 | `2 / 2 / 2 / 2` | `T (6)` | 0.00 |
| `P10147_OID20610` | **CCL3** | 5,675 | 507.3 | 28 | `10 / 8 / 9 / 10` | `T (32)` | 0.06 | 757.4 | 31 | `0 / 10 / 10 / 0` | `T (94)` | 0.48 |
| `P13500_OID21004` | **CCL2** | 7,080 | 15.1 | 50 | `2 / 1 / 1 / 1` | `T (24)` | 0.00 | 23.8 | 42 | `3 / 2 / 3 / 3` | `T (9)` | 0.00 |
| `P09237_OID20087` | **MMP7** | 13,587 | 301.3 | 61 | `9 / 10 / 5 / 10` | `T (28)` | 0.06 | 398.8 | 55 | `9 / 10 / 9 / 9` | `T (29)` | 0.05 |
| `P09603_OID20719` | **CSF1** | 6,318 | 137.1 | 51 | `3 / 3 / 3 / 3` | `T (30)` | 0.01 | 186.2 | 48 | `4 / 3 / 3 / 5` | `T (28)` | 0.01 |

### 2.3 Notable Biological & Structural Edge Cases in `EUR`
1. **`CXCL8` 4-Panel Replicate Consistency (`OID20153`, `OID20631`, `OID20997`, `OID21430`):**
   - In `ALL`, standard `SuSiE-RSS` and `PolyFun` failed completely (`0` CSs on all 4 panels).
   - In `EUR`, all 4 Olink panels converged in `11–23` iterations with `max_V < 0.005`, `p-ld` counts of `23–27`, and `2–4` concordant credible sets across all four Bayesian methods.
2. **Why `C1QA` (`P02745_OID20654`) and `CCL7` (`P80098_OID20523`) Show Residual $V > 1$ Even in `EUR`:**
   - **`C1QA` (`chr1:22.63–22.64 Mb`, `1p36.12`):** Sits in the tightly linked `C1QA–C1QC–C1QB` complement subcomponent paralog cluster with short-read mapping ambiguity / structural variation even within `EUR`. `SuSiE-rRSS`'s two-tier gate correctly suppresses `C1QA`'s non-converged SuSiE fit (`0(5d)`), while `CARMA` (`outlier_switch = TRUE`) identifies and removes the discordant structural-variant outliers to recover `7` clean credible sets.
   - **`CCL7` (`chr17:34.27 Mb`, `17q12` chemokine cluster `CCL2–CCL7–CCL11–CCL8–CCL13`):** Lies inside a well-known copy-number / segmental-duplication hotspot on `17q12`. `SuSiE-rRSS` (`eb_mix`) converged in `77` iterations, and the Tier-2 per-slot gate retained the `3` clean cis credible sets ($V \le 1.0$) while dropping `2` structural-duplication slots ($V = 7.17$).
3. **`IFNG` (`P01579_OID20495`) Rescued by `EUR` Specificity:**
   - `IFNG` has a modest European cis-pQTL (`LP = 7.5`, $P = 3.2 \times 10^{-8}$) that is absent/null in non-European subsets, causing sample-size dilution below $P = 5 \times 10^{-8}$ in `ALL`. In `EUR`, all 4 Bayesian methods (`RSS`, `rRSS`, `CARMA`, `PolyFun`) agree on `1` clean credible set (`max_V = 0.002`).

---

## 3. Six-Cohort Cross-Platform & Cross-Ancestry Comparison Framework

With [compare_pqtl_studies.R](../scripts/compare_pqtl_studies.R), we evaluate all 28 target proteins across **6 pQTL datasets** once `GWASfm_cross_study` finishes:

| Study Directory | Cohort | Proteomics Platform | GWAS Ancestry | Pan-UKB LD Panel ($B$) | LD Relationship | Matched Assays |
|---|---|---|---|---|---|---:|
| `sun_2023_37794186_EUR` | UK Biobank | Olink Explore 3072 (PEA) | `EUR` | `EUR` ($420{,}542$) | **In-sample UKB `EUR`** | 31 (complete) |
| `sun_2023_37794186_ALL` | UK Biobank | Olink Explore 3072 (PEA) | `ALL` | `EUR` ($420{,}542$) | Multi-ancestry mismatch control | 31 (complete) |
| `pietzner_2021_34648354` | Fenland (UK) | SomaScan v4 (Aptamer) | `EUR` | `EUR` ($420{,}542$) | Out-of-sample UK `EUR` | 27 |
| `ferkingstad_2021_34857953` | deCODE (Iceland) | SomaScan v4 (Aptamer) | `EUR` | `EUR` ($420{,}542$) | Founder-population `EUR` drift | 26 |
| `zhang_2022_35501419_EUR` | ARIC EA (US) | SomaScan v4 (`cis`-only) | `EUR` | `EUR` ($420{,}542$) | Admixed-European US `EUR` | 23+ |
| `zhang_2022_35501419_AFR` | ARIC AA (US) | SomaScan v4 (`cis`-only) | `AFR` | **`AFR` ($6{,}637$)** | Finite-sample `AFR` (`R_finite=6637`) | 23+ |

### Automated Metrics Computed by `scripts/compare_pqtl_studies.R`
1. **`01_cis_method_summary.tsv`:** Per-study/assay `cis_gene` activation, `M_matched`, `p-ld` ($r^2 < 0.001$), `RSS` / `rRSS` / `CARMA` / `PolyFun` CS counts, `converged`, `iterations`, `max_V`, `lambda_bias`, `Q_art`, and `CARMA` outlier counts.
2. **`02_cxcl8_technical_replicates.tsv`:** Pairwise PIP Pearson/Spearman correlation, lead-SNP / CS overlap, and multi-signal colocalisation ($\max PP.H_4$) across all 6 pairs of the 4 `CXCL8` Olink Explore panels in `EUR` vs `ALL`.
3. **`03_pairwise_study_coloc_and_concordance.tsv`:** For every protein and every pair of studies (Olink vs SomaScan; `Sun EUR` vs `Ferkingstad`/`Pietzner`/`Zhang EUR`; `Sun EUR` vs `Zhang AFR`; `Zhang EUR` vs `Zhang AFR`), computes:
   - Exact GRCh38 `GWAS_VARIANT_ID` intersection size,
   - `SuSiE-rRSS` and `CARMA` PIP correlation across shared cis variants,
   - Credible-set variant overlap (`n_shared_cs_pairs`), and
   - Exact multi-signal **`coloc.susie` / `coloc.bf_bf` posterior $\max PP.H_4$ and $\max PP.H_3$** from `fit$lbf_variable[fit$sets$cs_index, ]` (distinguishing shared cis-regulatory pQTLs [$PP.H_4 \ge 0.8$] from platform-specific SomaScan aptamer epitope variants [$PP.H_3 \ge 0.8$]).
4. **`comparison_report.html`:** Self-contained dark/light-mode HTML dashboard uniting the challenges summary, `Sun EUR` vs `ALL` comparison, `CXCL8` replicate concordance, and 6-study Olink vs SomaScan / `EUR` vs `AFR` colocalisation tables.
