# Running Analyses & Opt-In `MultiSuSiE`

`scripts/finemap.sh` is the unified local and SLURM entrypoint for **GWASfinemapping**. By default, it runs the 4-stage single-GWAS workflow across up to 5 methods (`p-ld`, `susie-rss`, `susie-rrss`, `carma`, `polyfun`). Multi-cohort / multi-ancestry joint `cis` fine-mapping (`MultiSuSiE`) is **disabled by default** and activated explicitly when you pass `--multisusie-studies` or `--multisusie-manifest`.

---

## 1. Standard Single-Study Execution Modes

Supply `--outdir` (defaults to `$GWASFM_WORK_ROOT`) plus **one** input source: `--input-dir`, `--manifest`, or `--vcf`.

### Mode A: Automatic `GWASprocessing` Tree Discovery (`--input-dir`)

Point `--input-dir` at any of three levels in a standardised [GWASprocessing](https://github.com/mattlee821/GWASprocessing) output hierarchy — no file renaming or intermediate conversion is needed:

```bash
# 1. Single trait directory
bash scripts/finemap.sh \
  --input-dir /path/to/GWASprocessing/GWAS/locke_2015_25673413/standardise/GWAScatalog_GCST002783 \
  --ld-pop EUR

# 2. Whole study directory (discovers every <GWAS_ID>/ under standardise/)
sbatch scripts/finemap.sh \
  --input-dir /path/to/GWASprocessing/GWAS/sun_2023_37794186_EUR \
  --ld-pop EUR \
  --max-forks 16 \
  --slurm

# 3. Root GWAS directory filtered to one study
sbatch scripts/finemap.sh \
  --input-dir /path/to/GWASprocessing/GWAS \
  --only-study pietzner_2021_34648354 \
  --ld-pop EUR \
  --slurm
```

> [!NOTE]
> **Verbatim `<study_dir>` Mirroring:** When discovering ancestry-split study folders (e.g. `sun_2023_37794186_EUR`, `zhang_2022_35501419_AFR`), `prepare_manifest.py` preserves the exact `<study_dir>` folder name in `<outdir>/<study_dir>/<gwas_id>/` even if upstream `metadata.json` shares a generic `study_id`, preventing cross-ancestry output collisions.

### Mode B: Custom Batch Manifest (`--manifest`)

Use a tab-separated file with one row per GWAS when combining traits across studies or assigning per-row LD panels (`EUR` vs `AFR`) in a single batch:

```tsv
study_id	gwas_id	gwas_dir	ld_population	gene_id
sun_2023_37794186_EUR	P39900_OID21439	/path/to/GWAS/sun_2023_37794186_EUR/standardise/P39900_OID21439	EUR	P39900
zhang_2022_35501419_EUR	P39900_4496-60	/path/to/GWAS/zhang_2022_35501419_EUR/standardise/P39900_4496-60	EUR	P39900
zhang_2022_35501419_AFR	P39900_4496-60	/path/to/GWAS/zhang_2022_35501419_AFR/standardise/P39900_4496-60	AFR	P39900
```

```bash
sbatch scripts/finemap.sh \
  --manifest manifests/my_batch.tsv \
  --outdir /scratch/user/finemapping \
  --max-forks 16 \
  --slurm
```

| Manifest Column | Required | Description |
|---|---|---|
| `gwas_dir` or `vcf` | Yes (one) | Path to `standardise/<GWAS_ID>` directory (`gwas_dir`) or `.gwas.vcf.gz` file (`vcf`) |
| `study_id`, `gwas_id` | Optional with `gwas_dir` | Output path `<outdir>/<study_id>/<gwas_id>/` (inferred automatically if omitted) |
| `ld_population` | Optional | Per-row Pan-UKB LD panel override (`AFR`, `AMR`, `CSA`, `EAS`, `EUR`, `MID`; defaults to `--ld-pop`) |
| `gene_id` | Optional | UniProt accession (`P39900`) or HGNC symbol (`MMP12`) for `cis_gene` (`±1 Mb`) activation |
| `metadata`, `qc`, `trait`, `ancestry` | Optional | Explicit sidecar paths or metadata overrides |

### Mode C: Single GWAS-VCF (`--vcf`)

```bash
bash scripts/finemap.sh \
  --vcf /path/to/P39900_OID21439.gwas.vcf.gz \
  --study-id sun_2023_37794186_EUR \
  --gwas-id P39900_OID21439 \
  --gene-id P39900 \
  --ld-pop EUR
```

---

## 2. Opt-In Multi-Cohort / Multi-Ancestry `MultiSuSiE` Mode

### When to Use `MultiSuSiE` (and When Not To)

Standard inverse-variance weighted (IVW) meta-analysis across cohorts with heterogeneous coverage (e.g. `±1 Mb` WGS/imputed arrays vs `±500 kb` cis-only summary stats) or multiple ancestries (`EUR` + `AFR`) creates sharp per-SNP sample-size steps and violates any single-population LD matrix.

Opt-in **`MultiSuSiE`** (`--multisusie-studies` or `--multisusie-manifest`) jointly models $K \ge 2$ cohorts **measuring the same molecular reagent or trait** without pooling summary statistics:

- Each cohort $k$ retains its own sample size $N_k$, per-allele `(ES, SE, AF)`, and **ancestry-matched Pan-UKB LD matrix** (e.g. `Pan-UKB EUR` for European cohorts and `Pan-UKB AFR` for African-ancestry cohorts).
- Within each LD ancestry, variants require complete observation across cohorts sharing that LD matrix (preventing cross-SNP coverage steps); across distinct ancestries, the **union** of panel-matched variants is retained (setting missing ancestry columns to `XTY = 0, X_l2 = 0` so ancestry-specific variants are preserved).
- Uses an unconditionally symmetric positive-definite 64-bit Woodbury LBF solver and the same two-tier credible-set gate (`converged` + per-slot $V \le 1.0$) as `susie-rrss`.

> [!IMPORTANT]
> **Best-Practice Cohort Selection Rules (from our [28-Protein Benchmark](sun_eur_vs_all_comparison_report.md)):**
> 1. **Combine cohorts within the same assay platform/reagent** (e.g. SomaScan SeqId `4496-60` across `Fenland EUR + ARIC EUR + ARIC AFR`), **never** across distinct platforms (Olink vs SomaScan), so cross-platform colocalisation remains an independent validation test of reagent binding specificity.
> 2. **Exclude founder-isolate cohorts lacking in-sample LD** (such as Icelandic deCODE `Ferkingstad` mapped against UK `Pan-UKB EUR`) from joint `MultiSuSiE` groups on high-$|Z|$ loci: run founder cohorts separately via single-study `susie-rrss` + `carma` instead.

### Option 1: Auto-Pair Shared `<gwas_id>` Folders Across Study Directories (`--multisusie-studies`)

When two or more `GWASprocessing` study directories share standardised `<gwas_id>` folder names (such as SomaScan `<UniProt>_<SeqId>` identifiers), pass them as a comma-separated list. Each study's ancestry (`EUR`, `AFR`, etc.) is auto-detected from its `metadata.json` or directory suffix (`_EUR`, `_AFR`):

```bash
# Dry-run preflight to inspect discovered multi-cohort groups
bash scripts/finemap.sh \
  --multisusie-studies /path/to/pietzner_2021_34648354,/path/to/zhang_2022_35501419_EUR,/path/to/zhang_2022_35501419_AFR \
  --group-study-id somascan_multisusie_3cohort \
  --dry-run

# Submit joint 3-cohort EUR+AFR cis fine-mapping on SLURM
# (--cis-significance-lp 0 evaluates all target cis windows even when single cohorts are sub-threshold!)
sbatch scripts/finemap.sh \
  --multisusie-studies /path/to/pietzner_2021_34648354,/path/to/zhang_2022_35501419_EUR,/path/to/zhang_2022_35501419_AFR \
  --group-study-id somascan_multisusie_3cohort \
  --cis-significance-lp 0 \
  --only-gwas P10145_3447-64,Q03403_9191-8,Q99988_4374-45 \
  --outdir /scratch/user/finemapping \
  --max-forks 8 \
  --slurm
```

### Option 2: Explicit Grouping Manifest (`--multisusie-manifest`)

When `<gwas_id>` names differ across cohorts or you want per-cohort `ld_population` / `cohort_id` labels, pass a TSV manifest where rows sharing `(group_study_id, group_id)` are jointly fine-mapped:

```tsv
group_study_id	group_id	cohort_id	gwas_dir	ld_population	gene_id
somascan_multisusie_3cohort	P39900_4496-60	pietzner_EUR	/path/to/pietzner_2021_34648354/standardise/P39900_4496-60	EUR	P39900
somascan_multisusie_3cohort	P39900_4496-60	zhang_EUR	/path/to/zhang_2022_35501419_EUR/standardise/P39900_4496-60	EUR	P39900
somascan_multisusie_3cohort	P39900_4496-60	zhang_AFR	/path/to/zhang_2022_35501419_AFR/standardise/P39900_4496-60	AFR	P39900
```

```bash
sbatch scripts/finemap.sh \
  --multisusie-manifest manifests/somascan_3cohort_groups.tsv \
  --cis-significance-lp 0 \
  --outdir /scratch/user/finemapping \
  --slurm
```

Each completed group writes `<outdir>/<group_study_id>/<group_id>/methods/multisusie.rda` (plus `{metadata,provenance,qc}.json` and `report.html`, including a per-cohort kriging & `s_hat` diagnostic table), which drops directly into `coloc::coloc.susie()` and `scripts/compare_pqtl_studies.R` alongside single-study runs.

---

## 2b. Kriging Diagnostics vs Opt-In Outlier Removal: When to Use for `X`, Why Off by Default for `Y`

Every `susie-rrss` locus ($M \le 20{,}000$) and every `MultiSuSiE` cohort automatically computes **non-destructive leave-one-out kriging & `DENTIST` diagnostics** (`s_hat`, `KRIGING_EXP_Z`, `KRIGING_STD_RESID`, `KRIGING_LOGLR`, `KRIGING_OUTLIER`, `kriging_cs_outlier_count`) in $O(M^2)$ time by reusing the existing correlation eigendecomposition — **without removing any variants by default (`--remove-kriging-outliers false`)**.

### Case `X` — When You SHOULD Enable `--remove-kriging-outliers true`

Enable active kriging outlier removal when fine-mapping **legacy or heterogeneous complex-trait summary statistics at modest $|Z|$** where `report.html` Panel 3 (`Kriging Expected Z vs Observed Z`) or `qc.json` (`kriging_cs_outlier_regions > 0`) reveals isolated sign-reversed variants caused by upstream harmonisation artifacts (such as rsID-only strand flips or fixed-effect IVW meta-analyses with unrecorded cohort allele swaps):

```bash
# Case X: Re-run a legacy complex-trait meta-analysis dropping flipped kriging outliers
# (Produces a distinct settings_hash so default runs in other studies are unaffected)
bash scripts/finemap.sh \
  --input-dir /path/to/GWASprocessing/GWAS/locke_2015_25673413/standardise/GWAScatalog_GCST002783 \
  --ld-pop EUR \
  --remove-kriging-outliers true
```

### Case `Y` — Why You Should NOT Enable `--remove-kriging-outliers` by Default (Especially for Proteins)

Keep `--remove-kriging-outliers false` (the default) for **standardised biobank GWAS-VCFs, multi-causal loci, founder/admixed cohorts, and all protein `cis-pQTL` / `cis-eQTL` studies**:

1. **Multi-causal `cis` loci (`MMP12`, `GDF15`):** An independent secondary causal variant with $z_2 = -6.5$ in modest positive LD ($r_{12} = +0.25$) with an extreme $+55$ primary lead has a positive leave-one-out expectation ($\hat{\mu}_2 > 0$) opposite to its true observed score ($z_2 < 0$), triggering `KRIGING_LOGLR > 2` even with exact alleles and in-sample LD. Active removal deletes genuine secondary causal signals before `SuSiE-rRSS` can jointly condition on the primary lead.
2. **Coding PAVs & founder drift (`deCODE`):** Low-LD coding missense variants and founder-drifted haplotypes exhibit large conditional residuals; dropping them deletes the true causal variant rather than a harmonisation error.
3. **Full diagnostic visibility with zero risk:** Leaving `--remove-kriging-outliers false` still flags any credible set containing a suspect variant via `kriging_cs_outlier_count`, `ld_mismatch_warning`, `qc.json`, and red `×` markers in `report.html` Panel 3. See the [Protein cis-pQTL Guide](protein-pqtl-example.md) and [LD-Mismatch & Kriging Diagnostics](diagnostics-and-kriging.md) for full details.

---

## 3. Complete CLI Options Reference

### Method, Threshold & Window Controls

| Flag | Default | Description |
|---|---|---|
| `--methods LIST` | `"p-ld susie-rss susie-rrss carma polyfun"` | Space- or comma-separated single-study methods to execute |
| `--ld-pop POP` | `EUR` | Default Pan-UKB LD panel (`AFR`, `AMR`, `CSA`, `EAS`, `EUR`, `MID`) |
| `--allow-ancestry-mismatch` | `false` | Permit running a single-ancestry GWAS against a different single-ancestry LD panel (multi-ancestry `ALL` GWAS are allowed automatically and flagged in `qc.json`) |
| `--min-maf FLOAT` | `0.01` | Minimum `FORMAT/AF` minor allele frequency $\min(\text{AF}, 1-\text{AF})$ in screening & locus prep (`0` disables) |
| `--significance-lp FLOAT` | `7.30103` | $-\log_{10} P$ activation threshold ($P \le 5 \times 10^{-8}$) |
| `--cis-significance-lp FLOAT` | Inherits `--significance-lp` | Optional $-\log_{10} P$ override for `cis_gene` windows (`0` forces cis evaluation regardless of peak $P$, useful in `MultiSuSiE`) |
| `--region-flank-bp INT` | `1500000` | Half-width around greedy lead variants (`±1.5 Mb` = 3 Mb nominal window; merged up to 9.5 Mb) |
| `--exclude-mhc` | `false` | Skip the extended MHC (`chr6:25.5–34.0 Mb` GRCh38) instead of fine-mapping it as 1 merged region |
| `--max-loci INT` | `0` (all) | Cap analyzed regions per GWAS (`cis_gene` ordered first, then descending lead `LP`) |
| `--max-iter INT` | `100` | Max IBSS iterations for SuSiE methods on `lead_locus` and `mhc` regions |
| `--cis-max-iter INT` | `200` | Max IBSS iterations for SuSiE methods on `cis_gene` windows (and `MultiSuSiE`) |
| `--susie-refine` | `false` | Enable `susieR` multi-restart CS refinement (`refine = TRUE`) |
| `--remove-kriging-outliers` | `false` | Opt-in active pre-fit removal of kriging sign-flip outliers (`KRIGING_LOGLR > 2` and `|z_adj| > 2`) in `susie-rrss` and `multisusie` (diagnostics are **always** recorded regardless of this flag) |
| `--output-format FMT` | `rda` | Saved method file format: `rda`, `rds`, or `both` |
| `--compress COMP` | `xz` | Compression for `.rda`/`.rds`: `xz` (~35% smaller) or `gzip` (faster write) |

### Execution, SLURM Scaling & Idempotency Controls

| Flag | Default | Description |
|---|---|---|
| `--slurm` / `--profile slurm` | Auto under `sbatch`, else `local` | Dispatch Nextflow worker tasks (`SCREEN_GWAS`, `FINEMAP_LOCUS`, `FINALISE_GWAS`, `MULTISUSIE_GROUP`) as SLURM jobs |
| `--slurm-partition PART` | `$GWASFM_SLURM_PARTITION` | SLURM partition(s) for worker jobs |
| `--slurm-exclude NODES` | `$GWASFM_SLURM_EXCLUDE` | Comma-separated SLURM nodes to exclude |
| `--max-forks INT` | `$GWASFM_MAX_FORKS` or `2` | Maximum concurrent Nextflow tasks / SLURM queue depth |
| `--temp-dir DIR` | `$GWASFM_TEMP_ROOT` | Fast scratch root for Nextflow workdirs, run logs, and `TMPDIR` |
| `--bm-path DIR` | Live S3 | Optional local mirror/cache of Pan-UKB `.ldadj.bm` BlockMatrices |
| `--container ENGINE` | Host tools | Execute worker tasks inside `apptainer` or `docker` |
| `--row-limit INT` | `0` (all) | Process at most `INT` incomplete GWAS / groups in this invocation |
| `--report true\|false` | `true` | Render self-contained `report.html` during finalisation (can defer and run `render_report.R` later) |
| `--keep-intermediates` | `false` | Retain per-locus lifted VCFs, `matched_variants.tsv`, and `matrix.npy` under `<outdir>/<study_id>/<gwas_id>/loci/` |
| `--dry-run` | `false` | Validate inputs, state ledger, and reference files without starting Nextflow |
| `--force` | `false` | Re-run rows even when `.finemap_state.tsv` records `COMPLETE` for the current `settings_hash` |
| `--resume` | Off | Pass `-resume` to Nextflow task cache |

---

## 4. SLURM Resource Scaling & Node-Local Scratch Hygiene

Under `--slurm` (`conf/slurm.config` + `finemapping.nf`):

- **Standard loci (`lead_locus`, `cis_gene`):** Request `4` CPUs and `24 GB` RAM on Attempt 1 (scaling automatically on retry up to `32 GB` / `48 GB`).
- **Extended MHC (`chr6:25.5–34.0 Mb`, $M \approx 28{,}000\text{–}45{,}000$ variants):** Automatically detected before submission and scheduled with **`16` CPUs, `160 GB` RAM, and `24h` walltime** on Attempt 1, paired with `OPENBLAS_NUM_THREADS=16` and R `options(matprod = "blas")`.
- **Zero home-directory LD bloat:** `FINEMAP_LOCUS` writes `matrix.npy`, regional `.parquet` slices, and R/DuckDB temp files to compute-node local `/tmp` (`$SLURM_TMPDIR`) and wipes them in a `finally:` block as soon as all requested methods complete for that locus.

---

## 5. Idempotent State Tracking & Deferred Reports

Every single-study row and `MultiSuSiE` group computes a deterministic 16-character `settings_hash` over input VCF size(s), LD panel(s), target gene, active methods, output format, MAF & significance gates, window flanks, MHC policy, `L`, and `max_iter`, recorded in `<outdir>/.finemap_state.tsv`.

- **Instant skip:** Re-running the exact same `scripts/finemap.sh` command skips completed traits in `<1 second` without booting the Nextflow JVM.
- **Automatic invalidation:** Changing `--methods`, `--ld-pop`, `--min-maf`, or passing `--force` automatically re-queues only the affected rows.
- **Deferred HTML rendering:** If a 10,000-protein batch was run with `--report false`, render any individual GWAS report on demand from its saved `.rda`/`.rds` and `.json` files:

```bash
source scripts/source.sh
Rscript scripts/render_report.R --gwas-dir /scratch/user/finemapping/sun_2023_37794186_EUR/P39900_OID21439
```
