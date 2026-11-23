# Project Handover & Persistent Memory (`GWASfinemapping`)

**Last Updated:** `2026-11-23`  
**Repository:** [https://github.com/mattlee821/GWASfinemapping](https://github.com/mattlee821/GWASfinemapping) (`v0.1.1` released on `main`)  
**Documentation & Benchmark Site:** [https://mattlee821.github.io/GWASfinemapping/](https://mattlee821.github.io/GWASfinemapping/)  
**Published Container Image:** `ghcr.io/mattlee821/gwasfinemapping:latest` (`:0.1.1`, `linux/amd64`)  
**Prior Antigravity Session Transcript & Artifacts (Pre-Rename `finemapping/`):**
- Conversation ID: `5a6f5d75-6853-4976-bb66-0bd2d4b37c6f` ([session link](conversation://5a6f5d75-6853-4976-bb66-0bd2d4b37c6f))
- Compact transcript: `~/.gemini/antigravity/brain/5a6f5d75-6853-4976-bb66-0bd2d4b37c6f/.system_generated/logs/transcript.jsonl`
- Artifacts directory: `~/.gemini/antigravity/brain/5a6f5d75-6853-4976-bb66-0bd2d4b37c6f/`

---

## 1. Standing Rules for Agents Working on `GWASfinemapping`

1. **Git hygiene:** Never run `git commit`, `git push`, or `git merge` without Matt's explicit go-ahead. Author commits as `mattlee821 <30179533+mattlee821@users.noreply.github.com>` with **no** `Co-Authored-By` trailer. Before pushing to GitHub, scan tracked files for private cohort identifiers, hostnames, or tokens (`0` allowed).
2. **Local Mac vs HPC workloads:** Never run real GWAS fine-mapping workloads or S3 BlockMatrix downloads on the local Mac (which is inside a cloud-synced OneDrive folder). Run only fast unit/package checks locally (`<60s`); run real GWAS workloads on the SLURM cluster (`~/GWASfinemapping` with outputs/scratch on `/scratch/`).

---

## 2. Current State (`v0.1.1` — Complete & Production-Ready)

- **Default Single-Study Pipeline (`scripts/finemap.sh` + `finemapping.nf`):**
  - 4-stage cascading DAG (`PREPARE_MANIFEST` → `SCREEN_GWAS` → `FINEMAP_LOCUS` → `FINALISE_GWAS` → `UPDATE_STATE`) across 5 default methods (`p-ld`, `susie-rss`, `susie-rrss`, `carma`, `polyfun`) against the 6 public Pan-UKB LD panels (`AFR`, `AMR`, `CSA`, `EAS`, `EUR`, `MID`).
  - Storage-light output contract (~25–40 MB per active GWAS): `xz`-compressed `coloc::coloc.susie()`-ready `.rda`/`.rds` files, `{metadata,provenance,qc}.json`, self-contained `report.html`, and batch `manifest.html` + `.finemap_state.tsv` (`settings_hash` idempotency). Node-local `$SLURM_TMPDIR` / `/tmp` LD matrices (`matrix.npy`) are wiped immediately after each locus completes.
  - **Zero-extra-$O(M^3)$ Kriging & `DENTIST` LD-Mismatch Diagnostics (`R/finemap_susie.R`, `src/fast_matrix.c`, `scripts/render_report.R`):** `susie-rrss` reuses its existing `eig <- eigen(R)` ($M \le 20{,}000$) via C `gwasfm_eigen_diag_quadform` ($O(M^2)$ time, `0` extra $M \times M$ matrices) to record null-MLE `s_hat`, `KRIGING_EXP_Z`, `KRIGING_STD_RESID` (`DENTIST`), `KRIGING_LOGLR`, `KRIGING_OUTLIER`, and `kriging_cs_outlier_count` non-destructively by default (`--remove-kriging-outliers false`), with opt-in active outlier removal (`--remove-kriging-outliers true`) for legacy/meta-analysis GWAS (`docs/diagnostics-and-kriging.md`, `docs/protein-pqtl-example.md`).
- **Opt-In Multi-Cohort / Multi-Ancestry Mode (`MultiSuSiE` — Off by Default):**
  - Activated explicitly via `--multisusie-studies <dir1,dir2,...>` or `--multisusie-manifest <tsv>` (`PREPARE_MULTISUSIE` → `MULTISUSIE_GROUP` → `UPDATE_STATE` in `scripts/run_multisusie.py` and `scripts/save_multisusie_rda.R`).
  - Uses a 64-bit unconditionally SPD Woodbury LBF solver $B_j = I_K + Q_j^{1/2} A Q_j^{1/2}$ ($\lambda_{\min}(B_j) \ge 1$) with a $\sigma^2_k \ge 10^{-2}$ floor, within-ancestry variant intersection + cross-ancestry variant union, per-cohort `(ancestry, obs_idx)`-cached kriging & `DENTIST` summaries (`kriging_by_cohort` in `qc.json` and `report.html`), and two-tier CS validity gating (`converged` + per-slot $V_l \le 1.0$).
- **CI / CD & Documentation:**
  - `.github/workflows/docs.yml`: Strict MkDocs Material build (`mkdocs build --strict`) deployed to `https://mattlee821.github.io/GWASfinemapping/`.
  - `.github/workflows/container.yml`: Builds `Dockerfile`, runs full `bash scripts/test.sh && bash scripts/benchmark.sh --offline` inside the container, and publishes `ghcr.io/mattlee821/gwasfinemapping:latest` (and semver tags).

---

## 3. Key Empirical Findings & Safeguards Solved (`28` CRC Proteins / `216` Fits)

Full tables and vector figures are published in `docs/sun_eur_vs_all_comparison_report.md` (and cross-referenced for downstream pQTL classification in `../pQTLmapping/docs/gwasfinemapping_empirical_findings.md`):

1. **Ancestry & Founder-LD Regimes:** In-sample/matched cohorts achieved `96.8%–100%` `SuSiE-rRSS` convergence (`Sun EUR`: `30/31`; `Pietzner EUR`: `14/14`; `Zhang EUR`: `12/12`; `Zhang AFR` on `Pan-UKB AFR` [$B = 6{,}637$]: `10/10` in median `6.5` iterations). Multi-ancestry `Sun ALL` and founder-drifted `deCODE` (`Ferkingstad`) collapsed standard `SuSiE-RSS` to `43%`, rescued to `77%–79%` by always-on `eb_mix` + two-tier `V <= 1.0` gating and `100%` by `CARMA`.
2. **Within-SomaScan `MultiSuSiE` (`Pietzner EUR + Zhang EUR + Zhang AFR`):** Achieved `26/26 (100%)` convergence, `0` $V > 1$ blow-ups, `+79%` total CSs (`77` vs `43`), `+114%` 1-SNP CSs (`30` vs `14`), rescued sub-threshold `CXCL8` (`PPH4 = 0.8891` vs `Sun EUR`), boosted `TFF2` colocalisation (`0.5385 → 0.9998`), and produced `0/11` false positives on cis-null Somamers. Do not pool founder-isolate cohorts (`deCODE`) without in-sample LD inside `MultiSuSiE` on high-$|Z|$ loci.
3. **Edge-Case Fixes Built In:** `--min-maf 0.01` rare zero-$r^2$ filter; pre-`bcftools norm` filtering of cross-contig / `_random` liftover targets (for deCODE 49x WGS pericentromeric loci); `<= 9.6 Mb` GRCh37 pericentromeric band clamping (`REG4`, `CA9`); automatic 16-CPU / 160-GB SLURM MHC scaling with `options(matprod = "blas")`; and exact asymptotic `diff_arg > 35` overflow protection in `scripts/compare_pqtl_studies.R` when $\text{LBF}_1 + \text{LBF}_2 > 709.8$.
