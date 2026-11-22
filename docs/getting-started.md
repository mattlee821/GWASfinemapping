# Getting Started

This guide walks through installing **GWASfinemapping**, configuring storage paths for local workstations or HPC clusters, building the one-time reference bundle (~3.8 GB), and verifying the installation with the automated test and benchmark suites.

---

## System Requirements

- **OS & Shell:** macOS (`arm64` or `x86_64`; run `xcode-select --install` if Command Line Tools are missing) or Linux (`x86_64` / `aarch64`), with `bash`, `curl`, and `git`.
- **Network Access:** Outbound HTTPS access to GitHub, Ensembl, UCSC, and anonymous read access to the public Pan-UKB AWS S3 bucket (`s3://pan-ukb-us-east-1`, no AWS credentials or requester-pays account needed).
- **Disk Space:** ~20 GB free total:
  - ~8 GB for the isolated conda/Python + R runtime (`GWASFM_TOOLS_ROOT`),
  - ~3.8 GB for the installed reference data (`GWASFM_REFERENCE_ROOT`),
  - ~5–10 GB of fast temporary scratch space (`GWASFM_TEMP_ROOT` / node-local `/tmp`) while active loci are running.

---

## 1. Clone & Configure Storage Locations (`.env`)

```bash
git clone https://github.com/mattlee821/GWASfinemapping.git
cd GWASfinemapping
cp -n .env.example .env
```

By default, heavy runtime environments, reference indexes, Nextflow work directories, and temporary LD matrices live **outside** the repository (`~/.cache/gwasfinemapping`), protecting cloud-synced folders (OneDrive, iCloud, Dropbox) from multi-gigabyte file churn. On an HPC cluster, edit `.env` so compute nodes write work and scratch files to a shared high-speed `/scratch` volume:

| `.env` Variable | Purpose | Recommended HPC Setting |
|---|---|---|
| `GWASFM_WORK_ROOT` | Published fine-mapping outputs (`<study_id>/<gwas_id>/`), `manifest.html`, and `.finemap_state.tsv` | `/scratch/<user>/GWASfinemapping/finemapping` |
| `GWASFM_TOOLS_ROOT` | Isolated micromamba `env/`, compiled `renv-library/`, and pulled `.sif` container | `~/.cache/gwasfinemapping` |
| `GWASFM_REFERENCE_ROOT` | Installed GRCh37 FASTA, liftover chain, UniProt gene map, PolyFun priors, and 6 Pan-UKB Parquet indexes (~3.8 GB) | `/scratch/<user>/GWASfinemapping/reference` or `~/.cache/gwasfinemapping/reference` |
| `GWASFM_TEMP_ROOT` | Unified parent directory for `GWASFM_CACHE_DIR`, `GWASFM_NXF_WORK`, and `TMPDIR` (overridable per run via `--temp-dir`) | `/scratch/<user>/GWASfinemapping` |
| `GWASFM_SLURM_PARTITION` | Default SLURM `--partition` for worker jobs under `--slurm` | Cluster compute partition(s) |
| `GWASFM_SLURM_ACCOUNT` | Optional SLURM `--account` | Project allocation (if required) |
| `GWASFM_SLURM_EXCLUDE` | Optional comma-separated faulty SLURM nodes to exclude | Leave empty unless needed |
| `GWASFM_MAX_FORKS` | Default concurrent Nextflow tasks / SLURM queue depth | `8`–`32` on HPC; `2`–`4` on laptop |

> [!WARNING]
> Never place `GWASFM_REFERENCE_ROOT`, `GWASFM_CACHE_DIR`, or `GWASFM_NXF_WORK` inside a cloud-synced directory (OneDrive, Dropbox, iCloud) or an unmounted login-only path invisible to SLURM compute nodes.

---

## 2. Install Runtime Tools (`Local` vs `Container`)

Both execution routes use the exact same repository, reference data, manifests, and `scripts/finemap.sh` commands; choose whichever fits your machine:

=== "Local Conda + R Runtime (Default)"

    Installs OpenJDK 21, Nextflow, Python 3.12 (`hail`, `ldcov[s3]`, `lazybgen`, `pyarrow`, `pyliftover`, `MultiSuSiE`), `bcftools`, `gsl`, R 4.6, pinned `susieR` 2.0.0, pinned `CARMA`, and the compiled local `gwasfinemapping` R/C package with **no root/admin privileges required**:

    ```bash
    bash scripts/source.sh
    ```

    Targeted sub-steps can also be run individually when updating:

    ```bash
    bash scripts/source.sh --env-only        # Create .env from .env.example if absent
    bash scripts/source.sh --tools-only      # Install/update micromamba + Python packages (incl. MultiSuSiE)
    bash scripts/source.sh --r-packages      # Restore renv.lock, CARMA, and compile gwasfinemapping
    bash scripts/source.sh --force           # Rebuild conda environment from scratch
    source scripts/source.sh                 # Activate tools & R library in your interactive shell
    ```

=== "Apptainer / Singularity (HPC Cluster)"

    Pulls the pre-built image `ghcr.io/mattlee821/gwasfinemapping:latest` into `GWASFM_TOOLS_ROOT` and installs Nextflow on the host to orchestrate containerised tasks:

    ```bash
    # Load cluster apptainer module first if applicable: module load apptainer
    bash scripts/source.sh --container apptainer
    ```

=== "Docker"

    Pulls `ghcr.io/mattlee821/gwasfinemapping:latest` locally:

    ```bash
    bash scripts/source.sh --container docker
    ```

---

## 3. Install & Verify Reference Data (~3.8 GB)

A one-time setup downloads and checksums the UCSC `hg38ToHg19` chain, Ensembl GRCh37 FASTA + `.fai`, GRCh38 UniProt `cis_gene` coordinates (`11,105` unambiguous genes), `19.48M` UK Biobank PolyFun functional priors, and all six whole-population Pan-UKB GRCh37 variant indexes (`AFR`, `AMR`, `CSA`, `EAS`, `EUR`, `MID`):

=== "Parallel SLURM Build (~26 minutes)"

    ```bash
    # 1. Install fast genome, annotation, and PolyFun references on the login node (~2 min)
    bash scripts/install_reference_data.sh --only genome,annotation,polyfun

    # 2. Build all 6 Pan-UKB Parquet variant indexes in parallel on SLURM (~26 min walltime, ~2.2 GB RAM each)
    for p in AFR AMR CSA EAS EUR MID; do
      sbatch --job-name="gwasfm_${p}" scripts/install_reference_data.sbatch --only panukb --populations "${p}"
    done

    # 3. Verify integrity of all reference files once jobs finish
    bash scripts/install_reference_data.sh check
    ```

=== "Sequential Single Machine (~2.5 hours)"

    ```bash
    bash scripts/install_reference_data.sh install
    bash scripts/install_reference_data.sh check
    ```

See [Reference Data](reference-data.md) and [Pan-UKB LD Panels](panukb-ld.md) for checksum manifests and exact population sample sizes.

---

## 4. Verify Installation (`test.sh` & `benchmark.sh`)

### Automated 4-Stage Test Suite (~3–4 minutes, 100% Offline)

```bash
bash scripts/test.sh
```

Every stage must pass (`=== ALL TESTS PASSED ===`):

1. **`[1/4]` Python Unit Suite (`tests/python/`):** Liftover & `_random` contig guards, allele-swap sign invariance, Kanai 3 Mb window merging, single-region MHC policy, PyArrow regional cuts, manifest preflight, and offline `MultiSuSiE` SPD Woodbury kernel + `.rda`/`.rds` round-trips.
2. **`[2/4]` R `testthat` Suite (`tests/testthat/`):** Input validation, `p-ld` Selected-Push Sieve, `susie-rss` / `susie-rrss` two-tier validity gating, `carma`, and `polyfun`.
3. **`[3/4]` `R CMD check --no-manual`:** Full CRAN-style C/R package check (`0 errors | 0 warnings | 0 notes`).
4. **`[4/4]` End-to-End Offline Nextflow Integration:** Batch `--manifest` execution across all 5 methods, idempotent `.finemap_state.tsv` skipping, `--exclude-mhc`, deferred `render_report.R`, per-locus `/tmp` cleanup, and `coloc::coloc.susie()` variant-ID alignment.

### Reproducibility Benchmarks

```bash
# Offline zero-network synthetic benchmark (~1 min)
bash scripts/benchmark.sh --offline

# Live S3 real-data benchmark (Locke 2015 BMI + UKB-PPP NPPB P16860_OID20049, --max-loci 2)
bash scripts/benchmark.sh --gwas-root /path/to/GWASprocessing/GWAS
```
