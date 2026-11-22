#!/usr/bin/env python3
"""Multi-cohort / multi-ancestry SomaScan v4 cis-pQTL fine-mapping via MultiSuSiE.

Evaluates two SomaScan-only joint fine-mapping arms across all 26 SomaScan target aptamers
(25 unique CRC proteins, with no P <= 5e-8 pre-screening gate):
  1. somascan_multisusie_3cohort : Pietzner (UK EUR) + Zhang (ARIC EUR) + Zhang (ARIC AFR)
  2. somascan_multisusie_4cohort : Ferkingstad (Iceland EUR) + Pietzner + Zhang EUR + Zhang AFR

Each cohort is paired with its own sample size N_k, per-allele (ES, SE, AF) summary statistics,
and ancestry-matched Pan-UKB LD matrix (Pan-UKB EUR for EUR cohorts; Pan-UKB AFR for Zhang AFR),
writing coloc.susie-ready `.rda` and `qc.json` outputs compatible with `compare_pqtl_studies.R`.
"""

from __future__ import annotations

import argparse
import csv
import gc
import json
import math
import subprocess
import sys
import tempfile
import time
import types
from pathlib import Path
from typing import Any

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from prep_locus import (  # noqa: E402
    DEFAULT_REF_ROOT,
    _apply_flips_inplace,
    _assemble_bounded,
    _clamp_to_grch37_ld_band,
    _normalize_ld_matrix_inplace,
    _prune_degenerate_panukb_ld_rows,
    default_panukb_dir,
    open_blockmatrix_reader,
)
from run_finemapping import (  # noqa: E402
    chromosome_sort_key,
    ensure_regional_index,
    filter_records_by_min_maf,
    format_vcf_region,
    inspect_vcf_tags,
    match_panel_variants,
    normalize_chromosome,
    normalize_lifted_vcf,
    query_vcf,
    r_runtime_environment,
    read_gene_annotation,
    resolve_panel_matches,
    run_checked,
    write_json,
    write_lifted_vcf,
    write_tsv,
)

TARGET_PROTEINS: dict[str, str] = {
    "Q9BYZ8": "REG4", "Q03403": "TFF2", "O95994": "AGR2", "P08727": "KRT19",
    "Q07654": "TFF3", "Q99988": "GDF15", "P39900": "MMP12", "O95407": "TNFRSF6B",
    "P04155": "TFF1", "Q16790": "CA9", "P50120": "RBP2", "Q8IZP9": "ADGRG2",
    "P10145": "CXCL8", "Q9P0M4": "IL17C", "O15444": "CCL25", "O14836": "TNFRSF13B",
    "O00220": "TNFRSF10A", "P80098": "CCL7", "P01589": "IL2RA", "P02745": "C1QA",
    "Q9ULX7": "CA14", "O43508": "TNFSF12", "P01579": "IFNG", "P78556": "CCL20",
    "P10147": "CCL3", "P13500": "CCL2", "P09237": "MMP7", "P09603": "CSF1",
}

SOMASCAN_COHORTS: list[tuple[str, str, int]] = [
    ("ferkingstad_2021_34857953", "EUR", 35559),
    ("pietzner_2021_34648354", "EUR", 10708),
    ("zhang_2022_35501419_EUR", "EUR", 7213),
    ("zhang_2022_35501419_AFR", "AFR", 1871),
]


def load_multisusie_backend() -> Any:
    """Import MultiSuSiE with fast batched LAPACK compute_lbf kernels (no numba/LLVM required)."""
    try:
        import numba  # noqa: F401
    except ImportError:
        shim = types.ModuleType("numba")
        shim.jit = lambda *a, **kw: (lambda f: f)  # type: ignore[attr-defined]
        sys.modules["numba"] = shim

    try:
        import MultiSuSiE.susiepy_ss as mss
    except ImportError as exc:
        raise ImportError(
            "MultiSuSiE is not installed in the active Python environment. Install via:\n"
            "  python3 -m pip install --no-deps git+https://github.com/jordanero/MultiSuSiE.git"
        ) from exc

    def _fast_lbf_no_moments(V, XTY, X_l2_arr, rho, residual_variance, float_type=np.float64):
        YT_invD_Z = XTY / residual_variance[None, :]
        A = rho * np.sqrt(np.outer(V, V))
        inv_A = np.linalg.inv(A)
        sign_a, logdetA = np.linalg.slogdet(A)
        if sign_a <= 0:
            raise np.linalg.LinAlgError("Prior covariance A is not positive definite")
        Q_diag = (X_l2_arr / residual_variance[:, None]).T
        K = A.shape[0]
        M = np.broadcast_to(inv_A, (XTY.shape[0], K, K)).copy()
        diag_idx = np.arange(K)
        M[:, diag_idx, diag_idx] += Q_diag
        sign_m, logdet_M = np.linalg.slogdet(M)
        if np.any(sign_m <= 0):
            raise np.linalg.LinAlgError("Posterior precision M is not positive definite")
        sol = np.linalg.solve(M, YT_invD_Z[..., None])[..., 0]
        lbf = 0.5 * np.sum(YT_invD_Z * sol, axis=1) - 0.5 * (logdetA + logdet_M)
        return lbf.astype(float_type, copy=False)

    def _fast_lbf_and_moments(V, XTY, X_l2_arr, rho, residual_variance, float_type=np.float64):
        YT_invD_Z = XTY / residual_variance[None, :]
        A = rho * np.sqrt(np.outer(V, V))
        inv_A = np.linalg.inv(A)
        sign_a, logdetA = np.linalg.slogdet(A)
        if sign_a <= 0:
            raise np.linalg.LinAlgError("Prior covariance A is not positive definite")
        Q_diag = (X_l2_arr / residual_variance[:, None]).T
        K = A.shape[0]
        M = np.broadcast_to(inv_A, (XTY.shape[0], K, K)).copy()
        diag_idx = np.arange(K)
        M[:, diag_idx, diag_idx] += Q_diag
        sign_m, logdet_M = np.linalg.slogdet(M)
        if np.any(sign_m <= 0):
            raise np.linalg.LinAlgError("Posterior precision M is not positive definite")
        sol = np.linalg.solve(M, YT_invD_Z[..., None])[..., 0]
        lbf = 0.5 * np.sum(YT_invD_Z * sol, axis=1) - 0.5 * (logdetA + logdet_M)
        AQ = A[None, :, :] * Q_diag[:, None, :]
        post_mean = (YT_invD_Z @ A.T) - np.einsum("pij,pj->pi", AQ, sol)
        Minv_AQT = np.linalg.solve(M, np.swapaxes(AQ, 1, 2))
        post_cov = A[None, :, :] - np.matmul(AQ, A) + np.matmul(AQ, Minv_AQT)
        pm = post_mean.T
        pm2 = np.maximum(np.moveaxis(post_cov, 0, 2) + pm[:, None, :] * pm[None, :, :], 0.0)
        return (
            lbf.astype(float_type, copy=False),
            pm.astype(float_type, copy=False),
            pm2.astype(float_type, copy=False),
        )

    mss.compute_lbf_no_moments = _fast_lbf_no_moments
    mss.compute_lbf_and_moments = _fast_lbf_and_moments
    mss.compute_lbf_and_moments_safe = _fast_lbf_and_moments
    return mss


def build_rho_matrix(ancestries: list[str], rho_within: float = 0.95, rho_across: float = 0.75) -> np.ndarray:
    K = len(ancestries)
    rho = np.empty((K, K), dtype=np.float64)
    for i in range(K):
        for j in range(K):
            if i == j:
                rho[i, j] = 1.0
            elif ancestries[i] == ancestries[j]:
                rho[i, j] = rho_within
            else:
                rho[i, j] = rho_across
    return rho


def fit_multisusie_locus(
    cohort_inputs: list[dict[str, Any]],
    R_by_ancestry: dict[str, np.ndarray],
    L: int = 10,
    max_iter: int = 200,
    coverage: float = 0.95,
    min_abs_corr: float = 0.5,
    single_pop_mac_thresh: float = 20.0,
    rho_within: float = 0.95,
    rho_across: float = 0.75,
) -> dict[str, Any]:
    """Run MultiSuSiE in float64 with exact lbf_variable recovery and two-tier CS gating."""
    mss = load_multisusie_backend()
    K = len(cohort_inputs)
    if K == 0:
        raise ValueError("cohort_inputs cannot be empty")
    P = len(cohort_inputs[0]["b"])
    L_eff = min(int(L), P)

    ancestries = [str(c["ancestry"]) for c in cohort_inputs]
    pop_sizes = np.array([float(c["n"]) for c in cohort_inputs], dtype=np.float64)
    YTY_list = [float(n - 1.0) for n in pop_sizes]
    rho = build_rho_matrix(ancestries, rho_within=rho_within, rho_across=rho_across)

    R_clean_list: list[np.ndarray] = []
    XTX_list: list[np.ndarray] = []
    XTY_list: list[np.ndarray] = []

    for k, c in enumerate(cohort_inputs):
        anc = ancestries[k]
        b = np.asarray(c["b"], dtype=np.float64).copy()
        s = np.asarray(c["s"], dtype=np.float64).copy()
        maf = np.asarray(c["maf"], dtype=np.float64).copy()
        n_k = pop_sizes[k]

        valid_mask = np.isfinite(b) & np.isfinite(s) & (s > 0.0) & np.isfinite(maf) & (maf > 0.0) & (maf < 1.0)
        if single_pop_mac_thresh > 0:
            mac = 2.0 * n_k * np.minimum(maf, 1.0 - maf)
            valid_mask &= mac >= single_pop_mac_thresh

        b_use = np.where(valid_mask, b, np.nan)
        s_use = np.where(valid_mask, s, np.nan)

        R_corr = np.asarray(R_by_ancestry[anc], dtype=np.float64).copy()
        censor_idx = ~valid_mask
        if np.any(censor_idx):
            R_corr[censor_idx, :] = 0.0
            R_corr[:, censor_idx] = 0.0
        R_clean_list.append(R_corr)

        R_work = R_corr.copy()
        XTY = mss.recover_XTX_and_XTY(b=b_use, s=s_use, R=R_work, YTY=YTY_list[k], n=n_k)
        if np.any(censor_idx):
            R_work[censor_idx, :] = 0.0
            R_work[:, censor_idx] = 0.0
            XTY[censor_idx] = 0.0
        XTX_list.append(R_work)
        XTY_list.append(XTY)

    fit = mss.susie_multi_ss(
        XTX_list=[m.copy() for m in XTX_list],
        XTY_list=[v.copy() for v in XTY_list],
        YTY_list=YTY_list,
        rho=rho,
        population_sizes=pop_sizes,
        L=L_eff,
        scaled_prior_variance=0.2,
        standardize=False,
        pop_spec_standardization=False,
        estimate_residual_variance=True,
        estimate_prior_variance=True,
        estimate_prior_method="early_EM",
        pop_spec_effect_priors=True,
        iter_before_zeroing_effects=5,
        prior_tol=1e-9,
        max_iter=int(max_iter),
        tol=1e-3,
        verbose=False,
        R_list=R_clean_list,
        coverage=float(coverage),
        min_abs_corr=float(min_abs_corr),
        float_type=np.float64,
        low_memory_mode=False,
    )

    X_l2_arr = np.array([np.diag(xtx) for xtx in XTX_list], dtype=np.float64)
    lbf_variable = np.zeros((L_eff, P), dtype=np.float64)
    std_V_matrix = np.zeros((K, L_eff), dtype=np.float64)

    for l_idx in range(L_eff):
        V_l = np.asarray(fit.V[:, l_idx], dtype=np.float64)
        for k in range(K):
            var_x_k = X_l2_arr[k] / max(pop_sizes[k] - 1.0, 1.0)
            std_V_matrix[k, l_idx] = float(V_l[k] * np.sum(fit.alpha[l_idx] * var_x_k))
        if np.any(V_l > 1e-9):
            XTY_resid_l = [
                XTY_list[k] - fit.Xr_list[k] + XTX_list[k].dot(fit.alpha[l_idx] * fit.mu[k, l_idx])
                for k in range(K)
            ]
            lbf_variable[l_idx, :] = mss.compute_lbf(
                V=V_l,
                XTY_list=XTY_resid_l,
                XTX_list=XTX_list,
                X_l2_arr=X_l2_arr,
                rho=rho,
                residual_variance=np.asarray(fit.sigma2, dtype=np.float64),
                return_moments=False,
                float_type=np.float64,
            )

    is_converged = bool(fit.converged)
    raw_indices_list, raw_purity, raw_cov, raw_mask = fit.sets[:4]
    raw_passed = [i for i, ok in enumerate(raw_mask) if bool(ok) and len(raw_indices_list[i]) > 0]

    max_std_V_all = float(np.max(std_V_matrix)) if std_V_matrix.size else 0.0
    credible_sets: list[dict[str, Any]] = []
    dropped_cs = 0
    dropped_v_gt_1 = 0
    suppress_cs = False
    suppress_reason: str | None = None

    if not is_converged:
        dropped_cs = len(raw_passed)
        suppress_cs = dropped_cs > 0
        suppress_reason = "non_converged" if suppress_cs else None
    else:
        abs_meta_R = np.zeros((P, P), dtype=np.float64)
        for R_c in R_clean_list:
            np.maximum(abs_meta_R, np.abs(R_c), out=abs_meta_R)

        for eff_i in raw_passed:
            max_v_slot = float(np.max(std_V_matrix[:, eff_i]))
            if max_v_slot > 1.0:
                dropped_cs += 1
                dropped_v_gt_1 += 1
                continue
            members_0 = np.asarray(raw_indices_list[eff_i], dtype=int)
            # Order credible-set members by descending alpha[eff_i, :]
            members_0 = members_0[np.argsort(-fit.alpha[eff_i, members_0])]
            if len(members_0) == 1:
                min_c = mean_c = med_c = 1.0
            else:
                sub = abs_meta_R[np.ix_(members_0, members_0)]
                tri = sub[np.triu_indices(len(members_0), k=1)]
                min_c = float(np.min(tri))
                mean_c = float(np.mean(tri))
                med_c = float(np.median(tri))
            credible_sets.append({
                "credible_set_id": f"L{eff_i + 1}",
                "effect_index": int(eff_i + 1),
                "indices_1based": (members_0 + 1).tolist(),
                "size": int(len(members_0)),
                "coverage": float(raw_cov[eff_i]) if raw_cov is not None else float(np.sum(fit.alpha[eff_i, members_0])),
                "min_abs_corr": round(min_c, 6),
                "mean_abs_corr": round(mean_c, 6),
                "median_abs_corr": round(med_c, 6),
                "max_std_prior_variance_V": round(max_v_slot, 6),
                "lbf": round(float(fit.lbf[eff_i]), 6),
            })

        if dropped_v_gt_1 > 0:
            if not credible_sets:
                suppress_cs = True
                suppress_reason = "all_cs_dropped_prior_variance_gt_1"
            else:
                suppress_reason = "partial_cs_dropped_prior_variance_gt_1"

    max_cs_V = max((cs["max_std_prior_variance_V"] for cs in credible_sets), default=0.0)
    return {
        "pip": np.asarray(fit.pip, dtype=np.float64),
        "alpha": np.asarray(fit.alpha, dtype=np.float64),
        "lbf_variable": lbf_variable,
        "lbf": np.asarray(fit.lbf, dtype=np.float64),
        "V": np.asarray(fit.V, dtype=np.float64),
        "std_V": std_V_matrix,
        "sigma2": np.asarray(fit.sigma2, dtype=np.float64),
        "converged": is_converged,
        "iterations": int(fit.niter),
        "credible_sets": credible_sets,
        "diagnostics": {
            "converged": is_converged,
            "iterations": int(fit.niter),
            "cohorts_k": K,
            "cohort_ids": [c["study_id"] for c in cohort_inputs],
            "cohort_ancestries": ancestries,
            "cohort_sample_sizes": [int(round(x)) for x in pop_sizes],
            "max_prior_variance_V": round(max_std_V_all, 6),
            "max_cs_prior_variance_V": round(max_cs_V, 6),
            "implausible_prior_variance_flag": bool(max_std_V_all > 1.0),
            "credible_sets_suppressed": suppress_cs,
            "credible_sets_suppressed_reason": suppress_reason,
            "unconverged_cs_dropped": int(dropped_cs),
            "dropped_implausible_v_cs": int(dropped_v_gt_1),
            "sigma2_by_cohort": [round(float(x), 6) for x in fit.sigma2],
        },
    }


def _find_gwas_files(study_root: Path, gwas_id: str) -> tuple[Path, Path | None] | None:
    gdir = study_root / "standardise" / gwas_id
    if not gdir.is_dir():
        return None
    vcfs = sorted(gdir.glob("*.vcf.gz"))
    tbis = list(gdir.glob("*.vcf.gz.tbi")) + list(gdir.glob("*.vcf.gz.csi"))
    if not vcfs or not tbis:
        return None
    metas = sorted(gdir.glob("*.metadata.json"))
    return vcfs[0], (metas[0] if metas else None)


def _resolve_cohort_n(records: list[dict[str, Any]], meta_path: Path | None, fallback_n: int) -> float:
    ss_vals = [float(r["SS"]) for r in records if r.get("SS") is not None and float(r["SS"]) > 0]
    if ss_vals:
        return float(np.median(ss_vals))
    if meta_path and meta_path.is_file():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        for k in ("sample_size", "n_total", "effective_n"):
            if meta.get(k):
                return float(meta[k])
    return float(fallback_n)


def _extract_population_ld(
    population: str,
    chrom: str,
    normalized_rows: list[dict[str, Any]],
    cache_dir: Path,
    panukb_dir: Path,
    bm_disk_cache: str | None,
    fetch_workers: int = 8,
) -> dict[tuple[str, int, str, str], tuple[int, bool, dict[str, Any]]] | tuple[np.ndarray, dict[str, int]]:
    """Match lifted variants against Pan-UKB `population` and extract sign-aligned correlation matrix."""
    if not normalized_rows:
        return np.zeros((0, 0), dtype=np.float64), {}
    target_start = max(1, min(int(r["LD_POS"]) for r in normalized_rows) - 100)
    target_end = max(int(r["LD_POS"]) for r in normalized_rows) + 100
    try:
        index_path = ensure_regional_index(
            population, chrom, target_start, target_end, cache_dir, SCRIPT_DIR, panukb_dir
        )
    except ValueError as exc:
        if "has no variants in" in str(exc):
            return np.zeros((0, 0), dtype=np.float64), {}
        raise

    matches = match_panel_variants(index_path, normalized_rows)
    matched, _ = resolve_panel_matches(normalized_rows, matches)
    if not matched:
        return np.zeros((0, 0), dtype=np.float64), {}

    matched.sort(key=lambda r: (chromosome_sort_key(r["LD_CHROM"]), int(r["LD_POS"]), r["LD_REF"], r["LD_ALT"]))
    global_idxs = [int(r["PANEL_IDX"]) for r in matched]
    bm_uri = f"s3://pan-ukb-us-east-1/ld_release/UKBB.{population}.ldadj.bm"
    reader = open_blockmatrix_reader(bm_uri, block_cache=2, disk_cache_root=bm_disk_cache)
    matrix = _assemble_bounded(reader, global_idxs, max_workers=fetch_workers)
    del reader
    gc.collect()

    matrix, matched, _ = _prune_degenerate_panukb_ld_rows(matrix, matched)
    if not matched:
        return np.zeros((0, 0), dtype=np.float64), {}

    flip_positions = [i for i, r in enumerate(matched) if r["PANEL_ALLELE_SWAPPED"]]
    matrix = _apply_flips_inplace(matrix, flip_positions)
    ids = [f"{r['LD_CHROM']}:{int(r['LD_POS'])}:{r['LD_REF']}:{r['LD_ALT']}" for r in matched]
    _normalize_ld_matrix_inplace(matrix, ids)
    pos_by_gwas_id = {str(r["GWAS_VARIANT_ID"]): idx for idx, r in enumerate(matched)}
    return matrix, pos_by_gwas_id


def save_multisusie_study_output(
    arm_study_id: str,
    gwas_id: str,
    gene: dict[str, Any],
    region_id: str,
    source_region_grch38: str,
    variant_rows: list[dict[str, Any]],
    fit_out: dict[str, Any],
    runtime_seconds: float,
    outdir: Path,
    r_env: dict[str, str],
) -> None:
    gwas_out = outdir / arm_study_id / gwas_id
    methods_dir = gwas_out / "methods"
    methods_dir.mkdir(parents=True, exist_ok=True)

    cs_lookup = [None] * len(variant_rows)
    for cs in fit_out["credible_sets"]:
        for idx1 in cs["indices_1based"]:
            cs_lookup[idx1 - 1] = cs["credible_set_id"]

    table_rows: list[dict[str, Any]] = []
    for idx, vr in enumerate(variant_rows):
        table_rows.append({
            "VARIANT_ID": vr["VARIANT_ID"],
            "GWAS_VARIANT_ID": vr["GWAS_VARIANT_ID"],
            "CHROM": vr["LD_CHROM"],
            "POS": vr["LD_POS"],
            "REF": vr["LD_REF"],
            "ALT": vr["LD_ALT"],
            "GWAS_CHROM": vr["GWAS_CHROM"],
            "GWAS_POS": vr["GWAS_POS"],
            "GWAS_REF": vr["GWAS_REF"],
            "GWAS_ALT": vr["GWAS_ALT"],
            "MAX_LP": round(float(vr["MAX_LP"]), 5),
            "PIP": float(fit_out["pip"][idx]),
            "credible_set_id": cs_lookup[idx],
        })

    cols = [
        "VARIANT_ID", "GWAS_VARIANT_ID", "CHROM", "POS", "REF", "ALT",
        "GWAS_CHROM", "GWAS_POS", "GWAS_REF", "GWAS_ALT", "MAX_LP", "PIP", "credible_set_id",
    ]
    L_eff, P = fit_out["alpha"].shape
    K = fit_out["V"].shape[0]

    with tempfile.TemporaryDirectory(prefix="gwasfm-msusie-rda-") as tmp:
        tdir = Path(tmp)
        write_tsv(tdir / "variants.tsv", table_rows, cols)
        np.ascontiguousarray(fit_out["pip"], dtype="<f8").tofile(tdir / "pip.bin")
        np.ascontiguousarray(fit_out["alpha"], dtype="<f8").tofile(tdir / "alpha.bin")
        np.ascontiguousarray(fit_out["lbf_variable"], dtype="<f8").tofile(tdir / "lbf_variable.bin")
        np.ascontiguousarray(fit_out["V"], dtype="<f8").tofile(tdir / "V.bin")

        fit_meta = {
            "region_id": region_id,
            "L": int(L_eff),
            "K": int(K),
            "converged": fit_out["converged"],
            "iterations": fit_out["iterations"],
            "requested_coverage": 0.95,
            "lbf": [float(x) for x in fit_out["lbf"]],
            "credible_sets": fit_out["credible_sets"],
            "diagnostics": fit_out["diagnostics"],
            "provenance": {
                "method_id": "multisusie",
                "method_version": "1.0",
                "cohorts": fit_out["diagnostics"]["cohort_ids"],
                "ancestries": fit_out["diagnostics"]["cohort_ancestries"],
                "sample_sizes": fit_out["diagnostics"]["cohort_sample_sizes"],
            },
        }
        write_json(tdir / "fit_meta.json", fit_meta)
        run_checked(
            [
                "Rscript", str(SCRIPT_DIR / "save_multisusie_rda.R"),
                "--work-dir", str(tdir),
                "--output", str(methods_dir / "multisusie.rda"),
                "--compress", "xz",
            ],
            env=r_env,
            cwd=PROJECT_ROOT,
        )

    lead_row = max(variant_rows, key=lambda r: (float(r["MAX_LP"]), -int(r["GWAS_POS"])))
    cs_sizes = [int(cs["size"]) for cs in fit_out["credible_sets"]]
    diag = fit_out["diagnostics"]
    ancs = sorted(set(diag["cohort_ancestries"]))
    anc_label = "+".join(ancs)

    method_qc = {
        "method_id": "multisusie",
        "credible_sets": len(fit_out["credible_sets"]),
        "credible_set_sizes": cs_sizes,
        "median_cs_size": float(np.median(cs_sizes)) if cs_sizes else None,
        "single_variant_cs_count": sum(1 for s in cs_sizes if s == 1),
        "cs_le_5_count": sum(1 for s in cs_sizes if s <= 5),
        "max_pip": round(float(np.max(fit_out["pip"])) if P > 0 else 0.0, 6),
        **diag,
    }
    overall_status = (
        "PASS"
        if fit_out["converged"] and not diag["credible_sets_suppressed"] and not diag["implausible_prior_variance_flag"]
        else "WARN_DIAGNOSTIC_FLAGS"
    )
    qc_payload = {
        "study_id": arm_study_id,
        "gwas_id": gwas_id,
        "overall_status": overall_status,
        "analyzed_regions_count": 1,
        "cis_gene": {
            "status": "analyzed_multisusie_cis",
            "uniprot_id": gene["uniprot_id"],
            "gene_symbol": gene["gene_symbol"],
            "region_id": region_id,
            "max_lp": round(float(lead_row["MAX_LP"]), 4),
        },
        "loci": [
            {
                "region_id": region_id,
                "region_type": "cis_gene",
                "source_region_grch38": source_region_grch38,
                "status": "analyzed",
                "variants_ld_matched": P,
                "variants_removed_min_maf": 0,
                "lead_gwas_variant_id": lead_row["GWAS_VARIANT_ID"],
                "lead_lp": round(float(lead_row["MAX_LP"]), 4),
                "methods": {
                    "multisusie": method_qc,
                },
            }
        ],
    }
    prov_payload = {
        "study_id": arm_study_id,
        "gwas_id": gwas_id,
        "gwas_ancestry": anc_label,
        "ld_population": anc_label,
        "gene_annotation": gene,
        "cohorts_combined": diag["cohort_ids"],
        "cohort_sample_sizes": diag["cohort_sample_sizes"],
        "total_sample_size": int(sum(diag["cohort_sample_sizes"])),
        "fine_mapping_methods": ["multisusie"],
        "runtime_seconds": round(runtime_seconds, 2),
    }
    write_json(gwas_out / "qc.json", qc_payload)
    write_json(gwas_out / "provenance.json", prov_payload)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gwas-root", required=True, help="Path to GWASprocessing/GWAS directory")
    p.add_argument("--outdir", required=True, help="Path to finemapping output root")
    p.add_argument("--chain", default=str(DEFAULT_REF_ROOT / "genome/hg38ToHg19.over.chain.gz"))
    p.add_argument("--grch37-fasta", default=str(DEFAULT_REF_ROOT / "genome/GRCh37.fa"))
    p.add_argument("--gene-annotation", default=str(DEFAULT_REF_ROOT / "annotation/gene_regions_GRCh38.tsv"))
    p.add_argument("--panukb-dir", default=None)
    p.add_argument("--cache-dir", default=None, help="Directory for regional Pan-UKB index slices")
    p.add_argument("--bm-disk-cache", default=None, help="Shared Pan-UKB BlockMatrix disk cache")
    p.add_argument("--min-maf", type=float, default=0.01)
    p.add_argument("--max-iter", type=int, default=200)
    p.add_argument("--L", type=int, default=10)
    p.add_argument("--fetch-workers", type=int, default=8)
    p.add_argument("--only-assay", default=None, help="Optional comma-separated gwas_id filter")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    load_multisusie_backend()

    gwas_root = Path(args.gwas_root).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    chain = Path(args.chain).expanduser().resolve()
    fasta = Path(args.grch37_fasta).expanduser().resolve()
    gene_annot = Path(args.gene_annotation).expanduser().resolve()
    panukb_dir = Path(args.panukb_dir).expanduser().resolve() if args.panukb_dir else default_panukb_dir()
    cache_dir = (
        Path(args.cache_dir).expanduser().resolve()
        if args.cache_dir
        else (outdir / ".cache_multisusie").resolve()
    )
    cache_dir.mkdir(parents=True, exist_ok=True)
    r_env = r_runtime_environment(PROJECT_ROOT)

    only_set = (
        {x.strip() for x in args.only_assay.split(",") if x.strip()}
        if args.only_assay
        else None
    )

    # Discover all SomaScan aptamer IDs across the 4 SomaScan cohorts for the 28 target proteins
    aptamer_cohorts: dict[str, dict[str, tuple[Path, Path | None, str, int]]] = {}
    for study_id, anc, default_n in SOMASCAN_COHORTS:
        sroot = gwas_root / study_id
        std = sroot / "standardise"
        if not std.is_dir():
            continue
        for u in TARGET_PROTEINS:
            for gdir in sorted(p for p in std.glob(f"{u}_*") if p.is_dir()):
                gid = gdir.name
                if only_set and gid not in only_set:
                    continue
                found = _find_gwas_files(sroot, gid)
                if found:
                    aptamer_cohorts.setdefault(gid, {})[study_id] = (found[0], found[1], anc, default_n)

    print(f"[MultiSuSiE] Discovered {len(aptamer_cohorts)} SomaScan target aptamer groups across {gwas_root}")
    from pyliftover import LiftOver
    lifter = LiftOver(str(chain))

    summary_log: list[dict[str, Any]] = []
    for idx_a, gwas_id in enumerate(sorted(aptamer_cohorts), start=1):
        t_assay0 = time.monotonic()
        uniprot = gwas_id.split("_", 1)[0]
        gene = read_gene_annotation(gene_annot, uniprot, strict=True)
        assert gene is not None
        chrom = normalize_chromosome(gene["chrom"])
        nom_start = max(1, int(gene["start"]) - 1_000_000)
        nom_end = int(gene["end"]) + 1_000_000

        c_map = aptamer_cohorts[gwas_id]
        cohort_records: dict[str, dict[tuple[int, str, str], dict[str, Any]]] = {}
        cohort_bounds: dict[str, tuple[int, int]] = {}
        cohort_n: dict[str, float] = {}
        canonical_alleles: dict[tuple[int, frozenset[str]], tuple[str, str]] = {}

        for study_id, (vcf_p, meta_p, anc, def_n) in c_map.items():
            vtags = inspect_vcf_tags(vcf_p)
            raw = query_vcf(vcf_p, region=format_vcf_region(vtags, chrom, nom_start, nom_end), vcf_tags=vtags)
            recs, _ = filter_records_by_min_maf(raw, float(args.min_maf))
            if not recs:
                continue
            cohort_n[study_id] = _resolve_cohort_n(recs, meta_p, def_n)
            positions = [int(r["POS"]) for r in recs]
            cohort_bounds[study_id] = (min(positions), max(positions))
            rec_by_var: dict[tuple[int, str, str], dict[str, Any]] = {}
            for r in recs:
                if r.get("ES") is None or r.get("SE") is None or float(r["SE"]) <= 0:
                    continue
                pos_i = int(r["POS"])
                ref_u, alt_u = str(r["REF"]).upper(), str(r["ALT"]).upper()
                aset = frozenset((ref_u, alt_u))
                if len(aset) != 2:
                    continue
                c_ref, c_alt = canonical_alleles.setdefault((pos_i, aset), (ref_u, alt_u))
                if (ref_u, alt_u) == (c_ref, c_alt):
                    es_val = float(r["ES"])
                    af_val = float(r["AF"]) if r.get("AF") is not None else 0.5
                elif (ref_u, alt_u) == (c_alt, c_ref):
                    es_val = -float(r["ES"])
                    af_val = (1.0 - float(r["AF"])) if r.get("AF") is not None else 0.5
                else:
                    continue
                rec_by_var[(pos_i, c_ref, c_alt)] = {
                    "CHROM": chrom,
                    "POS": pos_i,
                    "REF": c_ref,
                    "ALT": c_alt,
                    "ES": es_val,
                    "SE": float(r["SE"]),
                    "LP": float(r["LP"] or 0.0),
                    "SS": cohort_n[study_id],
                    "AF": af_val,
                }
            if rec_by_var:
                cohort_records[study_id] = rec_by_var

        if not cohort_records:
            print(f"  [{idx_a:02d}/{len(aptamer_cohorts)}] {gwas_id} ({gene['gene_symbol']}): no cis variants")
            continue

        # Shared cis interval covered by all available cohorts for this aptamer (avoids ±500 kb edge steps)
        shared_start = max(b[0] for b in cohort_bounds.values())
        shared_end = min(b[1] for b in cohort_bounds.values())
        if shared_end <= shared_start:
            shared_start, shared_end = nom_start, nom_end

        union_keys = sorted({
            k for rmap in cohort_records.values() for k in rmap if shared_start <= k[0] <= shared_end
        })
        union_source = [
            {
                "CHROM": chrom,
                "POS": pos_i,
                "REF": ref_u,
                "ALT": alt_u,
                "ES": 0.0,
                "SE": 1.0,
                "LP": max(
                    (rmap[(pos_i, ref_u, alt_u)]["LP"] for rmap in cohort_records.values() if (pos_i, ref_u, alt_u) in rmap),
                    default=0.0,
                ),
                "SS": 10000.0,
            }
            for (pos_i, ref_u, alt_u) in union_keys
        ]

        with tempfile.TemporaryDirectory(prefix=f"gwasfm-mslift-{uniprot}-") as ltmp:
            lvcf = Path(ltmp) / "lifted.vcf"
            lifted_by_id, _ = write_lifted_vcf(union_source, lifter, lvcf, expected_chrom=chrom)
            norm_rows = normalize_lifted_vcf(lvcf, fasta, None) if lifted_by_id else []

        for row in norm_rows:
            src = lifted_by_id[row["source_variant_key"]]
            row["CHROM"] = src["CHROM"]
            row["POS"] = src["POS"]
            row["REF"] = src["REF"]
            row["ALT"] = src["ALT"]
            row["LP"] = src["LP"]
            row["GWAS_VARIANT_ID"] = f"{src['CHROM']}:{src['POS']}:{src['REF']}:{src['ALT']}"

        norm_rows = [r for r in norm_rows if r["LD_CHROM"] == chrom]
        norm_rows, _, _ = _clamp_to_grch37_ld_band(
            norm_rows, {"region_type": "cis_gene"}, gene=gene, max_span_bp=9_600_000
        )
        if not norm_rows:
            print(f"  [{idx_a:02d}/{len(aptamer_cohorts)}] {gwas_id}: 0 lifted cis variants")
            continue

        R_eur, map_eur = _extract_population_ld(
            "EUR", chrom, norm_rows, cache_dir, panukb_dir, args.bm_disk_cache, args.fetch_workers
        )
        has_afr = "zhang_2022_35501419_AFR" in cohort_records
        if has_afr:
            R_afr, map_afr = _extract_population_ld(
                "AFR", chrom, norm_rows, cache_dir, panukb_dir, args.bm_disk_cache, args.fetch_workers
            )
        else:
            R_afr, map_afr = np.zeros((0, 0), dtype=np.float64), {}

        norm_by_gid = {str(r["GWAS_VARIANT_ID"]): r for r in norm_rows}
        region_id = f"cis_{gene['gene_symbol']}_chr{chrom}_{shared_start}_{shared_end}"
        src_reg = f"{chrom}:{shared_start}-{shared_end}"

        arms = [
            (
                "somascan_multisusie_3cohort",
                [c for c in ("pietzner_2021_34648354", "zhang_2022_35501419_EUR", "zhang_2022_35501419_AFR") if c in cohort_records],
            ),
            (
                "somascan_multisusie_4cohort",
                [c for c in ("ferkingstad_2021_34857953", "pietzner_2021_34648354", "zhang_2022_35501419_EUR", "zhang_2022_35501419_AFR") if c in cohort_records],
            ),
        ]

        for arm_id, arm_cohorts in arms:
            if not arm_cohorts:
                continue
            t_arm0 = time.monotonic()
            eur_cohorts = [c for c in arm_cohorts if not c.endswith("_AFR")]
            afr_cohorts = [c for c in arm_cohorts if c.endswith("_AFR")]

            selected_gids: list[str] = []
            for r in sorted(norm_rows, key=lambda x: (int(x["LD_POS"]), x["LD_REF"], x["LD_ALT"])):
                gid_v = str(r["GWAS_VARIANT_ID"])
                key3 = (int(r["POS"]), str(r["REF"]), str(r["ALT"]))
                in_all_eur = bool(eur_cohorts) and (gid_v in map_eur) and all(
                    key3 in cohort_records[c] for c in eur_cohorts
                )
                in_afr_only = (
                    bool(afr_cohorts)
                    and (gid_v in map_afr)
                    and all(key3 in cohort_records[c] for c in afr_cohorts)
                    and not any(key3 in cohort_records[c] for c in eur_cohorts)
                )
                if in_all_eur or in_afr_only:
                    selected_gids.append(gid_v)

            if len(selected_gids) < 5:
                continue

            P_arm = len(selected_gids)
            R_arm_by_anc: dict[str, np.ndarray] = {}
            if eur_cohorts:
                R_sub_e = np.zeros((P_arm, P_arm), dtype=np.float64)
                e_valid = [i for i, g in enumerate(selected_gids) if g in map_eur]
                e_src = [map_eur[selected_gids[i]] for i in e_valid]
                if e_valid:
                    R_sub_e[np.ix_(e_valid, e_valid)] = R_eur[np.ix_(e_src, e_src)]
                R_arm_by_anc["EUR"] = R_sub_e
            if afr_cohorts:
                R_sub_a = np.zeros((P_arm, P_arm), dtype=np.float64)
                a_valid = [i for i, g in enumerate(selected_gids) if g in map_afr]
                a_src = [map_afr[selected_gids[i]] for i in a_valid]
                if a_valid:
                    R_sub_a[np.ix_(a_valid, a_valid)] = R_afr[np.ix_(a_src, a_src)]
                R_arm_by_anc["AFR"] = R_sub_a

            c_inputs: list[dict[str, Any]] = []
            for cid in arm_cohorts:
                anc = "AFR" if cid.endswith("_AFR") else "EUR"
                pmap = map_afr if anc == "AFR" else map_eur
                b_arr = np.full(P_arm, np.nan, dtype=np.float64)
                s_arr = np.full(P_arm, np.nan, dtype=np.float64)
                maf_arr = np.full(P_arm, np.nan, dtype=np.float64)
                rmap = cohort_records[cid]
                for i, gid_v in enumerate(selected_gids):
                    if gid_v not in pmap:
                        continue
                    nr = norm_by_gid[gid_v]
                    key3 = (int(nr["POS"]), str(nr["REF"]), str(nr["ALT"]))
                    rec = rmap.get(key3)
                    if rec is not None:
                        b_arr[i] = rec["ES"]
                        s_arr[i] = rec["SE"]
                        maf_arr[i] = min(rec["AF"], 1.0 - rec["AF"])
                c_inputs.append({
                    "study_id": cid,
                    "ancestry": anc,
                    "n": cohort_n[cid],
                    "b": b_arr,
                    "s": s_arr,
                    "maf": maf_arr,
                })

            variant_meta_rows = []
            for i, gid_v in enumerate(selected_gids):
                nr = norm_by_gid[gid_v]
                key3 = (int(nr["POS"]), str(nr["REF"]), str(nr["ALT"]))
                max_lp_v = max(
                    (cohort_records[c][key3]["LP"] for c in arm_cohorts if key3 in cohort_records[c]),
                    default=0.0,
                )
                variant_meta_rows.append({
                    "VARIANT_ID": f"{nr['LD_CHROM']}:{int(nr['LD_POS'])}:{nr['LD_REF']}:{nr['LD_ALT']}",
                    "GWAS_VARIANT_ID": gid_v,
                    "LD_CHROM": nr["LD_CHROM"],
                    "LD_POS": int(nr["LD_POS"]),
                    "LD_REF": nr["LD_REF"],
                    "LD_ALT": nr["LD_ALT"],
                    "GWAS_CHROM": nr["CHROM"],
                    "GWAS_POS": int(nr["POS"]),
                    "GWAS_REF": nr["REF"],
                    "GWAS_ALT": nr["ALT"],
                    "MAX_LP": max_lp_v,
                })

            fit_res = fit_multisusie_locus(
                cohort_inputs=c_inputs,
                R_by_ancestry=R_arm_by_anc,
                L=int(args.L),
                max_iter=int(args.max_iter),
                coverage=0.95,
                min_abs_corr=0.5,
            )
            arm_sec = time.monotonic() - t_arm0
            save_multisusie_study_output(
                arm_study_id=arm_id,
                gwas_id=gwas_id,
                gene=gene,
                region_id=region_id,
                source_region_grch38=src_reg,
                variant_rows=variant_meta_rows,
                fit_out=fit_res,
                runtime_seconds=arm_sec,
                outdir=outdir,
                r_env=r_env,
            )
            cs_n = len(fit_res["credible_sets"])
            cs_sz = [c["size"] for c in fit_res["credible_sets"]]
            print(
                f"  [{idx_a:02d}/{len(aptamer_cohorts)}] {arm_id:28s} {gwas_id:17s} ({gene['gene_symbol']:9s}) | "
                f"K={len(arm_cohorts)} M={P_arm:5d} conv={fit_res['converged']}({fit_res['iterations']}it) "
                f"CS={cs_n} sizes={cs_sz} ({arm_sec:.1f}s)"
            )
            summary_log.append({
                "arm": arm_id,
                "gwas_id": gwas_id,
                "gene_symbol": gene["gene_symbol"],
                "K": len(arm_cohorts),
                "M": P_arm,
                "converged": fit_res["converged"],
                "iterations": fit_res["iterations"],
                "credible_sets": cs_n,
                "median_cs_size": float(np.median(cs_sz)) if cs_sz else None,
            })

        del R_eur, R_afr
        gc.collect()
        _ = t_assay0

    print(f"[MultiSuSiE] Completed {len(summary_log)} multi-cohort fits -> {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
