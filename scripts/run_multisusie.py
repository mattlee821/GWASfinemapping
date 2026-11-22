#!/usr/bin/env python3
"""Opt-in multi-cohort / multi-ancestry cis fine-mapping via MultiSuSiE.

Provides two Nextflow/CLI subcommands (paralleling prepare_manifest.py + locus_worker.py):
  - `prepare`   : Discover & validate multi-cohort groups from `--multisusie-manifest <tsv>`
                  or `--multisusie-studies <dir1,dir2,...>`, check `.finemap_state.tsv`
                  idempotency, and emit per-group JSON task specs (`rows/row_*.json`).
  - `run-group` : Execute 64-bit MultiSuSiE (unconditionally SPD Woodbury LBF kernel +
                  sigma2 >= 1e-2 floor + two-tier CS validity gate) for one target group
                  and write `{metadata,provenance,qc}.json`, `methods/multisusie.{rda,rds}`,
                  and `report.html`.
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import html
import json
import re
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
from prepare_manifest import (  # noqa: E402
    atomic_write_json,
    load_state,
)
from run_finemapping import (  # noqa: E402
    REFERENCE_SAMPLES,
    chromosome_sort_key,
    ensure_regional_index,
    extract_ld_from_local_panel,
    filter_records_by_min_maf,
    format_vcf_region,
    inspect_vcf_tags,
    match_panel_variants,
    matrix_to_npy,
    normalize_chromosome,
    normalize_lifted_vcf,
    query_vcf,
    r_runtime_environment,
    resolve_panel_matches,
    resolve_target_gene,
    run_checked,
    write_json,
    write_lifted_vcf,
    write_tsv,
)


def parse_bool(val: Any, default: bool = False) -> bool:
    if val is None:
        return default
    if isinstance(val, bool):
        return val
    return str(val).strip().lower() in ("1", "true", "yes", "y")


def load_multisusie_backend() -> Any:
    """Import MultiSuSiE with unconditionally SPD Woodbury LBF kernels and sigma2 >= 1e-2 floor."""
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
            "MultiSuSiE is not installed in the active Python environment. Run `bash scripts/source.sh` or:\n"
            "  python3 -m pip install --no-deps https://github.com/jordanero/MultiSuSiE/archive/refs/heads/main.zip"
        ) from exc

    def _fast_lbf_no_moments(V, XTY, X_l2_arr, rho, residual_variance, float_type=np.float64):
        rv = np.maximum(np.asarray(residual_variance, dtype=np.float64), 1e-2)
        V64 = np.maximum(np.asarray(V, dtype=np.float64), 0.0)
        P, K = XTY.shape
        z = np.asarray(XTY, dtype=np.float64) / rv[None, :]
        q = np.maximum((np.asarray(X_l2_arr, dtype=np.float64) / rv[:, None]).T, 0.0)
        s = np.sqrt(q)
        A = np.asarray(rho, dtype=np.float64) * np.sqrt(np.outer(V64, V64))
        B = (s[:, :, None] * s[:, None, :]) * A[None, :, :]
        diag_idx = np.arange(K)
        B[:, diag_idx, diag_idx] += 1.0
        _, logdet_B = np.linalg.slogdet(B)
        w = z @ A.T
        v = np.linalg.solve(B, (s * w)[..., None])[..., 0]
        u = w - ((s * v) @ A.T)
        lbf = 0.5 * np.sum(z * u, axis=1) - 0.5 * logdet_B
        return lbf.astype(float_type, copy=False)

    def _fast_lbf_and_moments(V, XTY, X_l2_arr, rho, residual_variance, float_type=np.float64):
        rv = np.maximum(np.asarray(residual_variance, dtype=np.float64), 1e-2)
        V64 = np.maximum(np.asarray(V, dtype=np.float64), 0.0)
        P, K = XTY.shape
        z = np.asarray(XTY, dtype=np.float64) / rv[None, :]
        q = np.maximum((np.asarray(X_l2_arr, dtype=np.float64) / rv[:, None]).T, 0.0)
        s = np.sqrt(q)
        A = np.asarray(rho, dtype=np.float64) * np.sqrt(np.outer(V64, V64))
        B = (s[:, :, None] * s[:, None, :]) * A[None, :, :]
        diag_idx = np.arange(K)
        B[:, diag_idx, diag_idx] += 1.0
        _, logdet_B = np.linalg.slogdet(B)
        w = z @ A.T
        v = np.linalg.solve(B, (s * w)[..., None])[..., 0]
        u = w - ((s * v) @ A.T)
        lbf = 0.5 * np.sum(z * u, axis=1) - 0.5 * logdet_B
        AS = A[None, :, :] * s[:, None, :]
        Binv_AST = np.linalg.solve(B, np.swapaxes(AS, 1, 2))
        post_cov = A[None, :, :] - np.matmul(AS, Binv_AST)
        pm = post_mean = u.T
        pm2 = np.maximum(np.moveaxis(post_cov, 0, 2) + pm[:, None, :] * pm[None, :, :], 0.0)
        _ = post_mean
        return (
            lbf.astype(float_type, copy=False),
            pm.astype(float_type, copy=False),
            pm2.astype(float_type, copy=False),
        )

    def _safe_estimate_residual_variance(
        XTX_list, XTY_list, YTY_list, s_obj, X_l2_arr, population_sizes, float_type=np.float64
    ):
        sigma2_arr = np.asarray(s_obj.ER2, dtype=np.float64) / np.asarray(population_sizes, dtype=np.float64)
        return np.maximum(sigma2_arr, 1e-2).astype(float_type, copy=False)

    mss.compute_lbf_no_moments = _fast_lbf_no_moments
    mss.compute_lbf_and_moments = _fast_lbf_and_moments
    mss.compute_lbf_and_moments_safe = _fast_lbf_and_moments
    mss.estimate_residual_variance_func = _safe_estimate_residual_variance
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

    ancestries = [str(c.get("ld_population") or c["ancestry"]) for c in cohort_inputs]
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


def _render_multisusie_report_html(
    arm_study_id: str,
    gwas_id: str,
    gene: dict[str, Any],
    region_id: str,
    variant_rows: list[dict[str, Any]],
    fit_out: dict[str, Any],
    runtime_seconds: float,
) -> str:
    diag = fit_out["diagnostics"]
    cs_rows = []
    for cs in fit_out["credible_sets"]:
        lead_i = cs["indices_1based"][0] - 1
        lead_v = variant_rows[lead_i]["GWAS_VARIANT_ID"]
        lead_pip = fit_out["pip"][lead_i]
        cs_rows.append(
            f"<tr><td><code>{html.escape(cs['credible_set_id'])}</code></td>"
            f"<td class='num'>{cs['size']}</td>"
            f"<td class='num'>{cs['coverage']:.3f}</td>"
            f"<td class='num'>{cs['min_abs_corr']:.3f}</td>"
            f"<td class='num'>{cs['max_std_prior_variance_V']:.4f}</td>"
            f"<td class='num'>{cs['lbf']:.2f}</td>"
            f"<td><code>{html.escape(lead_v)}</code> (PIP={lead_pip:.4f})</td></tr>"
        )
    cs_body = "\n".join(cs_rows) or "<tr><td colspan='7' class='muted'>No 95% credible sets retained.</td></tr>"

    top_idx = np.argsort(-fit_out["pip"])[:15]
    v_rows = []
    for i in top_idx:
        vr = variant_rows[int(i)]
        v_rows.append(
            f"<tr><td><code>{html.escape(vr['GWAS_VARIANT_ID'])}</code></td>"
            f"<td><code>{html.escape(vr['VARIANT_ID'])}</code></td>"
            f"<td class='num'>{float(vr['MAX_LP']):.2f}</td>"
            f"<td class='num'>{float(fit_out['pip'][i]):.4f}</td></tr>"
        )
    cohorts_str = ", ".join(
        f"{cid} ({anc}, N={n:,})"
        for cid, anc, n in zip(diag["cohort_ids"], diag["cohort_ancestries"], diag["cohort_sample_sizes"])
    )
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"/>
<title>MultiSuSiE Report: {html.escape(arm_study_id)}/{html.escape(gwas_id)}</title>
<style>
:root {{ --bg:#f8fafc; --card:#fff; --fg:#0f172a; --muted:#64748b; --border:#e2e8f0; }}
@media (prefers-color-scheme: dark) {{ :root {{ --bg:#0f172a; --card:#1e293b; --fg:#f1f5f9; --muted:#94a3b8; --border:#334155; }} }}
body {{ margin:0; padding:24px; background:var(--bg); color:var(--fg); font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; font-size:14px; }}
.container {{ max-width:1100px; margin:0 auto; }}
.panel {{ background:var(--card); border:1px solid var(--border); border-radius:8px; padding:16px; margin-bottom:16px; }}
table {{ width:100%; border-collapse:collapse; font-size:13px; }}
th, td {{ padding:7px 10px; border-bottom:1px solid var(--border); text-align:left; }}
.num {{ text-align:right; font-variant-numeric:tabular-nums; }}
.muted {{ color:var(--muted); }}
</style></head><body><div class="container">
<h1>MultiSuSiE Joint Fine-Mapping: <code>{html.escape(gwas_id)}</code> ({html.escape(gene['gene_symbol'])})</h1>
<p class="muted">Group: <code>{html.escape(arm_study_id)}</code> &middot; Region: <code>{html.escape(region_id)}</code> &middot; Runtime: {runtime_seconds:.1f}s</p>
<div class="panel">
  <strong>Cohorts ({diag['cohorts_k']}):</strong> {html.escape(cohorts_str)}<br/>
  <strong>Variants matched:</strong> {len(variant_rows):,} &middot;
  <strong>Converged:</strong> {fit_out['converged']} ({fit_out['iterations']} iterations) &middot;
  <strong>Max standardized V:</strong> {diag['max_prior_variance_V']:.4f} &middot;
  <strong>Dropped CSs:</strong> {diag['unconverged_cs_dropped']}
</div>
<h2>95% Credible Sets ({len(fit_out['credible_sets'])})</h2>
<div class="panel"><table><thead><tr><th>CS</th><th class="num">Size</th><th class="num">Coverage</th><th class="num">Min |r|</th><th class="num">Max Std V</th><th class="num">LBF</th><th>Top Variant</th></tr></thead><tbody>{cs_body}</tbody></table></div>
<h2>Top 15 Variants by PIP</h2>
<div class="panel"><table><thead><tr><th>GRCh38 Variant</th><th>GRCh37 Panel ID</th><th class="num">Max LP</th><th class="num">MultiSuSiE PIP</th></tr></thead><tbody>{"".join(v_rows)}</tbody></table></div>
</div></body></html>
"""


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
    output_format: str = "rda",
    compress: str = "xz",
    write_report: bool = True,
) -> Path:
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
                "--format", output_format,
                "--compress", compress,
            ],
            env=r_env,
            cwd=PROJECT_ROOT,
        )

    lead_row = max(variant_rows, key=lambda r: (float(r["MAX_LP"]), -int(r["GWAS_POS"])))
    cs_sizes = [int(cs["size"]) for cs in fit_out["credible_sets"]]
    diag = fit_out["diagnostics"]
    ancs = sorted(set(diag["cohort_ancestries"]))
    anc_label = "+".join(ancs)
    n_cs = len(fit_out["credible_sets"])

    method_qc = {
        "method_id": "multisusie",
        "credible_sets": n_cs,
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
        "genomewide_significant_variants": sum(1 for r in variant_rows if float(r["MAX_LP"]) >= 7.30103),
        "analyzed_regions_count": 1,
        "total_credible_sets_multisusie": n_cs,
        "nonconverged_regions": [] if fit_out["converged"] else [region_id],
        "mismatch_flagged_regions": [region_id] if diag["implausible_prior_variance_flag"] else [],
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
        "trait": f"{gene['gene_symbol']} ({gene['uniprot_id']}) multi-cohort cis-pQTL",
        "gwas_ancestry": anc_label,
        "ld_population": anc_label,
        "gene_annotation": gene,
        "cohorts_combined": diag["cohort_ids"],
        "cohort_sample_sizes": diag["cohort_sample_sizes"],
        "total_sample_size": int(sum(diag["cohort_sample_sizes"])),
        "fine_mapping_methods": ["multisusie"],
        "runtime_seconds": round(runtime_seconds, 2),
    }
    meta_payload = {
        "study_id": arm_study_id,
        "output_id": gwas_id,
        "uniprot_id": gene["uniprot_id"],
        "gene_symbol": gene["gene_symbol"],
        "phenotype": f"{gene['gene_symbol']} ({gwas_id}) MultiSuSiE ({len(diag['cohort_ids'])} cohorts)",
        "ancestry": anc_label,
        "sample_size": int(sum(diag["cohort_sample_sizes"])),
        "cohorts": diag["cohort_ids"],
    }
    write_json(gwas_out / "qc.json", qc_payload)
    write_json(gwas_out / "provenance.json", prov_payload)
    write_json(gwas_out / "metadata.json", meta_payload)
    if write_report:
        (gwas_out / "report.html").write_text(
            _render_multisusie_report_html(
                arm_study_id, gwas_id, gene, region_id, variant_rows, fit_out, runtime_seconds
            ),
            encoding="utf-8",
        )
    return gwas_out


def _infer_cohort_ancestry_and_ld(study_name: str, meta: dict[str, Any], default_pop: str) -> tuple[str, str]:
    anc = str(meta.get("ancestry") or "").strip().upper()
    suffix_match = re.search(r"_(AFR|AMR|CSA|EAS|EUR|MID|ALL)$", study_name.upper())
    if suffix_match:
        anc = suffix_match.group(1)
    if not anc:
        anc = default_pop.upper()
    ld_pop = anc if anc in REFERENCE_SAMPLES else default_pop.upper()
    return anc, ld_pop


def _resolve_gwas_entry(raw_dir_or_vcf: str, default_study_id: str, default_pop: str) -> dict[str, Any] | None:
    p = Path(raw_dir_or_vcf).expanduser().resolve()
    if p.is_file() and str(p).endswith(".vcf.gz"):
        gdir = p.parent
        vcfs = [p]
    elif p.is_dir():
        gdir = p
        vcfs = sorted(gdir.glob("*.vcf.gz"))
    else:
        return None
    tbis = list(gdir.glob("*.vcf.gz.tbi")) + list(gdir.glob("*.vcf.gz.csi"))
    if not vcfs or not tbis:
        return None
    metas = sorted(gdir.glob("*.metadata.json"))
    meta = json.loads(metas[0].read_text(encoding="utf-8")) if metas else {}
    anc, ld_pop = _infer_cohort_ancestry_and_ld(default_study_id, meta, default_pop)
    return {
        "cohort_study_id": default_study_id,
        "gwas_id": gdir.name,
        "vcf": str(vcfs[0]),
        "metadata_path": str(metas[0]) if metas else "",
        "metadata": meta,
        "ancestry": anc,
        "ld_population": ld_pop,
    }


def _extract_ld_for_group(
    population: str,
    chrom: str,
    normalized_rows: list[dict[str, Any]],
    cache_dir: Path,
    panukb_dir: Path | None,
    local_ld_panel: Path | None,
    bm_path: str | None,
    bm_disk_cache: str | None,
    fetch_workers: int = 4,
) -> tuple[np.ndarray, dict[str, int]]:
    if not normalized_rows:
        return np.zeros((0, 0), dtype=np.float64), {}

    target_start = max(1, min(int(r["LD_POS"]) for r in normalized_rows) - 100)
    target_end = max(int(r["LD_POS"]) for r in normalized_rows) + 100
    pdir = panukb_dir if panukb_dir else default_panukb_dir()
    try:
        index_path = ensure_regional_index(
            population, chrom, target_start, target_end, cache_dir, SCRIPT_DIR, pdir
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
    if local_ld_panel is not None and local_ld_panel.is_file():
        with tempfile.TemporaryDirectory(prefix="gwasfm-mslocal-") as tmp:
            tdir = Path(tmp)
            ld_prefix = tdir / f"local_{population.lower()}"
            extract_ld_from_local_panel(local_ld_panel, matched, ld_prefix)
            npy_p = tdir / "matrix.npy"
            matrix_to_npy(Path(f"{ld_prefix}.npz"), Path(f"{ld_prefix}.variants.tsv"), npy_p)
            mat = np.array(np.load(npy_p), dtype=np.float64, copy=True)
            pos_map = {str(r["GWAS_VARIANT_ID"]): idx for idx, r in enumerate(matched)}
            return mat, pos_map

    global_idxs = [int(r["PANEL_IDX"]) for r in matched]
    bm_uri = bm_path if bm_path else f"s3://pan-ukb-us-east-1/ld_release/UKBB.{population}.ldadj.bm"
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


def cmd_prepare(args: argparse.Namespace) -> int:
    outdir = Path(args.outdir).expanduser().resolve()
    state_file = Path(args.state_file).expanduser().resolve()
    rows_dir = Path(args.rows_dir).expanduser().resolve()
    rows_dir.mkdir(parents=True, exist_ok=True)
    gene_annot = Path(args.gene_annotation).expanduser().resolve() if args.gene_annotation else None
    dry_run = parse_bool(args.dry_run, False)
    force = parse_bool(args.force, False)
    report = parse_bool(args.report, True)

    only_filter = (
        {t.strip() for t in re.split(r"[\s,]+", args.only_gwas) if t.strip()}
        if args.only_gwas
        else None
    )

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    invalid: list[dict[str, Any]] = []
    source_desc = ""

    if args.multisusie_manifest:
        mpath = Path(args.multisusie_manifest).expanduser().resolve()
        source_desc = str(mpath)
        with mpath.open(encoding="utf-8", newline="") as fh:
            for raw in csv.DictReader(fh, delimiter="\t"):
                gid = (raw.get("gwas_id") or "").strip()
                loc_str = (raw.get("gwas_dir") or raw.get("vcf") or "").strip()
                if not loc_str:
                    continue
                c_sid = (raw.get("cohort_study_id") or raw.get("study_id") or Path(loc_str).parent.parent.name).strip()
                entry = _resolve_gwas_entry(loc_str, c_sid, args.population)
                if entry is None:
                    invalid.append({"study_id": c_sid, "gwas_id": gid or loc_str, "errors": ["Missing VCF or .tbi/.csi index"]})
                    continue
                gid_final = gid or entry["gwas_id"]
                u_tok = gid_final.split("_", 1)[0]
                if only_filter and gid_final not in only_filter and u_tok not in only_filter:
                    continue
                if raw.get("ld_population"):
                    entry["ld_population"] = raw["ld_population"].strip().upper()
                if raw.get("ancestry"):
                    entry["ancestry"] = raw["ancestry"].strip().upper()
                if raw.get("gene_id"):
                    entry["explicit_gene_id"] = raw["gene_id"].strip()
                if raw.get("sample_size"):
                    entry["sample_size"] = float(raw["sample_size"])
                grp_sid = (raw.get("group_study_id") or args.group_study_id or "multisusie").strip()
                grouped.setdefault((grp_sid, gid_final), []).append(entry)
    elif args.multisusie_studies:
        study_paths = [
            Path(s.strip()).expanduser().resolve()
            for s in re.split(r"[\s,]+", args.multisusie_studies)
            if s.strip()
        ]
        source_desc = ",".join(str(p) for p in study_paths)
        default_grp = args.group_study_id or f"multisusie_{len(study_paths)}cohort"
        for sp in study_paths:
            std_dir = sp if sp.name == "standardise" else (sp / "standardise" if (sp / "standardise").is_dir() else sp)
            study_name = sp.parent.name if sp.name == "standardise" else sp.name
            if not std_dir.is_dir():
                invalid.append({"study_id": study_name, "gwas_id": "*", "errors": [f"Directory not found: {std_dir}"]})
                continue
            for gdir in sorted(d for d in std_dir.iterdir() if d.is_dir() and not d.name.startswith(".")):
                gid_cand = gdir.name
                u_tok = gid_cand.split("_", 1)[0]
                if only_filter and gid_cand not in only_filter and u_tok not in only_filter:
                    continue
                entry = _resolve_gwas_entry(str(gdir), study_name, args.population)
                if entry is not None:
                    grouped.setdefault((default_grp, gid_cand), []).append(entry)
    else:
        raise SystemExit("Supply --multisusie-manifest or --multisusie-studies.")

    state = load_state(state_file)
    eligible: list[dict[str, Any]] = []
    skipped_complete = 0

    for (grp_sid, gid), cohorts in sorted(grouped.items()):
        if len(cohorts) < 2 and not parse_bool(args.allow_single_cohort, False):
            continue
        explicit_gid = next((c.get("explicit_gene_id") for c in cohorts if c.get("explicit_gene_id")), None)
        gene, g_reason = resolve_target_gene(gene_annot, gid, cohorts[0].get("metadata") or {}, explicit_gid)
        if gene is None:
            invalid.append({"study_id": grp_sid, "gwas_id": gid, "errors": [f"Could not resolve target cis gene ({g_reason})"]})
            continue

        vcf_sig = [
            {
                "cohort": c["cohort_study_id"],
                "vcf": c["vcf"],
                "size": Path(c["vcf"]).stat().st_size,
                "ld": c["ld_population"],
                "anc": c["ancestry"],
            }
            for c in cohorts
        ]
        h_payload = json.dumps(
            {
                "mode": "multisusie_v1",
                "cohorts": vcf_sig,
                "gene": gene["uniprot_id"],
                "min_maf": float(args.min_maf),
                "gene_window_bp": int(args.gene_window_bp),
                "max_causal": int(args.max_causal),
                "cis_max_iter": int(args.cis_max_iter),
                "cis_significance_lp": float(args.cis_significance_lp or 0.0),
                "output_format": args.output_format,
            },
            sort_keys=True,
        )
        shash = hashlib.sha256(h_payload.encode("utf-8")).hexdigest()[:16]
        rkey = f"{grp_sid}/{gid}"
        gwas_out = outdir / grp_sid / gid
        prev = state.get(rkey)
        expected_ext = ".rds" if args.output_format == "rds" else ".rda"
        if (
            not force
            and prev
            and prev.get("settings_hash") == shash
            and (gwas_out / "qc.json").is_file()
            and (gwas_out / "methods" / f"multisusie{expected_ext}").is_file()
        ):
            skipped_complete += 1
            continue

        eligible.append({
            "row_key": rkey,
            "study_id": grp_sid,
            "gwas_id": gid,
            "gene": gene,
            "cohorts": cohorts,
            "settings_hash": shash,
        })

    row_limit = int(args.row_limit or 0)
    if row_limit > 0 and len(eligible) > row_limit:
        to_run = eligible[:row_limit]
        deferred = len(eligible) - row_limit
    else:
        to_run = eligible
        deferred = 0

    manifest_out = Path(args.output).expanduser().resolve()
    manifest_out.parent.mkdir(parents=True, exist_ok=True)
    with manifest_out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["row_id", "study_id", "gwas_id", "n_cohorts", "gene_symbol", "settings_hash"])
        for idx, item in enumerate(to_run, start=1):
            slug = re.sub(r"[^A-Za-z0-9._-]+", "_", f"{item['study_id']}__{item['gwas_id']}")[:80]
            rid = f"row_{idx:05d}_{slug}"
            w.writerow([rid, item["study_id"], item["gwas_id"], len(item["cohorts"]), item["gene"]["gene_symbol"], item["settings_hash"]])
            if not dry_run:
                payload = {
                    **item,
                    "row_id": rid,
                    "outdir": str(outdir),
                    "chain": str(Path(args.chain).expanduser().resolve()),
                    "grch37_fasta": str(Path(args.grch37_fasta).expanduser().resolve()),
                    "panukb_dir": str(Path(args.panukb_dir).expanduser().resolve()) if args.panukb_dir else "",
                    "local_ld_panel": str(Path(args.local_ld_panel).expanduser().resolve()) if args.local_ld_panel else "",
                    "bm_path": str(Path(args.bm_path).expanduser().resolve()) if args.bm_path else "",
                    "cache_dir": str(Path(args.cache_dir).expanduser().resolve()),
                    "output_format": args.output_format,
                    "compress": args.compress,
                    "min_maf": float(args.min_maf),
                    "gene_window_bp": int(args.gene_window_bp),
                    "max_causal": int(args.max_causal),
                    "cis_max_iter": int(args.cis_max_iter),
                    "cis_significance_lp": float(args.cis_significance_lp or 0.0),
                    "report": report,
                }
                atomic_write_json(rows_dir / f"{rid}.json", payload)

    summary = {
        "mode": "multisusie",
        "source": source_desc,
        "outdir": str(outdir),
        "state_file": str(state_file),
        "dry_run": dry_run,
        "rows_discovered": len(grouped),
        "rows_skipped_complete": skipped_complete,
        "rows_eligible": len(eligible),
        "rows_deferred_due_to_limit": deferred,
        "rows_to_run": len(to_run),
        "invalid": invalid,
    }
    atomic_write_json(Path(args.summary).expanduser().resolve(), summary)
    return 0


def cmd_run_group(args: argparse.Namespace) -> int:
    t0 = time.monotonic()
    row = json.loads(Path(args.row_json).read_text(encoding="utf-8"))
    rec_out = Path(args.output_record).resolve()
    try:
        load_multisusie_backend()
        outdir = Path(row["outdir"]).resolve()
        chain = Path(row["chain"]).resolve()
        fasta = Path(row["grch37_fasta"]).resolve()
        panukb_dir = Path(row["panukb_dir"]).resolve() if row.get("panukb_dir") else default_panukb_dir()
        local_ld = Path(row["local_ld_panel"]).resolve() if row.get("local_ld_panel") else None
        bm_path = row.get("bm_path") or None
        cache_dir = Path(row["cache_dir"]).resolve() / "panukb_multisusie"
        cache_dir.mkdir(parents=True, exist_ok=True)

        gene = row["gene"]
        chrom = normalize_chromosome(gene["chrom"])
        win_bp = int(row.get("gene_window_bp", 1_000_000))
        nom_start = max(1, int(gene["start"]) - win_bp)
        nom_end = int(gene["end"]) + win_bp
        min_maf = float(row.get("min_maf", 0.01))
        min_cis_lp = float(row.get("cis_significance_lp") or 0.0)

        cohort_records: list[dict[str, Any]] = []
        cohort_bounds: list[tuple[int, int]] = []
        canonical_alleles: dict[tuple[int, frozenset[str]], tuple[str, str]] = {}
        observed_max_lp = 0.0

        for c in row["cohorts"]:
            vcf_p = Path(c["vcf"]).resolve()
            vtags = inspect_vcf_tags(vcf_p)
            raw = query_vcf(vcf_p, region=format_vcf_region(vtags, chrom, nom_start, nom_end), vcf_tags=vtags)
            recs, _ = filter_records_by_min_maf(raw, min_maf)
            if not recs:
                continue
            ss_vals = [float(r["SS"]) for r in recs if r.get("SS") is not None and float(r["SS"]) > 0]
            meta = c.get("metadata") or {}
            n_val = (
                float(c["sample_size"])
                if c.get("sample_size")
                else (
                    float(np.median(ss_vals))
                    if ss_vals
                    else float(meta.get("sample_size") or meta.get("n_total") or 10000.0)
                )
            )
            positions = [int(r["POS"]) for r in recs]
            cohort_bounds.append((min(positions), max(positions)))
            rmap: dict[tuple[int, str, str], dict[str, Any]] = {}
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
                    es_v = float(r["ES"])
                    af_v = float(r["AF"]) if r.get("AF") is not None else 0.5
                elif (ref_u, alt_u) == (c_alt, c_ref):
                    es_v = -float(r["ES"])
                    af_v = (1.0 - float(r["AF"])) if r.get("AF") is not None else 0.5
                else:
                    continue
                lp_v = float(r["LP"] or 0.0)
                if lp_v > observed_max_lp:
                    observed_max_lp = lp_v
                rmap[(pos_i, c_ref, c_alt)] = {
                    "ES": es_v,
                    "SE": float(r["SE"]),
                    "LP": lp_v,
                    "AF": af_v,
                }
            if rmap:
                cohort_records.append({
                    "study_id": c["cohort_study_id"],
                    "ancestry": c["ancestry"],
                    "ld_population": c["ld_population"],
                    "n": n_val,
                    "rmap": rmap,
                })

        if not cohort_records or observed_max_lp < min_cis_lp:
            gwas_out = outdir / row["study_id"] / row["gwas_id"]
            gwas_out.mkdir(parents=True, exist_ok=True)
            write_json(gwas_out / "qc.json", {
                "study_id": row["study_id"],
                "gwas_id": row["gwas_id"],
                "overall_status": "NO_SIGNIFICANT_LOCI",
                "analyzed_regions_count": 0,
                "total_credible_sets_multisusie": 0,
                "cis_gene": {"status": "below_cis_significance_gate", "max_lp": round(observed_max_lp, 4)},
                "loci": [],
            })
            write_json(gwas_out / "provenance.json", {
                "study_id": row["study_id"],
                "gwas_id": row["gwas_id"],
                "ld_population": "+".join(sorted({c["ld_population"] for c in row["cohorts"]})),
                "fine_mapping_methods": ["multisusie"],
                "runtime_seconds": round(time.monotonic() - t0, 2),
            })
            write_json(gwas_out / "metadata.json", {
                "study_id": row["study_id"],
                "output_id": row["gwas_id"],
                "phenotype": f"{gene['gene_symbol']} ({row['gwas_id']}) MultiSuSiE",
            })
            atomic_write_json(rec_out, {
                "ok": True,
                "row_key": row["row_key"],
                "study_id": row["study_id"],
                "gwas_id": row["gwas_id"],
                "settings_hash": row["settings_hash"],
                "overall_status": "NO_SIGNIFICANT_LOCI",
                "analyzed_regions": 0,
                "cs_rss": 0,
                "cs_rrss": 0,
                "output_dir": str(gwas_out),
            })
            return 0

        shared_start = max(b[0] for b in cohort_bounds)
        shared_end = min(b[1] for b in cohort_bounds)
        if shared_end <= shared_start:
            shared_start, shared_end = nom_start, nom_end

        union_keys = sorted({
            k for cr in cohort_records for k in cr["rmap"] if shared_start <= k[0] <= shared_end
        })
        union_source = [
            {
                "CHROM": chrom,
                "POS": pos_i,
                "REF": ref_u,
                "ALT": alt_u,
                "ES": 0.0,
                "SE": 1.0,
                "LP": max((cr["rmap"][(pos_i, ref_u, alt_u)]["LP"] for cr in cohort_records if (pos_i, ref_u, alt_u) in cr["rmap"]), default=0.0),
                "SS": 10000.0,
            }
            for (pos_i, ref_u, alt_u) in union_keys
        ]

        from pyliftover import LiftOver
        lifter = LiftOver(str(chain))
        with tempfile.TemporaryDirectory(prefix="gwasfm-mslift-") as ltmp:
            lvcf = Path(ltmp) / "lifted.vcf"
            lifted_by_id, _ = write_lifted_vcf(union_source, lifter, lvcf, expected_chrom=chrom)
            norm_rows = normalize_lifted_vcf(lvcf, fasta, None) if lifted_by_id else []

        for r in norm_rows:
            src = lifted_by_id[r["source_variant_key"]]
            r["CHROM"] = src["CHROM"]
            r["POS"] = src["POS"]
            r["REF"] = src["REF"]
            r["ALT"] = src["ALT"]
            r["LP"] = src["LP"]
            r["GWAS_VARIANT_ID"] = f"{src['CHROM']}:{src['POS']}:{src['REF']}:{src['ALT']}"

        norm_rows = [r for r in norm_rows if r["LD_CHROM"] == chrom]
        norm_rows, _, _ = _clamp_to_grch37_ld_band(
            norm_rows, {"region_type": "cis_gene"}, gene=gene, max_span_bp=9_600_000
        )

        ld_pops = sorted({cr["ld_population"] for cr in cohort_records})
        R_full_by_pop: dict[str, np.ndarray] = {}
        map_by_pop: dict[str, dict[str, int]] = {}
        for pop in ld_pops:
            R_p, m_p = _extract_ld_for_group(
                pop, chrom, norm_rows, cache_dir, panukb_dir, local_ld, bm_path,
                args.bm_disk_cache or None, int(args.threads),
            )
            R_full_by_pop[pop] = R_p
            map_by_pop[pop] = m_p

        # Within each LD population, require complete-case observation across cohorts sharing that LD population
        # so within-ancestry sample size is uniform, while retaining the union across distinct LD populations!
        cohorts_by_pop: dict[str, list[dict[str, Any]]] = {}
        for cr in cohort_records:
            cohorts_by_pop.setdefault(cr["ld_population"], []).append(cr)

        selected_rows: list[dict[str, Any]] = []
        for r in sorted(norm_rows, key=lambda x: (int(x["LD_POS"]), x["LD_REF"], x["LD_ALT"])):
            gid_v = str(r["GWAS_VARIANT_ID"])
            key3 = (int(r["POS"]), str(r["REF"]), str(r["ALT"]))
            keep_any_pop = False
            conflict_partial_pop = False
            for pop, p_cohorts in cohorts_by_pop.items():
                n_obs = sum(1 for cr in p_cohorts if key3 in cr["rmap"])
                if n_obs == len(p_cohorts) and gid_v in map_by_pop[pop]:
                    keep_any_pop = True
                elif 0 < n_obs < len(p_cohorts):
                    conflict_partial_pop = True
            if keep_any_pop and not conflict_partial_pop:
                selected_rows.append(r)

        P_sel = len(selected_rows)
        if P_sel < 5:
            raise RuntimeError(f"Insufficient multi-cohort LD-matched cis variants ({P_sel}) for {row['row_key']}")

        R_sub_by_pop: dict[str, np.ndarray] = {}
        for pop in ld_pops:
            R_sub = np.zeros((P_sel, P_sel), dtype=np.float64)
            pmap = map_by_pop[pop]
            valid_i = [i for i, r in enumerate(selected_rows) if str(r["GWAS_VARIANT_ID"]) in pmap]
            src_i = [pmap[str(selected_rows[i]["GWAS_VARIANT_ID"])] for i in valid_i]
            if valid_i:
                R_sub[np.ix_(valid_i, valid_i)] = R_full_by_pop[pop][np.ix_(src_i, src_i)]
            R_sub_by_pop[pop] = R_sub

        c_inputs: list[dict[str, Any]] = []
        for cr in cohort_records:
            pop = cr["ld_population"]
            pmap = map_by_pop[pop]
            b_arr = np.full(P_sel, np.nan, dtype=np.float64)
            s_arr = np.full(P_sel, np.nan, dtype=np.float64)
            maf_arr = np.full(P_sel, np.nan, dtype=np.float64)
            for i, r in enumerate(selected_rows):
                gid_v = str(r["GWAS_VARIANT_ID"])
                if gid_v not in pmap:
                    continue
                rec = cr["rmap"].get((int(r["POS"]), str(r["REF"]), str(r["ALT"])))
                if rec is not None:
                    b_arr[i] = rec["ES"]
                    s_arr[i] = rec["SE"]
                    maf_arr[i] = min(rec["AF"], 1.0 - rec["AF"])
            c_inputs.append({
                "study_id": cr["study_id"],
                "ancestry": cr["ancestry"],
                "ld_population": pop,
                "n": cr["n"],
                "b": b_arr,
                "s": s_arr,
                "maf": maf_arr,
            })

        variant_meta = [
            {
                "VARIANT_ID": f"{r['LD_CHROM']}:{int(r['LD_POS'])}:{r['LD_REF']}:{r['LD_ALT']}",
                "GWAS_VARIANT_ID": str(r["GWAS_VARIANT_ID"]),
                "LD_CHROM": r["LD_CHROM"],
                "LD_POS": int(r["LD_POS"]),
                "LD_REF": r["LD_REF"],
                "LD_ALT": r["LD_ALT"],
                "GWAS_CHROM": r["CHROM"],
                "GWAS_POS": int(r["POS"]),
                "GWAS_REF": r["REF"],
                "GWAS_ALT": r["ALT"],
                "MAX_LP": max(
                    (cr["rmap"][(int(r["POS"]), str(r["REF"]), str(r["ALT"]))]["LP"]
                     for cr in cohort_records
                     if (int(r["POS"]), str(r["REF"]), str(r["ALT"])) in cr["rmap"]),
                    default=0.0,
                ),
            }
            for r in selected_rows
        ]

        fit_res = fit_multisusie_locus(
            cohort_inputs=c_inputs,
            R_by_ancestry=R_sub_by_pop,
            L=int(row.get("max_causal", 10)),
            max_iter=int(row.get("cis_max_iter", 200)),
            coverage=0.95,
            min_abs_corr=0.5,
        )
        region_id = f"cis_{gene['gene_symbol']}_chr{chrom}_{shared_start}_{shared_end}"
        r_env = r_runtime_environment(PROJECT_ROOT)
        elapsed = time.monotonic() - t0
        gwas_out = save_multisusie_study_output(
            arm_study_id=row["study_id"],
            gwas_id=row["gwas_id"],
            gene=gene,
            region_id=region_id,
            source_region_grch38=f"{chrom}:{shared_start}-{shared_end}",
            variant_rows=variant_meta,
            fit_out=fit_res,
            runtime_seconds=elapsed,
            outdir=outdir,
            r_env=r_env,
            output_format=row.get("output_format", "rda"),
            compress=row.get("compress", "xz"),
            write_report=bool(row.get("report", True)),
        )
        qc_written = json.loads((gwas_out / "qc.json").read_text(encoding="utf-8"))
        atomic_write_json(rec_out, {
            "ok": True,
            "row_key": row["row_key"],
            "study_id": row["study_id"],
            "gwas_id": row["gwas_id"],
            "settings_hash": row["settings_hash"],
            "overall_status": qc_written["overall_status"],
            "analyzed_regions": 1,
            "cs_rss": 0,
            "cs_rrss": len(fit_res["credible_sets"]),
            "output_dir": str(gwas_out),
        })
        return 0
    except Exception as exc:  # noqa: BLE001
        atomic_write_json(rec_out, {
            "ok": False,
            "row_key": row.get("row_key", ""),
            "study_id": row.get("study_id", ""),
            "gwas_id": row.get("gwas_id", ""),
            "settings_hash": row.get("settings_hash", ""),
            "error": str(exc),
        })
        return 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    prep = sub.add_parser("prepare", help="Prepare multi-cohort MultiSuSiE manifest & task JSONs")
    prep.add_argument("--multisusie-manifest", default=None)
    prep.add_argument("--multisusie-studies", default=None)
    prep.add_argument("--group-study-id", default=None)
    prep.add_argument("--only-gwas", default=None)
    prep.add_argument("--allow-single-cohort", default="false")
    prep.add_argument("--outdir", required=True)
    prep.add_argument("--state-file", required=True)
    prep.add_argument("--output", required=True)
    prep.add_argument("--summary", required=True)
    prep.add_argument("--rows-dir", required=True)
    prep.add_argument("--chain", default=str(DEFAULT_REF_ROOT / "genome/hg38ToHg19.over.chain.gz"))
    prep.add_argument("--grch37-fasta", default=str(DEFAULT_REF_ROOT / "genome/GRCh37.fa"))
    prep.add_argument("--gene-annotation", default=str(DEFAULT_REF_ROOT / "annotation/gene_regions_GRCh38.tsv"))
    prep.add_argument("--panukb-dir", default=None)
    prep.add_argument("--local-ld-panel", default=None)
    prep.add_argument("--bm-path", default=None)
    prep.add_argument("--cache-dir", required=True)
    prep.add_argument("--population", default="EUR")
    prep.add_argument("--output-format", default="rda")
    prep.add_argument("--compress", default="xz")
    prep.add_argument("--min-maf", type=float, default=0.01)
    prep.add_argument("--gene-window-bp", type=int, default=1_000_000)
    prep.add_argument("--max-causal", type=int, default=10)
    prep.add_argument("--cis-max-iter", type=int, default=200)
    prep.add_argument("--cis-significance-lp", type=float, default=0.0)
    prep.add_argument("--row-limit", type=int, default=0)
    prep.add_argument("--force", default="false")
    prep.add_argument("--dry-run", default="false")
    prep.add_argument("--report", default="true")

    rg = sub.add_parser("run-group", help="Execute MultiSuSiE for one prepared group JSON")
    rg.add_argument("--row-json", required=True)
    rg.add_argument("--output-record", required=True)
    rg.add_argument("--bm-disk-cache", default="")
    rg.add_argument("--threads", type=int, default=4)
    return p


def main() -> int:
    args = build_parser().parse_args()
    if args.cmd == "prepare":
        return cmd_prepare(args)
    return cmd_run_group(args)


if __name__ == "__main__":
    raise SystemExit(main())
