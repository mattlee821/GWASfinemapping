import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from run_finemapping import r_runtime_environment  # noqa: E402
from run_multisusie_somascan import (  # noqa: E402
    fit_multisusie_locus,
    load_multisusie_backend,
    save_multisusie_study_output,
)


def _marginal_sumstats(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    n = X.shape[0]
    Xc = X - X.mean(0)
    yc = y - y.mean()
    xtx = (Xc * Xc).sum(0)
    b = (Xc.T @ yc) / xtx
    resid = yc[:, None] - Xc * b[None, :]
    s2 = (resid * resid).sum(0) / (n - 2)
    se = np.sqrt(s2 / xtx)
    R = np.corrcoef(Xc, rowvar=False)
    return b, se, R, n


class MultiSuSiEIntegrationTests(unittest.TestCase):
    def setUp(self):
        try:
            load_multisusie_backend()
        except ImportError:
            self.skipTest("MultiSuSiE package is not installed in this environment")

    def test_cross_ancestry_narrowing_censored_variants_and_coloc_rda(self):
        rng = np.random.default_rng(123)
        P = 30
        X_eur1 = rng.standard_normal((1800, P))
        X_eur1[:, 6] = 0.98 * X_eur1[:, 5] + np.sqrt(1 - 0.98**2) * rng.standard_normal(1800)
        X_eur2 = rng.standard_normal((1200, P))
        X_eur2[:, 6] = 0.98 * X_eur2[:, 5] + np.sqrt(1 - 0.98**2) * rng.standard_normal(1200)
        X_afr = rng.standard_normal((900, P))
        X_afr[:, 6] = 0.20 * X_afr[:, 5] + np.sqrt(1 - 0.20**2) * rng.standard_normal(900)

        y1 = 0.4 * X_eur1[:, 5] + rng.standard_normal(1800)
        y2 = 0.4 * X_eur2[:, 5] + rng.standard_normal(1200)
        ya = 0.4 * X_afr[:, 5] + rng.standard_normal(900)
        y1 /= y1.std(ddof=1)
        y2 /= y2.std(ddof=1)
        ya /= ya.std(ddof=1)

        b1, s1, R_eur, n1 = _marginal_sumstats(X_eur1, y1)
        b2, s2, _, n2 = _marginal_sumstats(X_eur2, y2)
        ba, sa, R_afr, na = _marginal_sumstats(X_afr, ya)

        maf_eur = np.full(P, 0.25)
        maf_afr = np.full(P, 0.25)
        # Variant 28 is EUR-only (monomorphic/absent in AFR); Variant 29 is AFR-only
        ba[28] = sa[28] = maf_afr[28] = np.nan
        b1[29] = s1[29] = b2[29] = s2[29] = maf_eur[29] = np.nan

        out = fit_multisusie_locus(
            cohort_inputs=[
                {"study_id": "pietzner_2021_34648354", "ancestry": "EUR", "n": n1, "b": b1, "s": s1, "maf": maf_eur},
                {"study_id": "zhang_2022_35501419_EUR", "ancestry": "EUR", "n": n2, "b": b2, "s": s2, "maf": maf_eur},
                {"study_id": "zhang_2022_35501419_AFR", "ancestry": "AFR", "n": na, "b": ba, "s": sa, "maf": maf_afr},
            ],
            R_by_ancestry={"EUR": R_eur, "AFR": R_afr},
            L=5,
            max_iter=200,
        )
        self.assertTrue(out["converged"])
        self.assertEqual(len(out["credible_sets"]), 1)
        cs0 = out["credible_sets"][0]
        self.assertEqual(cs0["indices_1based"], [6])  # 0-based index 5 -> 1-based index 6!
        self.assertFalse(np.isnan(cs0["min_abs_corr"]))

        # Verify exact lbf_variable -> alpha softmax identity on active component
        eff_0 = cs0["effect_index"] - 1
        lv = out["lbf_variable"][eff_0]
        w = np.exp(lv - np.max(lv))
        alpha_check = w / np.sum(w)
        self.assertLess(float(np.max(np.abs(alpha_check - out["alpha"][eff_0]))), 1e-10)

        # Verify .rda serialization and compare_pqtl_studies.R compatibility
        project_root = Path(__file__).resolve().parents[2]
        r_env = r_runtime_environment(project_root)
        rscript = shutil.which("Rscript", path=r_env.get("PATH"))
        if not rscript:
            return

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fmdir = root / "finemapping"
            var_rows = [
                {
                    "VARIANT_ID": f"21:{1000 + i}:A:G",
                    "GWAS_VARIANT_ID": f"21:{2000 + i}:A:G",
                    "LD_CHROM": "21", "LD_POS": 1000 + i, "LD_REF": "A", "LD_ALT": "G",
                    "GWAS_CHROM": "21", "GWAS_POS": 2000 + i, "GWAS_REF": "A", "GWAS_ALT": "G",
                    "MAX_LP": 50.0 if i == 5 else 0.5,
                }
                for i in range(P)
            ]
            gene = {"uniprot_id": "P04155", "gene_symbol": "TFF1", "chrom": "21", "start": 1500, "end": 2500}
            save_multisusie_study_output(
                arm_study_id="somascan_multisusie_3cohort",
                gwas_id="P04155_9185-15",
                gene=gene,
                region_id="cis_TFF1_chr21_1000_3000",
                source_region_grch38="21:1000-3000",
                variant_rows=var_rows,
                fit_out=out,
                runtime_seconds=1.2,
                outdir=fmdir,
                r_env=r_env,
            )
            save_multisusie_study_output(
                arm_study_id="somascan_multisusie_4cohort",
                gwas_id="P04155_9185-15",
                gene=gene,
                region_id="cis_TFF1_chr21_1000_3000",
                source_region_grch38="21:1000-3000",
                variant_rows=var_rows,
                fit_out=out,
                runtime_seconds=1.4,
                outdir=fmdir,
                r_env=r_env,
            )
            cmp_out = root / "cmp"
            res = subprocess.run(
                [
                    rscript,
                    str(project_root / "scripts" / "compare_pqtl_studies.R"),
                    "--finemap-dir", str(fmdir),
                    "--outdir", str(cmp_out),
                ],
                cwd=project_root,
                env=r_env,
                capture_output=True,
                text=True,
            )
            self.assertEqual(res.returncode, 0, res.stderr)
            import csv
            with (cmp_out / "02_pairwise_coloc_and_concordance.tsv").open(encoding="utf-8") as fh:
                prows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(len(prows), 1)
            self.assertGreater(float(prows[0]["coloc_max_PPH4"]), 0.999)
            self.assertEqual(prows[0]["coloc_top_pair"], "L1×L1 (21:2005:A:G)")
            with (cmp_out / "00_study_level_summary.tsv").open(encoding="utf-8") as fh:
                srows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(int(srows[0]["single_snp_cs"]), 1)
            self.assertEqual(float(srows[0]["median_cs_size"]), 1.0)

    def test_severe_ld_mismatch_negative_er2_never_crashes_spd_solver(self):
        P = 15
        R_eur = np.eye(P)
        R_eur[0, 1] = R_eur[1, 0] = 0.99
        R_afr = np.eye(P)
        # Contradictory extreme effects on tightly correlated SNPs (drives ER2 < 0 without sigma2 floor)
        b_mis = np.zeros(P); b_mis[0] = 2.5; b_mis[1] = -2.5
        s_mis = np.full(P, 0.02)
        b_afr = np.zeros(P); b_afr[0] = 0.3
        s_afr = np.full(P, 0.05)
        maf = np.full(P, 0.25)

        out = fit_multisusie_locus(
            cohort_inputs=[
                {"study_id": "ferkingstad_2021_34857953", "ancestry": "EUR", "n": 35559, "b": b_mis, "s": s_mis, "maf": maf},
                {"study_id": "zhang_2022_35501419_AFR", "ancestry": "AFR", "n": 1871, "b": b_afr, "s": s_afr, "maf": maf},
            ],
            R_by_ancestry={"EUR": R_eur, "AFR": R_afr},
            L=5,
            max_iter=50,
        )
        self.assertTrue(np.all(np.isfinite(out["pip"])))
        self.assertTrue(np.all(np.isfinite(out["lbf_variable"])))
        self.assertTrue(np.all(out["sigma2"] >= 1e-2))


if __name__ == "__main__":
    unittest.main()
