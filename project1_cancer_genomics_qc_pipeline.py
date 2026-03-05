"""
Cancer Genomics Data Quality Pipeline — OICR / GSI
Author: Nicholas Steven
Target Role: Bioinformatician II / Research Data Analyst — OICR
Repo: github.com/nicholasstevenr/OICR-health-data-project

Sequencing QC metric ingestion, sample-level flagging, somatic variant
confidence filtering, batch effect detection, and cohort QC summary.
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
import json, os, glob
import warnings
warnings.filterwarnings("ignore")

# ── OICR / ICGC QC Thresholds ────────────────────────────────────────────────

QC_THRESHOLDS = {
    "mean_coverage_min_wgs":  30.0,   # × for WGS
    "mean_coverage_min_wxs":  80.0,   # × for WXS
    "duplication_rate_max":    0.30,   # >30% = flag
    "contamination_max":       0.03,   # FREEMIX >0.03 = flag
    "uniformity_min":          0.80,   # % bases ≥20× >80%
    "on_target_rate_min":      0.85,   # WXS only
}

VARIANT_HARD_FILTER = {
    "filter_pass":   True,
    "alt_depth_min": 5,
    "vaf_min":       0.05,
}


# ── 1. Parse QC Metrics ───────────────────────────────────────────────────────

def parse_qc_metrics(qc_dir: str) -> pd.DataFrame:
    """
    Parse Picard CollectWgsMetrics / CollectHsMetrics JSON outputs per sample.
    Falls back to synthetic demo data if no files found.
    """
    records = []
    json_files = glob.glob(os.path.join(qc_dir, "*.qc.json")) if os.path.isdir(qc_dir) else []

    if not json_files:
        return _synthetic_qc_metrics()

    for fpath in json_files:
        with open(fpath) as f:
            data = json.load(f)
        sample_id = os.path.basename(fpath).replace(".qc.json","")
        records.append({
            "sample_id":         sample_id,
            "sequencing_type":   data.get("sequencing_type","WGS"),
            "mean_coverage":     data.get("MEAN_COVERAGE", np.nan),
            "pct_bases_20x":     data.get("PCT_20X", np.nan),
            "duplication_rate":  data.get("PERCENT_DUPLICATION", np.nan),
            "freemix":           data.get("FREEMIX", np.nan),
            "on_target_rate":    data.get("PCT_SELECTED_BASES", np.nan),
            "insert_size_median":data.get("MEDIAN_INSERT_SIZE", np.nan),
            "sequencing_run":    data.get("run_id","unknown"),
        })
    return pd.DataFrame(records)


def _synthetic_qc_metrics(n: int = 100) -> pd.DataFrame:
    np.random.seed(42)
    return pd.DataFrame({
        "sample_id":         [f"OICR_{i:04d}" for i in range(n)],
        "sequencing_type":   np.random.choice(["WGS","WXS"], n, p=[0.7,0.3]),
        "mean_coverage":     np.random.normal(45, 12, n).clip(10, 120),
        "pct_bases_20x":     np.random.normal(0.88, 0.08, n).clip(0.3, 1.0),
        "duplication_rate":  np.random.normal(0.18, 0.09, n).clip(0.02, 0.65),
        "freemix":           np.random.exponential(0.01, n).clip(0, 0.15),
        "on_target_rate":    np.random.normal(0.92, 0.05, n).clip(0.5, 1.0),
        "insert_size_median":np.random.normal(360, 50, n).clip(100, 600),
        "sequencing_run":    np.random.choice([f"RUN_{r:03d}" for r in range(1,10)], n),
    })


# ── 2. Sample-Level QC Flagging ───────────────────────────────────────────────

def flag_samples(qc: pd.DataFrame) -> pd.DataFrame:
    qc = qc.copy()
    cov_min = qc["sequencing_type"].map(
        {"WGS": QC_THRESHOLDS["mean_coverage_min_wgs"],
         "WXS": QC_THRESHOLDS["mean_coverage_min_wxs"]}
    ).fillna(QC_THRESHOLDS["mean_coverage_min_wgs"])

    qc["flag_low_coverage"]    = qc["mean_coverage"] < cov_min
    qc["flag_high_duplication"]= qc["duplication_rate"] > QC_THRESHOLDS["duplication_rate_max"]
    qc["flag_contamination"]   = qc["freemix"] > QC_THRESHOLDS["contamination_max"]
    qc["flag_low_uniformity"]  = qc["pct_bases_20x"] < QC_THRESHOLDS["uniformity_min"]

    flag_cols = [c for c in qc.columns if c.startswith("flag_")]
    qc["any_flag"] = qc[flag_cols].any(axis=1)
    qc["n_flags"]  = qc[flag_cols].sum(axis=1)
    qc["qc_status"]= qc["any_flag"].map({True: "FAIL", False: "PASS"})

    print(f"\n── Sample QC ──")
    print(f"  PASS: {(qc['qc_status']=='PASS').sum()}  |  FAIL: {(qc['qc_status']=='FAIL').sum()}")
    for col in flag_cols:
        print(f"  {col}: {qc[col].sum()} samples")
    return qc


# ── 3. Variant Confidence Classifier ─────────────────────────────────────────

def variant_classifier(variants_path: str = None) -> dict:
    """
    Train logistic regression classifier on variant QC features to
    distinguish true somatic variants from technical artifacts.
    Features: strand_bias_fs, mapping_quality_rms, base_quality_rms,
              alt_depth, vaf, read_position_rank_sum.
    """
    if variants_path and os.path.exists(variants_path):
        df = pd.read_csv(variants_path)
    else:
        # Synthetic variant data
        np.random.seed(99)
        n = 5000
        df = pd.DataFrame({
            "strand_bias_fs":         np.random.exponential(2, n),
            "mapping_quality_rms":    np.random.normal(55, 10, n).clip(0, 70),
            "base_quality_rms":       np.random.normal(35, 6, n).clip(0, 45),
            "alt_depth":              np.random.poisson(15, n),
            "vaf":                    np.random.beta(2, 8, n),
            "read_pos_rank_sum":      np.random.normal(0, 2, n),
            # Label: 1=true somatic, 0=artifact (simulated)
            "true_somatic":           np.random.choice([0,1], n, p=[0.35, 0.65]),
        })
        # Hard-filter first
        df = df[(df["alt_depth"] >= VARIANT_HARD_FILTER["alt_depth_min"]) &
                (df["vaf"] >= VARIANT_HARD_FILTER["vaf_min"])].copy()

    features = ["strand_bias_fs","mapping_quality_rms","base_quality_rms",
                "alt_depth","vaf","read_pos_rank_sum"]
    feat_avail = [f for f in features if f in df.columns]
    model_df = df[feat_avail + ["true_somatic"]].dropna()

    X = StandardScaler().fit_transform(model_df[feat_avail].values)
    y = model_df["true_somatic"].values
    lr = LogisticRegression(C=1.0, max_iter=300, random_state=42)
    cv_auc = cross_val_score(lr, X, y, cv=5, scoring="roc_auc")
    lr.fit(X, y)

    print(f"\n── Variant Classifier ──")
    print(f"  5-fold CV AUC: {cv_auc.mean():.3f} ± {cv_auc.std():.3f}")
    print(f"  Hard-filter pass: {len(model_df):,} variants")
    return {"cv_auc": round(cv_auc.mean(), 3), "n_variants_post_hardfilter": len(model_df)}


# ── 4. Batch Effect Detection ─────────────────────────────────────────────────

def batch_effect_detection(qc: pd.DataFrame) -> pd.DataFrame:
    """Flag sequencing runs where mean coverage is >2 SD below the overall mean."""
    overall_mean = qc["mean_coverage"].mean()
    overall_std  = qc["mean_coverage"].std()

    batch_summary = (
        qc.groupby("sequencing_run")
        .agg(n=("mean_coverage","count"),
             batch_mean_coverage=("mean_coverage","mean"),
             batch_std=("mean_coverage","std"))
        .round(2).reset_index()
    )
    batch_summary["z_vs_overall"] = (
        (batch_summary["batch_mean_coverage"] - overall_mean) / overall_std
    ).round(3)
    batch_summary["flag_batch_effect"] = batch_summary["z_vs_overall"] < -2.0

    flagged = batch_summary[batch_summary["flag_batch_effect"]]
    print(f"\n── Batch Effect Detection ──")
    print(f"  Flagged runs: {len(flagged)}")
    if len(flagged):
        print(flagged[["sequencing_run","n","batch_mean_coverage","z_vs_overall"]].to_string(index=False))
    return batch_summary


# ── Export ────────────────────────────────────────────────────────────────────

def export_all(results: dict, outdir: str = "output") -> None:
    os.makedirs(outdir, exist_ok=True)
    for name, obj in results.items():
        if isinstance(obj, pd.DataFrame) and len(obj):
            obj.to_csv(f"{outdir}/{name}.csv", index=False)
            print(f"  Exported → output/{name}.csv")


if __name__ == "__main__":
    qc_raw   = parse_qc_metrics("data/qc_json_outputs/")
    qc_flagged = flag_samples(qc_raw)
    variant_res = variant_classifier()
    batch_df  = batch_effect_detection(qc_flagged)

    export_all({"sample_qc_flagged": qc_flagged, "batch_effect_summary": batch_df})
    print(f"\n  Variant classifier AUC: {variant_res['cv_auc']}")
