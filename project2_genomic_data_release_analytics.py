"""
Genomic Data Release Analytics — ICGC Harmonization & Access Tracking
Author: Nicholas Steven
Target Role: Bioinformatician II / Research Data Analyst — OICR
Repo: github.com/nicholasstevenr/OICR-health-data-project

ICGC schema validation, consent tier classification, cross-institution
harmonization, release readiness scorecard, and data access tracking.
"""

import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")

# ── ICGC Required Metadata Fields ────────────────────────────────────────────

ICGC_REQUIRED_FIELDS = [
    "sample_id", "donor_id", "project_code", "cancer_type_icd10",
    "tumour_histological_type_icd_o_3", "specimen_type",
    "sample_type", "sequencing_strategy", "analysis_type",
]

ICGC_VOCAB = {
    "specimen_type":       ["Primary tumour","Metastatic tumour","Recurrent tumour","Normal"],
    "sequencing_strategy": ["WGS","WXS","RNA-Seq","Targeted-Seq","ChIP-Seq"],
    "analysis_type":       ["germline_variation","somatic_mutations","gene_expression"],
}

CONSENT_TIERS = {"open": 1, "controlled": 2, "embargoed": 3}
HARMONIZATION_INCONSISTENCY_THRESHOLD = 0.05   # 5%


# ── Load ──────────────────────────────────────────────────────────────────────

def load(metadata_path: str, consent_path: str, access_path: str = None):
    meta    = pd.read_csv(metadata_path) if metadata_path and pd.io.common.file_exists(metadata_path) \
              else _synthetic_metadata()
    consent = pd.read_csv(consent_path)  if consent_path  and pd.io.common.file_exists(consent_path)  \
              else _synthetic_consent(meta)
    access  = pd.read_csv(access_path, parse_dates=["application_date","approval_date"]) \
              if access_path and pd.io.common.file_exists(access_path) else pd.DataFrame()
    print(f"Samples: {len(meta):,}  |  Consent records: {len(consent):,}")
    return meta, consent, access


def _synthetic_metadata(n: int = 200) -> pd.DataFrame:
    np.random.seed(7)
    df = pd.DataFrame({
        "sample_id":     [f"OICR_{i:04d}" for i in range(n)],
        "donor_id":      [f"D_{i:04d}" for i in range(n)],
        "project_code":  np.random.choice(["PACA-CA","BRCA-CA","LICA-CA","COCA-CA"], n),
        "institution":   np.random.choice(["UHN","PMH","McMaster","Ottawa"], n),
        "cancer_type_icd10": np.random.choice(["C25","C50","C18","C22"], n),
        "tumour_histological_type_icd_o_3": np.random.choice(["8140/3","8500/3","8480/3",None], n, p=[0.3,0.3,0.3,0.1]),
        "specimen_type": np.random.choice(["Primary tumour","Normal","Metastatic tumour",None], n, p=[0.5,0.3,0.15,0.05]),
        "sample_type":   np.random.choice(["FF","FFPE"], n),
        "sequencing_strategy": np.random.choice(["WGS","WXS"], n, p=[0.7,0.3]),
        "analysis_type": np.random.choice(["somatic_mutations","germline_variation"], n),
        "tnm_stage":     np.random.choice(["I","II","III","IV","stage1","stage2"], n),  # deliberate inconsistency
        "release_tier_requested": np.random.choice(["open","controlled"], n, p=[0.3,0.7]),
    })
    return df


def _synthetic_consent(meta: pd.DataFrame) -> pd.DataFrame:
    n = len(meta)
    return pd.DataFrame({
        "donor_id":          meta["donor_id"],
        "consent_tier_max":  np.where(np.random.random(n) > 0.05, "controlled", "open"),
    })


# ── 1. ICGC Schema Validation ─────────────────────────────────────────────────

def icgc_schema_validation(meta: pd.DataFrame) -> pd.DataFrame:
    results = []
    for _, row in meta.iterrows():
        flags = []
        # Required fields
        for field in ICGC_REQUIRED_FIELDS:
            if field not in meta.columns or pd.isna(row.get(field)):
                flags.append(f"missing_{field}")
        # Vocabulary compliance
        for field, valid_vals in ICGC_VOCAB.items():
            if field in meta.columns and pd.notna(row.get(field)):
                if row[field] not in valid_vals:
                    flags.append(f"invalid_vocab_{field}:{row[field]}")
        results.append({
            "sample_id":      row["sample_id"],
            "n_schema_flags": len(flags),
            "schema_flags":   "; ".join(flags) if flags else "",
            "schema_pass":    len(flags) == 0,
        })
    df = pd.DataFrame(results)
    print(f"\n── ICGC Schema Validation ──")
    print(f"  Pass: {df['schema_pass'].sum()}  |  Fail: {(~df['schema_pass']).sum()}")
    return df


# ── 2. Consent Tier Classification ───────────────────────────────────────────

def consent_verification(meta: pd.DataFrame, consent: pd.DataFrame) -> pd.DataFrame:
    merged = meta[["sample_id","donor_id","release_tier_requested"]].merge(
        consent[["donor_id","consent_tier_max"]], on="donor_id", how="left"
    )
    merged["requested_level"] = merged["release_tier_requested"].map(CONSENT_TIERS).fillna(0)
    merged["consented_level"] = merged["consent_tier_max"].map(CONSENT_TIERS).fillna(0)
    merged["consent_breach"]  = merged["requested_level"] > merged["consented_level"]

    print(f"\n── Consent Verification ──")
    print(f"  Consent breaches detected: {merged['consent_breach'].sum()}")
    return merged


# ── 3. Cross-Institution Harmonization ───────────────────────────────────────

def harmonization_check(meta: pd.DataFrame, target_vars: list = None) -> pd.DataFrame:
    """
    Check consistency of clinical variable coding across institutions.
    Flags variables with >5% coding inconsistency (non-standard values).
    """
    if target_vars is None:
        target_vars = ["tnm_stage","specimen_type","cancer_type_icd10"]

    STANDARD_VALUES = {
        "tnm_stage":  ["I","II","III","IV"],
        "specimen_type": ICGC_VOCAB["specimen_type"],
    }

    results = []
    for var in target_vars:
        if var not in meta.columns:
            continue
        if var in STANDARD_VALUES:
            non_standard = ~meta[var].isin(STANDARD_VALUES[var]) & meta[var].notna()
            inconsistency_rate = non_standard.mean()
        else:
            # Mode-based: flag values deviating from most common within each cancer type
            inconsistency_rate = 0.0
            if "institution" in meta.columns:
                for inst, grp in meta.groupby("institution"):
                    val_counts = grp[var].value_counts(normalize=True)
                    if len(val_counts) > 1:
                        inconsistency_rate += (1 - val_counts.iloc[0]) * len(grp) / len(meta)

        results.append({
            "variable":            var,
            "inconsistency_rate":  round(inconsistency_rate, 4),
            "flag":                inconsistency_rate > HARMONIZATION_INCONSISTENCY_THRESHOLD,
        })
    df = pd.DataFrame(results)
    print(f"\n── Harmonization Check ──")
    print(df.to_string(index=False))
    return df


# ── 4. Release Readiness Scorecard ───────────────────────────────────────────

def release_readiness(schema_df: pd.DataFrame, consent_df: pd.DataFrame,
                       qc_df: pd.DataFrame = None) -> pd.DataFrame:
    scorecard = schema_df[["sample_id","schema_pass"]].merge(
        consent_df[["sample_id","consent_breach"]], on="sample_id", how="left"
    )
    if qc_df is not None and "qc_status" in qc_df.columns:
        scorecard = scorecard.merge(
            qc_df[["sample_id","qc_status"]], on="sample_id", how="left"
        )
        scorecard["qc_pass"] = scorecard["qc_status"] == "PASS"
    else:
        scorecard["qc_pass"] = True

    scorecard["release_ready"] = (
        scorecard["schema_pass"] &
        ~scorecard["consent_breach"].fillna(False) &
        scorecard["qc_pass"]
    )
    n_ready = scorecard["release_ready"].sum()
    n_total = len(scorecard)
    print(f"\n── Release Readiness: {n_ready}/{n_total} samples ready for release ({n_ready/n_total*100:.1f}%) ──")
    return scorecard


# ── 5. Data Access Tracking ───────────────────────────────────────────────────

def access_tracking(access: pd.DataFrame) -> pd.DataFrame:
    if access.empty:
        return pd.DataFrame()
    access["days_to_approval"] = (access["approval_date"] - access["application_date"]).dt.days
    summary = (
        access.groupby("data_tier")
        .agg(n_requests=("days_to_approval","count"),
             median_days=("days_to_approval","median"),
             p90_days=("days_to_approval",lambda x: x.quantile(0.90)))
        .round(1).reset_index()
    )
    print(f"\n── Access Tracking ──")
    print(summary.to_string(index=False))
    return summary


# ── Export ────────────────────────────────────────────────────────────────────

def export_all(results: dict, outdir: str = "output") -> None:
    import os; os.makedirs(outdir, exist_ok=True)
    for name, df in results.items():
        if isinstance(df, pd.DataFrame) and len(df):
            df.to_csv(f"{outdir}/{name}.csv", index=False)
            print(f"  Exported → output/{name}.csv")


if __name__ == "__main__":
    meta, consent, access = load(
        "data/oicr_sample_metadata.csv",
        "data/oicr_consent_records.csv",
        "data/oicr_data_access_requests.csv",
    )
    schema_df   = icgc_schema_validation(meta)
    consent_df  = consent_verification(meta, consent)
    harmonic_df = harmonization_check(meta)
    readiness   = release_readiness(schema_df, consent_df)
    access_df   = access_tracking(access)

    export_all({"schema_validation": schema_df, "consent_verification": consent_df,
                "harmonization_check": harmonic_df, "release_readiness": readiness,
                "access_tracking": access_df})
