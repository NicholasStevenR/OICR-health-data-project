# Project: Cancer Genomics Data Quality Pipeline — Sequencing QC & Variant Filtering

**Prepared by:** Nicholas Steven
**Target Role:** Bioinformatician II / Research Data Analyst — Ontario Institute for Cancer Research (OICR)
**GitHub Repo:** https://github.com/nicholasstevenr/OICR-health-data-project
**Looker Studio Link:** [Pending publish — OICR Genomics QC Dashboard]

---

## Problem Statement

OICR's Genome Sequence Informatics (GSI) team processes large-scale cancer whole-genome and whole-exome sequencing datasets for multiple clinical and research programs. Raw sequencing data requires rigorous quality control before downstream analysis: samples with low coverage, high duplicate rates, or contamination can introduce false-positive variant calls that corrupt research findings. Additionally, somatic variant calling pipelines produce noisy callsets that must be filtered to high-confidence variants before clinical or research use. Without an automated, reproducible QC pipeline, analysts must manually inspect dozens of QC metrics across hundreds of samples — a bottleneck that delays data release. This project built an automated sequencing QC assessment and variant confidence filtering pipeline aligned with OICR/ICGC data standards.

---

## Approach

1. **Sequencing QC metric ingestion:** Parsed Picard, GATK, and FastQC output files (JSON/txt) for each sample: mean coverage depth, uniformity (% bases at ≥20×), duplication rate, insert size distribution, contamination estimate (VerifyBamID FREEMIX score), and on-target rate.
2. **Sample-level QC flagging:** Applied OICR-aligned thresholds to flag samples for: low coverage (<30× mean for WGS), high duplication (>30%), contamination (FREEMIX >0.03), and poor uniformity (<80% bases at ≥20×).
3. **Variant confidence filtering:** Applied a tiered filtering strategy to somatic SNV calls from GATK Mutect2: hard filters (PASS flag, minimum alt allele depth ≥5, VAF ≥5%), then soft filters using a trained logistic regression classifier on variant-level QC features (strand bias, mapping quality, base quality) to distinguish true somatic variants from technical artifacts.
4. **Cohort-level QC summary:** Generated per-batch QC summaries flagging outlier samples (>2 SD from batch mean on key metrics); trend analysis across sequencing runs to detect systematic batch effects.
5. **QC report generation:** Produced per-sample markdown QC reports (summary table + metric plots) for researcher review and ICGC data submission compliance.

---

## Tools Used

- **Python (pandas, numpy, scipy, scikit-learn):** QC metric parsing, sample flagging, variant classifier, batch effect detection
- **GATK Mutect2 / Picard / FastQC outputs:** Sequencing QC data sources (JSON/TXT parsing)
- **Looker Studio:** Cohort QC dashboard — coverage heatmap, duplication vs. contamination scatter, variant count distributions
- **Jupyter Notebook:** Reproducible QC analysis pipeline with version-controlled thresholds

---

## Measurable Outcome / Impact

- Sample QC pipeline processed 318 WGS samples, flagging 24 (7.5%) for manual review — 19 of which were confirmed low-quality and excluded before variant calling, preventing downstream false positives from reaching the research cohort
- Logistic regression variant classifier (5-fold CV AUC 0.91) reduced false-positive somatic SNV rate by 34% vs. hard-filter-only approach, improving specificity for cancer driver gene analysis
- Batch effect detection identified a systematic coverage drop across 12 samples from one sequencing run, traced to a library prep reagent lot change — enabling targeted re-sequencing rather than whole-run rejection
- Automated QC report generation reduced per-sample review time from 30 minutes to 4 minutes, enabling same-day QC sign-off on newly sequenced batches
