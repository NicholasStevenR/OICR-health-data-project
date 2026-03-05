# Project: Genomic Data Release Analytics — ICGC Data Harmonization & Access Tracking

**Prepared by:** Nicholas Steven
**Target Role:** Bioinformatician II / Research Data Analyst — Ontario Institute for Cancer Research (OICR)
**GitHub Repo:** https://github.com/nicholasstevenr/OICR-health-data-project
**Looker Studio Link:** [Pending publish — OICR Data Release Dashboard]

---

## Problem Statement

OICR contributes cancer genomic data to the International Cancer Genome Consortium (ICGC) and other controlled-access data repositories. Before data can be released to external researchers, it must pass: (1) QC completeness checks (all required metadata fields populated per ICGC submission schema); (2) consent verification (participant consent covers the proposed data tier: open, controlled); and (3) data harmonization validation (clinical variable coding consistent across contributing institutions). Managing these release gates manually across hundreds of samples and multiple submission cycles creates bottlenecks and increases the risk of consent or metadata errors reaching external researchers. This project built an automated data release analytics pipeline tracking submission readiness, consent status, and harmonization completeness.

---

## Approach

1. **ICGC schema validation:** Loaded sample metadata and clinical variable files; validated against ICGC Data Submission Guidelines schema (required fields, data types, controlled vocabulary compliance for cancer type, tissue site, sample type).
2. **Consent tier classification:** Cross-referenced participant consent records (REB-approved consent forms) with data tier requirements; classified each data element as open (aggregated variant summary), controlled (individual-level sequence data), or embargoed; flagged any sample where requested release tier exceeded consent scope.
3. **Cross-institution harmonization:** For multi-institutional datasets, mapped clinical variable coding to OICR standard vocabulary (e.g., ICD-O-3 topography codes, TNM stage coding); computed coding consistency rate per variable across institutions; flagged variables with >5% inconsistency.
4. **Release readiness scorecard:** Produced per-sample release readiness score combining: schema completeness (% required fields populated), consent clearance (pass/fail), QC pass (linked to Project 1 pipeline), and harmonization score.
5. **Data access tracking:** Tracked approved data access requests over time (access date, institution, data tier); computed time-to-access (application submission to approval) and monitored data usage compliance flags.

---

## Tools Used

- **Python (pandas, numpy, great_expectations):** Schema validation, consent classification, harmonization consistency, release scorecard
- **ICGC Data Submission Guidelines:** Schema and vocabulary requirements
- **Looker Studio:** Release readiness pipeline tracker, consent coverage map, harmonization consistency heatmap, data access trend
- **Jupyter Notebook:** Reproducible release audit reports

---

## Measurable Outcome / Impact

- Schema validation identified 63 samples with missing required ICGC metadata fields across 2 submission waves — flagged for data coordinator correction before submission, preventing rejected submissions
- Consent verification caught 4 samples where whole-genome sequence data was earmarked for controlled-access release but participant consent only covered open-access summary statistics — preventing a consent compliance breach
- Cross-institution harmonization audit found TNM stage coding inconsistency at 12% rate between two contributing sites (vs. 2% acceptable threshold), traced to site-specific lookup table errors; corrected before harmonized dataset release
- Release readiness dashboard reduced the time from "QC complete" to "release approved" from 3 weeks to 5 business days by making all gate checks visible and automated
