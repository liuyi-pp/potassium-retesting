# Reproducibility and provenance

## Frozen analysis scope

Version 1.0 was written after feasibility work. It is not prospective preregistration or independent confirmation. The source implements an all-enrolled descriptive analysis with no added variables, clinical effect models or P values. The existing observation-duration description was requested after manuscript review and is identified as an addendum.

Required files below the authorized MIMIC-IV version 3.1 root:

```text
hosp/patients.csv.gz
hosp/admissions.csv.gz
icu/icustays.csv.gz
hosp/transfers.csv.gz
hosp/d_labitems.csv.gz
icu/d_items.csv.gz
hosp/labevents.csv.gz
icu/chartevents.csv.gz
SHA256SUMS.txt
```

`scripts/run_analysis.py` loads the data once for main/sensitivity/tables stages; the historical duration script reloads the required source tables. The `main` stage creates the provenance checks used by later stages. `all` runs main, sensitivity, tables and duration in order. Existing target outputs cause the public runner to stop rather than overwrite them. Stage-specific invocations must use outputs made with these same distributed modules; original historical hashes and public-copy hashes intentionally differ after path packaging.

The primary engine changes only its machine-specific input root to an environment configuration. All its function ASTs are identical to the historical source. The dependency module contains only the unchanged `classify_comments` and `unique_intervals` functions; unused exploratory analyses were omitted. Sensitivity, table and observation-duration scripts, protocol and execution specifications are byte-identical copies. `provenance/source_manifest.json` records both historical and public hashes, without private absolute paths.

## Local validation for this release

Compilation, import checks, artificial boundary tests, public-file inspection and aggregate figure reproduction are performed for release. Passing artificial tests does not validate the original cohort, access permissions, ethics, references or clinical interpretation. The restricted patient dataset is not reopened for this packaging task. Full extraction and bootstrap analysis are not rerun during this release and must not be reported as freshly reproduced.

The aggregate figure inputs remain byte-identical to their original publication-figure sources. They include the original flow counts and 16 overall/stratum-window summaries. Their denominators and Wilson intervals are checked independently by the plotting implementation. They cannot reproduce S5 patient-cluster resampling or reconstruct records.

## Data and availability boundaries

Source field names such as `subject_id` and `hadm_id` occur in code because they are database schema names. No actual values from those fields, patient dates, narrative specimen comments or fitted model objects are shipped. All examples in `tests/` are explicitly artificial. Any new restricted or derived patient datasets remain subject to the original data-use agreement; do not upload them to this public repository.

The source bundle and software tests do not resolve manuscript declarations still awaiting author confirmation. No institution, all-author approval, personal access credential or completed human adjudication is inferred from the presence of public code.
