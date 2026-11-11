# Potassium retesting after a first recorded hyperkalemic result in intensive care

Research code accompanying the descriptive MIMIC-IV version 3.1 study of documented potassium retesting after a first qualifying recorded chemistry potassium result. This repository preserves the version 1.0 analysis definitions. The manuscript is being prepared for submission; no paper DOI or software archival DOI has been assigned.

## What is shared

| Component | Contents | Verification boundary |
| --- | --- | --- |
| Clinical analysis source | Cohort selection, first recorded high result, ICU continuity, primary outcomes, five planned sensitivity analyses, baseline and 2/4/6/24-hour windows | Source functions preserved; local packaging did not rerun the full restricted-data extraction |
| Observation-duration addendum | Existing 24-hour capped observation-duration description | Post-review reporting addition, not prospectively prespecified |
| Figure reproduction | The original plotting implementation and three frozen aggregate source files | Recreates flow and retesting figures without patient data; checks counts and Wilson intervals |
| Offline rule tests | Artificial records authored for this release | Tests software boundary behavior, not clinical adjudication or patient-level cohort validation |

No patient-level records, identifiers, individual timestamps, source comments, restricted datasets, credentials, manuscript Word files or author-contact forms are included. The three figure-source files contain group-level counts, proportions and intervals only. They cannot reconstruct cohort extraction, patient trajectories or bootstrap resampling.

## Reproduce the figures without MIMIC access

Use Python 3.12 and a separate figure environment:

```bash
python -m venv .venv-figures
# Activate the environment using your operating system's command.
python -m pip install -r requirements-figure-validation.txt
python scripts/validate_public_bundle.py
python scripts/reproduce_figures.py --source-dir data/figure_source --output-dir reproduced_figures
python scripts/validate_figures.py --figure-dir reproduced_figures
```

Figures are exported as PNG, PDF, SVG and TIFF. Point positions and intervals use the frozen summaries without fitting a model. Font availability can change layout; Arial has a DejaVu Sans fallback. Analysis and plotting originally used different NumPy versions; keep their environments separate.

## Review or run the analysis source

MIMIC-IV must be obtained independently by appropriately credentialed researchers from [PhysioNet](https://physionet.org/content/mimiciv/3.1/). The repository grants no access to the data.

```bash
python -m venv .venv-analysis
# Activate this environment instead of the figure environment.
python -m pip install -r requirements-analysis.txt
python -m unittest discover -s tests -v
python scripts/run_analysis.py --help
python scripts/run_analysis.py --data-root /path/to/authorized/mimic-iv-3.1 --stage all
```

The source directory must include the eight compressed tables listed in `docs/REPRODUCIBILITY.md` and the official `SHA256SUMS.txt`. Use an independent working copy with no pre-existing analysis outputs. This is a large in-memory CSV workflow, scanning the original lab and chart tables; it is not a lightweight download-and-run demo. No patient data are fetched by these commands. Outputs remain local and are excluded from Git by default.

The executable source has cohort-size and hash checks specific to the original version 3.1 release. They are intentional safeguards, not a method for adapting to another release. This public package has not been independently validated as an end-to-end reproduction from a fresh restricted database. See [reproducibility and provenance](docs/REPRODUCIBILITY.md).

## Interpretation and reuse

The primary denominator retains all enrolled patients. First-repeat high-potassium detection is conditional on sampling and first-value interpretability. A later clean value does not replace an uninterpretable first candidate. The study describes recorded processes; it does not establish guideline adherence, nursing omission, treatment effectiveness or causality.

No open-source license has been selected for this release. Public availability is not a grant of rights to redistribute restricted MIMIC data. Cite the repository and fixed release when referring to the code. The [citation file](CITATION.cff) contains the author order confirmed for this manuscript.

Release v1.0.1 fixes Git line-ending preservation so distributed file hashes match the original aggregate inputs. Use v1.0.1; the initial upload is superseded.
