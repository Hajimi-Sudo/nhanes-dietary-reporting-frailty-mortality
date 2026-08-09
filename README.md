# NHANES Dietary Reporting, Frailty, and Mortality

Recommended GitHub repository name: `nhanes-dietary-reporting-frailty-mortality`

This repository contains the reproducible analysis code for the manuscript:

> An Estimand-Aware Survey-Weighted Framework for Linking Dietary Reporting, Frailty, and Mortality in Older Adults: An NHANES Analysis

The workflow uses NHANES 2017--2018 files and the corresponding linked mortality file. It constructs the prespecified Frailty Index, performs survey-weighted logistic and Cox models, runs sensitivity analyses, and generates manuscript tables and figures.

## Repository Boundary

Raw participant-level data are not stored in this repository. The download script retrieves the files directly from the official source URLs listed in `data/CATALOG.csv`. Generated participant-level intermediate data, result files, logs, and figures are ignored by Git unless explicitly selected for release.

The analysis is observational. The code and manuscript do not support causal or protective-effect claims.

## Requirements

- R 4.1.2 or newer with the packages listed in `requirements-r.txt`.
- Python 3.10 or newer with the packages listed in `requirements-python.txt`.
- PowerShell is optional; the cross-platform Python downloader is recommended.

## Reproduction

From the repository root, run:

```bash
bash run_pipeline.sh
```

The pipeline downloads the source files, audits the mapped variables, builds the participant-level cohort, runs the primary and sensitivity analyses, performs diagnostics, and generates the reporting tables and figures.

Individual stages can be run with:

```bash
python scripts/acquisition/download_nhanes_2017_2018.py .
Rscript scripts/cleaning/frailty_component_coverage_audit.R .
Rscript scripts/cleaning/variable_audit.R .
Rscript scripts/analysis/pilot_build_cohort.R .
Rscript scripts/analysis/run_full_analysis.R .
Rscript scripts/analysis/run_sensitivity_analysis.R .
Rscript scripts/analysis/run_selection_and_quantile_sensitivity.R .
Rscript scripts/analysis/run_mortality_no_fi_sensitivity.R .
Rscript scripts/analysis/run_model_diagnostics.R .
Rscript scripts/analysis/run_result_visualization_summaries.R .
python scripts/reporting/generate_manuscript_outputs.py
```

The scripts accept the repository root as their only positional argument so that they can run on a server without hard-coded local paths.

Server-specific environment and run-record guidance is in [`SERVER_RUN.md`](SERVER_RUN.md).

## Data Sources

NHANES documentation: <https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2017>

Linked mortality documentation: <https://www.cdc.gov/nchs/data-linkage/mortality-public.htm>

## Citation

The repository URL and a versioned archive DOI should be added here after the GitHub repository and Zenodo release are created.
