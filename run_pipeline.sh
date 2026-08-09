#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python scripts/acquisition/download_nhanes_2017_2018.py "$REPO_ROOT"
Rscript scripts/cleaning/frailty_component_coverage_audit.R "$REPO_ROOT"
Rscript scripts/cleaning/variable_audit.R "$REPO_ROOT"
Rscript scripts/analysis/pilot_build_cohort.R "$REPO_ROOT"
Rscript scripts/analysis/run_full_analysis.R "$REPO_ROOT"
Rscript scripts/analysis/run_sensitivity_analysis.R "$REPO_ROOT"
Rscript scripts/analysis/run_selection_and_quantile_sensitivity.R "$REPO_ROOT"
Rscript scripts/analysis/run_mortality_no_fi_sensitivity.R "$REPO_ROOT"
Rscript scripts/analysis/run_model_diagnostics.R "$REPO_ROOT"
Rscript scripts/analysis/run_result_visualization_summaries.R "$REPO_ROOT"
python scripts/reporting/generate_manuscript_outputs.py
