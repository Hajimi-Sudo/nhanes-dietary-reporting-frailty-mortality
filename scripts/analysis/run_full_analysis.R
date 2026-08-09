#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survey)
  library(survival)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript run_full_analysis.R /path/to/project")

project_root <- normalizePath(args[[1]], mustWork = TRUE)
cohort_path <- file.path(project_root, "data", "processed", "cohort.rds")
result_root <- file.path(project_root, "results")
log_root <- file.path(project_root, "logs")
if (!file.exists(cohort_path)) stop("Pilot cohort is missing: ", cohort_path)
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)

options(survey.lonely.psu = "adjust")
cohort <- readRDS(cohort_path)
cohort$diet_variety_5 <- cohort$DR1TNUMF / 5
cohort$age_10 <- cohort$RIDAGEYR / 10

required <- c(
  "SEQN", "FI_primary", "frail_primary", "DR1TNUMF", "diet_variety_5", "age_10",
  "RIAGENDR", "RIDRETH3", "INDFMPIR", "WTDRD1", "WTMEC2YR", "SDMVPSU", "SDMVSTRA",
  "mortstat", "permth_exm", "mortality_exm_eligible"
)
missing_columns <- setdiff(required, names(cohort))
if (length(missing_columns) > 0) stop("Missing cohort columns: ", paste(missing_columns, collapse = ", "))

complete_for_frailty <- with(cohort, !is.na(frail_primary) & !is.na(FI_primary) &
  !is.na(diet_variety_5) & !is.na(age_10) & !is.na(RIAGENDR) & !is.na(RIDRETH3) &
  !is.na(INDFMPIR) & !is.na(WTDRD1) & WTDRD1 > 0 & !is.na(SDMVPSU) & !is.na(SDMVSTRA))
analysis <- cohort[complete_for_frailty, , drop = FALSE]
analysis$sex <- factor(analysis$RIAGENDR)
analysis$race_ethnicity <- factor(analysis$RIDRETH3)

if (nrow(analysis) < 1500) stop("FULL_ANALYSIS_GATE_FAIL: primary analytic sample below 1500")
design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTDRD1,
  nest = TRUE,
  data = analysis
)

ci_from <- function(estimates, standard_errors, z = qnorm(0.975)) {
  cbind(ci_low = estimates - z * standard_errors, ci_high = estimates + z * standard_errors)
}

mean_terms <- c("DR1TNUMF", "FI_primary", "frail_primary")
mean_formula <- as.formula(paste("~", paste(mean_terms, collapse = " + ")))
mean_fit <- svymean(mean_formula, design, na.rm = TRUE)
mean_est <- as.numeric(coef(mean_fit))
mean_se <- as.numeric(SE(mean_fit))
descriptive <- data.frame(
  measure = mean_terms,
  estimate = mean_est,
  standard_error = mean_se,
  ci_from(mean_est, mean_se),
  row.names = NULL
)
write_csv(descriptive, file.path(result_root, "F1_weighted_descriptive.csv"))

frailty_model <- svyglm(
  frail_primary ~ diet_variety_5 + age_10 + sex + race_ethnicity + INDFMPIR,
  design = design,
  family = quasibinomial()
)

tidy_glm <- function(model, model_name) {
  table <- coef(summary(model))
  estimate <- as.numeric(table[, "Estimate"])
  standard_error <- as.numeric(table[, "Std. Error"])
  statistic <- as.numeric(table[, grep("value$", colnames(table))[1]])
  p_value <- as.numeric(table[, grep("^Pr", colnames(table))[1]])
  ci <- ci_from(estimate, standard_error)
  data.frame(
    model = model_name,
    term = rownames(table),
    estimate = estimate,
    odds_ratio = exp(estimate),
    standard_error = standard_error,
    statistic = statistic,
    p_value = p_value,
    ci_low = ci[, "ci_low"],
    ci_high = ci[, "ci_high"],
    odds_ratio_ci_low = exp(ci[, "ci_low"]),
    odds_ratio_ci_high = exp(ci[, "ci_high"]),
    row.names = NULL
  )
}
F2 <- tidy_glm(frailty_model, "F2_svy_logistic_frailty")
write_csv(F2, file.path(result_root, "F2_svy_logistic_frailty.csv"))

zero_followup_excluded_n <- with(analysis, mortality_exm_eligible & !is.na(permth_exm) & permth_exm <= 0)
mortality_complete <- with(analysis, mortality_exm_eligible & !is.na(mortstat) & !is.na(permth_exm) & permth_exm > 0)
mortality_analysis <- analysis[mortality_complete, , drop = FALSE]
if (nrow(mortality_analysis) < 1500 || sum(mortality_analysis$mortstat == 1) < 50) {
  stop("FULL_ANALYSIS_GATE_FAIL: mortality analysis sample or events below prespecified minimum")
}
mortality_design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTDRD1,
  nest = TRUE,
  data = mortality_analysis
)

cox_model <- svycoxph(
  Surv(permth_exm, mortstat) ~ FI_primary + diet_variety_5 + age_10 + sex + race_ethnicity + INDFMPIR,
  design = mortality_design
)

tidy_cox <- function(model, model_name) {
  table <- coef(summary(model))
  estimate <- as.numeric(table[, "coef"])
  se_column <- if ("robust se" %in% colnames(table)) "robust se" else "se(coef)"
  standard_error <- as.numeric(table[, se_column])
  z <- as.numeric(table[, "z"])
  p_value <- as.numeric(table[, grep("^Pr", colnames(table))[1]])
  ci <- ci_from(estimate, standard_error)
  data.frame(
    model = model_name,
    term = rownames(table),
    estimate = estimate,
    hazard_ratio = exp(estimate),
    standard_error = standard_error,
    z = z,
    p_value = p_value,
    ci_low = ci[, "ci_low"],
    ci_high = ci[, "ci_high"],
    hazard_ratio_ci_low = exp(ci[, "ci_low"]),
    hazard_ratio_ci_high = exp(ci[, "ci_high"]),
    row.names = NULL
  )
}
F3 <- tidy_cox(cox_model, "F3_svy_cox_mortality")
write_csv(F3, file.path(result_root, "F3_svy_cox_mortality.csv"))
F3_reporting <- F3[F3$term %in% c("FI_primary", "diet_variety_5"), , drop = FALSE]
F3_reporting$reporting_scale <- ifelse(F3_reporting$term == "FI_primary", "per 0.1 FI", "per 5 foods/beverages")
F3_reporting$hazard_ratio_reporting <- ifelse(
  F3_reporting$term == "FI_primary", exp(F3_reporting$estimate * 0.1), F3_reporting$hazard_ratio
)
F3_reporting$hazard_ratio_ci_low_reporting <- ifelse(
  F3_reporting$term == "FI_primary", exp(F3_reporting$ci_low * 0.1), F3_reporting$hazard_ratio_ci_low
)
F3_reporting$hazard_ratio_ci_high_reporting <- ifelse(
  F3_reporting$term == "FI_primary", exp(F3_reporting$ci_high * 0.1), F3_reporting$hazard_ratio_ci_high
)
write_csv(F3_reporting, file.path(result_root, "F3_reporting_scales.csv"))

contract <- data.frame(
  metric = c(
    "primary_model_n", "primary_model_frail_events", "mortality_model_n",
    "mortality_model_events", "mortality_followup_min_months", "mortality_followup_max_months",
    "zero_followup_excluded_n", "weight", "same_HUQ020_primary_score", "claim_boundary"
  ),
  value = c(
    nrow(analysis), sum(analysis$frail_primary == 1), nrow(mortality_analysis),
    sum(mortality_analysis$mortstat == 1), min(mortality_analysis$permth_exm),
    max(mortality_analysis$permth_exm), sum(zero_followup_excluded_n), "WTDRD1", "0.5", "association only; no causal claim"
  )
)
write_csv(contract, file.path(result_root, "analysis_contract.csv"))
writeLines(c(
  "FULL_ANALYSIS_PASS",
  paste0("primary_model_n=", nrow(analysis)),
  paste0("primary_model_frail_events=", sum(analysis$frail_primary == 1)),
  paste0("mortality_model_n=", nrow(mortality_analysis)),
  paste0("mortality_model_events=", sum(mortality_analysis$mortstat == 1)),
  paste0("result_root=", result_root)
), file.path(log_root, "full_analysis.log"))
cat("FULL_ANALYSIS_PASS\n")
cat("primary_model_n=", nrow(analysis), "\n", sep = "")
cat("mortality_model_n=", nrow(mortality_analysis), "\n", sep = "")
cat("mortality_model_events=", sum(mortality_analysis$mortstat == 1), "\n", sep = "")


