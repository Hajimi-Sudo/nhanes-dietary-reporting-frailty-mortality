#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survey)
  library(survival)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript run_mortality_no_fi_sensitivity.R /path/to/project")

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
cohort$sex <- factor(cohort$RIAGENDR)
cohort$race_ethnicity <- factor(cohort$RIDRETH3)

required <- c(
  "DR1TNUMF", "diet_variety_5", "age_10", "RIAGENDR", "RIDRETH3", "INDFMPIR",
  "WTDRD1", "SDMVPSU", "SDMVSTRA", "mortstat", "permth_exm",
  "mortality_exm_eligible"
)
missing_columns <- setdiff(required, names(cohort))
if (length(missing_columns) > 0) stop("Missing cohort columns: ", paste(missing_columns, collapse = ", "))

keep <- with(cohort,
  !is.na(DR1TNUMF) & !is.na(diet_variety_5) & !is.na(age_10) &
    !is.na(RIAGENDR) & !is.na(RIDRETH3) & !is.na(INDFMPIR) &
    !is.na(WTDRD1) & WTDRD1 > 0 & !is.na(SDMVPSU) & !is.na(SDMVSTRA) &
    !is.na(mortality_exm_eligible) & mortality_exm_eligible &
    !is.na(mortstat) & !is.na(permth_exm) & permth_exm > 0
)
analysis <- cohort[keep, , drop = FALSE]
if (nrow(analysis) < 1500 || sum(analysis$mortstat == 1) < 50) {
  stop("MORTALITY_NO_FI_GATE_FAIL: sample or event count below prespecified minimum")
}

design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTDRD1,
  nest = TRUE,
  data = analysis
)

model <- svycoxph(
  Surv(permth_exm, mortstat) ~ diet_variety_5 + age_10 + sex + race_ethnicity + INDFMPIR,
  design = design
)

summary_table <- coef(summary(model))
term <- "diet_variety_5"
estimate <- as.numeric(summary_table[term, "coef"])
se_column <- if ("robust se" %in% colnames(summary_table)) "robust se" else "se(coef)"
standard_error <- as.numeric(summary_table[term, se_column])
ci_low <- estimate - qnorm(0.975) * standard_error
ci_high <- estimate + qnorm(0.975) * standard_error

result <- data.frame(
  model = "F4_mortality_no_FI_PERMTH_EXM",
  term = term,
  n = nrow(analysis),
  events = sum(analysis$mortstat == 1),
  weight = "WTDRD1",
  followup = "PERMTH_EXM",
  adjustment = "age + sex + race/ethnicity + PIR; FI omitted",
  estimate = estimate,
  hazard_ratio = exp(estimate),
  standard_error = standard_error,
  p_value = as.numeric(summary_table[term, grep("^Pr", colnames(summary_table))[1]]),
  hazard_ratio_ci_low = exp(ci_low),
  hazard_ratio_ci_high = exp(ci_high),
  row.names = NULL
)
write_csv(result, file.path(result_root, "F4_mortality_no_FI.csv"))
writeLines(c(
  "MORTALITY_NO_FI_PASS",
  paste0("n=", nrow(analysis)),
  paste0("events=", sum(analysis$mortstat == 1)),
  paste0("hazard_ratio_per_5_items=", sprintf("%.10f", result$hazard_ratio)),
  paste0("ci_low=", sprintf("%.10f", result$hazard_ratio_ci_low)),
  paste0("ci_high=", sprintf("%.10f", result$hazard_ratio_ci_high)),
  paste0("p_value=", sprintf("%.10f", result$p_value)),
  "adjustment=age + sex + race/ethnicity + PIR; FI omitted",
  "weight=WTDRD1",
  "followup=PERMTH_EXM"
), file.path(log_root, "mortality_no_FI.log"))
cat("MORTALITY_NO_FI_PASS\n")
cat("n=", nrow(analysis), "\n", sep = "")
cat("events=", sum(analysis$mortstat == 1), "\n", sep = "")
cat("hazard_ratio_per_5_items=", sprintf("%.10f", result$hazard_ratio), "\n", sep = "")
cat("p_value=", sprintf("%.10f", result$p_value), "\n", sep = "")


