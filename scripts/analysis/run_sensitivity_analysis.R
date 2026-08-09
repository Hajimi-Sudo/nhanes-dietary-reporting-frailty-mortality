#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survey)
  library(survival)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript run_sensitivity_analysis.R /path/to/project")

project_root <- normalizePath(args[[1]], mustWork = TRUE)
cohort <- readRDS(file.path(project_root, "data", "processed", "cohort.rds"))
result_root <- file.path(project_root, "results")
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
options(survey.lonely.psu = "adjust")

cohort$diet_variety_5 <- cohort$DR1TNUMF / 5
cohort$age_10 <- cohort$RIDAGEYR / 10
cohort$energy_1000 <- cohort$DR1TKCAL / 1000
cohort$sex <- factor(cohort$RIAGENDR)
cohort$race_ethnicity <- factor(cohort$RIDRETH3)
cohort$frail_same0 <- ifelse(is.na(cohort$FI_same0), NA, as.integer(cohort$FI_same0 >= 0.21))
cohort$frail_same1 <- ifelse(is.na(cohort$FI_same1), NA, as.integer(cohort$FI_same1 >= 0.21))

base_complete <- with(cohort, !is.na(DR1TNUMF) & !is.na(diet_variety_5) & !is.na(age_10) &
  !is.na(RIAGENDR) & !is.na(RIDRETH3) & !is.na(INDFMPIR) & !is.na(SDMVPSU) & !is.na(SDMVSTRA))

tidy_glm_term <- function(model, term, model_name, outcome_name, weight_name, n, events) {
  table <- coef(summary(model))
  estimate <- as.numeric(table[term, "Estimate"])
  standard_error <- as.numeric(table[term, "Std. Error"])
  p_value <- as.numeric(table[term, grep("^Pr", colnames(table))[1]])
  ci_low <- estimate - qnorm(0.975) * standard_error
  ci_high <- estimate + qnorm(0.975) * standard_error
  data.frame(
    model = model_name, outcome = outcome_name, weight = weight_name,
    term = term, n = n, events = events,
    estimate = estimate, odds_ratio = exp(estimate), standard_error = standard_error,
    p_value = p_value, odds_ratio_ci_low = exp(ci_low), odds_ratio_ci_high = exp(ci_high)
  )
}

run_logistic <- function(outcome_name, weight_name, model_name) {
  keep <- base_complete & !is.na(cohort[[outcome_name]]) & !is.na(cohort[[weight_name]]) & cohort[[weight_name]] > 0
  data <- cohort[keep, , drop = FALSE]
  design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = as.formula(paste0("~", weight_name)), nest = TRUE, data = data)
  model <- svyglm(as.formula(paste(outcome_name, "~ diet_variety_5 + age_10 + sex + race_ethnicity + INDFMPIR")), design = design, family = quasibinomial())
  tidy_glm_term(model, "diet_variety_5", model_name, outcome_name, weight_name, nrow(data), sum(data[[outcome_name]] == 1))
}

logistic_results <- rbind(
  run_logistic("frail_same0", "WTDRD1", "F4_same0_WTDRD1"),
  run_logistic("frail_same1", "WTDRD1", "F4_same1_WTDRD1"),
  run_logistic("frail_primary", "WTMEC2YR", "F4_primary_WTMEC2YR")
)

run_logistic_energy <- function() {
  keep <- base_complete & !is.na(cohort$frail_primary) & !is.na(cohort$FI_primary) &
    !is.na(cohort$DR1TKCAL) & cohort$DR1TKCAL > 0 & !is.na(cohort$WTDRD1) & cohort$WTDRD1 > 0
  data <- cohort[keep, , drop = FALSE]
  design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTDRD1, nest = TRUE, data = data)
  model <- svyglm(frail_primary ~ diet_variety_5 + energy_1000 + age_10 + sex + race_ethnicity + INDFMPIR,
    design = design, family = quasibinomial())
  tidy_glm_term(model, "diet_variety_5", "F4_primary_energy_WTDRD1", "frail_primary", "WTDRD1", nrow(data), sum(data$frail_primary == 1))
}

logistic_results <- rbind(logistic_results, run_logistic_energy())

tidy_cox_term <- function(model, term, model_name, followup_name, n, events) {
  table <- coef(summary(model))
  estimate <- as.numeric(table[term, "coef"])
  se_column <- if ("robust se" %in% colnames(table)) "robust se" else "se(coef)"
  standard_error <- as.numeric(table[term, se_column])
  p_value <- as.numeric(table[term, grep("^Pr", colnames(table))[1]])
  ci_low <- estimate - qnorm(0.975) * standard_error
  ci_high <- estimate + qnorm(0.975) * standard_error
  data.frame(
    model = model_name, followup = followup_name, term = term,
    n = n, events = events,
    estimate = estimate, hazard_ratio = exp(estimate), standard_error = standard_error,
    p_value = p_value, hazard_ratio_ci_low = exp(ci_low), hazard_ratio_ci_high = exp(ci_high)
  )
}

run_cox <- function(followup_name, model_name) {
  keep <- base_complete & !is.na(cohort$FI_primary) & !is.na(cohort$mortstat) &
    !is.na(cohort[[followup_name]]) & !is.na(cohort$WTDRD1) & cohort$WTDRD1 > 0 &
    cohort$eligstat == 1 & cohort[[followup_name]] > 0
  data <- cohort[keep, , drop = FALSE]
  data$event <- data$mortstat
  data$followup <- data[[followup_name]]
  design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTDRD1, nest = TRUE, data = data)
  model <- svycoxph(Surv(followup, event) ~ FI_primary + diet_variety_5 + age_10 + sex + race_ethnicity + INDFMPIR, design = design)
  rbind(
    tidy_cox_term(model, "FI_primary", model_name, followup_name, nrow(data), sum(data$event == 1)),
    tidy_cox_term(model, "diet_variety_5", model_name, followup_name, nrow(data), sum(data$event == 1))
  )
}

cox_results <- run_cox("permth_int", "F4_primary_PERMTH_INT")

run_cox_energy <- function() {
  keep <- base_complete & !is.na(cohort$FI_primary) & !is.na(cohort$mortstat) &
    !is.na(cohort$permth_exm) & cohort$permth_exm > 0 & !is.na(cohort$DR1TKCAL) &
    cohort$DR1TKCAL > 0 & !is.na(cohort$WTDRD1) & cohort$WTDRD1 > 0 & cohort$eligstat == 1
  data <- cohort[keep, , drop = FALSE]
  data$event <- data$mortstat
  data$followup <- data$permth_exm
  design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTDRD1, nest = TRUE, data = data)
  model <- svycoxph(Surv(followup, event) ~ FI_primary + diet_variety_5 + energy_1000 + age_10 + sex + race_ethnicity + INDFMPIR,
    design = design)
  rbind(
    tidy_cox_term(model, "FI_primary", "F4_primary_energy_WTDRD1", "permth_exm", nrow(data), sum(data$event == 1)),
    tidy_cox_term(model, "diet_variety_5", "F4_primary_energy_WTDRD1", "permth_exm", nrow(data), sum(data$event == 1))
  )
}

cox_results <- rbind(cox_results, run_cox_energy())

run_cox_parsimonious <- function() {
  keep <- base_complete & !is.na(cohort$FI_primary) & !is.na(cohort$mortstat) &
    !is.na(cohort$permth_exm) & cohort$permth_exm > 0 & !is.na(cohort$WTDRD1) &
    cohort$WTDRD1 > 0 & cohort$eligstat == 1
  data <- cohort[keep, , drop = FALSE]
  data$event <- data$mortstat
  data$followup <- data$permth_exm
  design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTDRD1, nest = TRUE, data = data)
  model <- svycoxph(Surv(followup, event) ~ FI_primary + diet_variety_5 + age_10 + sex + INDFMPIR,
    design = design)
  model_name <- "F4_mortality_parsimonious_PERMTH_EXM"
  rbind(
    tidy_cox_term(model, "FI_primary", model_name, "permth_exm", nrow(data), sum(data$event == 1)),
    tidy_cox_term(model, "diet_variety_5", model_name, "permth_exm", nrow(data), sum(data$event == 1))
  )
}

cox_results <- rbind(cox_results, run_cox_parsimonious())
write_csv(logistic_results, file.path(result_root, "F4_sensitivity_logistic.csv"))
write_csv(cox_results, file.path(result_root, "F4_sensitivity_cox.csv"))
writeLines(c(
  "SENSITIVITY_ANALYSIS_PASS",
  paste0("logistic_rows=", nrow(logistic_results)),
  paste0("cox_rows=", nrow(cox_results)),
  paste0("output_root=", result_root)
), file.path(project_root, "logs", "sensitivity_analysis.log"))
cat("SENSITIVITY_ANALYSIS_PASS\n")
cat("logistic_rows=", nrow(logistic_results), "\n", sep = "")
cat("cox_rows=", nrow(cox_results), "\n", sep = "")


