#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survey)
  library(survival)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript run_selection_and_quantile_sensitivity.R /path/to/project")

project_root <- normalizePath(args[[1]], mustWork = TRUE)
cohort_path <- file.path(project_root, "data", "processed", "cohort.rds")
result_root <- file.path(project_root, "results")
log_root <- file.path(project_root, "logs")
if (!file.exists(cohort_path)) stop("Pilot cohort is missing: ", cohort_path)
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)

options(survey.lonely.psu = "adjust")
cohort <- readRDS(cohort_path)
cohort$age_10 <- cohort$RIDAGEYR / 10
cohort$sex <- factor(cohort$RIAGENDR)
cohort$race_ethnicity <- factor(cohort$RIDRETH3)

required <- c(
  "SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3", "INDFMPIR", "DR1TNUMF",
  "FI_primary", "frail_primary", "WTDRD1", "WTMEC2YR", "SDMVPSU", "SDMVSTRA",
  "mortstat", "permth_exm", "mortality_exm_eligible"
)
missing_columns <- setdiff(required, names(cohort))
if (length(missing_columns) > 0) stop("Missing cohort columns: ", paste(missing_columns, collapse = ", "))

frailty_keep <- with(cohort,
  !is.na(frail_primary) & !is.na(FI_primary) & !is.na(DR1TNUMF) &
    !is.na(RIDAGEYR) & !is.na(RIAGENDR) & !is.na(RIDRETH3) & !is.na(INDFMPIR) &
    !is.na(WTDRD1) & WTDRD1 > 0 & !is.na(SDMVPSU) & !is.na(SDMVSTRA)
)

# Selection audit: report observed values and missingness without adding p-values.
selection_numeric <- c("RIDAGEYR", "DR1TNUMF", "FI_primary", "INDFMPIR", "WTDRD1")
selection_categorical <- c("RIAGENDR", "RIDRETH3")
selection_rows <- list()
row_index <- 1L
groups <- c("primary_complete_case", "excluded_from_primary")

for (variable in selection_numeric) {
  values <- cohort[[variable]]
  observed <- !is.na(values)
  for (group_name in groups) {
    in_group <- if (group_name == "primary_complete_case") frailty_keep else !frailty_keep
    group_values <- values[in_group]
    group_observed <- !is.na(group_values)
    selection_rows[[row_index]] <- data.frame(
      variable = variable,
      level = "numeric",
      group = group_name,
      n_total = sum(in_group),
      n_observed = sum(group_observed),
      missing_percent = 100 * mean(!group_observed),
      mean_or_proportion = if (any(group_observed)) mean(group_values[group_observed]) else NA_real_,
      standard_deviation = if (sum(group_observed) > 1) sd(group_values[group_observed]) else NA_real_,
      row.names = NULL
    )
    row_index <- row_index + 1L
  }
}

for (variable in selection_categorical) {
  values <- cohort[[variable]]
  levels_present <- sort(unique(values[!is.na(values)]))
  for (level_value in levels_present) {
    for (group_name in groups) {
      in_group <- if (group_name == "primary_complete_case") frailty_keep else !frailty_keep
      group_values <- values[in_group]
      group_observed <- !is.na(group_values)
      selection_rows[[row_index]] <- data.frame(
        variable = variable,
        level = as.character(level_value),
        group = group_name,
        n_total = sum(in_group),
        n_observed = sum(group_observed),
        missing_percent = 100 * mean(!group_observed),
        mean_or_proportion = if (any(group_observed)) mean(group_values[group_observed] == level_value) else NA_real_,
        standard_deviation = NA_real_,
        row.names = NULL
      )
      row_index <- row_index + 1L
    }
  }
}
selection_comparison <- do.call(rbind, selection_rows)
write_csv(selection_comparison, file.path(result_root, "selection_comparison.csv"))

exposure_values <- cohort$DR1TNUMF[!is.na(cohort$DR1TNUMF)]
quartile_cutpoints <- as.numeric(quantile(exposure_values, probs = c(0, 0.25, 0.5, 0.75, 1), names = FALSE, type = 7))
if (length(unique(quartile_cutpoints)) != 5) stop("QUARTILE_GATE_FAIL: DR1TNUMF has non-unique quartile cutpoints")
quartile_labels <- c("Q1", "Q2", "Q3", "Q4")
cohort$exposure_q4 <- cut(
  cohort$DR1TNUMF,
  breaks = quartile_cutpoints,
  include.lowest = TRUE,
  labels = quartile_labels,
  right = TRUE
)
write_csv(data.frame(
  quantile = c("0%", "25%", "50%", "75%", "100%"),
  cutpoint = quartile_cutpoints,
  row.names = NULL
), file.path(result_root, "exposure_quartile_cutpoints.csv"))

tidy_quartile_glm <- function(model, data, model_name, design_df) {
  table <- coef(summary(model))
  terms <- grep("^exposure_q4", rownames(table), value = TRUE)
  critical_value <- if (is.finite(design_df) && design_df > 0) qt(0.975, df = design_df) else qnorm(0.975)
  rows <- lapply(terms, function(term) {
    estimate <- as.numeric(table[term, "Estimate"])
    standard_error <- as.numeric(table[term, "Std. Error"])
    ci_low <- estimate - critical_value * standard_error
    ci_high <- estimate + critical_value * standard_error
    data.frame(
      model = model_name,
      term = term,
      n = nrow(data),
      events = sum(data$frail_primary == 1),
      design_df = design_df,
      odds_ratio = exp(estimate),
      odds_ratio_ci_low = exp(ci_low),
      odds_ratio_ci_high = exp(ci_high),
      p_value = 2 * pt(-abs(estimate / standard_error), df = design_df),
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}

frailty_data <- cohort[frailty_keep & !is.na(cohort$exposure_q4), , drop = FALSE]
frailty_design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTDRD1,
  nest = TRUE,
  data = frailty_data
)
frailty_q4_model <- svyglm(
  frail_primary ~ exposure_q4 + age_10 + sex + race_ethnicity + INDFMPIR,
  design = frailty_design,
  family = quasibinomial()
)
frailty_q4_results <- tidy_quartile_glm(
  frailty_q4_model,
  frailty_data,
  "F4_exposure_quartile_frailty_WTDRD1",
  degf(frailty_design)
)
frailty_q4_results$outcome <- "frail_primary"
frailty_q4_results$weight <- "WTDRD1"
write_csv(frailty_q4_results, file.path(result_root, "F4_exposure_quartile_frailty.csv"))

tidy_quartile_cox <- function(model, data, model_name) {
  table <- coef(summary(model))
  terms <- grep("^exposure_q4", rownames(table), value = TRUE)
  rows <- lapply(terms, function(term) {
    estimate <- as.numeric(table[term, "coef"])
    se_column <- if ("robust se" %in% colnames(table)) "robust se" else "se(coef)"
    standard_error <- as.numeric(table[term, se_column])
    ci_low <- estimate - qnorm(0.975) * standard_error
    ci_high <- estimate + qnorm(0.975) * standard_error
    data.frame(
      model = model_name,
      term = term,
      n = nrow(data),
      events = sum(data$mortstat == 1),
      hazard_ratio = exp(estimate),
      hazard_ratio_ci_low = exp(ci_low),
      hazard_ratio_ci_high = exp(ci_high),
      p_value = as.numeric(table[term, grep("^Pr", colnames(table))[1]]),
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}

mortality_keep <- with(cohort,
  !is.na(DR1TNUMF) & !is.na(exposure_q4) & !is.na(RIDAGEYR) &
    !is.na(RIAGENDR) & !is.na(RIDRETH3) & !is.na(INDFMPIR) &
    !is.na(WTDRD1) & WTDRD1 > 0 & !is.na(SDMVPSU) & !is.na(SDMVSTRA) &
    !is.na(mortality_exm_eligible) & mortality_exm_eligible &
    !is.na(mortstat) & !is.na(permth_exm) & permth_exm > 0
)
mortality_data <- cohort[mortality_keep, , drop = FALSE]
mortality_design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTDRD1,
  nest = TRUE,
  data = mortality_data
)
mortality_q4_model <- svycoxph(
  Surv(permth_exm, mortstat) ~ exposure_q4 + age_10 + sex + race_ethnicity + INDFMPIR,
  design = mortality_design
)
mortality_q4_results <- tidy_quartile_cox(mortality_q4_model, mortality_data, "F4_exposure_quartile_mortality_no_FI")
mortality_q4_results$outcome <- "mortstat"
mortality_q4_results$weight <- "WTDRD1"
mortality_q4_results$followup <- "PERMTH_EXM"
mortality_q4_results$adjustment <- "age + sex + race/ethnicity + PIR; FI omitted"
write_csv(mortality_q4_results, file.path(result_root, "F4_exposure_quartile_mortality_no_FI.csv"))

quartile_count_rows <- lapply(c("frailty", "mortality_no_FI"), function(outcome_name) {
  data <- if (outcome_name == "frailty") frailty_data else mortality_data
  event_values <- if (outcome_name == "frailty") data$frail_primary == 1 else data$mortstat == 1
  levels_present <- quartile_labels[quartile_labels %in% as.character(data$exposure_q4)]
  data.frame(
    exposure_q4 = levels_present,
    n = vapply(levels_present, function(level) sum(as.character(data$exposure_q4) == level), integer(1)),
    events = vapply(levels_present, function(level) sum(as.character(data$exposure_q4) == level & event_values), integer(1)),
    outcome = outcome_name,
    row.names = NULL
  )
})
write_csv(do.call(rbind, quartile_count_rows), file.path(result_root, "exposure_quartile_counts.csv"))

writeLines(c(
  "SELECTION_AND_QUANTILE_SENSITIVITY_PASS",
  paste0("source_population_n=", nrow(cohort)),
  paste0("frailty_complete_case_n=", sum(frailty_keep)),
  paste0("frailty_quartile_model_n=", nrow(frailty_data)),
  paste0("frailty_quartile_events=", sum(frailty_data$frail_primary == 1)),
  paste0("mortality_no_FI_quartile_model_n=", nrow(mortality_data)),
  paste0("mortality_no_FI_quartile_events=", sum(mortality_data$mortstat == 1)),
  paste0("quartile_cutpoints=", paste(quartile_cutpoints, collapse = ",")),
  "selection_comparison=unweighted observed-value audit; no p-values",
  "frailty_quartile_weight=WTDRD1",
  "mortality_quartile_weight=WTDRD1",
  "mortality_quartile_followup=PERMTH_EXM",
  "mortality_quartile_adjustment=age + sex + race/ethnicity + PIR; FI omitted"
), file.path(log_root, "selection_and_quantile_sensitivity.log"))
cat("SELECTION_AND_QUANTILE_SENSITIVITY_PASS\n")
cat("source_population_n=", nrow(cohort), "\n", sep = "")
cat("frailty_quartile_model_n=", nrow(frailty_data), "\n", sep = "")
cat("mortality_no_FI_quartile_model_n=", nrow(mortality_data), "\n", sep = "")
cat("mortality_no_FI_quartile_events=", sum(mortality_data$mortstat == 1), "\n", sep = "")


