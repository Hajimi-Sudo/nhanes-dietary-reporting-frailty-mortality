#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survey)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript run_internal_validation.R /path/to/project")

project_root <- normalizePath(args[[1]], mustWork = TRUE)
cohort_path <- file.path(project_root, "02_data", "processed", "pilot_cohort_20260808.rds")
result_root <- file.path(project_root, "04_results", "internal_validation_20260812")
log_root <- file.path(project_root, "logs")
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(cohort_path)) stop("Missing cohort: ", cohort_path)

options(survey.lonely.psu = "adjust")
cohort <- readRDS(cohort_path)
cohort$diet_variety_5 <- cohort$DR1TNUMF / 5
cohort$age_10 <- cohort$RIDAGEYR / 10
required <- c("FI_primary", "frail_primary", "DR1TNUMF", "diet_variety_5", "age_10",
  "RIAGENDR", "RIDRETH3", "INDFMPIR", "WTDRD1", "SDMVPSU", "SDMVSTRA")
missing_columns <- setdiff(required, names(cohort))
if (length(missing_columns) > 0) stop("Missing columns: ", paste(missing_columns, collapse = ", "))

complete <- with(cohort, !is.na(frail_primary) & !is.na(FI_primary) &
  !is.na(diet_variety_5) & !is.na(age_10) & !is.na(RIAGENDR) & !is.na(RIDRETH3) &
  !is.na(INDFMPIR) & !is.na(WTDRD1) & WTDRD1 > 0 & !is.na(SDMVPSU) & !is.na(SDMVSTRA))
analysis <- cohort[complete, , drop = FALSE]
analysis$sex <- factor(analysis$RIAGENDR)
analysis$race_ethnicity <- factor(analysis$RIDRETH3)
analysis$psu_key <- interaction(analysis$SDMVSTRA, analysis$SDMVPSU, drop = TRUE, lex.order = TRUE)
analysis$stratum_key <- analysis$SDMVSTRA

weighted_auc <- function(y, p, w) {
  keep <- is.finite(y) & is.finite(p) & is.finite(w) & w > 0
  y <- as.numeric(y[keep]); p <- as.numeric(p[keep]); w <- as.numeric(w[keep])
  total_case <- sum(w[y == 1]); total_control <- sum(w[y == 0])
  if (total_case <= 0 || total_control <= 0) return(NA_real_)
  ord <- order(p)
  y <- y[ord]; p <- p[ord]; w <- w[ord]
  groups <- split(seq_along(p), p)
  controls_before <- 0
  concordance <- 0
  for (idx in groups) {
    case_inc <- sum(w[idx][y[idx] == 1])
    control_inc <- sum(w[idx][y[idx] == 0])
    concordance <- concordance + case_inc * (controls_before + control_inc / 2)
    controls_before <- controls_before + control_inc
  }
  as.numeric(concordance / (total_case * total_control))
}

weighted_brier <- function(y, p, w) {
  keep <- is.finite(y) & is.finite(p) & is.finite(w) & w > 0
  sum(w[keep] * (y[keep] - p[keep])^2) / sum(w[keep])
}

fit_validation_metrics <- function(train, validation, repeat_id) {
  train_design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTDRD1,
    nest = TRUE, data = train)
  validation_design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTDRD1,
    nest = TRUE, data = validation)
  formula_primary <- frail_primary ~ diet_variety_5 + age_10 + sex + race_ethnicity + INDFMPIR
  fit <- tryCatch(svyglm(formula_primary, design = train_design, family = quasibinomial()),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  predicted <- as.numeric(predict(fit, newdata = validation, type = "response"))
  predicted <- pmin(pmax(predicted, 1e-8), 1 - 1e-8)
  linear_predictor <- qlogis(predicted)
  validation$linear_predictor <- linear_predictor
  validation_design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTDRD1,
    nest = TRUE, data = validation)
  calibration <- tryCatch(svyglm(frail_primary ~ linear_predictor, design = validation_design,
    family = quasibinomial()), error = function(e) NULL)
  validation_fit <- tryCatch(svyglm(formula_primary, design = validation_design,
    family = quasibinomial()), error = function(e) NULL)
  calibration_table <- if (is.null(calibration)) NULL else coef(summary(calibration))
  validation_table <- if (is.null(validation_fit)) NULL else coef(summary(validation_fit))
  calibration_intercept <- if (!is.null(calibration) && "(Intercept)" %in% rownames(calibration_table)) calibration_table["(Intercept)", "Estimate"] else NA_real_
  calibration_slope <- if (!is.null(calibration) && "linear_predictor" %in% rownames(calibration_table)) calibration_table["linear_predictor", "Estimate"] else NA_real_
  validation_beta <- if (!is.null(validation_fit) && "diet_variety_5" %in% rownames(validation_table)) validation_table["diet_variety_5", "Estimate"] else NA_real_
  validation_se <- if (!is.null(validation_fit) && "diet_variety_5" %in% rownames(validation_table)) validation_table["diet_variety_5", "Std. Error"] else NA_real_
  data.frame(
    repeat_id = repeat_id,
    train_n = nrow(train),
    validation_n = nrow(validation),
    train_frailty_events = sum(train$frail_primary == 1),
    validation_frailty_events = sum(validation$frail_primary == 1),
    validation_weighted_auc = weighted_auc(validation$frail_primary, predicted, validation$WTDRD1),
    validation_weighted_brier = weighted_brier(validation$frail_primary, predicted, validation$WTDRD1),
    calibration_intercept = calibration_intercept,
    calibration_slope = calibration_slope,
    validation_or_per_five = exp(validation_beta),
    validation_or_low = exp(validation_beta - qnorm(.975) * validation_se),
    validation_or_high = exp(validation_beta + qnorm(.975) * validation_se),
    row.names = NULL
  )
}

set.seed(20260812)
strata_values <- sort(unique(analysis$SDMVSTRA))
repeats <- vector("list", 50)
for (r in seq_len(50)) {
  holdout_keys <- vapply(strata_values, function(s) {
    candidates <- unique(as.character(analysis$psu_key[analysis$SDMVSTRA == s]))
    sample(candidates, size = 1)
  }, character(1))
  validation_rows <- as.character(analysis$psu_key) %in% holdout_keys
  train <- analysis[!validation_rows, , drop = FALSE]
  validation <- analysis[validation_rows, , drop = FALSE]
  repeats[[r]] <- tryCatch(fit_validation_metrics(train, validation, r), error = function(e) NULL)
}
results <- do.call(rbind, repeats[!vapply(repeats, is.null, logical(1))])
if (is.null(results) || nrow(results) < 40) stop("VALIDATION_GATE_FAIL: fewer than 40 successful repeats")
write_csv(results, file.path(result_root, "internal_validation_repeats.csv"))

metric_names <- c("validation_weighted_auc", "validation_weighted_brier", "calibration_intercept",
  "calibration_slope", "validation_or_per_five")
summary_rows <- lapply(metric_names, function(metric) {
  x <- results[[metric]]
  data.frame(metric = metric, successful_repeats = sum(is.finite(x)),
    mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE),
    q025 = as.numeric(quantile(x, .025, na.rm = TRUE)),
    median = median(x, na.rm = TRUE), q975 = as.numeric(quantile(x, .975, na.rm = TRUE)),
    row.names = NULL)
})
summary_table <- do.call(rbind, summary_rows)
summary_table$share_below_null <- vapply(summary_table$metric, function(metric) {
  if (metric != "validation_or_per_five") return(NA_real_)
  mean(results[[metric]] < 1, na.rm = TRUE)
}, numeric(1))
write_csv(summary_table, file.path(result_root, "internal_validation_summary.csv"))

contract <- data.frame(
  metric = c("complete_case_n", "strata_n", "psu_per_stratum", "repeats_requested", "repeats_successful",
    "split_unit", "weight", "validation_scope"),
  value = c(nrow(analysis), length(strata_values), "2", 50, nrow(results),
    "one PSU held out within every SDMVSTRA", "WTDRD1",
    "internal transportability and calibration check; not causal validation"),
  row.names = NULL
)
write_csv(contract, file.path(result_root, "internal_validation_contract.csv"))
writeLines(c(
  "INTERNAL_VALIDATION_PASS",
  paste0("complete_case_n=", nrow(analysis)),
  paste0("successful_repeats=", nrow(results)),
  paste0("validation_result_root=", result_root),
  "split_unit=one PSU held out within every SDMVSTRA",
  "interpretation=internal transportability and calibration check; not causal validation"
), file.path(log_root, "internal_validation_20260812.log"))
cat("INTERNAL_VALIDATION_PASS\n")
cat("complete_case_n=", nrow(analysis), "\n", sep = "")
cat("successful_repeats=", nrow(results), "\n", sep = "")
print(summary_table)
