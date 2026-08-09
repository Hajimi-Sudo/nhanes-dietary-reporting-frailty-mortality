#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survey)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript run_result_visualization_summaries.R /path/to/project")

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
  "RIDAGEYR", "RIAGENDR", "RIDRETH3", "INDFMPIR", "DR1TNUMF", "FI_primary",
  "frail_primary", "WTDRD1", "SDMVPSU", "SDMVSTRA", "mortstat", "permth_exm",
  "mortality_exm_eligible"
)
missing_columns <- setdiff(required, names(cohort))
if (length(missing_columns) > 0) stop("Missing columns: ", paste(missing_columns, collapse = ", "))

exposure_values <- cohort$DR1TNUMF[!is.na(cohort$DR1TNUMF)]
cutpoints <- as.numeric(quantile(exposure_values, probs = c(0, 0.25, 0.5, 0.75, 1), names = FALSE, type = 7))
if (length(unique(cutpoints)) != 5) stop("VISUALIZATION_GATE_FAIL: non-unique exposure quartile cutpoints")
quartile_labels <- c("Q1", "Q2", "Q3", "Q4")
cohort$exposure_q4 <- cut(
  cohort$DR1TNUMF,
  breaks = cutpoints,
  include.lowest = TRUE,
  labels = quartile_labels,
  right = TRUE
)

frailty_keep <- with(cohort,
  !is.na(frail_primary) & !is.na(FI_primary) & !is.na(DR1TNUMF) &
    !is.na(RIDAGEYR) & !is.na(RIAGENDR) & !is.na(RIDRETH3) & !is.na(INDFMPIR) &
    !is.na(WTDRD1) & WTDRD1 > 0 & !is.na(SDMVPSU) & !is.na(SDMVSTRA)
)
frailty_data <- cohort[frailty_keep, , drop = FALSE]
frailty_design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTDRD1,
  nest = TRUE,
  data = frailty_data
)

extract_svyby <- function(design, formula, outcome_name, n_values, event_values) {
  summary <- svyby(
    formula,
    ~exposure_q4,
    design,
    svymean,
    na.rm = TRUE,
    vartype = c("se", "ci"),
    keep.var = TRUE
  )
  estimate_name <- setdiff(names(summary), c("exposure_q4", "se", "ci_l", "ci_u"))[1]
  data.frame(
    outcome = outcome_name,
    exposure_q4 = as.character(summary$exposure_q4),
    weighted_estimate = as.numeric(summary[[estimate_name]]),
    standard_error = as.numeric(summary$se),
    ci_low = as.numeric(summary$ci_l),
    ci_high = as.numeric(summary$ci_u),
    n = as.integer(n_values[match(as.character(summary$exposure_q4), quartile_labels)]),
    events = as.integer(event_values[match(as.character(summary$exposure_q4), quartile_labels)]),
    row.names = NULL
  )
}

frailty_counts <- as.integer(table(factor(frailty_data$exposure_q4, levels = quartile_labels)))
names(frailty_counts) <- quartile_labels
frailty_events <- tapply(frailty_data$frail_primary == 1, frailty_data$exposure_q4, sum)
frailty_events <- frailty_events[quartile_labels]
frailty_summary <- extract_svyby(
  frailty_design,
  ~frail_primary,
  "frailty_prevalence",
  as.integer(frailty_counts),
  as.integer(frailty_events)
)
fi_summary <- extract_svyby(
  frailty_design,
  ~FI_primary,
  "fi_mean",
  as.integer(frailty_counts),
  as.integer(frailty_events)
)
exposure_summary <- extract_svyby(
  frailty_design,
  ~DR1TNUMF,
  "exposure_mean",
  as.integer(frailty_counts),
  as.integer(frailty_events)
)

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
mortality_counts <- as.integer(table(factor(mortality_data$exposure_q4, levels = quartile_labels)))
names(mortality_counts) <- quartile_labels
mortality_events <- tapply(mortality_data$mortstat == 1, mortality_data$exposure_q4, sum)
mortality_events <- mortality_events[quartile_labels]
mortality_summary <- extract_svyby(
  mortality_design,
  ~mortstat,
  "linked_mortality_proportion",
  as.integer(mortality_counts),
  as.integer(mortality_events)
)

output <- rbind(frailty_summary, fi_summary, exposure_summary, mortality_summary)
write_csv(output, file.path(result_root, "result_visualization_quartiles.csv"))
writeLines(c(
  "RESULT_VISUALIZATION_SUMMARIES_PASS",
  paste0("frailty_n=", nrow(frailty_data)),
  paste0("mortality_n=", nrow(mortality_data)),
  paste0("mortality_events=", sum(mortality_data$mortstat == 1)),
  paste0("quartile_cutpoints=", paste(cutpoints, collapse = ",")),
  "frailty_weight=WTDRD1",
  "mortality_weight=WTDRD1",
  "mortality_followup=PERMTH_EXM",
  "summary_type=survey-weighted descriptive outcome proportions with 95% confidence intervals"
), file.path(log_root, "result_visualization_summaries.log"))
cat("RESULT_VISUALIZATION_SUMMARIES_PASS\n")
cat("frailty_n=", nrow(frailty_data), "\n", sep = "")
cat("mortality_n=", nrow(mortality_data), "\n", sep = "")
cat("mortality_events=", sum(mortality_data$mortstat == 1), "\n", sep = "")


