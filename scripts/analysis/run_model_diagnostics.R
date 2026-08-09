#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survival)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript run_model_diagnostics.R /path/to/project")

project_root <- normalizePath(args[[1]], mustWork = TRUE)
cohort <- readRDS(file.path(project_root, "data", "processed", "cohort.rds"))
result_root <- file.path(project_root, "results")
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)

cohort$diet_variety_5 <- cohort$DR1TNUMF / 5
cohort$age_10 <- cohort$RIDAGEYR / 10
cohort$sex <- factor(cohort$RIAGENDR)
cohort$race_ethnicity <- factor(cohort$RIDRETH3)
candidate <- with(cohort, !is.na(FI_primary) & !is.na(DR1TNUMF) & !is.na(age_10) &
  !is.na(RIAGENDR) & !is.na(RIDRETH3) & !is.na(INDFMPIR) &
  !is.na(permth_exm) & !is.na(mortstat) & eligstat == 1 &
  !is.na(WTDRD1) & WTDRD1 > 0)
zero_followup_excluded_n <- sum(candidate & cohort$permth_exm <= 0, na.rm = TRUE)
keep <- candidate & cohort$permth_exm > 0
data <- cohort[keep, , drop = FALSE]

diagnostic_model <- coxph(
  Surv(permth_exm, mortstat) ~ FI_primary + diet_variety_5 + age_10 + sex + race_ethnicity + INDFMPIR,
  data = data,
  ties = "efron"
)
ph <- cox.zph(diagnostic_model)
ph_table <- as.data.frame(ph$table)
ph_table$term <- rownames(ph_table)
rownames(ph_table) <- NULL
ph_table <- ph_table[, c("term", "chisq", "df", "p")]
names(ph_table) <- c("term", "chisq", "df", "p_value")

cell_counts <- rbind(
  transform(as.data.frame(table(data$sex, data$mortstat)), dimension = "sex"),
  transform(as.data.frame(table(data$race_ethnicity, data$mortstat)), dimension = "race_ethnicity")
)
names(cell_counts)[seq_len(2)] <- c("category", "mortstat")

summary <- data.frame(
  metric = c("n", "events", "zero_followup_excluded_n", "weight_min", "weight_max", "ph_global_p"),
  value = c(
    nrow(data), sum(data$mortstat == 1), zero_followup_excluded_n,
    min(data$WTDRD1), max(data$WTDRD1), ph_table$p_value[ph_table$term == "GLOBAL"]
  )
)
write_csv(summary, file.path(result_root, "model_diagnostics_summary.csv"))
write_csv(ph_table, file.path(result_root, "cox_ph_diagnostics.csv"))
write_csv(cell_counts, file.path(result_root, "mortality_cell_counts.csv"))
writeLines(c(
  "MODEL_DIAGNOSTICS_PASS",
  paste0("n=", nrow(data)),
  paste0("events=", sum(data$mortstat == 1)),
  paste0("zero_followup_excluded_n=", zero_followup_excluded_n),
  paste0("ph_global_p=", ph_table$p_value[ph_table$term == "GLOBAL"])
), file.path(project_root, "logs", "model_diagnostics.log"))
cat("MODEL_DIAGNOSTICS_PASS\n")
cat("n=", nrow(data), "\n", sep = "")
cat("events=", sum(data$mortstat == 1), "\n", sep = "")
cat("zero_followup_excluded_n=", zero_followup_excluded_n, "\n", sep = "")


