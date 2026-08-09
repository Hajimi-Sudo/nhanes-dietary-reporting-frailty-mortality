#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(readr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript variable_audit.R /path/to/project")
}

project_root <- normalizePath(args[[1]], mustWork = TRUE)
raw_root <- file.path(project_root, "data", "raw", "nhanes_2017_2018")
log_root <- file.path(project_root, "logs")
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)

read_xpt <- function(name) {
  path <- file.path(raw_root, paste0(name, ".xpt"))
  if (!file.exists(path)) stop("Missing XPT file: ", path)
  haven::read_xpt(path)
}

datasets <- list(
  DEMO = read_xpt("DEMO_J"),
  DR1TOT = read_xpt("DR1TOT_J"),
  BMX = read_xpt("BMX_J"),
  PFQ = read_xpt("PFQ_J"),
  MCQ = read_xpt("MCQ_J"),
  PAQ = read_xpt("PAQ_J"),
  SMQ = read_xpt("SMQ_J"),
  ALQ = read_xpt("ALQ_J")
)

key <- datasets$DEMO %>%
  filter(RIDAGEYR >= 60) %>%
  select(SEQN, RIDAGEYR)

candidate_vars <- list(
  DEMO = c("SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3", "INDFMPIR", "WTMEC2YR", "WTINT2YR", "SDMVPSU", "SDMVSTRA"),
  DR1TOT = c("SEQN", "DR1TNUMF", "DR1TKCAL", "WTDRD1", "DR1DRSTZ"),
  BMX = c("SEQN", "BMXBMI", "BMXWT", "BMXHT", "BMXWAIST"),
  PFQ = c("SEQN", "PFQ049", "PFQ051", "PFQ054", "PFQ057", "PFQ059", "PFQ061A", "PFQ061B", "PFQ061C", "PFQ063A", "PFQ063B", "PFQ090"),
  MCQ = c("SEQN", "MCQ160A", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F", "MCQ160G", "MCQ160K", "MCQ160M", "MCQ160N", "MCQ160O"),
  PAQ = c("SEQN", "PAQ605", "PAQ610", "PAD615", "PAQ620", "PAQ625", "PAD630", "PAQ635", "PAQ640", "PAD645", "PAQ650", "PAQ655", "PAD660", "PAQ665", "PAQ670", "PAD675", "PAD680"),
  SMQ = c("SEQN", "SMQ020", "SMQ040"),
  ALQ = c("SEQN", "ALQ111", "ALQ121", "ALQ130", "ALQ142", "ALQ151", "ALQ170")
)

safe_unique_nonmissing <- function(x) {
  x <- as.character(x[!is.na(x)])
  length(unique(x))
}

audit <- bind_rows(lapply(names(candidate_vars), function(dataset_name) {
  data <- datasets[[dataset_name]]
  vars <- intersect(candidate_vars[[dataset_name]], names(data))
  missing_vars <- setdiff(candidate_vars[[dataset_name]], names(data))
  data <- data %>% inner_join(key, by = "SEQN")
  present_rows <- lapply(vars, function(var) {
    x <- data[[var]]
    tibble(
      dataset = dataset_name,
      variable = var,
      exists = TRUE,
      older_n = nrow(data),
      nonmissing_n = sum(!is.na(x)),
      missing_n = sum(is.na(x)),
      unique_nonmissing = safe_unique_nonmissing(x),
      label = as.character(attr(datasets[[dataset_name]][[var]], "label") %||% "")
    )
  })
  missing_rows <- lapply(missing_vars, function(var) {
    tibble(
      dataset = dataset_name,
      variable = var,
      exists = FALSE,
      older_n = nrow(data),
      nonmissing_n = NA_integer_,
      missing_n = NA_integer_,
      unique_nonmissing = NA_integer_,
      label = ""
    )
  })
  bind_rows(c(present_rows, missing_rows))
}))

write_csv(audit, file.path(log_root, "variable_audit.csv"))

mortality_path <- file.path(raw_root, "NHANES_2017_2018_MORT_2019_PUBLIC.dat")
mortality_lines <- readLines(mortality_path, warn = FALSE)
mortality_seqn <- suppressWarnings(as.integer(substr(mortality_lines, 1, 6)))
mortality_check <- tibble(
  metric = c("mortality_rows", "mortality_unique_seqn", "older_adults", "older_adults_linked"),
  value = c(
    length(mortality_lines),
    dplyr::n_distinct(mortality_seqn, na.rm = TRUE),
    nrow(key),
    sum(as.integer(key$SEQN) %in% mortality_seqn, na.rm = TRUE)
  )
)
write_csv(mortality_check, file.path(log_root, "mortality_linkage_audit.csv"))

cat("VARIABLE_AUDIT_OK\n")
cat("older_adults=", nrow(key), "\n", sep = "")
cat("audit_file=", file.path(log_root, "variable_audit.csv"), "\n", sep = "")
cat("mortality_file=", file.path(log_root, "mortality_linkage_audit.csv"), "\n", sep = "")


