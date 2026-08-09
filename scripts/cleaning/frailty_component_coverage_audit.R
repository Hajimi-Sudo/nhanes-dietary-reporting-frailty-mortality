#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(haven)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript frailty_component_coverage_audit.R /path/to/project")

project_root <- normalizePath(args[[1]], mustWork = TRUE)
raw_root <- file.path(project_root, "data", "raw", "nhanes_2017_2018")
metadata_root <- file.path(project_root, "data", "metadata")
log_root <- file.path(project_root, "logs")
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)

mapping <- read_csv(
  file.path(metadata_root, "FRAILTY_INDEX_49_MAPPING.csv"),
  show_col_types = FALSE,
  progress = FALSE
)
if (nrow(mapping) != 49 || !all(mapping$item == seq_len(49))) {
  stop("The frailty mapping must contain items 1 through 49 in order")
}

needed_files <- unique(c("DEMO_J", mapping$nhanes_file))
datasets <- setNames(lapply(needed_files, function(name) {
  path <- file.path(raw_root, paste0(name, ".xpt"))
  if (!file.exists(path)) stop("Missing XPT file: ", path)
  haven::read_xpt(path)
}), needed_files)

demo <- datasets[["DEMO_J"]]
if (is.null(demo)) stop("DEMO_J is required")
key <- demo[!is.na(demo$RIDAGEYR) & demo$RIDAGEYR >= 60, c("SEQN", "RIAGENDR")]
wide <- key
for (name in needed_files) {
  vars <- unique(c(mapping$variable[mapping$nhanes_file == name], if (name == "RXQ_RX_J") "RXDUSE" else character(0)))
  data <- datasets[[name]][, intersect(c("SEQN", vars), names(datasets[[name]])), drop = FALSE]
  data <- data[!duplicated(data$SEQN), , drop = FALSE]
  wide <- merge(wide, data, by = "SEQN", all.x = TRUE, sort = FALSE)
}

numeric_values <- function(x) suppressWarnings(as.numeric(x))
valid_yes_no <- function(x) !is.na(x) & x %in% c(1, 2)
valid_difficulty <- function(x) !is.na(x) & x %in% c(1, 2, 3, 4)
valid_depression <- function(x) !is.na(x) & x %in% c(0, 1, 2, 3)

valid_item <- function(item, data) {
  x <- numeric_values(data[[mapping$variable[item]]])
  if (item == 1 || item %in% 2:17) return(if (item == 1) valid_yes_no(x) else valid_difficulty(x))
  if (item %in% 18:24) return(valid_depression(x))
  if (item %in% c(25:34, 36, 40)) return(valid_yes_no(x))
  if (item == 35) return(!is.na(x) & x %in% c(1, 2, 3))
  if (item == 37) return(!is.na(x) & x %in% 1:5)
  if (item == 38) return(!is.na(x) & x %in% 1:5)
  if (item == 39) return(!is.na(x) & x %in% 1:3)
  if (item == 41) return(!is.na(x) & x >= 0 & x <= 90)
  if (item == 42) {
    use <- numeric_values(data$RXDUSE)
    return((!is.na(use) & use == 2) | (!is.na(use) & use == 1 & !is.na(x) & x >= 1 & x <= 50))
  }
  if (item == 43) return(!is.na(x) & is.finite(x) & x > 0 & x < 100)
  if (item %in% 44:49) return(!is.na(x) & is.finite(x))
  stop("No validity rule for item ", item)
}

valid_matrix <- sapply(seq_len(nrow(mapping)), valid_item, data = wide)
if (nrow(valid_matrix) != nrow(wide)) valid_matrix <- t(valid_matrix)
colnames(valid_matrix) <- paste0("item_", mapping$item)
wide$valid_component_n <- rowSums(valid_matrix)

same_value_n <- sum(numeric_values(wide$HUQ020) == 2, na.rm = TRUE)
coverage <- data.frame(
  item = mapping$item,
  domain = mapping$domain,
  component = mapping$component,
  nhanes_file = mapping$nhanes_file,
  variable = mapping$variable,
  valid_n = colSums(valid_matrix),
  missing_or_invalid_n = nrow(wide) - colSums(valid_matrix)
)
summary <- data.frame(
  metric = c(
    "older_adults",
    "complete_49_items_n",
    "at_least_40_items_n",
    "at_least_40_items_percent",
    "HUQ020_same_response_n",
    "HUQ020_same_response_note"
  ),
  value = c(
    nrow(wide),
    sum(wide$valid_component_n == 49),
    sum(wide$valid_component_n >= 40),
    round(100 * mean(wide$valid_component_n >= 40), 2),
    same_value_n,
    "HUQ020 code 2 is scored as 0.5 in the primary FI; 0 and 1 are sensitivity analyses"
  )
)

write_csv(coverage, file.path(log_root, "frailty_component_coverage_audit.csv"))
write_csv(summary, file.path(log_root, "frailty_coverage_summary.csv"))
cat("FRAILTY_COVERAGE_AUDIT_OK\n")
cat("older_adults=", nrow(wide), "\n", sep = "")
cat("at_least_40_items_n=", sum(wide$valid_component_n >= 40), "\n", sep = "")
cat("HUQ020_same_response_n=", same_value_n, "\n", sep = "")


