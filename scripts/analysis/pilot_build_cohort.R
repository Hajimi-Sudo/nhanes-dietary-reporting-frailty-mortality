#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(haven)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript pilot_build_cohort.R /path/to/project")

project_root <- normalizePath(args[[1]], mustWork = TRUE)
raw_root <- file.path(project_root, "data", "raw", "nhanes_2017_2018")
mapping_path <- file.path(project_root, "data", "metadata", "FRAILTY_INDEX_49_MAPPING.csv")
result_root <- file.path(project_root, "results")
processed_root <- file.path(project_root, "data", "processed")
log_root <- file.path(project_root, "logs")
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)

mapping <- read_csv(mapping_path, show_col_types = FALSE, progress = FALSE)
if (nrow(mapping) != 49 || !all(mapping$item == seq_len(49))) stop("Invalid 49-item mapping")

read_xpt <- function(name) {
  path <- file.path(raw_root, paste0(name, ".xpt"))
  if (!file.exists(path)) stop("Missing XPT file: ", path)
  haven::read_xpt(path)
}

collapse_by_seqn <- function(data, variables) {
  variables <- intersect(unique(variables), names(data))
  groups <- split(seq_len(nrow(data)), as.character(data$SEQN))
  out <- data.frame(SEQN = as.numeric(names(groups)), stringsAsFactors = FALSE)
  for (variable in variables) {
    out[[variable]] <- vapply(groups, function(index) {
      values <- data[[variable]][index]
      values <- values[!is.na(values)]
      if (length(values) == 0) NA_real_ else as.numeric(values[[1]])
    }, numeric(1))
  }
  out
}

needed_files <- unique(c("DEMO_J", "DR1TOT_J", mapping$nhanes_file))
datasets <- setNames(lapply(needed_files, read_xpt), needed_files)

demo <- datasets[["DEMO_J"]]
older <- demo[!is.na(demo$RIDAGEYR) & as.numeric(demo$RIDAGEYR) >= 60, c(
  "SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3", "INDFMPIR",
  "WTMEC2YR", "WTINT2YR", "SDMVPSU", "SDMVSTRA"
)]
older <- collapse_by_seqn(older, setdiff(names(older), "SEQN"))

extra_variables <- list(
  DR1TOT_J = c("DR1TNUMF", "DR1TKCAL", "WTDRD1", "DR1DRSTZ"),
  RXQ_RX_J = c("RXDUSE")
)
wide <- older
for (file_name in setdiff(needed_files, "DEMO_J")) {
  variables <- c(mapping$variable[mapping$nhanes_file == file_name], extra_variables[[file_name]])
  collapsed <- collapse_by_seqn(datasets[[file_name]], variables)
  wide <- merge(wide, collapsed, by = "SEQN", all.x = TRUE, sort = FALSE)
}

numeric_values <- function(x) suppressWarnings(as.numeric(x))
map_values <- function(x, source, target) {
  x <- numeric_values(x)
  result <- rep(NA_real_, length(x))
  for (index in seq_along(source)) result[!is.na(x) & x == source[index]] <- target[index]
  result
}
score_yes_no <- function(x) map_values(x, c(1, 2), c(1, 0))
score_difficulty <- function(x) map_values(x, c(1, 2, 3, 4), c(0, 1 / 3, 2 / 3, 1))
score_depression <- function(x) map_values(x, c(0, 1, 2, 3), c(0, 1 / 3, 2 / 3, 1))

score_item <- function(item, data, same_score) {
  x <- data[[mapping$variable[item]]]
  if (item == 1 || item %in% 25:34 || item == 36 || item == 40) return(score_yes_no(x))
  if (item %in% 2:17) return(score_difficulty(x))
  if (item %in% 18:24) return(score_depression(x))
  if (item == 35) return(map_values(x, c(1, 2, 3), c(1, 0, 0.5)))
  if (item == 37) return(map_values(x, 1:5, c(1, 0.75, 0.5, 0.25, 0)))
  if (item == 38) return(map_values(x, 1:5, c(0, 0, 0, 1, 1)))
  if (item == 39) return(map_values(x, c(1, 2, 3), c(0, same_score, 1)))
  if (item == 41) {
    x <- numeric_values(x)
    return(ifelse(is.na(x) | x < 0 | x > 90, NA_real_, ifelse(x == 0, 0, ifelse(x <= 4, 0.5, 1))))
  }
  if (item == 42) {
    x <- numeric_values(x); use <- numeric_values(data$RXDUSE)
    result <- ifelse(!is.na(use) & use == 2, 0, NA_real_)
    eligible_count <- !is.na(use) & use == 1 & !is.na(x) & x >= 1 & x <= 50
    result[eligible_count] <- ifelse(x[eligible_count] <= 4, 0.5, 1)
    return(result)
  }
  if (item == 43) {
    x <- numeric_values(x)
    return(ifelse(is.na(x) | x <= 0 | x >= 100, NA_real_, ifelse(x < 18.5 | x >= 30, 1, ifelse(x < 25, 0, 0.5))))
  }
  if (item == 44) {
    x <- numeric_values(x)
    return(ifelse(is.na(x), NA_real_, ifelse(x <= 5.7, 0, 1)))
  }
  if (item == 45) {
    x <- numeric_values(x); sex <- numeric_values(data$RIAGENDR)
    valid <- ifelse(sex == 1, x >= 4.7 & x < 6.1, ifelse(sex == 2, x >= 4.2 & x < 5.4, NA))
    return(ifelse(is.na(x) | is.na(valid), NA_real_, ifelse(valid, 0, 1)))
  }
  if (item == 46) {
    x <- numeric_values(x); sex <- numeric_values(data$RIAGENDR)
    valid <- ifelse(sex == 1, x >= 13.5 & x < 18, ifelse(sex == 2, x >= 12 & x < 16, NA))
    return(ifelse(is.na(x) | is.na(valid), NA_real_, ifelse(valid, 0, 1)))
  }
  if (item == 47) {
    x <- numeric_values(x)
    return(ifelse(is.na(x), NA_real_, ifelse(x >= 11.6 & x < 14.6, 0, 1)))
  }
  if (item == 48) {
    x <- numeric_values(x)
    return(ifelse(is.na(x), NA_real_, ifelse(x >= 20 & x < 40, 0, 1)))
  }
  if (item == 49) {
    x <- numeric_values(x)
    return(ifelse(is.na(x), NA_real_, ifelse(x >= 40 & x < 80, 0, 1)))
  }
  stop("No scoring rule for item ", item)
}

build_fi <- function(data, same_score) {
  scores <- sapply(seq_len(nrow(mapping)), score_item, data = data, same_score = same_score)
  if (nrow(scores) != nrow(data)) scores <- t(scores)
  valid_n <- rowSums(!is.na(scores))
  fi <- ifelse(valid_n >= 40, rowMeans(scores, na.rm = TRUE), NA_real_)
  list(scores = scores, valid_n = valid_n, fi = fi)
}

primary <- build_fi(wide, 0.5)
same_zero <- build_fi(wide, 0)
same_one <- build_fi(wide, 1)
wide$frailty_valid_n <- primary$valid_n
wide$FI_primary <- primary$fi
wide$FI_same0 <- same_zero$fi
wide$FI_same1 <- same_one$fi
wide$frail_primary <- ifelse(is.na(wide$FI_primary), NA, as.integer(wide$FI_primary >= 0.21))
wide$DR1TNUMF <- numeric_values(wide$DR1TNUMF)

mortality_path <- file.path(raw_root, "NHANES_2017_2018_MORT_2019_PUBLIC.dat")
mortality_lines <- readLines(mortality_path, warn = FALSE)
mortality <- read_fwf(
  file = mortality_path,
  col_types = "iiiiiiii",
  col_positions = fwf_cols(
    seqn = c(1, 6),
    eligstat = c(15, 15),
    mortstat = c(16, 16),
    ucod_leading = c(17, 19),
    diabetes = c(20, 20),
    hyperten = c(21, 21),
    permth_int = c(43, 45),
    permth_exm = c(46, 48)
  ),
  na = c("", ".")
)
mortality <- mortality[!duplicated(mortality$seqn), ]
mortality_seqn <- mortality$seqn
wide <- merge(wide, mortality, by.x = "SEQN", by.y = "seqn", all.x = TRUE, sort = FALSE)
wide$mortality_linked <- wide$SEQN %in% mortality_seqn
wide$mortality_int_eligible <- wide$eligstat == 1 & !is.na(wide$mortstat) & !is.na(wide$permth_int)
wide$mortality_exm_eligible <- wide$eligstat == 1 & !is.na(wide$mortstat) & !is.na(wide$permth_exm)

fi_overlap <- !is.na(wide$FI_primary) & !is.na(wide$DR1TNUMF)
fi_mortality_int_overlap <- fi_overlap & wide$mortality_int_eligible
fi_mortality_exm_overlap <- fi_overlap & wide$mortality_exm_eligible
summary <- data.frame(
  metric = c(
    "older_adults", "fi_at_least_40_n", "fi_complete_49_n", "fi_exposure_overlap_n",
    "mortality_linked_n", "mortality_linkage_percent", "mortality_deceased_n",
    "mortality_mortstat_valid_n", "mortality_permth_int_valid_n", "mortality_permth_exm_valid_n",
    "mortality_ucod_valid_n", "mortality_int_eligible_n", "mortality_exm_eligible_n",
    "fi_mortality_int_overlap_n", "fi_mortality_exm_overlap_n", "survey_weight_wtdrd1_n",
    "survey_weight_wtmec2yr_n", "frail_primary_n", "fi_primary_min", "fi_primary_max",
    "HUQ020_same_n", "same0_overlap_n", "same1_overlap_n"
  ),
  value = c(
    nrow(wide), sum(primary$valid_n >= 40), sum(primary$valid_n == 49), sum(fi_overlap),
    sum(wide$mortality_linked), round(100 * mean(wide$mortality_linked), 2),
    sum(wide$mortstat == 1, na.rm = TRUE), sum(!is.na(wide$mortstat)),
    sum(!is.na(wide$permth_int)), sum(!is.na(wide$permth_exm)), sum(!is.na(wide$ucod_leading)),
    sum(wide$mortality_int_eligible), sum(wide$mortality_exm_eligible),
    sum(fi_mortality_int_overlap), sum(fi_mortality_exm_overlap),
    sum(fi_overlap & !is.na(wide$WTDRD1) & wide$WTDRD1 > 0),
    sum(fi_overlap & !is.na(wide$WTMEC2YR) & wide$WTMEC2YR > 0),
    sum(wide$frail_primary == 1, na.rm = TRUE), min(wide$FI_primary, na.rm = TRUE),
    max(wide$FI_primary, na.rm = TRUE), sum(numeric_values(wide$HUQ020) == 2, na.rm = TRUE),
    sum(!is.na(wide$FI_same0) & !is.na(wide$DR1TNUMF)),
    sum(!is.na(wide$FI_same1) & !is.na(wide$DR1TNUMF))
  )
)

component_summary <- data.frame(
  item = mapping$item,
  variable = mapping$variable,
  valid_n = colSums(!is.na(primary$scores)),
  missing_n = nrow(wide) - colSums(!is.na(primary$scores))
)

if (sum(primary$valid_n >= 40) < 1800 || sum(fi_overlap) < 1700) {
  stop("PILOT_SIGNAL_GATE_FAIL: cohort size below prespecified threshold")
}
if (any(wide$FI_primary < 0 | wide$FI_primary > 1, na.rm = TRUE)) {
  stop("PILOT_SIGNAL_GATE_FAIL: FI outside [0,1]")
}
if (!all(wide$mortality_linked)) stop("PILOT_SIGNAL_GATE_FAIL: mortality linkage incomplete")
if (sum(wide$eligstat == 1 & is.na(wide$mortstat)) > 0 || sum(wide$eligstat == 1 & is.na(wide$permth_int)) > 0) {
  stop("PILOT_SIGNAL_GATE_FAIL: eligible interview mortality fields are incomplete")
}
if (sum(wide$mortality_exm_eligible) < 2000) {
  stop("PILOT_SIGNAL_GATE_FAIL: too few eligible exam follow-up records")
}

write_csv(summary, file.path(result_root, "pilot_summary.csv"))
write_csv(component_summary, file.path(result_root, "pilot_component_summary.csv"))
saveRDS(wide, file.path(processed_root, "cohort.rds"))
writeLines(c(
  "PILOT_PASS",
  paste0("older_adults=", nrow(wide)),
  paste0("fi_at_least_40_n=", sum(primary$valid_n >= 40)),
  paste0("fi_exposure_overlap_n=", sum(fi_overlap)),
  paste0("mortality_linked_n=", sum(wide$mortality_linked)),
  paste0("output_root=", result_root)
), file.path(log_root, "pilot_build.log"))
cat("PILOT_PASS\n")
cat("older_adults=", nrow(wide), "\n", sep = "")
cat("fi_at_least_40_n=", sum(primary$valid_n >= 40), "\n", sep = "")
cat("fi_exposure_overlap_n=", sum(fi_overlap), "\n", sep = "")
cat("mortality_linked_n=", sum(wide$mortality_linked), "\n", sep = "")


