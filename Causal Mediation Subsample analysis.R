# This script performs a case-sampling sensitivity analysis for mediation.
# For each disease, up to 1,000 unique case sets are generated using stratified
# sampling. The 10 sets with the smallest imbalance across all eight adjustment
# variables are retained for mediation analysis.

library(mediation)

set.seed(123)

args <- commandArgs(trailingOnly = TRUE)
exposure_index <- as.integer(args[1])
outcome_index <- as.integer(args[2])

METADATA_FILE <- "data/metadata.csv"
ABUNDANCE_FILE <- "data/motu_abundance.csv"
MAASLIN_FILE <- "data/maaslin_results.csv"
OUTPUT_DIR <- "results/sampling_sensitivity"

TARGET_CASES <- 134
MAX_CANDIDATE_SETS <- 1000
N_SELECTED_SETS <- 10
N_SIMS <- 1000
SAMPLING_SEED <- 123
MEDIATION_SEED <- 123

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(
  METADATA_FILE,
  check.names = FALSE
)

adjustment_positions <- c(3:5, 51:55)
adjustment_original_names <- colnames(metadata)[adjustment_positions]
colnames(metadata)[adjustment_positions] <- paste0("ADJ", 1:8)
colnames(metadata)[2] <- "SampleID"

ppi_position <- grep(
  "proton|pump|ppi",
  adjustment_original_names,
  ignore.case = TRUE
)[1]

antibiotic_position <- grep(
  "antibiotic",
  adjustment_original_names,
  ignore.case = TRUE
)[1]

PPI_VARIABLE <- paste0("ADJ", ppi_position)
ANTIBIOTIC_VARIABLE <- paste0("ADJ", antibiotic_position)

exposure_names <- colnames(metadata)[6:33]
outcome_names <- colnames(metadata)[35:50]

exposure <- exposure_names[exposure_index]
outcome <- outcome_names[outcome_index]

metadata <- metadata[!is.na(metadata$SampleID), ]
rownames(metadata) <- as.character(metadata$SampleID)

abundance <- read.csv(
  ABUNDANCE_FILE,
  row.names = 1,
  check.names = FALSE
)

abundance <- abundance[rownames(metadata), , drop = FALSE]
abundance[abundance == 0] <- min(abundance[abundance > 0]) / 2
abundance <- log10(abundance)

maaslin <- read.csv(
  MAASLIN_FILE,
  check.names = FALSE
)

exposure_features <- maaslin[
  maaslin$metadata == exposure &
    maaslin$qval < 0.1,
  ,
  drop = FALSE
]

outcome_features <- maaslin[
  maaslin$metadata == outcome &
    maaslin$qval < 0.1,
  ,
  drop = FALSE
]

selected_features <- intersect(
  exposure_features$feature,
  outcome_features$feature
)

selected_features <- colnames(abundance)[
  colnames(abundance) %in% selected_features
]

abundance <- abundance[, selected_features, drop = FALSE]

analysis_variables <- c(
  exposure,
  outcome,
  paste0("ADJ", 1:8)
)

metadata[analysis_variables] <- lapply(
  metadata[analysis_variables],
  as.numeric
)

complete_samples <- complete.cases(
  metadata[, analysis_variables, drop = FALSE]
)

metadata <- metadata[complete_samples, , drop = FALSE]
abundance <- abundance[rownames(metadata), , drop = FALSE]

metadata <- metadata[
  metadata[[outcome]] %in% c(0, 1),
  ,
  drop = FALSE
]

abundance <- abundance[rownames(metadata), , drop = FALSE]

case_data <- metadata[metadata[[outcome]] == 1, , drop = FALSE]
control_ids <- rownames(metadata)[metadata[[outcome]] == 0]

if (nrow(case_data) <= TARGET_CASES) {
  quit(save = "no", status = 0)
}

case_data$SexGroup <- as.integer(case_data$ADJ2)
case_data$AgeGroup <- as.integer(case_data$ADJ1 >= 65)
case_data$BMIGroup <- as.integer(case_data$ADJ3 >= 25)
case_data$PPIGroup <- as.integer(
  case_data[[PPI_VARIABLE]] > 0
)

case_data$AntibioticGroup <- as.integer(
  case_data[[ANTIBIOTIC_VARIABLE]] > 0
)

case_data$Stratum <- interaction(
  case_data$SexGroup,
  case_data$AgeGroup,
  case_data$BMIGroup,
  case_data$PPIGroup,
  case_data$AntibioticGroup,
  drop = TRUE
)

allocate_quota <- function(counts, target_n) {
  raw_quota <- target_n * counts / sum(counts)
  quota <- floor(raw_quota)
  remaining <- target_n - sum(quota)

  if (remaining > 0) {
    add_order <- order(
      raw_quota - quota,
      counts,
      decreasing = TRUE
    )

    quota[add_order[seq_len(remaining)]] <-
      quota[add_order[seq_len(remaining)]] + 1
  }

  quota
}

stratum_counts <- table(case_data$Stratum)
stratum_quota <- allocate_quota(
  stratum_counts,
  TARGET_CASES
)

draw_case_set <- function() {
  unlist(
    lapply(names(stratum_quota), function(stratum) {
      sample(
        rownames(case_data)[case_data$Stratum == stratum],
        size = stratum_quota[[stratum]],
        replace = FALSE
      )
    }),
    use.names = FALSE
  )
}

score_case_set <- function(selected_ids) {
  full_cases <- case_data[, paste0("ADJ", 1:8), drop = FALSE]
  selected_cases <- case_data[
    selected_ids,
    paste0("ADJ", 1:8),
    drop = FALSE
  ]

  standardized_differences <- vapply(
    paste0("ADJ", 1:8),
    function(variable) {
      denominator <- sd(full_cases[[variable]])

      if (is.na(denominator) || denominator == 0) {
        return(0)
      }

      abs(
        mean(selected_cases[[variable]]) -
          mean(full_cases[[variable]])
      ) / denominator
    },
    numeric(1)
  )

  c(
    mean_absolute_difference = mean(standardized_differences),
    maximum_absolute_difference = max(standardized_differences),
    standardized_differences
  )
}

set.seed(SAMPLING_SEED)

candidate_sets <- list()
candidate_signatures <- character(0)

attempts <- 0
maximum_attempts <- MAX_CANDIDATE_SETS * 100

while (
  length(candidate_sets) < MAX_CANDIDATE_SETS &&
    attempts < maximum_attempts
) {
  attempts <- attempts + 1

  selected_ids <- draw_case_set()
  signature <- paste(sort(selected_ids), collapse = "|")

  if (!(signature %in% candidate_signatures)) {
    candidate_sets[[length(candidate_sets) + 1]] <- selected_ids
    candidate_signatures <- c(candidate_signatures, signature)
  }
}

candidate_scores <- do.call(
  rbind,
  lapply(candidate_sets, score_case_set)
)

candidate_scores <- as.data.frame(candidate_scores)
candidate_scores$CandidateSet <- seq_len(nrow(candidate_scores))

candidate_scores <- candidate_scores[
  order(
    candidate_scores$mean_absolute_difference,
    candidate_scores$maximum_absolute_difference
  ),
  ,
  drop = FALSE
]

selected_set_indices <- head(
  candidate_scores$CandidateSet,
  N_SELECTED_SETS
)

sampling_sets <- candidate_sets[selected_set_indices]

selected_sampling_scores <- candidate_scores[
  seq_len(min(N_SELECTED_SETS, nrow(candidate_scores))),
  ,
  drop = FALSE
]

selected_sampling_scores$SelectedSet <- seq_len(
  nrow(selected_sampling_scores)
)

sampling_plan <- do.call(
  rbind,
  lapply(seq_along(sampling_sets), function(set_index) {
    data.frame(
      Disease = outcome,
      Set = set_index,
      SampleID = sampling_sets[[set_index]],
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  selected_sampling_scores,
  file.path(
    OUTPUT_DIR,
    paste0(
      "sampling_balance_",
      exposure_index,
      "_",
      outcome_index,
      ".csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  sampling_plan,
  file.path(
    OUTPUT_DIR,
    paste0(
      "sampling_plan_",
      exposure_index,
      "_",
      outcome_index,
      ".csv"
    )
  ),
  row.names = FALSE
)

extract_mediation_result <- function(result) {
  if (inherits(result, "try-error")) {
    return(
      data.frame(
        ACME = NA,
        ACME_p = NA,
        ACME_ci_low = NA,
        ACME_ci_high = NA,
        ADE = NA,
        ADE_p = NA,
        ADE_ci_low = NA,
        ADE_ci_high = NA,
        Proportion_mediated = NA,
        Proportion_mediated_p = NA,
        Total_effect = NA,
        Total_effect_p = NA
      )
    )
  }

  data.frame(
    ACME = result$d.avg,
    ACME_p = result$d.avg.p,
    ACME_ci_low = result$d.avg.ci[1],
    ACME_ci_high = result$d.avg.ci[2],
    ADE = result$z.avg,
    ADE_p = result$z.avg.p,
    ADE_ci_low = result$z.avg.ci[1],
    ADE_ci_high = result$z.avg.ci[2],
    Proportion_mediated = result$n.avg,
    Proportion_mediated_p = result$n.avg.p,
    Total_effect = result$tau.coef,
    Total_effect_p = result$tau.p
  )
}

get_maaslin_value <- function(data, feature, column) {
  value <- data[data$feature == feature, column]
  if (length(value) == 0) NA_real_ else as.numeric(value[1])
}

all_results <- lapply(seq_along(sampling_sets), function(set_index) {
  selected_ids <- c(
    control_ids,
    sampling_sets[[set_index]]
  )

  selected_ids <- rownames(metadata)[
    rownames(metadata) %in% selected_ids
  ]

  metadata_set <- metadata[selected_ids, , drop = FALSE]
  abundance_set <- abundance[selected_ids, , drop = FALSE]

  set_results <- lapply(colnames(abundance_set), function(feature) {
    analysis_data <- data.frame(
      mediator = abundance_set[[feature]],
      exposure = metadata_set[[exposure]],
      outcome = metadata_set[[outcome]],
      metadata_set[, paste0("ADJ", 1:8), drop = FALSE]
    )

    analysis_data <- na.omit(analysis_data)

    mediator_model <- lm(
      mediator ~
        exposure +
        ADJ1 + ADJ2 + ADJ3 + ADJ4 +
        ADJ5 + ADJ6 + ADJ7 + ADJ8,
      data = analysis_data
    )

    outcome_model <- glm(
      outcome ~
        exposure +
        mediator +
        ADJ1 + ADJ2 + ADJ3 + ADJ4 +
        ADJ5 + ADJ6 + ADJ7 + ADJ8,
      data = analysis_data,
      family = binomial()
    )

    set.seed(MEDIATION_SEED)

    mediation_result <- try(
      mediate(
        mediator_model,
        outcome_model,
        sims = N_SIMS,
        treat = "exposure",
        mediator = "mediator"
      ),
      silent = TRUE
    )

    result <- extract_mediation_result(mediation_result)

    data.frame(
      Set = set_index,
      CandidateSet = selected_set_indices[set_index],
      Mean_absolute_ADJ_difference =
        selected_sampling_scores$mean_absolute_difference[set_index],
      Maximum_absolute_ADJ_difference =
        selected_sampling_scores$maximum_absolute_difference[set_index],
      Treat = exposure,
      Disease = outcome,
      Mediator = feature,
      N = nrow(analysis_data),
      N_case = sum(analysis_data$outcome == 1),
      N_control = sum(analysis_data$outcome == 0),
      result,
      Exposure_coefficient = get_maaslin_value(
        exposure_features,
        feature,
        "coef"
      ),
      Exposure_p = get_maaslin_value(
        exposure_features,
        feature,
        "pval"
      ),
      Exposure_q = get_maaslin_value(
        exposure_features,
        feature,
        "qval"
      ),
      Disease_coefficient = get_maaslin_value(
        outcome_features,
        feature,
        "coef"
      ),
      Disease_p = get_maaslin_value(
        outcome_features,
        feature,
        "pval"
      ),
      Disease_q = get_maaslin_value(
        outcome_features,
        feature,
        "qval"
      )
    )
  })

  set_results <- do.call(rbind, set_results)
  set_results$ACME_q <- p.adjust(
    set_results$ACME_p,
    method = "fdr"
  )

  set_results
})

results <- do.call(rbind, all_results)

write.csv(
  results,
  file.path(
    OUTPUT_DIR,
    paste0(
      "mediation_sampling_sensitivity_",
      exposure_index,
      "_",
      outcome_index,
      ".csv"
    )
  ),
  row.names = FALSE
)
