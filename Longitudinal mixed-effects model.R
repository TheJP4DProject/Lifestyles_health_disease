# This script evaluates longitudinal changes in microbial features using
# baseline-adjusted linear mixed-effects models, Spearman correlation,
# and paired Wilcoxon tests.

library(openxlsx)
library(lme4)
library(lmerTest)

# Load metadata, abundance, and feature maps
METADATA_FILE <- "data/longitudinal_metadata.xlsx"
ABUNDANCE_FILE <- "data/motu_abundance.tsv"
FEATURE_MAP_FILE <- "data/feature_name_map.csv"
REFERENCE_ABUNDANCE_FILE <- "data/reference_motu_abundance.xlsx"
REFERENCE_METADATA_FILE <- "data/reference_metadata.csv"
OUTPUT_DIR <- "results/longitudinal_analysis"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
metadata <- read.xlsx(METADATA_FILE)

abundance <- read.csv(
  ABUNDANCE_FILE,
  sep = "\t",
  check.names = FALSE
)

colnames(abundance) <- c(colnames(abundance)[-1], "remove")
abundance <- as.data.frame(t(abundance[, -ncol(abundance)]))
abundance <- abundance / rowSums(abundance)

# Align samples/features and prepare longitudinal data
feature_map <- read.csv(
  FEATURE_MAP_FILE,
  row.names = 1,
  check.names = FALSE
)

rownames(feature_map) <- feature_map$new

reference_abundance <- read.xlsx(
  REFERENCE_ABUNDANCE_FILE,
  rowNames = TRUE,
  check.names = FALSE
)

reference_metadata <- read.csv(
  REFERENCE_METADATA_FILE,
  check.names = FALSE
)

rownames(reference_metadata) <- reference_metadata$SampleID

reference_abundance <- reference_abundance[
  rownames(reference_metadata),
  ,
  drop = FALSE
]

reference_abundance <- reference_abundance[
  ,
  colMeans(reference_abundance > 0) > 0.025,
  drop = FALSE
]

abundance <- abundance[
  ,
  feature_map[colnames(reference_abundance), "old"],
  drop = FALSE
]

abundance <- abundance[, -1, drop = FALSE]

metadata <- metadata[
  metadata$Sample_T1 %in% rownames(abundance) &
    metadata$Sample_T2 %in% rownames(abundance),
  ,
  drop = FALSE
]

abundance <- abundance[
  rownames(abundance) %in% c(metadata$Sample_T1, metadata$Sample_T2),
  ,
  drop = FALSE
]

raw_abundance <- abundance

abundance[abundance == 0] <- min(abundance[abundance > 0]) / 2
abundance <- scale(log10(abundance))

baseline_variables <- colnames(metadata)[12:39]
followup_variables <- colnames(metadata)[41:68]
change_variables <- colnames(metadata)[70:97]

rownames(feature_map) <- feature_map$old

safe_scale <- function(x) {
  x <- suppressWarnings(as.numeric(x))

  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }

  if (sd(x, na.rm = TRUE) == 0) {
    return(rep(NA_real_, length(x)))
  }

  as.numeric(scale(x))
}

# Run baseline-adjusted linear mixed-effects models
lm_results <- do.call(
  rbind,
  lapply(seq_along(baseline_variables), function(i) {
    do.call(
      rbind,
      lapply(seq_len(ncol(abundance)), function(j) {
        baseline_variable <- baseline_variables[i]
        followup_variable <- followup_variables[i]
        feature_id <- colnames(abundance)[j]

        subject_id <- metadata$ID
        if (all(is.na(subject_id))) {
          subject_id <- metadata$COMMONID
        }

        data_t1 <- data.frame(
          ID = subject_id,
          timepoint = "Baseline",
          species = abundance[metadata$Sample_T1, j],
          lifestyle = metadata[[baseline_variable]],
          baseline_species = abundance[metadata$Sample_T1, j],
          baseline_lifestyle = metadata[[baseline_variable]],
          stringsAsFactors = FALSE
        )

        data_t2 <- data.frame(
          ID = subject_id,
          timepoint = "Followup",
          species = abundance[metadata$Sample_T2, j],
          lifestyle = metadata[[followup_variable]],
          baseline_species = abundance[metadata$Sample_T1, j],
          baseline_lifestyle = metadata[[baseline_variable]],
          stringsAsFactors = FALSE
        )

        model_data <- rbind(data_t1, data_t2)

        model_data$timepoint <- factor(
          model_data$timepoint,
          levels = c("Baseline", "Followup")
        )

        model_data$species <- suppressWarnings(
          as.numeric(model_data$species)
        )
        model_data$lifestyle <- suppressWarnings(
          as.numeric(model_data$lifestyle)
        )
        model_data$baseline_species <- suppressWarnings(
          as.numeric(model_data$baseline_species)
        )
        model_data$baseline_lifestyle <- suppressWarnings(
          as.numeric(model_data$baseline_lifestyle)
        )

        model_data$lifestyle_z <- safe_scale(model_data$lifestyle)
        model_data$baseline_lifestyle_z <- safe_scale(
          model_data$baseline_lifestyle
        )

        model_data <- model_data[
          !is.na(model_data$ID) &
            !is.na(model_data$species) &
            !is.na(model_data$lifestyle_z) &
            !is.na(model_data$baseline_species) &
            !is.na(model_data$baseline_lifestyle_z) &
            !is.na(model_data$timepoint),
          ,
          drop = FALSE
        ]

        fit <- tryCatch(
          lmer(
            species ~
              lifestyle_z +
              baseline_species +
              baseline_lifestyle_z +
              timepoint +
              (1 | ID),
            data = model_data,
            REML = FALSE,
            control = lmerControl(
              optimizer = "bobyqa",
              optCtrl = list(maxfun = 100000),
              check.conv.singular = "ignore"
            )
          ),
          error = function(e) NULL
        )

        if (is.null(fit)) {
          return(
            data.frame(
              metadata_variable = baseline_variable,
              feature = feature_map[feature_id, "new"],
              feature_id = feature_id,
              coefficient = NA_real_,
              p_value = NA_real_
            )
          )
        }

        coefficient_table <- as.data.frame(summary(fit)$coefficients)

        if (!"lifestyle_z" %in% rownames(coefficient_table)) {
          return(
            data.frame(
              metadata_variable = baseline_variable,
              feature = feature_map[feature_id, "new"],
              feature_id = feature_id,
              coefficient = NA_real_,
              p_value = NA_real_
            )
          )
        }

        data.frame(
          metadata_variable = baseline_variable,
          feature = feature_map[feature_id, "new"],
          feature_id = feature_id,
          coefficient = coefficient_table[
            "lifestyle_z",
            "Estimate"
          ],
          p_value = coefficient_table[
            "lifestyle_z",
            "Pr(>|t|)"
          ]
        )
      })
    )
  })
)

write.csv(
  lm_results,
  file.path(OUTPUT_DIR, "linear_regression_results.csv"),
  row.names = FALSE
)

spearman_results <- do.call(
  rbind,
  lapply(seq_along(change_variables), function(i) {
    do.call(
      rbind,
      lapply(seq_len(ncol(abundance)), function(j) {
        delta_abundance <- abundance[metadata$Sample_T2, j] -
          abundance[metadata$Sample_T1, j]

        test <- cor.test(
          delta_abundance,
          metadata[[change_variables[i]]],
          method = "spearman"
        )

        data.frame(
          metadata_variable = baseline_variables[i],
          feature = feature_map[colnames(abundance)[j], "new"],
          feature_id = colnames(abundance)[j],
          coefficient = unname(test$estimate),
          p_value = test$p.value
        )
      })
    )
  })
)

write.csv(
  spearman_results,
  file.path(OUTPUT_DIR, "spearman_results.csv"),
  row.names = FALSE
)

# Run paired Wilcoxon tests and save all results
wilcoxon_results <- do.call(
  rbind,
  lapply(seq_along(change_variables), function(i) {
    do.call(
      rbind,
      lapply(seq_len(ncol(raw_abundance)), function(j) {
        keep <- metadata[[change_variables[i]]] > 0

        abundance_t1 <- raw_abundance[metadata$Sample_T1[keep], j]
        abundance_t2 <- raw_abundance[metadata$Sample_T2[keep], j]

        test <- wilcox.test(
          abundance_t2,
          abundance_t1,
          paired = TRUE,
          exact = TRUE
        )

        data.frame(
          metadata_variable = baseline_variables[i],
          feature = feature_map[colnames(raw_abundance)[j], "new"],
          feature_id = colnames(raw_abundance)[j],
          log2_fold_change = log2(
            mean(abundance_t2) / mean(abundance_t1)
          ),
          p_value = test$p.value
        )
      })
    )
  })
)

write.csv(
  wilcoxon_results,
  file.path(OUTPUT_DIR, "paired_wilcoxon_results.csv"),
  row.names = FALSE
)
