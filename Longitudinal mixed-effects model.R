# This script evaluates longitudinal changes in microbial features using
# standardized linear regression, Spearman correlation, and paired Wilcoxon tests.

library(openxlsx)

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
change_variables <- colnames(metadata)[70:97]

rownames(feature_map) <- feature_map$old

lm_results <- do.call(
  rbind,
  lapply(seq_along(change_variables), function(i) {
    do.call(
      rbind,
      lapply(seq_len(ncol(abundance)), function(j) {
        delta_abundance <- abundance[metadata$Sample_T2, j] -
          abundance[metadata$Sample_T1, j]

        delta_metadata <- metadata[[change_variables[i]]]

        fit <- summary(
          lm(
            scale(delta_abundance) ~ scale(delta_metadata)
          )
        )

        data.frame(
          metadata_variable = baseline_variables[i],
          feature = feature_map[colnames(abundance)[j], "new"],
          feature_id = colnames(abundance)[j],
          coefficient = fit$coefficients[2, "Estimate"],
          p_value = fit$coefficients[2, "Pr(>|t|)"]
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
