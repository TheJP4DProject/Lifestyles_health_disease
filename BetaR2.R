# This script tests associations of mOTU- and KO-based Bray-Curtis
# beta-diversity with lifestyle and disease variables using PERMANOVA.
#
# Each model is adjusted for predefined covariates and log10-transformed
# sequencing depth.

library(vegan)
library(openxlsx)
library(doParallel)
library(foreach)

# =============================================================================
# Parallel settings
# =============================================================================

N_CORES <- 10

registerDoParallel(
  cores = N_CORES
)

# =============================================================================
# File paths
# =============================================================================

METADATA_FILE <- "data/metadata.xlsx"
DEPTH_FILE <- "data/sequencing_depth.tsv"

MOTU_DISTANCE_FILE <- "data/bray_curtis_motu.rds"
KO_DISTANCE_FILE <- "data/bray_curtis_ko.rds"

OUTPUT_DIR <- "results/beta_diversity"

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# =============================================================================
# Metadata
# =============================================================================

metadata <- read.xlsx(
  METADATA_FILE,
  check.names = FALSE
)

metadata <- metadata[
  ,
  -c(2, 3),
  drop = FALSE
]

rownames(metadata) <- metadata[[1]]

# Predefined adjustment variables
adjustment_cols <- colnames(metadata)[
  2:11
]

# Lifestyle and disease variables
test_cols <- colnames(metadata)[
  12:72
]

# =============================================================================
# Sequencing depth
# =============================================================================

depth <- read.table(
  DEPTH_FILE,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Column 1 contains sample IDs and column 6 contains sequencing depth
depth <- depth[
  ,
  c(1, 6),
  drop = FALSE
]

colnames(depth) <- c(
  "SampleID",
  "Depth"
)

depth$Depth <- as.numeric(
  depth$Depth
)

# Sequencing depth was log10-transformed before inclusion in the model
depth$Depth_log10 <- ifelse(
  depth$Depth > 0,
  log10(depth$Depth),
  NA_real_
)

metadata$Depth_log10 <- depth$Depth_log10[
  match(
    rownames(metadata),
    depth$SampleID
  )
]

# =============================================================================
# Bray-Curtis distance matrices
# =============================================================================

bc_matrices <- list(
  mOTU = as.matrix(
    readRDS(MOTU_DISTANCE_FILE)
  ),
  KO = as.matrix(
    readRDS(KO_DISTANCE_FILE)
  )
)

# =============================================================================
# Analysis tasks
# =============================================================================

tasks <- expand.grid(
  DataType = names(bc_matrices),
  Variable = test_cols,
  stringsAsFactors = FALSE
)

# =============================================================================
# PERMANOVA
# =============================================================================

results <- foreach(
  i = seq_len(nrow(tasks)),
  .combine = rbind,
  .packages = "vegan"
) %dopar% {

  data_type <- tasks$DataType[i]
  variable <- tasks$Variable[i]

  distance_matrix <- bc_matrices[[data_type]]

  model_variables <- c(
    variable,
    adjustment_cols,
    "Depth_log10"
  )

  common_samples <- Reduce(
    intersect,
    list(
      rownames(metadata),
      rownames(distance_matrix),
      colnames(distance_matrix)
    )
  )

  analysis_metadata <- metadata[
    common_samples,
    model_variables,
    drop = FALSE
  ]

  # Samples with missing values in any model variable were excluded
  analysis_metadata <- analysis_metadata[
    complete.cases(analysis_metadata),
    ,
    drop = FALSE
  ]

  analysis_distance <- distance_matrix[
    rownames(analysis_metadata),
    rownames(analysis_metadata),
    drop = FALSE
  ]

  # Non-numeric variables were treated as categorical variables
  analysis_metadata[] <- lapply(
    analysis_metadata,
    function(x) {
      if (is.numeric(x)) {
        x
      } else {
        factor(x)
      }
    }
  )

  # Sequencing depth was retained as a continuous variable
  analysis_metadata$Depth_log10 <- as.numeric(
    analysis_metadata$Depth_log10
  )

  formula_text <- paste(
    "as.dist(analysis_distance) ~",
    paste(
      sprintf(
        "`%s`",
        model_variables
      ),
      collapse = " + "
    )
  )

  set.seed(123)

  fit <- adonis2(
    formula = as.formula(formula_text),
    data = analysis_metadata,
    permutations = 999,
    by = "term",
    parallel = 1
  )

  fit <- as.data.frame(
    fit,
    check.names = FALSE
  )

  fit$Term <- rownames(fit)

  # Retain the PERMANOVA result for the target metadata variable
  target_result <- fit[
    gsub(
      "`",
      "",
      fit$Term,
      fixed = TRUE
    ) == variable,
    ,
    drop = FALSE
  ]

  data.frame(
    DataType = data_type,
    Variable = variable,
    Adjustment = "Predefined covariates + log10 sequencing depth",
    N = nrow(analysis_metadata),
    Df = target_result$Df,
    SumOfSqs = target_result$SumOfSqs,
    R2 = target_result$R2,
    F = target_result$F,
    P = target_result[["Pr(>F)"]],
    check.names = FALSE
  )
}

# =============================================================================
# Multiple-testing correction
# =============================================================================

results$FDR <- ave(
  results$P,
  results$DataType,
  FUN = function(p) {
    p.adjust(
      p,
      method = "BH"
    )
  }
)

# =============================================================================
# Output
# =============================================================================

write.csv(
  results,
  file.path(
    OUTPUT_DIR,
    "permanova_depth_adjusted_results.csv"
  ),
  row.names = FALSE
)

stopImplicitCluster()
