# This script evaluates associations of species-level alpha-diversity measures
# and microbial load with lifestyle and disease variables using partial
# Spearman correlation.
#
# Two models are fitted:
#   Model 1: adjusted for eight predefined covariates
#   Model 2: adjusted for the same eight covariates plus log10 sequencing depth

library(vegan)
library(ppcor)

# =============================================================================
# File paths
# =============================================================================

ABUNDANCE_FILE <- "data/abundance/motu.rds"
METADATA_FILE <- "data/metadata.csv"
DEPTH_FILE <- "data/sequencing_depth.tsv"
MICROBIAL_LOAD_FILE <- "data/microbial_load.tsv"
OUTPUT_FILE <- "results/alpha_diversity/motu_alpha_associations.csv"

dir.create(
  dirname(OUTPUT_FILE),
  recursive = TRUE,
  showWarnings = FALSE
)

# =============================================================================
# Partial Spearman correlation
# =============================================================================

partial_spearman <- function(x, y, covariates) {

  analysis_data <- na.omit(
    data.frame(
      x = x,
      y = y,
      covariates,
      check.names = FALSE
    )
  )

  result <- pcor.test(
    x = analysis_data$x,
    y = analysis_data$y,
    z = analysis_data[, -(1:2), drop = FALSE],
    method = "spearman"
  )

  c(
    rho = unname(result$estimate),
    p = result$p.value,
    n = nrow(analysis_data)
  )
}

# =============================================================================
# Metadata
# =============================================================================

metadata <- read.csv(
  METADATA_FILE,
  check.names = FALSE
)

metadata <- metadata[
  ,
  -c(3, 4),
  drop = FALSE
]

metadata <- metadata[
  !is.na(metadata$SampleID),
  ,
  drop = FALSE
]

rownames(metadata) <- metadata$SampleID

# Eight predefined adjustment variables
adjustment_cols <- colnames(metadata)[
  c(3:5, 51:55)
]

# Lifestyle and disease variables
test_cols <- colnames(metadata)[
  6:50
]

# =============================================================================
# Sequencing depth
# =============================================================================

depth <- read.table(
  DEPTH_FILE,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE
)

# The first column contains sample IDs and the sixth column contains depth
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

depth$Depth_log10 <- log10(
  depth$Depth
)

metadata$Depth_log10 <- depth$Depth_log10[
  match(
    rownames(metadata),
    depth$SampleID
  )
]

# =============================================================================
# Species-level mOTU abundance
# =============================================================================

motu <- readRDS(
  ABUNDANCE_FILE
)

# Remove unassigned species
motu <- motu[
  ,
  colnames(motu) != "Unassigned species",
  drop = FALSE
]

# Match abundance data with metadata
samples <- intersect(
  rownames(metadata),
  rownames(motu)
)

motu <- motu[
  samples,
  ,
  drop = FALSE
]

metadata <- metadata[
  samples,
  ,
  drop = FALSE
]

# Retain species detected in more than 2.5% of samples
prevalence <- colMeans(
  motu > 0
)

motu <- motu[
  ,
  prevalence > 0.025,
  drop = FALSE
]

# =============================================================================
# Alpha-diversity measures
# =============================================================================

alpha <- data.frame(
  Shannon = diversity(
    motu,
    index = "shannon"
  ),
  Richness = rowSums(
    motu > 0
  ),
  row.names = rownames(motu)
)

# =============================================================================
# Microbial load
# =============================================================================

load_data <- read.table(
  MICROBIAL_LOAD_FILE,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

microbial_load <- as.numeric(
  load_data[, 1]
)

names(microbial_load) <- rownames(
  load_data
)

alpha$MicrobialLoad <- microbial_load[
  rownames(alpha)
]

# =============================================================================
# Association analysis
# =============================================================================

results <- lapply(
  test_cols,
  function(variable) {

    outcome <- metadata[
      rownames(alpha),
      variable
    ]

    # Model 1: eight adjustment variables
    covariates_8 <- metadata[
      rownames(alpha),
      adjustment_cols,
      drop = FALSE
    ]

    # Model 2: eight adjustment variables plus log10 sequencing depth
    covariates_8_depth <- metadata[
      rownames(alpha),
      c(adjustment_cols, "Depth_log10"),
      drop = FALSE
    ]

    shannon_8 <- partial_spearman(
      x = alpha$Shannon,
      y = outcome,
      covariates = covariates_8
    )

    richness_8 <- partial_spearman(
      x = alpha$Richness,
      y = outcome,
      covariates = covariates_8
    )

    load_8 <- partial_spearman(
      x = alpha$MicrobialLoad,
      y = outcome,
      covariates = covariates_8
    )

    shannon_8_depth <- partial_spearman(
      x = alpha$Shannon,
      y = outcome,
      covariates = covariates_8_depth
    )

    richness_8_depth <- partial_spearman(
      x = alpha$Richness,
      y = outcome,
      covariates = covariates_8_depth
    )

    load_8_depth <- partial_spearman(
      x = alpha$MicrobialLoad,
      y = outcome,
      covariates = covariates_8_depth
    )

    data.frame(
      Variable = variable,

      Spearman_Shannon_8cov =
        shannon_8["rho"],

      P_Shannon_8cov =
        shannon_8["p"],

      N_Shannon_8cov =
        shannon_8["n"],

      Spearman_Richness_8cov =
        richness_8["rho"],

      P_Richness_8cov =
        richness_8["p"],

      N_Richness_8cov =
        richness_8["n"],

      Spearman_MicrobialLoad_8cov =
        load_8["rho"],

      P_MicrobialLoad_8cov =
        load_8["p"],

      N_MicrobialLoad_8cov =
        load_8["n"],

      Spearman_Shannon_8cov_depth =
        shannon_8_depth["rho"],

      P_Shannon_8cov_depth =
        shannon_8_depth["p"],

      N_Shannon_8cov_depth =
        shannon_8_depth["n"],

      Spearman_Richness_8cov_depth =
        richness_8_depth["rho"],

      P_Richness_8cov_depth =
        richness_8_depth["p"],

      N_Richness_8cov_depth =
        richness_8_depth["n"],

      Spearman_MicrobialLoad_8cov_depth =
        load_8_depth["rho"],

      P_MicrobialLoad_8cov_depth =
        load_8_depth["p"],

      N_MicrobialLoad_8cov_depth =
        load_8_depth["n"],

      check.names = FALSE
    )
  }
)

results <- do.call(
  rbind,
  results
)

rownames(results) <- NULL

write.csv(
  results,
  OUTPUT_FILE,
  row.names = FALSE
)
