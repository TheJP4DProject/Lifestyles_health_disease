# This script evaluates associations of Shannon diversity, observed richness,
# and microbial load with lifestyle and disease variables using partial
# Spearman correlation adjusted for eight covariates, with and without depth.

library(vegan)
library(ppcor)

RDS_DIR <- "data/abundance"
METADATA_FILE <- "data/metadata.csv"
DEPTH_FILE <- "data/sequencing_depth.tsv"
MICROBIAL_LOAD_FILE <- "data/microbial_load.tsv"
OUTPUT_DIR <- "results/alpha_diversity"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

partial_spearman <- function(x, y, covariates) {
  dat <- na.omit(data.frame(x, y, covariates))
  fit <- pcor.test(
    dat$x,
    dat$y,
    dat[, -(1:2), drop = FALSE],
    method = "spearman"
  )

  c(
    rho = unname(fit$estimate),
    p = fit$p.value,
    n = nrow(dat)
  )
}

run_alpha <- function(data, output_file, microbial_load = NULL) {
  samples <- intersect(rownames(metadata), rownames(data))

  data <- data[samples, , drop = FALSE]
  meta <- metadata[samples, , drop = FALSE]

  prevalence <- colMeans(data > 0)
  data <- data[, prevalence > 0.025, drop = FALSE]

  alpha <- data.frame(
    Shannon = diversity(data, index = "shannon"),
    Richness = rowSums(data > 0),
    row.names = rownames(data)
  )

  if (!is.null(microbial_load)) {
    alpha$MicrobialLoad <- microbial_load[rownames(alpha)]
  }

  results <- lapply(test_cols, function(var) {
    covariates <- meta[rownames(alpha), fixkey, drop = FALSE]
    covariates_depth <- meta[
      rownames(alpha),
      c(fixkey, "Depth_log10"),
      drop = FALSE
    ]

    shannon <- partial_spearman(
      alpha$Shannon,
      meta[rownames(alpha), var],
      covariates
    )

    richness <- partial_spearman(
      alpha$Richness,
      meta[rownames(alpha), var],
      covariates
    )

    shannon_depth <- partial_spearman(
      alpha$Shannon,
      meta[rownames(alpha), var],
      covariates_depth
    )

    richness_depth <- partial_spearman(
      alpha$Richness,
      meta[rownames(alpha), var],
      covariates_depth
    )

    result <- data.frame(
      Variable = var,
      Spearman_Shannon = shannon["rho"],
      P_Shannon = shannon["p"],
      Spearman_Richness = richness["rho"],
      P_Richness = richness["p"],
      N = shannon["n"],
      Spearman_Shannon_8fix_depth = shannon_depth["rho"],
      P_Shannon_8fix_depth = shannon_depth["p"],
      Spearman_Richness_8fix_depth = richness_depth["rho"],
      P_Richness_8fix_depth = richness_depth["p"],
      N_8fix_depth = shannon_depth["n"]
    )

    if (!is.null(microbial_load)) {
      load <- partial_spearman(
        alpha$MicrobialLoad,
        meta[rownames(alpha), var],
        covariates
      )

      load_depth <- partial_spearman(
        alpha$MicrobialLoad,
        meta[rownames(alpha), var],
        covariates_depth
      )

      result$Spearman_ML <- load["rho"]
      result$P_ML <- load["p"]
      result$Spearman_ML_8fix_depth <- load_depth["rho"]
      result$P_ML_8fix_depth <- load_depth["p"]

      result <- result[, c(
        "Variable",
        "Spearman_Shannon", "P_Shannon",
        "Spearman_Richness", "P_Richness",
        "Spearman_ML", "P_ML", "N",
        "Spearman_Shannon_8fix_depth", "P_Shannon_8fix_depth",
        "Spearman_Richness_8fix_depth", "P_Richness_8fix_depth",
        "Spearman_ML_8fix_depth", "P_ML_8fix_depth",
        "N_8fix_depth"
      )]
    }

    result
  })

  write.csv(
    do.call(rbind, results),
    file.path(OUTPUT_DIR, output_file),
    row.names = FALSE
  )
}

metadata <- read.csv(METADATA_FILE, check.names = FALSE)
metadata <- metadata[, -c(3, 4)]
metadata <- metadata[!is.na(metadata$METAF), ]
rownames(metadata) <- metadata$METAF

fixkey <- colnames(metadata)[c(3:5, 51:55)]
test_cols <- colnames(metadata)[6:50]

depth <- read.table(
  DEPTH_FILE,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE
)[, c(1, 6)]

colnames(depth) <- c("SampleID", "Depth")
depth$Depth_log10 <- log10(as.numeric(depth$Depth))

metadata$Depth_log10 <- depth$Depth_log10[
  match(rownames(metadata), depth$SampleID)
]

mo <- readRDS(file.path(RDS_DIR, "mo.rds"))
ko <- readRDS(file.path(RDS_DIR, "ko.rds"))
motu <- readRDS(file.path(RDS_DIR, "motu.rds"))

motu <- motu[
  ,
  colnames(motu) != "Unassigned species",
  drop = FALSE
]

load_data <- read.table(
  MICROBIAL_LOAD_FILE,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

microbial_load <- as.numeric(load_data[, 1])
names(microbial_load) <- rownames(load_data)

run_alpha(mo, "mo_alpha_associations.csv")
run_alpha(ko, "ko_alpha_associations.csv")
run_alpha(motu, "motu_alpha_associations.csv", microbial_load)
