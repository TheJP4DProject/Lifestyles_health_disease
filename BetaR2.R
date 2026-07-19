# This script calculates Bray-Curtis distance matrices for MO, mOTU, and KO
# profiles, then tests their associations with metadata variables using adonis2.
# Models are adjusted for ten predefined covariates and log10 sequencing depth.

library(vegan)
library(openxlsx)
library(parallelDist)
library(doParallel)
library(foreach)

registerDoParallel(cores = 10)

METADATA_FILE <- "data/metadata.xlsx"
DEPTH_FILE <- "data/sequencing_depth.tsv"

MO_FILE <- "data/mo_abundance.tsv"
MOTU_FILE <- "data/motu_abundance.xlsx"
KO_FILE <- "data/ko_abundance.xlsx"

OUTPUT_DIR <- "results/beta_diversity"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

metadata <- read.xlsx(METADATA_FILE, check.names = FALSE)
metadata <- metadata[, -c(2, 3)]
rownames(metadata) <- metadata[[1]]

fixkey <- colnames(metadata)[2:11]
test_cols <- colnames(metadata)[12:72]

depth <- read.table(
  DEPTH_FILE,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE
)[, c(1, 6)]

colnames(depth) <- c("SampleID", "Depth")

metadata$Depth_log10 <- log10(
  depth$Depth[match(rownames(metadata), depth$SampleID)]
)

calculate_bc <- function(data, prevalence = 0.025) {
  common <- intersect(rownames(data), rownames(metadata))
  data <- data[common, , drop = FALSE]
  data <- data[, colMeans(data > 0) > prevalence, drop = FALSE]

  bc <- as.matrix(
    parDist(
      as.matrix(data),
      method = "bray",
      threads = 10
    )
  )

  rownames(bc) <- colnames(bc) <- rownames(data)
  bc
}

mo <- read.table(
  MO_FILE,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

colnames(mo) <- c(colnames(mo)[-1], "remove")
mo <- t(mo[, -ncol(mo)])

motu <- read.xlsx(
  MOTU_FILE,
  rowNames = TRUE,
  check.names = FALSE
)

motu <- motu[
  ,
  colnames(motu) != "Unassigned species",
  drop = FALSE
]

ko <- read.xlsx(
  KO_FILE,
  rowNames = TRUE,
  check.names = FALSE
)

ko <- t(ko)

bc_matrices <- list(
  MO = calculate_bc(mo),
  motu = calculate_bc(motu),
  KO = calculate_bc(ko)
)

saveRDS(
  bc_matrices$MO,
  file.path(OUTPUT_DIR, "bc_MO.rds")
)

saveRDS(
  bc_matrices$motu,
  file.path(OUTPUT_DIR, "bc_motu.rds")
)

saveRDS(
  bc_matrices$KO,
  file.path(OUTPUT_DIR, "bc_KO.rds")
)

tasks <- expand.grid(
  data_type = names(bc_matrices),
  metadata_var = test_cols,
  stringsAsFactors = FALSE
)

results <- foreach(
  i = seq_len(nrow(tasks)),
  .combine = rbind,
  .packages = "vegan"
) %dopar% {
  type <- tasks$data_type[i]
  var <- tasks$metadata_var[i]
  bc <- bc_matrices[[type]]

  samples <- intersect(rownames(metadata), rownames(bc))

  meta <- metadata[
    samples,
    c(var, fixkey, "Depth_log10"),
    drop = FALSE
  ]

  meta <- meta[complete.cases(meta), , drop = FALSE]
  bc <- bc[rownames(meta), rownames(meta), drop = FALSE]

  meta[] <- lapply(
    meta,
    function(x) if (is.character(x)) factor(x) else x
  )

  formula_text <- paste(
    "as.dist(bc) ~",
    paste(
      sprintf("`%s`", c(var, fixkey, "Depth_log10")),
      collapse = " + "
    )
  )

  set.seed(123)

  fit <- adonis2(
    as.formula(formula_text),
    data = meta,
    permutations = 999,
    by = "term"
  )

  result <- as.data.frame(fit)
  result$term <- rownames(result)
  result$data_type <- type
  result$metadata_var <- var
  result$n_samples <- nrow(meta)

  result <- result[
    gsub("`", "", result$term) == var,
    ,
    drop = FALSE
  ]

  result[, c(
    "data_type",
    "metadata_var",
    "term",
    "n_samples",
    setdiff(
      colnames(result),
      c("data_type", "metadata_var", "term", "n_samples")
    )
  )]
}

results$FDR <- ave(
  results[["Pr(>F)"]],
  results$data_type,
  FUN = function(p) p.adjust(p, method = "BH")
)

write.csv(
  results,
  file.path(OUTPUT_DIR, "adonis2_results.csv"),
  row.names = FALSE
)

stopImplicitCluster()
