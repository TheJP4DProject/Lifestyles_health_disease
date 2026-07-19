# This script evaluates a species-only microbiome model using stratified
# 5-fold outer cross-validation and inner cross-validation for LASSO tuning.
# The prevalence filter is defined once on the full dataset and reused across
# all outer folds.

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
  library(caret)
  library(dplyr)
  library(openxlsx)
  library(future)
  library(future.apply)
})

SEED <- 123
OUTER_FOLDS <- 5
INNER_FOLDS <- 10
N_CORES <- 10
MIN_PREVALENCE <- 0.025

METADATA_FILE <- "data/metadata.csv"
ABUNDANCE_FILE <- "data/motu_abundance.xlsx"
OUTPUT_DIR <- "results/made_auc"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(
  METADATA_FILE,
  check.names = FALSE
)

metadata <- metadata[!is.na(metadata$SampleID), ]
rownames(metadata) <- metadata$SampleID

outcomes <- colnames(metadata)[34:50]

motu <- read.xlsx(
  ABUNDANCE_FILE,
  rowNames = TRUE,
  check.names = FALSE
)

motu <- motu[rownames(metadata), , drop = FALSE]

selected_features <- colnames(motu)[
  colMeans(motu > 0, na.rm = TRUE) > MIN_PREVALENCE
]

clr_transform <- function(x, pseudocount = NULL) {
  x <- as.matrix(x)

  if (is.null(pseudocount)) {
    pseudocount <- min(x[x > 0], na.rm = TRUE) / 2
  }

  x[x == 0] <- pseudocount
  log_x <- log(x)

  list(
    data = sweep(log_x, 1, rowMeans(log_x), "-"),
    pseudocount = pseudocount
  )
}

scale_train_test <- function(train, test) {
  center <- colMeans(train)
  scale <- apply(train, 2, sd)
  scale[scale == 0 | is.na(scale)] <- 1

  list(
    train = sweep(sweep(train, 2, center, "-"), 2, scale, "/"),
    test = sweep(sweep(test, 2, center, "-"), 2, scale, "/")
  )
}

make_folds <- function(y, k, seed) {
  set.seed(seed)

  folds <- createFolds(
    factor(y, levels = c(0, 1)),
    k = k,
    returnTrain = FALSE
  )

  fold_id <- integer(length(y))

  for (i in seq_along(folds)) {
    fold_id[folds[[i]]] <- i
  }

  fold_id
}

calculate_auc <- function(y, prediction) {
  as.numeric(
    auc(
      roc(
        y,
        prediction,
        quiet = TRUE,
        direction = "<"
      )
    )
  )
}

run_outcome <- function(outcome) {
  y <- as.numeric(metadata[[outcome]])
  sample_ids <- rownames(metadata)[!is.na(y)]
  y <- y[!is.na(y)]
  names(y) <- sample_ids

  outer_fold_id <- make_folds(
    y,
    OUTER_FOLDS,
    SEED
  )

  names(outer_fold_id) <- sample_ids

  train_y_all <- numeric(0)
  train_pred_all <- numeric(0)
  test_y_all <- numeric(0)
  test_pred_all <- numeric(0)

  fold_train_auc <- numeric(0)
  fold_test_auc <- numeric(0)
  lambda_values <- numeric(0)
  selected_feature_counts <- numeric(0)

  for (fold in seq_len(OUTER_FOLDS)) {
    test_ids <- names(outer_fold_id)[outer_fold_id == fold]
    train_ids <- names(outer_fold_id)[outer_fold_id != fold]

    train_raw <- motu[train_ids, selected_features, drop = FALSE]
    test_raw <- motu[test_ids, selected_features, drop = FALSE]

    train_clr <- clr_transform(train_raw)
    test_clr <- clr_transform(
      test_raw,
      train_clr$pseudocount
    )

    scaled <- scale_train_test(
      train_clr$data,
      test_clr$data
    )

    x_train <- scaled$train
    x_test <- scaled$test

    y_train <- y[rownames(x_train)]
    y_test <- y[rownames(x_test)]

    inner_fold_id <- make_folds(
      y_train,
      min(
        INNER_FOLDS,
        sum(y_train == 0),
        sum(y_train == 1)
      ),
      SEED
    )

    set.seed(SEED)

    cv_fit <- cv.glmnet(
      x = x_train,
      y = y_train,
      family = "binomial",
      type.measure = "deviance",
      foldid = inner_fold_id,
      standardize = FALSE
    )

    lambda <- cv_fit$lambda.1se

    train_prediction <- as.numeric(
      predict(
        cv_fit$glmnet.fit,
        newx = x_train,
        s = lambda,
        type = "link"
      )
    )

    test_prediction <- as.numeric(
      predict(
        cv_fit$glmnet.fit,
        newx = x_test,
        s = lambda,
        type = "link"
      )
    )

    coefficients <- as.matrix(
      coef(
        cv_fit$glmnet.fit,
        s = lambda
      )
    )

    selected_feature_count <- sum(
      coefficients[rownames(coefficients) != "(Intercept)", 1] != 0
    )

    train_y_all <- c(train_y_all, y_train)
    train_pred_all <- c(train_pred_all, train_prediction)
    test_y_all <- c(test_y_all, y_test)
    test_pred_all <- c(test_pred_all, test_prediction)

    fold_train_auc <- c(
      fold_train_auc,
      calculate_auc(y_train, train_prediction)
    )

    fold_test_auc <- c(
      fold_test_auc,
      calculate_auc(y_test, test_prediction)
    )

    lambda_values <- c(lambda_values, lambda)
    selected_feature_counts <- c(
      selected_feature_counts,
      selected_feature_count
    )
  }

  data.frame(
    outcome = outcome,
    model = "species_only",
    N = length(y),
    n_case = sum(y == 1),
    n_control = sum(y == 0),
    n_features = length(selected_features),
    pooled_train_auc = calculate_auc(
      train_y_all,
      train_pred_all
    ),
    pooled_test_auc = calculate_auc(
      test_y_all,
      test_pred_all
    ),
    mean_fold_train_auc = mean(fold_train_auc),
    sd_fold_train_auc = sd(fold_train_auc),
    mean_fold_test_auc = mean(fold_test_auc),
    sd_fold_test_auc = sd(fold_test_auc),
    mean_lambda_1se = mean(lambda_values),
    sd_lambda_1se = sd(lambda_values),
    mean_selected_features = mean(selected_feature_counts),
    sd_selected_features = sd(selected_feature_counts),
    stringsAsFactors = FALSE
  )
}

plan(
  multisession,
  workers = N_CORES
)

results <- future_lapply(
  outcomes,
  run_outcome,
  future.seed = TRUE
)

results <- bind_rows(results)

write.csv(
  results,
  file.path(
    OUTPUT_DIR,
    "species_only_5fold_oof_auc.csv"
  ),
  row.names = FALSE
)

write.xlsx(
  results,
  file.path(
    OUTPUT_DIR,
    "species_only_5fold_oof_auc.xlsx"
  ),
  overwrite = TRUE
)
