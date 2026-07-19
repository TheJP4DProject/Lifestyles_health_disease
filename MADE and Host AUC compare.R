# This script compares Host-only, MADE-only, and Host+MADE models
# using 10-fold out-of-fold prediction and paired DeLong tests.

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
OUTER_FOLDS <- 10
INNER_FOLDS <- 10
N_CORES <- 10
MIN_PREVALENCE <- 0.025

METADATA_FILE <- "data/metadata.csv"
ABUNDANCE_FILE <- "data/motu_abundance.xlsx"
OUTPUT_DIR <- "results/host_vs_made"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(
  METADATA_FILE,
  check.names = FALSE
)

metadata <- metadata[!is.na(metadata$SampleID), ]
rownames(metadata) <- metadata$SampleID

outcomes <- colnames(metadata)[34:50]
host_variables <- c("Age", "Sex", "BMI")

motu <- read.xlsx(
  ABUNDANCE_FILE,
  rowNames = TRUE,
  check.names = FALSE
)

motu <- motu[rownames(metadata), , drop = FALSE]

selected_features <- colnames(motu)[
  colMeans(motu > 0, na.rm = TRUE) > MIN_PREVALENCE
]

clr_transform <- function(x) {
  x <- as.matrix(x)
  x[x == 0] <- min(x[x > 0], na.rm = TRUE) / 2
  log_x <- log(x)
  sweep(log_x, 1, rowMeans(log_x), "-")
}

motu_clr <- clr_transform(motu)
motu_clr <- motu_clr[, selected_features, drop = FALSE]

if (ncol(motu_clr) > 0) {
  motu_clr <- motu_clr[, -1, drop = FALSE]
}

X_made <- scale(motu_clr)

host_complete <- complete.cases(
  metadata[, host_variables, drop = FALSE]
)

host_design_complete <- model.matrix(
  ~ .,
  data = metadata[
    host_complete,
    host_variables,
    drop = FALSE
  ]
)[, -1, drop = FALSE]

X_host <- matrix(
  NA_real_,
  nrow = nrow(metadata),
  ncol = ncol(host_design_complete),
  dimnames = list(
    rownames(metadata),
    colnames(host_design_complete)
  )
)

X_host[rownames(host_design_complete), ] <-
  host_design_complete

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

scale_train_test <- function(train, test) {
  center <- colMeans(train)
  scale_value <- apply(train, 2, sd)
  scale_value[is.na(scale_value) | scale_value == 0] <- 1

  list(
    train = sweep(
      sweep(train, 2, center, "-"),
      2,
      scale_value,
      "/"
    ),
    test = sweep(
      sweep(test, 2, center, "-"),
      2,
      scale_value,
      "/"
    )
  )
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

fit_logistic_model <- function(
    x_train,
    x_test,
    y_train
) {
  x_train <- as.data.frame(x_train)
  x_test <- as.data.frame(x_test)

  colnames(x_train) <- paste0("V", seq_len(ncol(x_train)))
  colnames(x_test) <- colnames(x_train)

  fit <- glm(
    y_train ~ .,
    data = x_train,
    family = binomial()
  )

  list(
    train = as.numeric(
      predict(fit, newdata = x_train, type = "link")
    ),
    test = as.numeric(
      predict(fit, newdata = x_test, type = "link")
    )
  )
}

run_outcome <- function(outcome) {
  y_all <- as.numeric(metadata[[outcome]])
  names(y_all) <- rownames(metadata)

  sample_ids <- rownames(metadata)[
    !is.na(y_all) &
      complete.cases(X_made)
  ]

  y <- y_all[sample_ids]
  fold_id <- make_folds(
    y,
    OUTER_FOLDS,
    SEED
  )

  names(fold_id) <- sample_ids

  made_predictions <- data.frame()
  comparison_predictions <- data.frame()

  for (fold in seq_len(OUTER_FOLDS)) {
    train_ids <- names(fold_id)[fold_id != fold]
    test_ids <- names(fold_id)[fold_id == fold]

    y_train <- y[train_ids]
    y_test <- y[test_ids]

    x_train_made <- X_made[train_ids, , drop = FALSE]
    x_test_made <- X_made[test_ids, , drop = FALSE]

    inner_k <- min(
      INNER_FOLDS,
      sum(y_train == 0),
      sum(y_train == 1)
    )

    inner_fold_id <- make_folds(
      y_train,
      inner_k,
      SEED + fold
    )

    set.seed(SEED + fold)

    cv_fit <- cv.glmnet(
      x = x_train_made,
      y = y_train,
      family = "binomial",
      type.measure = "deviance",
      foldid = inner_fold_id,
      standardize = FALSE
    )

    made_train_score <- as.numeric(
      predict(
        cv_fit,
        newx = x_train_made,
        s = "lambda.1se",
        type = "link"
      )
    )

    made_test_score <- as.numeric(
      predict(
        cv_fit,
        newx = x_test_made,
        s = "lambda.1se",
        type = "link"
      )
    )

    names(made_train_score) <- train_ids
    names(made_test_score) <- test_ids

    made_predictions <- bind_rows(
      made_predictions,
      data.frame(
        outcome = outcome,
        sample_id = test_ids,
        fold = fold,
        y = y_test,
        made_score = made_test_score
      )
    )

    host_train <- X_host[train_ids, , drop = FALSE]
    host_test <- X_host[test_ids, , drop = FALSE]

    train_ok <- complete.cases(host_train)
    test_ok <- complete.cases(host_test)

    train_ids_complete <- train_ids[train_ok]
    test_ids_complete <- test_ids[test_ok]

    y_train_complete <- y[train_ids_complete]
    y_test_complete <- y[test_ids_complete]

    scaled_host <- scale_train_test(
      host_train[train_ok, , drop = FALSE],
      host_test[test_ok, , drop = FALSE]
    )

    made_train_matrix <- matrix(
      made_train_score[train_ids_complete],
      ncol = 1
    )

    made_test_matrix <- matrix(
      made_test_score[test_ids_complete],
      ncol = 1
    )

    host_fit <- fit_logistic_model(
      scaled_host$train,
      scaled_host$test,
      y_train_complete
    )

    made_fit <- fit_logistic_model(
      made_train_matrix,
      made_test_matrix,
      y_train_complete
    )

    combined_fit <- fit_logistic_model(
      cbind(
        scaled_host$train,
        made_train_matrix
      ),
      cbind(
        scaled_host$test,
        made_test_matrix
      ),
      y_train_complete
    )

    comparison_predictions <- bind_rows(
      comparison_predictions,
      data.frame(
        outcome = outcome,
        sample_id = test_ids_complete,
        fold = fold,
        y = y_test_complete,
        pred_host_only = host_fit$test,
        pred_made_only = made_fit$test,
        pred_host_made = combined_fit$test
      )
    )
  }

  made_auc <- data.frame(
    outcome = outcome,
    model = "MADE only",
    N = nrow(made_predictions),
    pooled_test_auc = calculate_auc(
      made_predictions$y,
      made_predictions$made_score
    )
  )

  model_predictions <- list(
    "Host only" =
      comparison_predictions$pred_host_only,
    "MADE only" =
      comparison_predictions$pred_made_only,
    "Host + MADE" =
      comparison_predictions$pred_host_made
  )

  comparison_auc <- bind_rows(
    lapply(names(model_predictions), function(model_name) {
      data.frame(
        outcome = outcome,
        model = model_name,
        N = nrow(comparison_predictions),
        pooled_test_auc = calculate_auc(
          comparison_predictions$y,
          model_predictions[[model_name]]
        )
      )
    })
  )

  comparisons <- list(
    c("Host only", "MADE only"),
    c("Host only", "Host + MADE"),
    c("MADE only", "Host + MADE")
  )

  delong_results <- bind_rows(
    lapply(comparisons, function(models) {
      prediction_1 <- model_predictions[[models[1]]]
      prediction_2 <- model_predictions[[models[2]]]

      roc_1 <- roc(
        comparison_predictions$y,
        prediction_1,
        quiet = TRUE,
        direction = "<"
      )

      roc_2 <- roc(
        comparison_predictions$y,
        prediction_2,
        quiet = TRUE,
        direction = "<"
      )

      data.frame(
        outcome = outcome,
        comparison = paste(
          models,
          collapse = " vs "
        ),
        model_1 = models[1],
        model_2 = models[2],
        auc_model_1 = as.numeric(auc(roc_1)),
        auc_model_2 = as.numeric(auc(roc_2)),
        delta_auc_model_2_minus_model_1 =
          as.numeric(auc(roc_2)) -
          as.numeric(auc(roc_1)),
        p_value = roc.test(
          roc_1,
          roc_2,
          paired = TRUE,
          method = "delong"
        )$p.value,
        n_paired = nrow(comparison_predictions)
      )
    })
  )

  list(
    made_auc = made_auc,
    comparison_auc = comparison_auc,
    predictions = comparison_predictions,
    delong = delong_results
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

plan(sequential)

made_auc_results <- bind_rows(
  lapply(results, `[[`, "made_auc")
)

comparison_auc_results <- bind_rows(
  lapply(results, `[[`, "comparison_auc")
)

prediction_results <- bind_rows(
  lapply(results, `[[`, "predictions")
)

delong_results <- bind_rows(
  lapply(results, `[[`, "delong")
)

delong_results$fdr <- p.adjust(
  delong_results$p_value,
  method = "BH"
)

write.csv(
  made_auc_results,
  file.path(
    OUTPUT_DIR,
    "made_10fold_oof_auc.csv"
  ),
  row.names = FALSE
)

write.csv(
  comparison_auc_results,
  file.path(
    OUTPUT_DIR,
    "host_vs_made_10fold_oof_auc.csv"
  ),
  row.names = FALSE
)

write.csv(
  prediction_results,
  file.path(
    OUTPUT_DIR,
    "host_vs_made_10fold_oof_predictions.csv"
  ),
  row.names = FALSE
)

write.csv(
  delong_results,
  file.path(
    OUTPUT_DIR,
    "host_vs_made_delong_results.csv"
  ),
  row.names = FALSE
)
