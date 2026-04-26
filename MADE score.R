library(doParallel)
library(foreach)
library(glmnet)
library(pROC)
library(caret)
library(dplyr)
library(openxlsx)
library(compositions)

setwd(".")

metadata <- as.data.frame(read.csv("meta.atc.filter.csv"))
metadata <- metadata[!is.na(metadata$METAF),]

fixkey <- colnames(metadata)[c(3:5,51:55)]
diskey <- colnames(metadata)[6:50]

rownames(metadata) <- metadata$METAF

motu <- as.data.frame(read.xlsx("RAdata_motu3_7562samples_newsp_20250629.xlsx",
                                rowNames = TRUE,
                                check.names = FALSE,
                                sep.names = " "))
motu <- motu[rownames(metadata), , drop = FALSE]

modify_rownames <- function(df) {
  rownames(df) <- sapply(rownames(df), function(x) {
    v <- suppressWarnings(as.numeric(x))
    if (!is.na(v)) as.character(v) else x
  })
  df
}

clr_transform <- function(M) {
  M[M == 0] <- min(M[M > 0]) / 2
  logM <- log(M)
  sweep(logM, 1, rowMeans(logM), "-")
}

standardize_cols <- function(M) {
  if (ncol(M) == 0) return(M)
  M <- scale(M)
  colnames(M) <- colnames(M)
  rownames(M) <- rownames(M)
  M
}

save_dir <- "model_results_sp_only"
dir.create(save_dir, showWarnings = FALSE)
rds_dir <- file.path(save_dir, "rds")
dir.create(rds_dir, showWarnings = FALSE)

motu_check <- colnames(motu)[(colSums(motu != 0) / nrow(motu)) > 0.025]

motu_clr <- clr_transform(motu)
motu_clr <- motu_clr[, motu_check, drop = FALSE]

X_sp <- standardize_cols(motu_clr)

demo_data <- metadata[, c("Age", "Sex", "BMI",
                          "HMG.CoA.reductase.inhibitors",
                          "Angiotensin.II.antagonists.plain",
                          "Proton.pump.inhibitors",
                          "Any_antibiotics",
                          "Biguanides"), drop = FALSE]

demo_matrix <- model.matrix(~ Age + Sex + BMI +
                              HMG.CoA.reductase.inhibitors +
                              Angiotensin.II.antagonists.plain +
                              Proton.pump.inhibitors +
                              Any_antibiotics + Biguanides,
                            data = demo_data)[, -1, drop = FALSE]

demo_covars <- matrix(NA, nrow(metadata), ncol(demo_matrix),
                      dimnames = list(rownames(metadata), colnames(demo_matrix)))

demo_covars[complete.cases(demo_data), ] <- scale(demo_matrix)

lifestyle_vars <- colnames(metadata)[6:33]

life_data <- metadata[, lifestyle_vars, drop = FALSE]

life_matrix <- model.matrix(~ ., data = life_data)[, -1, drop = FALSE]

lifestyle_covars <- matrix(NA, nrow(metadata), ncol(life_matrix),
                           dimnames = list(rownames(metadata), colnames(life_matrix)))

lifestyle_covars[complete.cases(life_data), ] <- scale(life_matrix)

model_defs <- list(
  sp = list(sp = "o", lifestyle = "x", ASB = "x")
)

build_X_and_cate <- function(def, sample_keep) {

  X_parts <- list()
  cate_vec <- character(0)
  features <- character(0)

  if (def$sp != "x") {
    X_parts <- c(X_parts, list(X_sp[sample_keep, , drop = FALSE]))
    cate_vec <- c(cate_vec, rep("sp", ncol(X_sp)))
    features <- c(features, colnames(X_sp))
  }

  if (def$lifestyle != "x") {
    X_parts <- c(X_parts, list(lifestyle_covars[sample_keep, , drop = FALSE]))
    cate_vec <- c(cate_vec, rep("lifestyle", ncol(lifestyle_covars)))
    features <- c(features, colnames(lifestyle_covars))
  }

  if (def$ASB != "x") {
    X_parts <- c(X_parts, list(demo_covars[sample_keep, , drop = FALSE]))
    cate_vec <- c(cate_vec, rep("ASB", ncol(demo_covars)))
    features <- c(features, colnames(demo_covars))
  }

  if (length(X_parts) == 0) return(list(X = NULL, cate = NULL, features = NULL))

  X_model <- do.call(cbind, X_parts)
  ok <- complete.cases(X_model)
  X_model <- X_model[ok, , drop = FALSE]

  list(X = X_model,
       cate = cate_vec,
       features = features,
       row_keep = rownames(X_model))
}

build_pf_for_def <- function(def, X_model) {

  pf <- c()

  if (def$sp != "x") {
    pf <- c(pf, rep(1, ncol(X_sp[, colnames(X_sp) %in% colnames(X_model), drop = FALSE])))
  }

  if (def$lifestyle != "x") {
    pf <- c(pf, rep(1, ncol(lifestyle_covars)))
  }

  if (def$ASB != "x") {
    pf <- c(pf, rep(1, ncol(demo_covars)))
  }

  pf
}

outcome_cols <- 34:50
outcome_names <- colnames(metadata)[outcome_cols]

cl <- makeCluster(10)
registerDoParallel(cl)

GLOBAL_SEED <- 123

task_list <- list()

for (m in names(model_defs)) {
  for (outcome in outcome_names) {
    task_list[[paste0(outcome, "_", m)]] <- list(outcome = outcome, model = m)
  }
}

clusterExport(cl, c("metadata","X_sp","demo_covars","lifestyle_covars",
                    "model_defs","build_X_and_cate",
                    "build_pf_for_def","GLOBAL_SEED","rds_dir"))

results <- foreach(task = task_list,
                   .packages = c("glmnet","pROC","dplyr"),
                   .errorhandling = "pass") %dopar% {

  tryCatch({

    outcome <- task$outcome
    m <- task$model
    def <- model_defs[[m]]

    y <- metadata[, outcome]
    keep <- !is.na(y)

    b <- build_X_and_cate(def, keep)
    X <- b$X

    if (is.null(X) || ncol(X) == 0)
      return(list(status = "nopredictor"))

    cc <- complete.cases(X)
    X <- X[cc, , drop = FALSE]

    y <- as.numeric(y[keep])[cc]

    y_factor <- factor(y, levels = c(0, 1))

    pf <- rep(1, ncol(X))

    set.seed(GLOBAL_SEED)

    cvfit <- cv.glmnet(
      x = as.matrix(X),
      y = y_factor,
      family = "binomial",
      nfold = 10,
      penalty.factor = pf,
      standardize = FALSE
    )

    lambda <- cvfit$lambda.1se

    score <- predict(cvfit$glmnet.fit,
                     as.matrix(X),
                     s = lambda,
                     type = "link")

    score <- as.numeric(score)

    auc_val <- auc(roc(y, score, quiet = TRUE))

    cf <- coef(cvfit$glmnet.fit, s = lambda)
    cf <- as.numeric(cf)
    names(cf) <- rownames(coef(cvfit$glmnet.fit))

    nonzero <- names(cf)[cf != 0]
    nonzero <- nonzero[nonzero != "(Intercept)"]

    res <- list(
      outcome = outcome,
      model = m,
      auc = auc_val,
      score = score,
      features = nonzero,
      n = length(y)
    )

    saveRDS(res, file.path(rds_dir, paste0(outcome, "_", m, ".rds")))

    return(res)

  }, error = function(e) {
    list(status = "error")
  })
}

stopCluster(cl)

cat("DONE\n")