library(vegan)
library(openxlsx)
library(doParallel)
library(foreach)

registerDoParallel(cores = 10)

metadata <- read.xlsx("../meta2_full.xlsx")
metadata <- metadata[, c(-2, -3)]
rownames(metadata) <- metadata[, 1]

modify_rownames <- function(df) {
old_names <- rownames(df)
new_names <- sapply(old_names, function(name) {
num_name <- suppressWarnings(as.numeric(name))
if (!is.na(num_name)) {
return(as.character(num_name))
} else {
return(name)
}
})
rownames(df) <- new_names
return(df)
}
metadata <- modify_rownames(metadata)

fixkey <- colnames(metadata)[2:11]
diskey <- colnames(metadata)[12:72]

cat("预载距离矩阵...\n")
bc_files <- c(
motu = "bc_motusp_6803.rds",
KO = "bc_KO_6803.rds",

)

bc_matrices <- list()
for (type in names(bc_files)) {
bc_mat <- readRDS(file.path(".", bc_files[type]))
bc_matrices[[type]] <- bc_mat

}

tasks <- list()
for (type in names(bc_matrices)) {
for (meta_var in diskey) {
tasks[[length(tasks) + 1]] <- list(
data_type = type,
metadata_var = meta_var,
bc_matrix = bc_matrices[[type]]
)
}
}

output_dir <- "adonis2"

start_time <- Sys.time()

progress_file <- file.path(output_dir, "progress.txt")

null_results <- foreach(task = tasks, .packages = "vegan") %dopar% {
type <- task$data_type
meta_var <- task$metadata_var
bc_matrix <- task$bc_matrix

task_id <- paste(type, gsub("[^a-zA-Z0-9]", "*", meta_var), sep = "*")
output_file <- file.path(output_dir, paste0(task_id, ".csv"))

if (file.exists(output_file)) {
return(NULL)
}

if (!(meta_var %in% colnames(metadata))) {
return(NULL)
}

temp_df <- metadata

for (col in c(meta_var, fixkey)) {
if (col %in% colnames(temp_df) && !is.numeric(temp_df[[col]])) {
temp_df[[col]] <- as.factor(temp_df[[col]])
}
}

needed_cols <- c(meta_var, fixkey)
needed_cols <- needed_cols[needed_cols %in% colnames(temp_df)]
complete_cases <- complete.cases(temp_df[, needed_cols, drop = FALSE])
temp_df_subset <- temp_df[complete_cases, , drop = FALSE]

if (nrow(temp_df_subset) < 10) {
return(NULL)
}

common_samples <- intersect(rownames(temp_df_subset), rownames(bc_matrix))
if (length(common_samples) < 10) {
return(NULL)
}

temp_df_subset <- temp_df_subset[common_samples, , drop = FALSE]
dist_subset <- as.dist(bc_matrix[common_samples, common_samples])

formula_str <- paste("dist_subset ~", meta_var, "+", paste(fixkey, collapse = " + "))

result_df <- data.frame()
try({
set.seed(123)
adonis_result <- adonis2(
as.formula(formula_str),
data = temp_df_subset,
permutations = 999,
by = "term",
parallel = 1
)

if (nrow(adonis_result) > 0) {
  result_df <- as.data.frame(adonis_result)
  result_df$term <- rownames(adonis_result)
  result_df$metadata_var <- meta_var
  result_df$n_samples <- length(common_samples)
  result_df$model <- formula_str
  result_df$data_type <- type
}

})

if (nrow(result_df) > 0) {
write.csv(result_df, file = output_file, row.names = FALSE)
}

return(NULL)
}

end_time <- Sys.time()

stopImplicitCluster()

all_csv_files <- list.files(output_dir, pattern = "\.csv$", full.names = TRUE)

if (length(all_csv_files) > 0) {
all_results <- data.frame()

for (csv_file in all_csv_files) {
try({
df <- read.csv(csv_file)
all_results <- rbind(all_results, df)
})
}

if (nrow(all_results) > 0) {
write.csv(all_results, file = "adonis2_all_results.csv", row.names = FALSE)
saveRDS(all_results, file = "adonis2_all_results.rds", compress = 6)


for (type in unique(all_results$data_type)) {
  type_data <- all_results[all_results$data_type == type, ]
  n_vars <- length(unique(type_data$metadata_var))
  n_rows <- nrow(type_data)
  n_sig <- sum(type_data$`Pr(>F)` < 0.05, na.rm = TRUE)
}


}
}
