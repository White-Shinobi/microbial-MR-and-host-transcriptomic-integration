#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glmnet)
  library(e1071)
})

set.seed(20260719)

args <- commandArgs(trailingOnly = TRUE)
bootstrap_reps <- if (length(args)) as.integer(args[[1]]) else 200L
if (!is.finite(bootstrap_reps) || bootstrap_reps < 100L || bootstrap_reps > 500L) {
  stop("bootstrap_reps must be an integer between 100 and 500")
}

script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
out_dir <- file.path(root, "results/bulk/crc/machine_learning_stability_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = TRUE, row.names = FALSE, na = "")
}

patient_fold_vector <- function(patient_ids, k = 10L, seed = 20260719L) {
  set.seed(seed)
  ids <- unique(patient_ids)
  assignment <- sample(rep(seq_len(k), length.out = length(ids)))
  unname(setNames(assignment, ids)[patient_ids])
}

lasso_select <- function(x, y_binary, foldid) {
  cv <- cv.glmnet(
    x, y_binary, family = "binomial", alpha = 1, standardize = TRUE,
    foldid = foldid, type.measure = "deviance"
  )
  extract <- function(s) {
    beta <- as.matrix(coef(cv, s = s))
    intersect(colnames(x), rownames(beta)[beta[, 1] != 0])
  }
  extract_coef <- function(s) {
    beta <- as.matrix(coef(cv, s = s))
    unname(beta[colnames(x), 1])
  }
  list(
    selected_min = extract("lambda.min"),
    selected_1se = extract("lambda.1se"),
    coefficients_min = setNames(extract_coef("lambda.min"), colnames(x)),
    coefficients_1se = setNames(extract_coef("lambda.1se"), colnames(x)),
    lambda_min = unname(cv$lambda.min),
    lambda_1se = unname(cv$lambda.1se)
  )
}

svm_rfe_select <- function(x, y_factor, foldid) {
  current <- colnames(x)
  subsets <- list()
  performance <- data.frame()
  elimination <- data.frame()
  repeat {
    pred <- factor(rep(NA_character_, length(y_factor)), levels = levels(y_factor))
    for (fold in sort(unique(foldid))) {
      train <- foldid != fold
      test <- foldid == fold
      fit <- svm(
        x = x[train, current, drop = FALSE], y = y_factor[train],
        kernel = "linear", cost = 1, scale = TRUE, type = "C-classification"
      )
      pred[test] <- predict(fit, x[test, current, drop = FALSE])
    }
    subsets[[as.character(length(current))]] <- current
    performance <- rbind(
      performance,
      data.frame(
        n_features = length(current),
        cv_accuracy = mean(pred == y_factor),
        retained_genes = paste(current, collapse = ";"),
        stringsAsFactors = FALSE
      )
    )
    if (length(current) == 1L) break
    full <- svm(
      x = x[, current, drop = FALSE], y = y_factor,
      kernel = "linear", cost = 1, scale = TRUE, type = "C-classification"
    )
    weights <- as.numeric(t(full$coefs) %*% full$SV)
    names(weights) <- current
    removed <- names(which.min(weights^2))
    elimination <- rbind(
      elimination,
      data.frame(
        elimination_step = nrow(elimination) + 1L,
        removed_gene = removed,
        weight_squared = unname(weights[removed]^2),
        stringsAsFactors = FALSE
      )
    )
    current <- setdiff(current, removed)
  }
  best_accuracy <- max(performance$cv_accuracy)
  optimal_n <- min(
    performance$n_features[abs(performance$cv_accuracy - best_accuracy) < 1e-12]
  )
  performance$selected_subset <- performance$n_features == optimal_n
  list(
    selected = subsets[[as.character(optimal_n)]],
    optimal_n = optimal_n,
    max_accuracy = best_accuracy,
    performance = performance,
    elimination = elimination
  )
}

read_ml_inputs <- function(project_root) {
  candidate_file <- file.path(
    project_root, "results/bulk/crc/intersection/CRC_DEG_GMRG_candidate_genes.tsv"
  )
  discovery_dir <- file.path(project_root, "results/bulk/crc/discovery")
  ml_dir <- file.path(project_root, "results/bulk/crc/machine_learning")
  candidates <- read.delim(candidate_file, stringsAsFactors = FALSE, check.names = FALSE)
  expr <- readRDS(file.path(discovery_dir, "GSE44076_gene_expression_matrix.rds"))
  manifest <- read.delim(
    file.path(discovery_dir, "GSE44076_sample_manifest.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  membership <- read.delim(
    file.path(ml_dir, "LASSO_SVM_RFE_intersection.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  gate <- read.delim(
    file.path(ml_dir, "biomarker_gate_summary.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  lasso <- read.delim(
    file.path(ml_dir, "LASSO_selection.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  grouped_fold_file <- file.path(
    ml_dir, "GSE44076_fixed_patient_grouped_10fold_assignment.tsv"
  )
  sample_fold_file <- if (file.exists(grouped_fold_file)) {
    grouped_fold_file
  } else {
    file.path(ml_dir, "GSE44076_stratified_10fold_assignment.tsv")
  }
  sample_folds <- read.delim(
    sample_fold_file,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  list(
    candidates = candidates,
    expr = expr,
    manifest = manifest,
    membership = membership,
    gate = gate,
    lasso = lasso,
    sample_folds = sample_folds
  )
}

v2 <- read_ml_inputs(root)
candidates <- unique(v2$candidates$gene_symbol)
stopifnot(
  length(candidates) > 1L,
  all(candidates %in% rownames(v2$expr)),
  all(v2$manifest$sample_id %in% colnames(v2$expr))
)

x <- t(v2$expr[candidates, v2$manifest$sample_id, drop = FALSE])
y_factor <- factor(v2$manifest$group, levels = c("Normal", "Tumor"))
y_binary <- as.integer(y_factor == "Tumor")
patients <- v2$manifest$individual_id
stopifnot(
  length(unique(patients)) == 98L,
  all(table(patients) == 2L),
  all(tapply(as.character(y_factor), patients, function(z) {
    setequal(z, c("Normal", "Tumor"))
  }))
)

# Fixed patient-grouped 10-fold assignment: paired samples always stay together.
grouped_fold <- patient_fold_vector(patients, 10L, 20260719L)
fold_manifest <- data.frame(
  sample_id = rownames(x),
  individual_id = patients,
  group = as.character(y_factor),
  patient_grouped_fold = grouped_fold,
  stringsAsFactors = FALSE
)
write_tsv(
  fold_manifest,
  file.path(out_dir, "GSE44076_fixed_patient_grouped_10fold_assignment.tsv")
)

# Full-data patient-grouped LASSO at lambda.min and lambda.1se and SVM-RFE.
full_lasso <- lasso_select(x, y_binary, grouped_fold)
full_svm <- svm_rfe_select(x, y_factor, grouped_fold)
full_membership <- data.frame(
  gene_symbol = candidates,
  LASSO_lambda_min = candidates %in% full_lasso$selected_min,
  LASSO_lambda_1se = candidates %in% full_lasso$selected_1se,
  SVM_RFE = candidates %in% full_svm$selected,
  intersection_lambda_min = candidates %in%
    intersect(full_lasso$selected_min, full_svm$selected),
  intersection_lambda_1se = candidates %in%
    intersect(full_lasso$selected_1se, full_svm$selected),
  LASSO_only_lambda_min = candidates %in%
    setdiff(full_lasso$selected_min, full_svm$selected),
  LASSO_only_lambda_1se = candidates %in%
    setdiff(full_lasso$selected_1se, full_svm$selected),
  SVM_RFE_only_vs_lambda_min = candidates %in%
    setdiff(full_svm$selected, full_lasso$selected_min),
  SVM_RFE_only_vs_lambda_1se = candidates %in%
    setdiff(full_svm$selected, full_lasso$selected_1se),
  stringsAsFactors = FALSE
)
full_lasso_coefficients <- data.frame(
  gene_symbol = candidates,
  coefficient_lambda_min = unname(full_lasso$coefficients_min[candidates]),
  coefficient_lambda_1se = unname(full_lasso$coefficients_1se[candidates]),
  stringsAsFactors = FALSE
)
full_summary <- data.frame(
  metric = c(
    "lambda_min", "lambda_1se", "LASSO_lambda_min_genes",
    "LASSO_lambda_1se_genes", "SVM_RFE_optimal_feature_count",
    "SVM_RFE_max_CV_accuracy", "SVM_RFE_genes",
    "intersection_lambda_min_genes", "intersection_lambda_1se_genes"
  ),
  value = c(
    signif(full_lasso$lambda_min, 10),
    signif(full_lasso$lambda_1se, 10),
    paste(full_lasso$selected_min, collapse = ";"),
    paste(full_lasso$selected_1se, collapse = ";"),
    full_svm$optimal_n,
    signif(full_svm$max_accuracy, 10),
    paste(full_svm$selected, collapse = ";"),
    paste(intersect(full_lasso$selected_min, full_svm$selected), collapse = ";"),
    paste(intersect(full_lasso$selected_1se, full_svm$selected), collapse = ";")
  ),
  stringsAsFactors = FALSE
)
write_tsv(
  full_membership,
  file.path(out_dir, "V2_patient_grouped_full_data_method_membership.tsv")
)
write_tsv(
  full_summary,
  file.path(out_dir, "V2_patient_grouped_full_data_summary.tsv")
)
write_tsv(
  full_lasso_coefficients,
  file.path(out_dir, "V2_patient_grouped_LASSO_coefficients_min_1se.tsv")
)
write_tsv(
  full_svm$performance,
  file.path(out_dir, "V2_patient_grouped_SVM_RFE_performance.tsv")
)

# Apply the unchanged two-dataset direction-and-P-value gate separately to the
# lambda.min and lambda.1se intersections for this sensitivity report.
validation_expr <- readRDS(
  file.path(root, "results/bulk/crc/validation/GSE87211_gene_expression_matrix.rds")
)
validation_manifest <- read.delim(
  file.path(root, "results/bulk/crc/validation/GSE87211_sample_manifest.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
expression_stat <- function(gene, matrix, manifest) {
  values <- as.numeric(matrix[gene, manifest$sample_id])
  tumor <- values[manifest$group == "Tumor"]
  normal <- values[manifest$group == "Normal"]
  difference <- median(tumor) - median(normal)
  c(
    median_difference = difference,
    wilcoxon_p = wilcox.test(tumor, normal, exact = FALSE)$p.value,
    direction = ifelse(difference > 0, "Up", "Down")
  )
}
min_intersection <- intersect(full_lasso$selected_min, full_svm$selected)
one_se_intersection <- intersect(full_lasso$selected_1se, full_svm$selected)
gate_genes <- unique(c(min_intersection, one_se_intersection))
lambda_gate <- do.call(rbind, lapply(gate_genes, function(gene) {
  discovery <- expression_stat(gene, v2$expr, v2$manifest)
  validation <- expression_stat(gene, validation_expr, validation_manifest)
  consistent <- discovery[["direction"]] == validation[["direction"]]
  significant <- as.numeric(discovery[["wilcoxon_p"]]) < 0.05 &&
    as.numeric(validation[["wilcoxon_p"]]) < 0.05
  data.frame(
    gene_symbol = gene,
    in_lambda_min_intersection = gene %in% min_intersection,
    in_lambda_1se_intersection = gene %in% one_se_intersection,
    GSE44076_direction = discovery[["direction"]],
    GSE44076_median_difference = as.numeric(discovery[["median_difference"]]),
    GSE44076_wilcoxon_p = as.numeric(discovery[["wilcoxon_p"]]),
    GSE87211_direction = validation[["direction"]],
    GSE87211_median_difference = as.numeric(validation[["median_difference"]]),
    GSE87211_wilcoxon_p = as.numeric(validation[["wilcoxon_p"]]),
    consistent_direction = consistent,
    significant_both = significant,
    final_at_lambda_min = gene %in% min_intersection && consistent && significant,
    final_at_lambda_1se = gene %in% one_se_intersection && consistent && significant,
    stringsAsFactors = FALSE
  )
}))
write_tsv(
  lambda_gate,
  file.path(out_dir, "V2_patient_grouped_lambda_min_1se_expression_gate.tsv")
)

# Candidate-gene Spearman correlations. These are descriptive all-sample
# correlations and can reflect the shared tumor-normal contrast.
cor_matrix <- cor(x, method = "spearman", use = "pairwise.complete.obs")
write.table(
  cbind(gene_symbol = rownames(cor_matrix), as.data.frame(cor_matrix)),
  file.path(out_dir, "V2_candidate_gene_Spearman_correlation_matrix.tsv"),
  sep = "\t", quote = TRUE, row.names = FALSE, na = ""
)
cor_index <- which(upper.tri(cor_matrix), arr.ind = TRUE)
cor_long <- data.frame(
  gene_1 = rownames(cor_matrix)[cor_index[, 1]],
  gene_2 = colnames(cor_matrix)[cor_index[, 2]],
  spearman_rho = cor_matrix[cor_index],
  stringsAsFactors = FALSE
)
cor_long <- cor_long[order(-abs(cor_long$spearman_rho)), ]
write_tsv(
  cor_long,
  file.path(out_dir, "V2_candidate_gene_Spearman_correlations_ranked.tsv")
)

# Ten outer-fold training-set reruns. The fixed patient folds remaining in each
# training set are reused as the inner CV folds.
outer_membership <- data.frame()
for (outer in sort(unique(grouped_fold))) {
  train <- grouped_fold != outer
  inner_fold <- as.integer(factor(grouped_fold[train], levels = sort(unique(grouped_fold[train]))))
  fit_lasso <- lasso_select(x[train, , drop = FALSE], y_binary[train], inner_fold)
  fit_svm <- svm_rfe_select(x[train, , drop = FALSE], y_factor[train], inner_fold)
  outer_membership <- rbind(
    outer_membership,
    data.frame(
      outer_fold = outer,
      gene_symbol = candidates,
      LASSO_lambda_min = candidates %in% fit_lasso$selected_min,
      LASSO_lambda_1se = candidates %in% fit_lasso$selected_1se,
      SVM_RFE = candidates %in% fit_svm$selected,
      intersection_lambda_min = candidates %in%
        intersect(fit_lasso$selected_min, fit_svm$selected),
      intersection_lambda_1se = candidates %in%
        intersect(fit_lasso$selected_1se, fit_svm$selected),
      lambda_min = fit_lasso$lambda_min,
      lambda_1se = fit_lasso$lambda_1se,
      SVM_RFE_optimal_n = fit_svm$optimal_n,
      SVM_RFE_max_accuracy = fit_svm$max_accuracy,
      stringsAsFactors = FALSE
    )
  )
}
outer_frequency <- aggregate(
  cbind(
    LASSO_lambda_min, LASSO_lambda_1se, SVM_RFE,
    intersection_lambda_min, intersection_lambda_1se
  ) ~ gene_symbol,
  outer_membership,
  mean
)
names(outer_frequency)[-1] <- paste0(names(outer_frequency)[-1], "_frequency")
write_tsv(
  outer_membership,
  file.path(out_dir, "V2_patient_grouped_outer_fold_selection_membership.tsv")
)
write_tsv(
  outer_frequency,
  file.path(out_dir, "V2_patient_grouped_outer_fold_selection_frequency.tsv")
)

# Patient-pair bootstrap of the current V2 candidate-gene ML selection stage.
bootstrap_membership <- data.frame()
bootstrap_diagnostics <- data.frame()
unique_patients <- unique(patients)
for (b in seq_len(bootstrap_reps)) {
  set.seed(20260719L + b)
  sampled <- sample(unique_patients, length(unique_patients), replace = TRUE)
  rows <- unlist(lapply(sampled, function(id) which(patients == id)), use.names = FALSE)
  boot_pair <- rep(sprintf("B%03d", seq_along(sampled)), each = 2L)
  bx <- x[rows, , drop = FALSE]
  by_factor <- factor(as.character(y_factor[rows]), levels = levels(y_factor))
  by_binary <- as.integer(by_factor == "Tumor")
  boot_fold <- patient_fold_vector(boot_pair, 10L, 20300000L + b)
  fit_lasso <- lasso_select(bx, by_binary, boot_fold)
  fit_svm <- svm_rfe_select(bx, by_factor, boot_fold)
  bootstrap_membership <- rbind(
    bootstrap_membership,
    data.frame(
      bootstrap = b,
      gene_symbol = candidates,
      LASSO_lambda_min = candidates %in% fit_lasso$selected_min,
      LASSO_lambda_1se = candidates %in% fit_lasso$selected_1se,
      SVM_RFE = candidates %in% fit_svm$selected,
      intersection_lambda_min = candidates %in%
        intersect(fit_lasso$selected_min, fit_svm$selected),
      intersection_lambda_1se = candidates %in%
        intersect(fit_lasso$selected_1se, fit_svm$selected),
      stringsAsFactors = FALSE
    )
  )
  bootstrap_diagnostics <- rbind(
    bootstrap_diagnostics,
    data.frame(
      bootstrap = b,
      unique_original_patients = length(unique(sampled)),
      lambda_min = fit_lasso$lambda_min,
      lambda_1se = fit_lasso$lambda_1se,
      n_LASSO_lambda_min = length(fit_lasso$selected_min),
      n_LASSO_lambda_1se = length(fit_lasso$selected_1se),
      n_SVM_RFE = length(fit_svm$selected),
      n_intersection_lambda_min = length(intersect(
        fit_lasso$selected_min, fit_svm$selected
      )),
      n_intersection_lambda_1se = length(intersect(
        fit_lasso$selected_1se, fit_svm$selected
      )),
      SVM_RFE_max_accuracy = fit_svm$max_accuracy,
      stringsAsFactors = FALSE
    )
  )
  if (b %% 20L == 0L) {
    message("Completed patient-pair bootstrap ", b, "/", bootstrap_reps)
  }
}
bootstrap_frequency <- aggregate(
  cbind(
    LASSO_lambda_min, LASSO_lambda_1se, SVM_RFE,
    intersection_lambda_min, intersection_lambda_1se
  ) ~ gene_symbol,
  bootstrap_membership,
  mean
)
names(bootstrap_frequency)[-1] <- paste0(
  names(bootstrap_frequency)[-1], "_frequency"
)
write_tsv(
  bootstrap_membership,
  file.path(out_dir, "V2_200_patient_bootstrap_selection_membership.tsv")
)
write_tsv(
  bootstrap_frequency,
  file.path(out_dir, "V2_200_patient_bootstrap_selection_frequency.tsv")
)
write_tsv(
  bootstrap_diagnostics,
  file.path(out_dir, "V2_200_patient_bootstrap_diagnostics.tsv")
)

# The internal comparison with an earlier project version was omitted from
# this public release. The reported patient-grouped stability analyses above
# are unchanged.

audit_manifest <- data.frame(
  item = c(
    "bootstrap_repetitions", "bootstrap_unit", "candidate_set",
    "fixed_fold_unit", "fixed_fold_seed", "bootstrap_seed_base",
    "LASSO_rule_min", "LASSO_rule_1se", "SVM_RFE_rule",
    "correlation_rule"
  ),
  value = c(
    bootstrap_reps,
    "GSE44076 patient pair sampled with replacement",
    paste(candidates, collapse = ";"),
    "patient; tumor-normal pair always in same fold",
    20260719L,
    20300000L,
    "nonzero coefficient at lambda.min",
    "nonzero coefficient at lambda.1se",
    "linear SVM cost=1; smallest subset tied at maximum grouped-CV accuracy",
    "Spearman across all 196 GSE44076 samples"
  ),
  stringsAsFactors = FALSE
)
write_tsv(audit_manifest, file.path(out_dir, "audit_method_manifest.tsv"))

message("V2 ML stability audit completed with ", bootstrap_reps, " bootstraps.")
