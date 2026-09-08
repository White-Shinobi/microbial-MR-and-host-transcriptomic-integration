#!/usr/bin/env Rscript

# Preserve the validated V1 ML/plotting workflow while preventing paired
# GSE44076 samples from leaking across cross-validation folds.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
v1_script <- file.path(root, "legacy_sources", "run_ref1_ml_expression_validation.R")
code <- readLines(v1_script, warn = FALSE)

root_line <- grep("^root <- normalizePath", code)
stopifnot(length(root_line) == 1L)
code[root_line] <- sprintf("root <- %s", deparse(root))

fold_start <- grep("^# One reproducible stratified 10-fold", code)
fold_end <- grep('write_tsv\\(fold_table, file.path\\(ml_dir, "GSE44076_stratified_10fold_assignment.tsv"\\)\\)', code)
stopifnot(length(fold_start) == 1L, length(fold_end) == 1L, fold_end > fold_start)
grouped_fold_code <- c(
  "# Fixed patient-grouped 10-fold partition shared by LASSO and SVM-RFE.",
  "# Each tumor/normal pair is kept in the same fold.",
  "patient_ids <- unique(disc_manifest$individual_id)",
  "set.seed(20260719)",
  "patient_fold <- sample(rep(seq_len(10), length.out = length(patient_ids)))",
  "names(patient_fold) <- patient_ids",
  "foldid <- unname(patient_fold[disc_manifest$individual_id])",
  "stopifnot(!anyNA(foldid), all(tapply(foldid, disc_manifest$individual_id,",
  "                                    function(z) length(unique(z)) == 1L)))",
  "fold_table <- data.frame(",
  "  sample_id = rownames(x), individual_id = disc_manifest$individual_id,",
  "  group = y, fold = foldid, stringsAsFactors = FALSE",
  ")",
  'write_tsv(fold_table, file.path(ml_dir, "GSE44076_fixed_patient_grouped_10fold_assignment.tsv"))'
)
code <- c(
  code[seq_len(fold_start - 1L)],
  grouped_fold_code,
  code[(fold_end + 1L):length(code)]
)
code <- gsub(
  "E  Ten-fold LASSO cross-validation",
  "E  LASSO cross-validation (10 folds; patient-level grouping)",
  code,
  fixed = TRUE
)
code <- gsub(
  "F  Ten-fold SVM-RFE",
  "F  SVM-RFE (10-fold CV; patient-level grouping)",
  code,
  fixed = TRUE
)
code <- gsub(
  "One fixed stratified ten-fold partition",
  "One fixed patient-grouped ten-fold partition",
  code,
  fixed = TRUE
)

eval(parse(text = code, keep.source = TRUE), envir = .GlobalEnv)
