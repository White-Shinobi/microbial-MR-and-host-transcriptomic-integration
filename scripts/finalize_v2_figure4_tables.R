#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
tab <- file.path(root, "manuscript/tables/supplementary")
dir.create(tab, recursive = TRUE, showWarnings = FALSE)

# The reusable V1 DEG script also prepares an unused GSE39582 matrix. The V2
# expression gate is GSE87211, so make the stage summary reflect the dataset
# actually used for biomarker validation.
bulk_summary_path <- file.path(
  root, "results/bulk/crc/intersection/CRC_bulk_DEG_GMRG_summary.tsv"
)
bulk_summary <- fread(bulk_summary_path)
bulk_summary[metric == "validation_dataset", value := "GSE87211"]
bulk_summary[metric == "validation_tumor_samples", value := "203"]
bulk_summary[metric == "validation_normal_samples", value := "160"]
fwrite(bulk_summary, bulk_summary_path, sep = "\t")

copy_table <- function(source, name) {
  stopifnot(file.exists(source))
  ok <- file.copy(source, file.path(tab, name), overwrite = TRUE)
  stopifnot(ok)
}

copy_table(
  file.path(root, "results/bulk/crc/external_validation/two_cohort_fig3",
            "two_cohort_external_validation_performance.tsv"),
  "Table_S14_external_validation_cohort_audit.tsv"
)

outer <- fread(file.path(
  root, "results/bulk/crc/machine_learning_stability_audit",
  "V2_patient_grouped_outer_fold_selection_frequency.tsv"
))
bootstrap <- fread(file.path(
  root, "results/bulk/crc/machine_learning_stability_audit",
  "V2_200_patient_bootstrap_selection_frequency.tsv"
))
setnames(outer, setdiff(names(outer), "gene_symbol"),
         paste0("outer10_", setdiff(names(outer), "gene_symbol")))
setnames(bootstrap, setdiff(names(bootstrap), "gene_symbol"),
         paste0("bootstrap200_", setdiff(names(bootstrap), "gene_symbol")))
stability <- merge(outer, bootstrap, by = "gene_symbol", all = TRUE)
fwrite(stability, file.path(tab, "Table_S15_patient_grouped_ML_stability.tsv"),
       sep = "\t")

eqtl_audit <- fread(file.path(
  root, "results/coloc/v2_biomarker_loci",
  "crc_biomarker_coloc_default_prior.tsv"
))
eqtl_audit <- eqtl_audit[comparison != "microbial exposure vs CRC"]
fwrite(
  eqtl_audit,
  file.path(tab, "Table_S40_targeted_biomarker_colocalization.tsv"),
  sep = "\t"
)
