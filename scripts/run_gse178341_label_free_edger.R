#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep("^--file=", args_all)][1]
root <- normalizePath(
  file.path(dirname(sub("^--file=", "", script_arg)), ".."),
  mustWork = TRUE
)
atlas <- file.path(
  root, "results", "single_cell", "GSE178341", "label_free_atlas"
)
pb_dir <- file.path(atlas, "pseudobulk")
input_dir <- file.path(pb_dir, "edger_input")

biomarkers <- c("GTF2IRD1", "NRXN1", "FAM135B", "KIAA1671", "WFDC2")
cell_type_slugs <- c(
  Epithelial = "Epithelial",
  `T/NK/ILC` = "T_NK_ILC",
  Myeloid = "Myeloid",
  B = "B",
  Plasma = "Plasma",
  Stromal = "Stromal",
  Mast = "Mast"
)
results <- list()
model_audit <- list()

for (cell_type in names(cell_type_slugs)) {
  slug <- unname(cell_type_slugs[[cell_type]])
  counts_dt <- fread(
    cmd = paste(
      "gzip -dc",
      shQuote(file.path(input_dir, paste0(slug, "_counts.tsv.gz")))
    )
  )
  sample_dt <- fread(
    file.path(input_dir, paste0(slug, "_samples.tsv"))
  )
  sample_names <- sample_dt$sample_name
  stopifnot(identical(names(counts_dt)[-(1:3)], sample_names))
  count_matrix <- as.matrix(counts_dt[, ..sample_names])
  storage.mode(count_matrix) <- "integer"
  rownames(count_matrix) <- paste0(
    counts_dt$ensembl_id, "__", counts_dt$feature_index
  )

  patient <- factor(sample_dt$PID)
  tissue <- factor(sample_dt$SPECIMEN_TYPE, levels = c("N", "T"))
  design <- model.matrix(~ patient + tissue)
  tissue_coefficient <- which(colnames(design) == "tissueT")
  stopifnot(
    length(tissue_coefficient) == 1L,
    qr(design)$rank == ncol(design)
  )

  y_all <- DGEList(
    counts = count_matrix,
    genes = data.frame(
      feature_index = counts_dt$feature_index,
      ensembl_id = counts_dt$ensembl_id,
      gene_symbol = counts_dt$gene_symbol,
      stringsAsFactors = FALSE
    )
  )
  keep_standard <- filterByExpr(y_all, design = design)
  keep <- keep_standard | y_all$genes$gene_symbol %in% biomarkers
  y <- y_all[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design = design, robust = TRUE)
  fit <- glmQLFit(y, design = design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = tissue_coefficient)
  table_all <- topTags(qlf, n = Inf, sort.by = "none")$table
  table_all$feature_key <- rownames(table_all)
  table_all$gene_symbol <- y$genes$gene_symbol[
    match(table_all$feature_key, rownames(y$counts))
  ]

  for (gene in biomarkers) {
    hit <- table_all[
      table_all$gene_symbol == gene, , drop = FALSE
    ]
    if (nrow(hit) != 1L) {
      stop(
        paste(
          cell_type, gene,
          "expected once after filtering; found", nrow(hit)
        )
      )
    }
    target_key <- hit$feature_key
    results[[length(results) + 1L]] <- data.table(
      gene = gene,
      cell_type = cell_type,
      n_paired_patients = uniqueN(sample_dt$PID),
      edgeR_log2_fold_change_CRC_vs_normal = hit$logFC,
      edgeR_logCPM = hit$logCPM,
      edgeR_QLF = hit$F,
      edgeR_QL_P = hit$PValue,
      target_gene_passed_filterByExpr = keep_standard[
        match(target_key, rownames(y_all$counts))
      ]
    )
  }
  model_audit[[length(model_audit) + 1L]] <- data.table(
    cell_type = cell_type,
    samples = nrow(sample_dt),
    paired_patients = uniqueN(sample_dt$PID),
    design_columns = ncol(design),
    design_rank = qr(design)$rank,
    genes_before_filter = nrow(count_matrix),
    genes_after_filter_plus_targets = nrow(y),
    common_dispersion = y$common.dispersion
  )
}

result <- rbindlist(results)
stopifnot(
  nrow(result) == 35L,
  !anyDuplicated(result[, .(gene, cell_type)])
)
result[is.na(edgeR_QL_P), edgeR_QL_P := 1]
result[, edgeR_BH_FDR_35_tests := p.adjust(edgeR_QL_P, method = "BH")]

wilcox <- fread(
  file.path(
    pb_dir,
    "label_free_paired_wilcoxon_biomarker_results_35_tests.tsv"
  )
)
setnames(
  wilcox,
  "n_paired_patients",
  "wilcoxon_n_paired_patients"
)
result <- merge(
  result,
  wilcox,
  by = c("gene", "cell_type"),
  all.x = TRUE,
  sort = FALSE
)
stopifnot(
  result$n_paired_patients == result$wilcoxon_n_paired_patients
)
result[, edgeR_FDR_lt_0_05 := edgeR_BH_FDR_35_tests < 0.05]
result[, wilcoxon_FDR_lt_0_05 :=
  paired_wilcoxon_BH_FDR_35_tests < 0.05]

fwrite(
  result,
  file.path(pb_dir, "label_free_paired_edger_biomarker_results_35_tests.tsv"),
  sep = "\t"
)
fwrite(
  rbindlist(model_audit),
  file.path(pb_dir, "label_free_paired_edger_model_audit.tsv"),
  sep = "\t"
)
writeLines(
  capture.output(sessionInfo()),
  file.path(pb_dir, "label_free_paired_edger_sessionInfo.txt")
)
cat(
  "Completed label-free paired edgeR analysis for",
  nrow(result), "gene-cell-type tests\n"
)
