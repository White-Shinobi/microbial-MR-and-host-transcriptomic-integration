#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Biobase))

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
script_path <- sub("^--file=", "", script_arg)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

processed_dir <- file.path(root, "data", "processed", "bulk", "crc", "GSE87211")
result_dir <- file.path(root, "results", "bulk", "crc", "validation")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

eset <- readRDS(file.path(processed_dir, "GSE87211_ExpressionSet.rds"))
annotation <- read.delim(file.path(processed_dir, "GPL13497_platform_annotation.tsv"),
                         check.names = FALSE, stringsAsFactors = FALSE, quote = "\"")

expr <- exprs(eset)
map_index <- match(rownames(expr), annotation$ID)
symbols <- trimws(annotation$GENE_SYMBOL[map_index])
symbols <- sub("[;,/].*$", "", symbols)
symbols[symbols %in% c("", "NA", "---")] <- NA_character_

keep <- !is.na(symbols)
expr_mapped <- expr[keep, , drop = FALSE]
symbols <- symbols[keep]
probe_ids <- rownames(expr_mapped)
probe_iqr <- apply(expr_mapped, 1, IQR, na.rm = TRUE)

ord <- order(symbols, -probe_iqr, probe_ids)
selected <- ord[!duplicated(symbols[ord])]
gene_expr <- expr_mapped[selected, , drop = FALSE]
rownames(gene_expr) <- symbols[selected]

pd <- pData(eset)
group <- ifelse(grepl("rectal tumor", pd$source_name_ch1, ignore.case = TRUE), "Tumor",
                ifelse(grepl("mucosa", pd$source_name_ch1, ignore.case = TRUE), "Normal", NA_character_))
patient_raw <- if ("patient:ch1" %in% names(pd)) pd[["patient:ch1"]] else pd$characteristics_ch1
patient_id <- sub("^patient: ", "", patient_raw, ignore.case = TRUE)
patient_id <- sub("^(Con|Tum) Patient ", "Patient ", patient_id, ignore.case = TRUE)

manifest <- data.frame(
  sample_id = colnames(gene_expr),
  title = pd$title,
  source = pd$source_name_ch1,
  group = group,
  patient_id = patient_id,
  stringsAsFactors = FALSE
)

probe_map <- data.frame(
  gene_symbol = rownames(gene_expr),
  selected_probe = probe_ids[selected],
  probe_IQR = probe_iqr[selected],
  stringsAsFactors = FALSE
)

saveRDS(gene_expr, file.path(result_dir, "GSE87211_gene_expression_matrix.rds"))
write.table(manifest, file.path(result_dir, "GSE87211_sample_manifest.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)
write.table(probe_map, file.path(result_dir, "GSE87211_gene_probe_map.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)

candidates <- c("R3HDM1", "FAM135B", "MCM6", "LPP", "DOCK5", "PELI2", "ZRANB3", "NAV3", "WFDC2")
coverage <- data.frame(
  gene_symbol = candidates,
  present_in_gene_matrix = candidates %in% rownames(gene_expr),
  selected_probe = probe_map$selected_probe[match(candidates, probe_map$gene_symbol)],
  stringsAsFactors = FALSE
)
write.table(coverage, file.path(result_dir, "GSE87211_candidate_gene_matrix_coverage.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE, na = "")

summary <- data.frame(
  metric = c("samples_total", "tumor_samples", "normal_samples", "genes_in_matrix",
             "candidate_genes_present", "complete_matched_patient_pairs"),
  value = c(ncol(gene_expr), sum(group == "Tumor"), sum(group == "Normal"), nrow(gene_expr),
            sum(coverage$present_in_gene_matrix),
            length(intersect(manifest$patient_id[manifest$group == "Tumor"],
                             manifest$patient_id[manifest$group == "Normal"]))),
  stringsAsFactors = FALSE
)
write.table(summary, file.path(result_dir, "GSE87211_preparation_summary.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)

cat("Prepared GSE87211:", sum(group == "Tumor"), "tumors,", sum(group == "Normal"),
    "normals,", sum(coverage$present_in_gene_matrix), "of 9 candidates present.\n")
