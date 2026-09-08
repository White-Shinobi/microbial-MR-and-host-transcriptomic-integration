#!/usr/bin/env Rscript

# Run the read-only V1 bulk-DEG implementation against the V2 project.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
v1_script <- file.path(root, "legacy_sources", "run_crc_bulk_deg_ref1.R")
code <- readLines(v1_script, warn = FALSE)
root_line <- grep("^root <- normalizePath", code)
stopifnot(length(root_line) == 1L)
code[root_line] <- sprintf("root <- %s", deparse(root))

# Replace ref1's nominal-P DEG criterion with gene-level BH-FDR < 0.05 while
# retaining its absolute log2-fold-change threshold. BH correction is applied
# after collapsing the probe results to the 19,468 tested gene symbols.
direction_start <- grep("^gene_results\\$direction <-", code)
stopifnot(length(direction_start) == 1L)
direction_end <- direction_start + 1L
stopifnot(grepl("Not significant", code[direction_end], fixed = TRUE))
code <- c(
  code[seq_len(direction_start - 1L)],
  'gene_results$BH_FDR <- p.adjust(gene_results$P.Value, method = "BH")',
  'gene_results$direction <- ifelse(gene_results$logFC > 0.5 & gene_results$BH_FDR < 0.05, "Up",',
  '                                 ifelse(gene_results$logFC < -0.5 & gene_results$BH_FDR < 0.05, "Down", "Not significant"))',
  code[(direction_end + 1L):length(code)]
)

candidate_columns <- grep(
  'candidate <- merge\\(gmrgs, degs\\[, c\\("gene_symbol", "probe_id", "logFC", "P.Value", "adj.P.Val", "direction"\\)\\],',
  code
)
stopifnot(length(candidate_columns) == 1L)
code[candidate_columns] <- sub(
  '"P.Value", "adj.P.Val", "direction"',
  '"P.Value", "adj.P.Val", "BH_FDR", "direction"',
  code[candidate_columns],
  fixed = TRUE
)

fdr_axis <- grep("^plot_data\\$minus_log10_p <-", code)
fdr_aes <- grep("^p_volcano <- ggplot\\(plot_data, aes\\(logFC, minus_log10_p", code)
fdr_label <- grep('y = expression\\(-log\\[10\\]~italic\\(P\\)\\)', code)
stopifnot(length(fdr_axis) == 1L, length(fdr_aes) == 1L, length(fdr_label) == 1L)
code[fdr_axis] <- 'plot_data$minus_log10_fdr <- -log10(pmax(plot_data$BH_FDR, .Machine$double.xmin))'
code[fdr_aes] <- 'p_volcano <- ggplot(plot_data, aes(logFC, minus_log10_fdr, color = direction)) +'
code[fdr_label] <- sub(
  'y = expression(-log[10]~italic(P))',
  'y = expression(-log[10]~"BH-FDR")',
  code[fdr_label],
  fixed = TRUE
)

eval(parse(text = code, keep.source = TRUE), envir = .GlobalEnv)
