#!/usr/bin/env Rscript

# Execute the read-only V1 three-layer coloc implementation using a dynamic V2
# biomarker-locus manifest and the already extracted revised-MR regional files.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
v1_script <- file.path(root, "legacy_sources", "run_crc_biomarker_coloc.R")
code <- readLines(v1_script, warn = FALSE)

root_line <- grep("^root <- normalizePath", code)
stopifnot(length(root_line) == 1L)
code[root_line] <- sprintf("root <- %s", deparse(root))
code <- gsub("crc_biomarker_loci", "v2_biomarker_loci", code, fixed = TRUE)

loci_start <- grep("^loci <- data.table\\(", code)
loci_end <- grep("^loci\\[, `:=`\\(start =", code)
stopifnot(length(loci_start) == 1L, length(loci_end) == 1L, loci_end > loci_start)
code <- c(
  code[seq_len(loci_start - 1L)],
  'loci <- fread(file.path(region_dir, "biomarker_locus_manifest.tsv"))',
  code[(loci_end + 1L):length(code)]
)

type_line <- grep('^  exposure_type <- if \\(locus\\$accession ==', code)
stopifnot(length(type_line) == 1L)
code[type_line] <- '  exposure_type <- if (tolower(locus$model) == "logistic") "cc" else "quant"'

eval(parse(text = code, keep.source = TRUE), envir = .GlobalEnv)
