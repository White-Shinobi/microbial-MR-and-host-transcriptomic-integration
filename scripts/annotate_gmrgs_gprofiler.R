#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(httr)
  library(jsonlite)
})

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

file_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
script_path <- sub("^--file=", "", file_arg[[1]])
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

mr_dir <- Sys.getenv("GMRG_MR_DIR", unset = file.path(project_dir, "results", "mr", "crc", "strict_ref1"))
out_dir <- Sys.getenv("GMRG_OUT_DIR", unset = file.path(project_dir, "results", "gmrg"))
processed_dir <- Sys.getenv("GMRG_PROCESSED_DIR", unset = file.path(project_dir, "data", "processed", "gmrg"))
input_label <- Sys.getenv(
  "GMRG_INPUT_LABEL",
  unset = "Harmonized instruments from strict-ref1 IVW-selected microbial features"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

mr <- fread(file.path(mr_dir, "strict_ref1_mr_results.tsv"))
instruments <- fread(file.path(mr_dir, "strict_ref1_harmonized_instruments.tsv"))

candidates <- unique(mr[primary_ivw == TRUE & ref1_candidate == TRUE,
                        .(accession, family, feature_id, feature_name)])
candidate_instruments <- merge(
  instruments,
  candidates,
  by = c("accession", "family", "feature_id")
)

rsids <- sort(unique(candidate_instruments$SNP))
stopifnot(length(rsids) > 0L)

endpoint <- "https://biit.cs.ut.ee/gprofiler/api/snpense/snpense"
response <- POST(
  endpoint,
  add_headers(`Content-Type` = "application/json"),
  body = list(query = rsids, output = "json"),
  encode = "json",
  config(connecttimeout = 60),
  timeout(120)
)
stop_for_status(response)
response_text <- content(response, as = "text", encoding = "UTF-8")
writeLines(response_text, file.path(processed_dir, "gprofiler_snpense_response.json"), useBytes = TRUE)

parsed <- fromJSON(response_text, simplifyVector = FALSE)
records <- parsed$result

annotation_rows <- rbindlist(lapply(records, function(record) {
  genes <- unlist(record$gene_names, use.names = FALSE)
  ensgs <- unlist(record$ensgs, use.names = FALSE)
  if (!length(genes) || !length(ensgs)) {
    return(data.table(
      SNP = record$rs_id,
      chromosome = as.character(record$chromosome %||% NA_character_),
      start = as.integer(record$start %||% NA_integer_),
      end = as.integer(record$end %||% NA_integer_),
      strand = as.character(record$strand %||% NA_character_),
      ensembl_gene_id = NA_character_,
      gene_symbol = NA_character_,
      consequence = paste(names(record$variants), collapse = ";")
    ))
  }
  n <- max(length(genes), length(ensgs))
  data.table(
    SNP = record$rs_id,
    chromosome = as.character(record$chromosome %||% NA_character_),
    start = as.integer(record$start %||% NA_integer_),
    end = as.integer(record$end %||% NA_integer_),
    strand = as.character(record$strand %||% NA_character_),
    ensembl_gene_id = rep(ensgs, length.out = n),
    gene_symbol = rep(genes, length.out = n),
    consequence = paste(names(record$variants), collapse = ";")
  )
}), fill = TRUE)

missing_rsids <- setdiff(rsids, annotation_rows$SNP)
if (length(missing_rsids)) {
  annotation_rows <- rbind(
    annotation_rows,
    data.table(
      SNP = missing_rsids, chromosome = NA_character_, start = NA_integer_, end = NA_integer_,
      strand = NA_character_, ensembl_gene_id = NA_character_, gene_symbol = NA_character_,
      consequence = NA_character_
    ),
    fill = TRUE
  )
}
setorder(annotation_rows, chromosome, start, SNP, gene_symbol, na.last = TRUE)

feature_snp_gene <- merge(
  unique(candidate_instruments[, .(accession, family, feature_id, feature_name, SNP)]),
  annotation_rows,
  by = "SNP",
  all.x = TRUE,
  allow.cartesian = TRUE
)
setorder(feature_snp_gene, family, feature_name, SNP, gene_symbol)

gmrg <- unique(annotation_rows[!is.na(gene_symbol) & gene_symbol != "",
                               .(ensembl_gene_id, gene_symbol)])
setorder(gmrg, gene_symbol, ensembl_gene_id)

summary_table <- data.table(
  metric = c(
    "mr_selected_features", "candidate_instrument_rows", "unique_candidate_snps",
    "snps_with_gene_annotation", "snps_without_gene_annotation", "unique_gmrgs"
  ),
  value = c(
    nrow(candidates), nrow(candidate_instruments), length(rsids),
    uniqueN(annotation_rows[!is.na(gene_symbol), SNP]),
    length(setdiff(rsids, annotation_rows[!is.na(gene_symbol), unique(SNP)])),
    nrow(gmrg)
  )
)

fwrite(annotation_rows, file.path(out_dir, "ref1_gprofiler_snp_to_gene.tsv"), sep = "\t", na = "NA")
fwrite(feature_snp_gene, file.path(out_dir, "ref1_candidate_feature_snp_gene.tsv"), sep = "\t", na = "NA")
fwrite(gmrg, file.path(out_dir, "ref1_gmrg_unique_genes.tsv"), sep = "\t", na = "NA")
fwrite(summary_table, file.path(out_dir, "ref1_gmrg_summary.tsv"), sep = "\t")

metadata <- data.table(
  field = c("timestamp", "service", "endpoint", "input", "mapping_rule", "R_version"),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    "g:Profiler g:SNPense",
    endpoint,
    input_label,
    "All Ensembl gene IDs and gene symbols returned by g:SNPense; merge and deduplicate genes",
    R.version.string
  )
)
fwrite(metadata, file.path(out_dir, "ref1_gmrg_metadata.tsv"), sep = "\t")

print(summary_table)
