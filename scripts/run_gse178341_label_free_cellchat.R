#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(Matrix)
  library(data.table)
  library(reshape2)
})

options(error = function() { traceback(20); q(status = 1) })
set.seed(20260723)
args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep("^--file=", args_all)][1]
root <- normalizePath(
  file.path(dirname(sub("^--file=", "", script_arg)), ".."),
  mustWork = TRUE
)
out_dir <- file.path(
  root, "results", "single_cell", "GSE178341",
  "label_free_figure7"
)

genes <- fread(
  file.path(out_dir, "label_free_cellchat_gene_symbols.tsv")
)$gene_symbol
meta <- fread(
  file.path(out_dir, "label_free_cellchat_cell_metadata.tsv")
)
counts <- readMM(
  gzfile(
    file.path(
      out_dir,
      "label_free_cellchat_counts_genes_by_cells.mtx.gz"
    )
  )
)
counts <- as(counts, "dgCMatrix")
stopifnot(nrow(counts) == length(genes), ncol(counts) == nrow(meta))
rownames(counts) <- genes
colnames(counts) <- meta$cellID
rownames(meta) <- meta$cellID
obj <- CreateSeuratObject(counts = counts, meta.data = as.data.frame(meta))
obj <- NormalizeData(
  obj, normalization.method = "LogNormalize",
  scale.factor = 10000, verbose = FALSE
)

is_epithelial <- function(value) {
  startsWith(as.character(value), "E")
}

extract_all <- function(cc, tissue) {
  net <- as.data.table(
    reshape2::melt(cc@net$prob, value.name = "prob")
  )
  setnames(net, 1:3, c("source", "target", "interaction_name"))
  pdt <- as.data.table(
    reshape2::melt(cc@net$pval, value.name = "pval")
  )
  stopifnot(nrow(net) == nrow(pdt))
  net[, pval := pdt$pval]
  net <- net[is.finite(prob) & prob > 0]
  lr <- as.data.table(
    cc@LR$LRsig, keep.rownames = "interaction_name"
  )
  keep <- intersect(
    c(
      "interaction_name", "interaction_name_2", "pathway_name",
      "ligand", "receptor", "annotation", "evidence"
    ),
    names(lr)
  )
  net <- merge(
    net, lr[, ..keep], by = "interaction_name", all.x = TRUE
  )
  net[, `:=`(
    tissue = tissue,
    p_fdr_bh = p.adjust(pval, method = "BH")
  )]
  net[, fdr_significant := p_fdr_bh < 0.05]
  net[, involves_epithelial :=
    is_epithelial(source) | is_epithelial(target)]
  setorder(net, p_fdr_bh, -prob)
  net
}

run_tissue <- function(code) {
  cells <- meta[SPECIMEN_TYPE == code, cellID]
  data_input <- GetAssayData(
    obj, assay = "RNA", layer = "data"
  )[, cells, drop = FALSE]
  meta_input <- as.data.frame(meta[match(cells, cellID)])
  rownames(meta_input) <- meta_input$cellID
  cc <- createCellChat(
    object = data_input, meta = meta_input,
    group.by = "cellchat_identity"
  )
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(
    cc, type = "triMean", raw.use = TRUE,
    population.size = TRUE, nboot = 200,
    seed.use = 20260723
  )
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  list(object = cc, interactions = extract_all(cc, code))
}

normal_checkpoint <- file.path(
  out_dir, "label_free_cellchat_normal_checkpoint.rds"
)
tumor_checkpoint <- file.path(
  out_dir, "label_free_cellchat_crc_checkpoint.rds"
)
if (file.exists(normal_checkpoint)) {
  normal <- readRDS(normal_checkpoint)
} else {
  normal <- run_tissue("N")
  saveRDS(normal, normal_checkpoint, compress = FALSE)
}
if (file.exists(tumor_checkpoint)) {
  tumor <- readRDS(tumor_checkpoint)
} else {
  tumor <- run_tissue("T")
  saveRDS(tumor, tumor_checkpoint, compress = FALSE)
}

all_interactions <- rbindlist(
  list(normal$interactions, tumor$interactions), fill = TRUE
)
significant <- all_interactions[fdr_significant == TRUE]
fwrite(
  all_interactions,
  file.path(
    out_dir,
    "label_free_cellchat_all_positive_interactions.tsv.gz"
  ),
  sep = "\t"
)
fwrite(
  significant,
  file.path(
    out_dir,
    "label_free_cellchat_BH_FDR_significant_interactions.tsv.gz"
  ),
  sep = "\t"
)
audit <- meta[, .(
  sampled_cells = .N, unique_patients = uniqueN(PID)
), by = .(SPECIMEN_TYPE, cellchat_identity)]
audit[, evaluable_min_10_cells := sampled_cells >= 10]
fwrite(
  audit,
  file.path(out_dir, "label_free_cellchat_model_audit.tsv"),
  sep = "\t"
)
saveRDS(
  list(Normal = normal$object, CRC = tumor$object),
  file.path(out_dir, "label_free_GSE178341_CellChat_nboot200.rds"),
  compress = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "label_free_cellchat_sessionInfo.txt")
)
cat(
  "Normal significant:", normal$interactions[, sum(fdr_significant)],
  "CRC significant:", tumor$interactions[, sum(fdr_significant)], "\n"
)
