#!/usr/bin/env Rscript

# Execute a validated V1 implementation without changing V1, replacing only a
# hard-coded gene vector with the corresponding V2-selected vector.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1]] %in% c("prepare87211", "gsea", "immune")) {
  stop("Usage: run_v2_dynamic_v1_stage.R <prepare87211|gsea|immune>")
}
stage <- args[[1]]
all_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", all_args[grep("^--file=", all_args)][1])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
if (stage == "prepare87211") {
  candidate_path <- file.path(
    root, "results/bulk/crc/intersection/CRC_DEG_GMRG_candidate_genes.tsv"
  )
  if (!file.exists(candidate_path)) stop("Missing V2 candidate file: ", candidate_path)
  candidate_table <- read.delim(candidate_path, stringsAsFactors = FALSE, check.names = FALSE)
  biomarkers <- unique(candidate_table$gene_symbol)
} else {
  gate_path <- file.path(root, "results/bulk/crc/machine_learning/biomarker_gate_summary.tsv")
  if (!file.exists(gate_path)) stop("Missing V2 biomarker gate: ", gate_path)
  gate <- read.delim(gate_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("gene_symbol", "final_biomarker") %in% names(gate))) {
    stop("Unexpected biomarker gate columns")
  }
  flag <- tolower(as.character(gate$final_biomarker)) %in% c("true", "t", "1", "yes")
  biomarkers <- unique(gate$gene_symbol[flag & !is.na(gate$gene_symbol) & gate$gene_symbol != ""])
}
if (!length(biomarkers)) stop("No V2 gene available for stage ", stage)

legacy_scripts <- file.path(root, "legacy_sources")
source_path <- file.path(
  legacy_scripts,
  switch(
    stage,
    prepare87211 = "prepare_gse87211_validation.R",
    gsea = "run_ref1_biomarker_gsea.R",
    immune = "run_ref1_immune_infiltration.R"
  )
)
code <- readLines(source_path, warn = FALSE)
vector_name <- switch(stage, prepare87211 = "candidates", ridge = "genes", "biomarkers")
replacement <- paste0(vector_name, " <- c(",
                      paste(encodeString(biomarkers, quote = '"'), collapse = ", "), ")")
bio_line <- grep(paste0("^", vector_name, " <- c\\("), code)
if (length(bio_line) != 1L) stop("Could not uniquely locate V1 gene declaration")
code[bio_line] <- replacement
if (stage == "prepare87211") {
  code <- sub(
    'sum\\(coverage\\$present_in_gene_matrix\\), "of 9 candidates present\\.\\\\n"\\)',
    'sum(coverage$present_in_gene_matrix), "of", length(candidates), "candidates present.\\\\n")',
    code
  )
}
if (stage == "gsea") {
  figure_line <- grep("^figure_biomarkers <- c\\(", code)
  if (length(figure_line) != 1L) stop("Could not uniquely locate V1 figure biomarker declaration")
  code[figure_line] <- sub("^biomarkers", "figure_biomarkers", replacement)
  exclusion_line <- grep("ZRANB3 excluded from Figure 4 only", code)
  if (length(exclusion_line) == 1L) {
    code[exclusion_line] <- '            "No V2 biomarker was excluded from the Figure 4 display",'
  }
  # V1 used a four-gene, hard-coded chromosome ideogram. Rebuild this block
  # from the current V2 SNP-gene mapping and independently verified GRCh38 gene
  # intervals, and overlay the mapped instrument SNPs rather than reusing V1
  # loci.
  block_start <- grep("^  gene_loci <- data.frame\\(", code)
  block_end <- grep(
    '^  ggsave\\(file.path\\(figure_dir, "Figure_4_biomarker_KEGG_GSEA.pdf"\\),',
    code
  )
  if (length(block_start) != 1L || length(block_end) != 1L || block_end <= block_start) {
    stop("Could not locate V1 Figure 4 hard-coded locus block")
  }
  dynamic_layout <- c(
    '  mapping_path <- file.path(root, "results", "gmrg", "ref1_candidate_feature_snp_gene.tsv")',
    '  mr_path <- file.path(root, "manuscript", "tables", "main", "Table_1_primary_IVW_associations.tsv")',
    '  stopifnot(file.exists(mapping_path), file.exists(mr_path))',
    '  mapping <- data.table::fread(mapping_path)',
    '  mr_primary <- data.table::fread(mr_path)',
    '  mapped <- unique(mapping[gene_symbol %in% figure_biomarkers, .(',
    '    gene_symbol, SNP, accession, family, feature_name,',
    '    snp_chromosome = as.character(chromosome), snp_start = start, snp_end = end',
    '  )])',
    '  mapped <- merge(mapped, mr_primary[, .(accession, OR, q_value)], by = "accession", all.x = TRUE)',
    '  if (nrow(mapped) != length(figure_biomarkers) || anyDuplicated(mapped$gene_symbol) ||',
    '      !setequal(mapped$gene_symbol, figure_biomarkers) || any(!is.finite(mapped$OR))) {',
    '    stop("Figure 4F requires exactly one current instrument mapping per V2 biomarker")',
    '  }',
    '  gene_intervals <- data.table::data.table(',
    '    gene_symbol = c("GTF2IRD1", "NRXN1", "FAM135B", "KIAA1671", "WFDC2"),',
    '    ensembl_gene_id = c("ENSG00000006704", "ENSG00000179915", "ENSG00000147724", "ENSG00000197077", "ENSG00000101443"),',
    '    Chromosome = c("chr7", "chr2", "chr8", "chr22", "chr20"),',
    '    chromStart = c(74452463, 49918503, 138130023, 24952643, 45469166),',
    '    chromEnd = c(74602619, 51225575, 138497267, 25197448, 45481533),',
    '    Cytoband = c("7q11.23", "2p16.3", "8q24.23", "22q11.23", "20q13.12")',
    '  )',
    '  locus_source <- merge(gene_intervals, mapped, by = "gene_symbol", all.x = TRUE)',
    '  locus_source[, snp_Chromosome := paste0("chr", snp_chromosome)]',
    '  locus_source[, mr_direction := ifelse(OR > 1, "Higher CRC susceptibility", "Lower CRC susceptibility")]',
    '  locus_source[, exposure_family := ifelse(family == "module", "Metabolic module", "Species presence/absence")]',
    '  chrom_order <- setNames(seq_along(c(1:22, "X", "Y")), paste0("chr", c(1:22, "X", "Y")))',
    '  locus_source[, chromosome_order := unname(chrom_order[Chromosome])]',
    '  data.table::setorder(locus_source, chromosome_order, chromStart)',
    '  data.table::fwrite(locus_source[, chromosome_order := NULL],',
    '                     file.path(source_dir, "Figure_4F_biomarker_instrument_loci_source.tsv"), sep = "\\t")',
    '  gene_loci <- locus_source[, .(',
    '    Chromosome, chromStart, chromEnd,',
    '    Gene = gene_symbol,',
    '    PlotColor = ifelse(family == "module", "#6A51A3", "#008C95")',
    '  )]',
    '  data.table::fwrite(gene_loci, file.path(source_dir, "Figure_4_biomarker_loci_source.tsv"), sep = "\\t")',
    '  module_gene_loci <- gene_loci[PlotColor == "#6A51A3"]',
    '  species_gene_loci <- gene_loci[PlotColor == "#008C95"]',
    '  snp_loci <- locus_source[, .(',
    '    Chromosome = snp_Chromosome, chromStart = snp_start, chromEnd = snp_end,',
    '    Direction = ifelse(OR > 1, 1, -1)',
    '  )]',
    '  data(UCSC.HG38.Human.CytoBandIdeogram, package = "RCircos")',
    '  cyto_info <- UCSC.HG38.Human.CytoBandIdeogram',
    '  chromosome_panel <- ggplotify::as.ggplot(~{',
    '    RCircos::RCircos.Set.Core.Components(',
    '      cyto.info = cyto_info, chr.exclude = NULL, tracks.inside = 5, tracks.outside = 2',
    '    )',
    '    rcircos_parameters <- RCircos::RCircos.Get.Plot.Parameters()',
    '    rcircos_parameters$text.size <- 0.70',
    '    rcircos_parameters$point.size <- 2.0',
    '    rcircos_parameters$point.type <- 16',
    '    RCircos::RCircos.Reset.Plot.Parameters(rcircos_parameters)',
    '    RCircos::RCircos.Set.Plot.Area()',
    '    RCircos::RCircos.Chromosome.Ideogram.Plot()',
    '    RCircos::RCircos.Gene.Connector.Plot(',
    '      genomic.data = module_gene_loci[, 1:4], track.num = 1, side = "in"',
    '    )',
    '    RCircos::RCircos.Gene.Name.Plot(',
    '      gene.data = module_gene_loci, name.col = 4, track.num = 2, side = "in"',
    '    )',
    '    RCircos::RCircos.Gene.Connector.Plot(',
    '      genomic.data = species_gene_loci[, 1:4], track.num = 1, side = "out"',
    '    )',
    '    RCircos::RCircos.Gene.Name.Plot(',
    '      gene.data = species_gene_loci, name.col = 4, track.num = 2, side = "out"',
    '    )',
    '    RCircos::RCircos.Scatter.Plot(',
    '      scatter.data = snp_loci, data.col = 4, track.num = 5, side = "in",',
    '      by.fold = 0.5, min.value = -1, max.value = 1',
    '    )',
    '    text(0, 0, "GRCh38", font = 2, cex = 0.78)',
    '  }) + labs(',
    '    tag = letters[[length(gsea_panels) + 1L]],',
    '    title = "Biomarker and instrument loci",',
    '    caption = "Gene label: purple = metabolic module; teal = species presence/absence\\nInstrument SNP: red = OR > 1; blue = OR < 1"',
    '  ) + theme(',
    '    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),',
    '    plot.caption = element_text(size = 7.5, hjust = 0.5, lineheight = 0.95),',
    '    plot.margin = margin(5, 5, 8, 5)',
    '  )',
    '  gsea_grid <- wrap_plots(gsea_panels, ncol = 2)',
    '  combined <- gsea_grid / chromosome_panel +',
    '    plot_layout(heights = c(3, 1.90)) +',
    '    plot_annotation() &',
    '    theme(plot.tag = element_text(face = "bold", size = 14))'
  )
  code <- c(code[seq_len(block_start - 1L)], dynamic_layout, code[block_end:length(code)])
  code <- gsub(
    "combined, width = 16, height = 20",
    "combined, width = 14, height = 22",
    code,
    fixed = TRUE
  )
}
if (stage == "immune") {
  code <- gsub(
    "four-biomarker associations",
    "V2 biomarker associations",
    code,
    fixed = TRUE
  )
  code <- gsub(
    "four direction-consistent biomarkers",
    "V2 direction-consistent biomarkers",
    code,
    fixed = TRUE
  )
}

message("Running V1 ", stage, " algorithm with V2 biomarkers: ",
        paste(biomarkers, collapse = ", "))
eval(parse(text = code, keep.source = TRUE), envir = .GlobalEnv)
