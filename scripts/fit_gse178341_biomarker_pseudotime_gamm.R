#!/usr/bin/env Rscript

# Fit patient-aware smooth expression trajectories for the two focal
# epithelial biomarkers in the 36 matched GSE178341 patients.

suppressPackageStartupMessages(library(mgcv))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg))
root <- normalizePath(file.path(dirname(script_path), ".."))
out_dir <- file.path(
  root, "results/single_cell/GSE178341/label_free_figure7"
)
source_dir <- file.path(root, "manuscript/source_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

units <- read.delim(
  file.path(out_dir, "label_free_biomarker_patient_pseudotime_bins.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
units$eligible_for_plot <- tolower(
  as.character(units$eligible_for_plot)
) %in% c("true", "t", "1")

patient_tissues <- unique(units[, c("PID", "SPECIMEN_TYPE")])
matched <- intersect(
  patient_tissues$PID[patient_tissues$SPECIMEN_TYPE == "N"],
  patient_tissues$PID[patient_tissues$SPECIMEN_TYPE == "T"]
)
analysis <- units[
  units$eligible_for_plot & units$PID %in% matched,
  ,
  drop = FALSE
]
analysis$tissue <- factor(
  ifelse(analysis$SPECIMEN_TYPE == "N", "Normal", "CRC"),
  levels = c("Normal", "CRC")
)
analysis$tissue_difference <- ordered(
  analysis$tissue, levels = c("Normal", "CRC")
)
analysis$PID <- factor(analysis$PID, levels = sort(matched))
analysis$graph_pseudotime <- analysis$median_cell_pseudotime
analysis$expression <- analysis$log2_CPM_plus_0_5

prediction_rows <- list()
statistics_rows <- list()
model_dir <- file.path(out_dir, "biomarker_pseudotime_gamm_models")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

for (gene_name in c("GTF2IRD1", "KIAA1671")) {
  data_gene <- analysis[analysis$gene == gene_name, , drop = FALSE]
  data_gene <- droplevels(data_gene)

  full_model <- gam(
    expression ~
      tissue +
      s(graph_pseudotime, bs = "cr", k = 4) +
      s(
        graph_pseudotime,
        by = tissue_difference,
        bs = "cr",
        k = 4
      ) +
      s(graph_pseudotime, PID, bs = "fs", k = 3, m = 1),
    data = data_gene,
    method = "ML"
  )
  null_model <- gam(
    expression ~
      s(graph_pseudotime, bs = "cr", k = 4) +
      s(graph_pseudotime, PID, bs = "fs", k = 3, m = 1),
    data = data_gene,
    method = "ML"
  )
  comparison <- anova(null_model, full_model, test = "Chisq")
  p_column <- grep("^Pr\\(", colnames(comparison), value = TRUE)
  overall_p <- as.numeric(comparison[2, p_column[1]])
  deviance_difference <- as.numeric(comparison[2, "Deviance"])
  df_difference <- as.numeric(comparison[2, "Df"])

  smooth_table <- summary(full_model)$s.table
  difference_row <- grep(
    "tissue_difference",
    rownames(smooth_table),
    value = TRUE
  )
  difference_p <- if (length(difference_row) == 1) {
    as.numeric(smooth_table[difference_row, "p-value"])
  } else {
    NA_real_
  }
  parametric_table <- summary(full_model)$p.table
  tissue_p <- as.numeric(parametric_table["tissueCRC", "Pr(>|t|)"])

  statistics_rows[[gene_name]] <- data.frame(
    gene = gene_name,
    matched_patients = length(unique(data_gene$PID)),
    patient_tissue_bin_units = nrow(data_gene),
    normal_units = sum(data_gene$tissue == "Normal"),
    CRC_units = sum(data_gene$tissue == "CRC"),
    spline_basis_dimension = 4,
    overall_curve_LRT_deviance = deviance_difference,
    overall_curve_LRT_df = df_difference,
    overall_curve_P = overall_p,
    overall_tissue_level_P = tissue_p,
    tissue_difference_smooth_P = difference_p,
    observed_pseudotime_min = min(data_gene$graph_pseudotime),
    observed_pseudotime_max = max(data_gene$graph_pseudotime),
    full_model_AIC = AIC(full_model),
    null_model_AIC = AIC(null_model),
    stringsAsFactors = FALSE
  )

  prediction_grid <- expand.grid(
    graph_pseudotime = seq(0, 100, length.out = 201),
    tissue = factor(c("Normal", "CRC"), levels = c("Normal", "CRC")),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  prediction_grid$tissue_difference <- ordered(
    prediction_grid$tissue,
    levels = c("Normal", "CRC")
  )
  prediction_grid$PID <- factor(
    levels(data_gene$PID)[1],
    levels = levels(data_gene$PID)
  )
  predicted <- predict(
    full_model,
    newdata = prediction_grid,
    type = "response",
    se.fit = TRUE,
    exclude = "s(graph_pseudotime,PID)"
  )
  prediction_grid$gene <- gene_name
  prediction_grid$fit <- as.numeric(predicted$fit)
  prediction_grid$standard_error <- as.numeric(predicted$se.fit)
  prediction_grid$confidence_95_lower <- (
    prediction_grid$fit - 1.96 * prediction_grid$standard_error
  )
  prediction_grid$confidence_95_upper <- (
    prediction_grid$fit + 1.96 * prediction_grid$standard_error
  )
  prediction_rows[[gene_name]] <- prediction_grid[
    ,
    c(
      "gene", "tissue", "graph_pseudotime", "fit",
      "standard_error", "confidence_95_lower",
      "confidence_95_upper"
    )
  ]

  saveRDS(
    list(full = full_model, null = null_model),
    file.path(model_dir, paste0(gene_name, "_patient_GAMM.rds"))
  )
}

predictions <- do.call(rbind, prediction_rows)
statistics <- do.call(rbind, statistics_rows)
rownames(predictions) <- NULL
rownames(statistics) <- NULL

write.table(
  analysis,
  file.path(out_dir, "label_free_biomarker_pseudotime_GAMM_units.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  predictions,
  file.path(
    out_dir,
    "label_free_biomarker_pseudotime_GAMM_predictions.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  statistics,
  file.path(
    out_dir,
    "label_free_biomarker_pseudotime_GAMM_statistics.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  analysis,
  file.path(
    source_dir,
    "Figure_7H_GSE178341_biomarker_pseudotime_GAMM_units.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  predictions,
  file.path(
    source_dir,
    "Figure_7H_GSE178341_biomarker_pseudotime_GAMM_predictions.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  statistics,
  file.path(
    source_dir,
    "Figure_7H_GSE178341_biomarker_pseudotime_GAMM_statistics.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(statistics)
