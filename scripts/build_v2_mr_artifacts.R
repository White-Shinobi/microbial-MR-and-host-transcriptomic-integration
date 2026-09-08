#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)][1])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
mr_dir <- Sys.getenv("V2_MR_DIR", unset = file.path(root, "results/mr/crc/v2_primary"))
fig_main <- file.path(root, "manuscript/figures/main")
fig_supp <- file.path(root, "manuscript/figures/supplementary")
tab_main <- file.path(root, "manuscript/tables/main")
tab_supp <- file.path(root, "manuscript/tables/supplementary")
source_dir <- file.path(root, "manuscript/source_data")
for (d in c(fig_main, fig_supp, tab_main, tab_supp, source_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

read_required <- function(stem) {
  path <- file.path(mr_dir, stem)
  if (!file.exists(path)) stop("Missing required MR output: ", path)
  fread(path)
}
save_plot <- function(p, stem, width, height, supplementary = FALSE) {
  target <- if (supplementary) fig_supp else fig_main
  ggsave(file.path(target, paste0(stem, ".pdf")), p, width = width, height = height,
         units = "in", device = "pdf", bg = "white")
  ggsave(file.path(target, paste0(stem, ".png")), p, width = width, height = height,
         units = "in", dpi = 300, bg = "white")
}

mr <- read_required("v2_mr_results.tsv")
selection_manifest <- fread(file.path(root, "config/v2_ldsc_selected_exposures.tsv"))
instruments <- read_required("strict_ref1_harmonized_instruments.tsv")
eligibility <- read_required("strict_ref1_exposure_audit.tsv")
egger <- read_required("strict_ref1_egger_pleiotropy.tsv")
heterogeneity <- read_required("strict_ref1_heterogeneity.tsv")
presso <- read_required("strict_ref1_mr_presso_summary.tsv")
steiger <- read_required("strict_ref1_steiger.tsv")
loo <- read_required("strict_ref1_leave_one_out.tsv")
heterogeneity <- heterogeneity[grepl("Inverse variance weighted", method, fixed = TRUE),
                               .(accession, heterogeneity_Q = Q,
                                 heterogeneity_df = Q_df,
                                 heterogeneity_p = Q_pval)]
sensitivity <- Reduce(
  function(x, y) merge(x, y, by = "accession", all = TRUE),
  list(
    unique(mr[, .(accession, exposure)]),
    egger[, .(accession, egger_intercept, egger_intercept_se = se,
              egger_intercept_p = pval)],
    heterogeneity,
    presso[, .(accession, mr_presso_global_p = global_p,
               mr_presso_distortion_p = distortion_p,
               mr_presso_status = status)],
    steiger[, .(accession, steiger_correct_direction = correct_causal_direction,
                steiger_p = steiger_pval)]
  )
)
ivw <- mr[primary_ivw == TRUE]
manifest_columns <- intersect(
  c("accession", "cohort", "phenotype_class", "ldsc_h2", "ldsc_h2_se", "ldsc_h2_z"),
  names(selection_manifest)
)
ivw <- merge(
  ivw,
  unique(selection_manifest[, ..manifest_columns]),
  by = "accession",
  all.x = TRUE
)
ivw[, domain := fifelse(
  family == "module",
  "Metabolic module",
  fifelse(
    species_model_class == "presence_absence",
    "Species - presence/absence",
    "Species - abundance"
  )
)]
ivw[, q_value := as.numeric(v2_BH_FDR)]
ivw[, fdr_signal := as.logical(v2_fdr_signal)]
signals <- ivw[fdr_signal == TRUE]

ivw[, display_text := paste0(exposure, " [", cohort, "; ",
                             fifelse(family == "module", "module",
                                     gsub("_", "/", species_model_class)), "]")]
ivw[, display_name := factor(
  display_text,
  levels = rev(unique(display_text[order(domain, OR, accession)]))
)]
ivw[, significance := fifelse(fdr_signal, "BH-FDR q<0.05", "Not FDR-significant")]
forest_data <- ivw[fdr_signal == TRUE | pval < 0.05]
forest_data[, display_name := factor(
  display_text,
  levels = rev(unique(display_text[order(domain, OR, accession)]))
)]
p_forest <- ggplot(forest_data, aes(OR, display_name, color = significance)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.35) +
  geom_errorbar(aes(xmin = CI_95_low, xmax = CI_95_high), orientation = "y",
                width = 0.12, linewidth = 0.45, alpha = 0.85) +
  geom_point(size = 1.7) +
  facet_grid(domain ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
  scale_x_log10() +
  scale_color_manual(values = c("BH-FDR q<0.05" = "#B2182B",
                                "Not FDR-significant" = "#6B7280")) +
  labs(x = "Odds ratio for colorectal cancer (95% CI)", y = NULL, color = NULL) +
  theme_bw(base_size = 8.3) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"), axis.text.y = element_text(size = 6.2))
save_plot(p_forest, "Figure_1_primary_MR_forest", 10.4,
          max(7.5, 0.27 * nrow(forest_data) + 2.7))
fwrite(forest_data, file.path(source_dir, "Figure_1_primary_MR_forest_source.tsv"),
       sep = "\t")

if (nrow(signals)) {
  signal_accessions <- signals$accession
  inst <- instruments[accession %chin% signal_accessions & mr_keep == TRUE]
  inst[, exposure := signals$exposure[match(accession, signals$accession)]]
  method_lines <- mr[accession %chin% signal_accessions,
                     .(accession, exposure, method, b)]
  method_lines[, intercept := 0]
  if ("egger_intercept" %in% names(sensitivity)) {
    method_lines[grepl("Egger", method, ignore.case = TRUE),
                 intercept := sensitivity$egger_intercept[
                   match(accession, sensitivity$accession)]]
  }
  method_lines[, method := factor(
    method,
    levels = c("Inverse variance weighted (fixed effects)",
               "Inverse variance weighted (random effects)",
               "Weighted median", "MR Egger", "Simple mode", "Weighted mode")
  )]
  p_scatter <- ggplot(inst, aes(beta.exposure, beta.outcome)) +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
    geom_errorbar(aes(ymin = beta.outcome - 1.96 * se.outcome,
                      ymax = beta.outcome + 1.96 * se.outcome),
                  width = 0, alpha = 0.35, linewidth = 0.3) +
    geom_point(size = 1.5) +
    geom_abline(data = method_lines,
                aes(slope = b, intercept = intercept, color = method),
                linewidth = 0.65) +
    facet_wrap(~exposure, scales = "free", ncol = 2) +
    labs(x = "SNP effect on microbial exposure", y = "SNP effect on CRC", color = NULL) +
    theme_bw(base_size = 8.5) +
    theme(strip.text = element_text(face = "italic"), legend.position = "bottom")
  save_plot(p_scatter, "Figure_S1_MR_scatter", 9,
            max(5, 3.1 * ceiling(nrow(signals) / 2)), TRUE)
  fwrite(inst, file.path(source_dir, "Figure_S1_MR_scatter_instruments.tsv"), sep = "\t")
  fwrite(method_lines, file.path(source_dir, "Figure_S1_MR_scatter_slopes.tsv"), sep = "\t")

  single <- copy(inst)
  single[, `:=`(
    wald_b = beta.outcome / beta.exposure,
    wald_se = se.outcome / abs(beta.exposure)
  )]
  single[, `:=`(lo = wald_b - 1.96 * wald_se, hi = wald_b + 1.96 * wald_se)]
  p_single <- ggplot(single, aes(wald_b, reorder(SNP, wald_b))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.35) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                  width = 0.12, linewidth = 0.35) +
    geom_point(size = 1.25, color = "#2C5282") +
    facet_wrap(~exposure, scales = "free", ncol = 2) +
    labs(x = "SNP-specific log odds ratio", y = NULL) +
    theme_bw(base_size = 8)
  save_plot(p_single, "Figure_S2_SNP_level_forest", 10,
            max(6, 3.5 * ceiling(nrow(signals) / 2)), TRUE)

  single[, precision := 1 / wald_se]
  p_funnel <- ggplot(single, aes(wald_b, precision)) +
    geom_point(size = 1.5, alpha = 0.75, color = "#2C5282") +
    geom_vline(data = signals, aes(xintercept = b), linetype = "dashed",
               color = "#C53030", linewidth = 0.45) +
    facet_wrap(~exposure, scales = "free", ncol = 2) +
    labs(x = "SNP-specific log odds ratio", y = "Precision (1/SE)") +
    theme_bw(base_size = 8.5)
  save_plot(p_funnel, "Figure_S3_MR_funnel", 9,
            max(5, 3.1 * ceiling(nrow(signals) / 2)), TRUE)
  fwrite(single, file.path(source_dir, "Figures_S2_S3_SNP_level_source.tsv"), sep = "\t")

  loo_signal <- loo[accession %chin% signal_accessions]
  loo_signal[, `:=`(lo = b - 1.96 * se, hi = b + 1.96 * se)]
  p_loo <- ggplot(loo_signal, aes(b, reorder(SNP, b))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.35) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                  width = 0.12, linewidth = 0.35) +
    geom_point(size = 1.25, color = "#2F855A") +
    facet_wrap(~exposure, scales = "free", ncol = 2) +
    labs(x = "Leave-one-out IVW log odds ratio", y = "Omitted SNP") +
    theme_bw(base_size = 8)
  save_plot(p_loo, "Figure_S4_leave_one_out", 10,
            max(6, 3.5 * ceiling(nrow(signals) / 2)), TRUE)
  fwrite(loo_signal, file.path(source_dir, "Figure_S4_leave_one_out_source.tsv"), sep = "\t")
}

main <- ivw[, .(domain, cohort, phenotype_class, exposure, accession, nsnp,
                ivw_effect_model, OR, CI_95_low, CI_95_high, pval, q_value,
                species_model_stratified_BH_FDR, fdr_signal)]
setorder(main, domain, -fdr_signal, pval)
fwrite(main, file.path(tab_main, "Table_1_primary_IVW_associations.tsv"), sep = "\t")
fwrite(mr, file.path(tab_supp, "Table_S1_all_MR_estimators.tsv"), sep = "\t")
fwrite(eligibility, file.path(tab_supp, "Table_S2_instrument_eligibility_audit.tsv"), sep = "\t")
fwrite(sensitivity, file.path(tab_supp, "Table_S3_MR_sensitivity_all_exposures.tsv"), sep = "\t")
species_model_fdr <- ivw[
  family == "species",
  .(cohort, phenotype_class, exposure, accession, nsnp, OR, CI_95_low,
    CI_95_high, pval, species_family_BH_FDR = q_value,
    species_model_stratified_BH_FDR, fdr_signal)
]
setorder(species_model_fdr, phenotype_class, pval)
fwrite(
  species_model_fdr,
  file.path(tab_supp, "Table_S3A_species_model_stratified_FDR.tsv"),
  sep = "\t"
)

if (nrow(signals)) {
  signal_sensitivity <- merge(
    signals[, .(accession, exposure, family, cohort, phenotype_class, nsnp,
                ivw_effect_model, OR, CI_95_low, CI_95_high, pval, q_value)],
    sensitivity, by = c("accession", "exposure"), all.x = TRUE
  )
  setorder(signal_sensitivity, family, pval)
} else {
  signal_sensitivity <- data.table()
}
fwrite(signal_sensitivity, file.path(tab_main, "Table_2_MR_sensitivity.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(mr_dir, "v2_mr_artifacts_sessionInfo.txt"))
message("V2 MR manuscript artifacts complete; FDR signals: ", nrow(signals))
