#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)][1])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
out_dir <- file.path(root, "results/bulk/crc/immune_infiltration")
figure_dir <- file.path(root, "manuscript/figures/main")
source_dir <- file.path(root, "manuscript/source_data")
gate <- fread(file.path(root, "results/bulk/crc/machine_learning/biomarker_gate_summary.tsv"))
biomarkers <- gate[final_biomarker == TRUE, gene_symbol]

fractions <- fread(file.path(out_dir, "GSE44076_LM22_CIBERSORT_compatible_fractions.tsv"))
long <- fread(file.path(out_dir, "GSE44076_LM22_fractions_long.tsv"))
differential <- fread(file.path(out_dir, "GSE44076_LM22_differential_immune_cells.tsv"))
biomarker_cor <- fread(file.path(out_dir, "GSE44076_biomarker_DIC_spearman_correlations.tsv"))
cell_columns <- setdiff(
  names(fractions),
  c("sample_id", "group", "individual_id", "P_value", "Correlation", "RMSE", "Nu")
)
dics <- differential[ref1_DIC == TRUE][
  order(-abs(median_difference_tumor_minus_normal)), immune_cell
]

p_stars <- function(p) {
  fifelse(p < 0.001, "***", fifelse(p < 0.01, "**", fifelse(p < 0.05, "*", "")))
}
cell_palette <- setNames(hcl.colors(length(cell_columns), "Dark 3"), cell_columns)
sample_levels_a <- fractions[
  order(factor(group, levels = c("Normal", "Tumor")), individual_id), sample_id
]
plot_a_data <- copy(long)
plot_a_data[, sample_order := factor(sample_id, levels = sample_levels_a)]
plot_a_data[, immune_cell := factor(immune_cell, levels = cell_columns)]
plot_a <- ggplot(plot_a_data, aes(sample_order, fraction, fill = immune_cell)) +
  geom_col(width = 1, linewidth = 0) +
  facet_grid(
    ~group, scales = "free_x", space = "free_x",
    labeller = as_labeller(c(Normal = "Control", Tumor = "CRC"))
  ) +
  scale_fill_manual(values = cell_palette, drop = FALSE) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(
    x = NULL, y = "Estimated immune-cell fraction",
    fill = "Immune cell", tag = "a"
  ) +
  theme_classic(base_size = 8.5) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.key.height = grid::unit(2.6, "mm"),
    legend.key.width = grid::unit(2.6, "mm"),
    legend.text = element_text(size = 7),
    strip.background = element_rect(fill = "#F2F2F2", color = "#BDBDBD"),
    panel.spacing.x = grid::unit(2, "mm")
  ) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

plot_b_data <- long[immune_cell %in% dics]
plot_b_data[, immune_cell := factor(immune_cell, levels = dics)]
cell_p_values <- setNames(
  differential[match(levels(plot_b_data$immune_cell), immune_cell), wilcoxon_rank_sum_p],
  levels(plot_b_data$immune_cell)
)
cell_facet_labels <- setNames(
  paste0(names(cell_p_values), "  ", p_stars(cell_p_values)),
  names(cell_p_values)
)
plot_b <- ggplot(plot_b_data, aes(group, fraction, fill = group, color = group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.72, linewidth = 0.35) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white",
               color = "#333333", linewidth = 0.3) +
  facet_wrap(~immune_cell, scales = "free_y", ncol = 4,
             labeller = as_labeller(cell_facet_labels)) +
  scale_x_discrete(labels = c(Normal = "Control", Tumor = "CRC")) +
  scale_fill_manual(values = c(Normal = "#4C9BB0", Tumor = "#E07A5F")) +
  scale_color_manual(values = c(Normal = "#2B7184", Tumor = "#B95742")) +
  labs(x = NULL, y = "Estimated immune-cell fraction", tag = "b") +
  theme_classic(base_size = 8.5) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#F5F5F5", color = "#D0D0D0"),
    strip.text = element_text(face = "bold", size = 7.2),
    axis.text.x = element_text(size = 7)
  )

plot_c_data <- copy(long[immune_cell %in% dics])
plot_c_data[, z_fraction := as.numeric(scale(fraction)), by = immune_cell]
plot_c_data[, z_fraction := pmax(pmin(z_fraction, 2.5), -2.5)]
sample_levels <- fractions[
  order(factor(group, levels = c("Normal", "Tumor")), individual_id), sample_id
]
plot_c_data[, sample_order := factor(sample_id, levels = sample_levels)]
plot_c_data[, immune_cell := factor(immune_cell, levels = rev(dics))]
plot_c <- ggplot(plot_c_data, aes(sample_order, immune_cell, fill = z_fraction)) +
  geom_tile() +
  facet_grid(~group, scales = "free_x", space = "free_x",
             labeller = as_labeller(c(Normal = "Control", Tumor = "CRC"))) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish
  ) +
  labs(x = NULL, y = NULL, fill = "Row z-score", tag = "c") +
  theme_classic(base_size = 8.5) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "#F2F2F2", color = "#BDBDBD"),
    panel.spacing.x = grid::unit(2, "mm"), legend.position = "right"
  )

plot_d_data <- copy(biomarker_cor)
plot_d_data[, significance := p_stars(p_value)]
plot_d_data[, label := sprintf("%.3f\n%s", rho, significance)]
plot_d_data[, text_color := fifelse(abs(rho) >= 0.55, "white", "#202020")]
plot_d_data[, biomarker := factor(biomarker, levels = biomarkers)]
plot_d_data[, immune_cell := factor(immune_cell, levels = rev(dics))]
plot_d <- ggplot(plot_d_data, aes(biomarker, immune_cell, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = label, color = text_color), lineheight = 0.88, size = 2.65) +
  scale_color_identity() +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(
    x = NULL, y = NULL, fill = "Spearman rho",
    caption = "* P < 0.05    ** P < 0.01    *** P < 0.001", tag = "d"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 8), legend.position = "right",
    plot.caption = element_text(hjust = 1, size = 8)
  )

figure <- plot_a / plot_b / plot_c / plot_d +
  plot_layout(heights = c(1.35, 1.65, 0.95, 1.35)) +
  plot_annotation() &
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold")
  )
ggsave(file.path(figure_dir, "Figure_5_immune_infiltration.pdf"), figure,
       width = 14, height = 21, device = cairo_pdf)
ggsave(file.path(figure_dir, "Figure_5_immune_infiltration.png"), figure,
       width = 14, height = 21, dpi = 300)
fwrite(plot_c_data, file.path(source_dir, "Figure_5C_differential_immune_heatmap_source.tsv"),
       sep = "\t")
fwrite(plot_a_data, file.path(source_dir, "Figure_5A_sample_level_stacked_immune_fractions.tsv"),
       sep = "\t")
message("Replotted Figure 5 with ", length(biomarkers), " V2 biomarkers.")
