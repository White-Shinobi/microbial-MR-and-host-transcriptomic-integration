#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
res <- file.path(root, "results/bulk/crc/ridge_nomogram_external_validation")
out <- file.path(root, "manuscript/figures/main")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read_tsv <- function(x) read.delim(file.path(res, x), check.names = FALSE, stringsAsFactors = FALSE)

coef <- read_tsv("Ridge_nomogram_coefficients.tsv")
pts <- read_tsv("Ridge_nomogram_point_scale.tsv")
cal <- read_tsv("external_smooth_calibration_curves.tsv")
perf <- read_tsv("Ridge_nomogram_external_performance.tsv")
roc <- read_tsv("external_ROC_coordinates.tsv")
chall <- read_tsv("two_cohort_disease_state_individual_scores.tsv")
chall_summary <- read_tsv("two_cohort_disease_state_summary.tsv")
chall_pairwise <- read_tsv("two_cohort_disease_state_pairwise_tests.tsv")
train_pred <- read_tsv("GSE44076_training_predictions.tsv")
expression_long <- read.delim(
  file.path(root, "manuscript/source_data/Figures_2H_2I_expression_source.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)

navy <- "#17365D"; blue <- "#2F75B5"; teal <- "#00A6A6"; orange <- "#E07A38"
theme_pub <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.tag = element_text(face = "bold", size = 14),
        legend.title = element_blank(), legend.position = "inside")

# A: ridge nomogram with development-cohort distributions, following the
# distribution-enhanced nomogram presentation used in ref1 Figure 3A.
genes <- c("GTF2IRD1", "NRXN1", "FAM135B", "KIAA1671", "WFDC2")
ys <- setNames(6:2, genes)
training_expression <- expression_long[
  expression_long$dataset == "GSE44076" & expression_long$gene_symbol %in% genes,
  c("gene_symbol", "expression")
]
expression_scale <- do.call(rbind, lapply(genes, function(g) {
  x <- training_expression$expression[training_expression$gene_symbol == g]
  data.frame(gene_symbol = g, expression_mean = mean(x), expression_sd = sd(x),
             expression_min = min(x), expression_max = max(x),
             stringsAsFactors = FALSE)
}))
nom_axis_source <- merge(pts, expression_scale, by = "gene_symbol", sort = FALSE)
nom_axis_source$log2_normalized_expression <- with(
  nom_axis_source, expression_mean + z_score * expression_sd
)
nom_display_source <- do.call(rbind, lapply(genes, function(g) {
  d <- nom_axis_source[nom_axis_source$gene_symbol == g, ]
  raw_ticks <- seq(d$expression_min[1], d$expression_max[1], length.out = 5)
  z_ticks <- (raw_ticks - d$expression_mean[1]) / d$expression_sd[1]
  data.frame(label = g, y = unname(ys[g]),
             x = approx(d$z_score, d$points, xout = z_ticks, rule = 2)$y,
             tick = sprintf("%.1f", raw_ticks), z_score = z_ticks,
             log2_normalized_expression = raw_ticks)
}))
nom <- nom_display_source[, c("label", "y", "x", "tick")]
top_ticks <- data.frame(label = "Points", y = 7, x = seq(0, 100, 20), tick = seq(0, 100, 20))
max_total <- sum(tapply(pts$points, pts$gene_symbol, max))
total_values <- seq(0, floor(max_total / 50) * 50, 50)
total_ticks <- data.frame(label = "Total points", y = 1, x = total_values / max_total * 100, tick = total_values)
b <- setNames(coef$coefficient, coef$term)
max_span <- max(abs(b[genes]) * 5); scale_lp <- 100 / max_span
lp_min <- b["Intercept"] + sum(sapply(genes, function(g) b[g] * ifelse(b[g] > 0, -2.5, 2.5)))
pr <- c(.01, .10, .50, .90, .99); pr_x <- (qlogis(pr) - lp_min) * scale_lp / max_total * 100
keep <- pr_x >= 0 & pr_x <= 100
prob_ticks <- data.frame(label = "Predicted probability", y = 0, x = pr_x[keep], tick = sprintf("%.2f", pr[keep]))
all_ticks <- rbind(top_ticks, nom, total_ticks, prob_ticks)
axis_df <- do.call(rbind, lapply(split(all_ticks, paste(all_ticks$label, all_ticks$y)), function(d) {
  data.frame(label = d$label[1], y = d$y[1], xmin = min(d$x), xmax = max(d$x))
}))
distribution_polygon <- function(x, y, from = min(x), to = max(x), height = .20) {
  den <- density(x, from = from, to = to, n = 256, adjust = .85)
  h <- den$y / max(den$y) * height
  data.frame(x = c(den$x, rev(den$x)), y = c(y + h, rev(y - h)))
}
dist_poly <- data.frame()
dist_summary <- data.frame()
for (g in genes) {
  d <- pts[pts$gene_symbol == g, ]
  x <- approx(d$z_score, d$points, xout = train_pred[[g]], rule = 2)$y
  poly <- distribution_polygon(x, ys[g], min(x), max(x))
  poly$label <- g
  dist_poly <- rbind(dist_poly, poly)
  dist_summary <- rbind(dist_summary, data.frame(
    label = g, y = ys[g], q1 = quantile(x, .25), median = median(x), q3 = quantile(x, .75)))
}
total_x <- train_pred$nomogram_total_points / max_total * 100
poly <- distribution_polygon(total_x, 1, 0, 100)
poly$label <- "Total points"
dist_poly <- rbind(dist_poly, poly)
dist_summary <- rbind(dist_summary, data.frame(
  label = "Total points", y = 1, q1 = quantile(total_x, .25),
  median = median(total_x), q3 = quantile(total_x, .75)))

# A real development-cohort sample nearest to a predicted probability of 0.50
# illustrates how predictor values are projected to points, summed and then
# converted to a predicted probability. This is a reading aid, not a separate
# estimate or validation result.
example_sample <- train_pred[which.min(abs(train_pred$predicted_probability - .50)), ]
example_long <- do.call(rbind, lapply(genes, function(g) {
  d <- pts[pts$gene_symbol == g, ]
  scale_row <- expression_scale[expression_scale$gene_symbol == g, ]
  z <- example_sample[[g]]
  data.frame(
    sample_id = example_sample$sample_id,
    individual_id = example_sample$individual_id,
    observed_group = example_sample$group,
    gene_symbol = g,
    y = unname(ys[g]),
    z_score = z,
    log2_normalized_expression = scale_row$expression_mean + z * scale_row$expression_sd,
    gene_points = approx(d$z_score, d$points, xout = z, rule = 2)$y,
    nomogram_total_points = example_sample$nomogram_total_points,
    predicted_probability = example_sample$predicted_probability,
    stringsAsFactors = FALSE
  )
}))
example_total_x <- example_sample$nomogram_total_points / max_total * 100
red <- "#D73027"
pA <- ggplot() +
  geom_polygon(data = dist_poly, aes(x = x, y = y, group = label),
               fill = "#DCE9F5", color = "#3C78D8", linewidth = .35) +
  geom_segment(data = axis_df, aes(x = xmin, xend = xmax, y = y, yend = y), color = navy, linewidth = .55) +
  geom_segment(data = dist_summary, aes(x = q1, xend = q3, y = y, yend = y),
               color = "#222222", linewidth = 1.4) +
  geom_point(data = dist_summary, aes(x = median, y = y), shape = 21,
             fill = "#B7C9DA", color = "#222222", size = 1.7, stroke = .35) +
  geom_segment(data = example_long,
               aes(x = gene_points, xend = gene_points, y = y, yend = 7),
               color = red, linetype = "dotted", linewidth = .48, alpha = .72) +
  geom_point(data = example_long, aes(x = gene_points, y = y),
             shape = 21, fill = red, color = red, size = 2.0) +
  geom_point(data = example_long, aes(x = gene_points, y = 7),
             shape = 21, fill = red, color = red, size = 2.0) +
  geom_segment(data = all_ticks, aes(x = x, xend = x, y = y - .07, yend = y + .07), color = navy, linewidth = .3) +
  geom_text(data = subset(all_ticks, label != "Predicted probability"),
            aes(x = x, y = y - .16, label = tick), size = 2.15, vjust = 1) +
  geom_text(data = subset(all_ticks, label == "Predicted probability"),
            aes(x = x, y = y - .14, label = tick), size = 2.05, angle = 45, hjust = 1, vjust = 1) +
  geom_text(data = axis_df, aes(x = -5, y = y, label = label), hjust = 1, fontface = "bold", size = 2.7) +
  annotate("point", x = example_total_x, y = 1, shape = 23,
           size = 2.8, stroke = .45, fill = red, color = red) +
  annotate("segment", x = example_total_x, xend = example_total_x,
           y = .84, yend = .13, color = red, linewidth = .65,
           arrow = grid::arrow(length = grid::unit(.075, "inches"), type = "closed")) +
  annotate("point", x = example_total_x, y = 0, shape = 25,
           size = 2.7, stroke = .4, fill = red, color = red) +
  annotate("text", x = example_total_x - 1.8, y = .67,
           label = sprintf("%.1f", example_sample$nomogram_total_points),
           hjust = 1, vjust = .5, color = red, fontface = "bold", size = 2.3) +
  annotate("text", x = example_total_x + 1.5, y = .16,
           label = sprintf("p = %.3f", example_sample$predicted_probability),
           hjust = 0, vjust = .5, color = red, fontface = "bold", size = 2.3) +
  annotate("text", x = 0, y = -.55,
           label = "Axes: GSE44076 log2-normalized expression; ridge model: corresponding z scores",
           hjust = 0, size = 2.15, color = "#666666") +
  coord_cartesian(xlim = c(-25, 104), ylim = c(-.7, 7.35), clip = "off") +
  labs(title = "Five-gene ridge nomogram", tag = "A") + theme_void(base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 11), plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(8, 10, 3, 25))

# B: smooth parametric calibration curves with Hosmer-Lemeshow tests. The
# fitted curves use the already reported calibration intercept and slope.
cal_grid <- seq(.0005, .9995, length.out = 500)
external_datasets <- c("GSE41258", "GSE37364")
cal <- do.call(rbind, lapply(external_datasets, function(ds) {
  d <- perf[perf$dataset == ds, ]
  data.frame(dataset = ds, predicted_probability = cal_grid,
             actual_probability = plogis(d$calibration_intercept +
                                          d$calibration_slope * qlogis(cal_grid)))
}))
cal$dataset <- factor(cal$dataset, levels = external_datasets)
format_hl_p <- function(p) ifelse(p < .001, sprintf("%.2e", p), sprintf("%.3f", p))
hl_text <- paste(sprintf("%s: Hosmer–Lemeshow P = %s",
                         external_datasets,
                         format_hl_p(perf$Hosmer_Lemeshow_p[match(external_datasets, perf$dataset)])),
                 collapse = "\n")
pB <- ggplot(cal, aes(predicted_probability, actual_probability, color = dataset)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#777777") +
  geom_line(linewidth = .9) +
  scale_color_manual(values = c(blue, teal)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  annotate("label", x = .04, y = .96, label = hl_text, hjust = 0, vjust = 1,
           size = 2.55, fill = "white", color = "#222222") +
  labs(title = "External calibration: CRC versus normal + adenoma", x = "Nomogram-predicted probability",
       y = "Actual CRC probability", tag = "B") + theme_pub +
  theme(legend.position.inside = c(.76, .17))

# C: empirical ROC step curves. Repeated coordinates are reduced to the upper
# envelope at each false-positive rate without changing the empirical AUCs.
roc <- do.call(rbind, lapply(split(roc, roc$dataset), function(d) {
  upper <- aggregate(true_positive_rate ~ false_positive_rate, d, max)
  upper <- rbind(data.frame(false_positive_rate = 0, true_positive_rate = 0), upper,
                 data.frame(false_positive_rate = 1, true_positive_rate = 1))
  upper <- unique(upper[order(upper$false_positive_rate, upper$true_positive_rate), ])
  upper$dataset <- d$dataset[1]
  upper
}))
roc$dataset <- factor(roc$dataset, levels = external_datasets)
labels <- setNames(sprintf("%s: AUC %.3f (%.3f–%.3f)", perf$dataset, perf$AUC,
                           perf$AUC_CI_low, perf$AUC_CI_high), perf$dataset)
pC <- ggplot(roc, aes(false_positive_rate, true_positive_rate, color = dataset)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#999999") +
  geom_step(linewidth = .85, direction = "hv") +
  scale_color_manual(values = c(blue, teal), labels = labels) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "CRC versus composite controls", x = "1 − Specificity",
       y = "Sensitivity", tag = "C") + theme_pub +
  theme(legend.position.inside = c(.62, .15), legend.text = element_text(size = 7.5))

# D: normal--adenoma--CRC score distributions in both external cohorts.
chall$dataset <- factor(chall$dataset, levels = external_datasets)
chall$group <- factor(chall$group, levels = c("Normal", "Adenoma", "Tumor"),
                      labels = c("Normal", "Adenoma", "CRC"))
challenge_annotation <- do.call(rbind, lapply(external_datasets, function(ds) {
  counts <- table(chall$group[chall$dataset == ds])
  kw <- unique(chall_summary$Kruskal_Wallis_p[chall_summary$dataset == ds])[1]
  q <- chall_pairwise$BH_adjusted_p[chall_pairwise$dataset == ds & (
    (chall_pairwise$group_1 == "Tumor" & chall_pairwise$group_2 == "Adenoma") |
    (chall_pairwise$group_1 == "Adenoma" & chall_pairwise$group_2 == "Tumor"))][1]
  data.frame(dataset = factor(ds, levels = external_datasets), x = .55,
             y = max(chall$ridge_nomogram_total_points[chall$dataset == ds]),
             label = sprintf("n Normal/Adenoma/CRC = %d/%d/%d\nKruskal–Wallis P = %.2e\nCRC vs adenoma q = %.3g",
                             counts["Normal"], counts["Adenoma"], counts["CRC"], kw, q))
}))
pD <- ggplot(chall, aes(group, ridge_nomogram_total_points, fill = group, color = group)) +
  geom_violin(alpha = .2, width = .82, linewidth = .4, trim = FALSE) +
  geom_boxplot(width = .24, alpha = .85, outlier.shape = NA, color = "#444444", linewidth = .35) +
  geom_jitter(width = .10, size = 1, alpha = .6, show.legend = FALSE) +
  scale_fill_manual(values = c("#7F8C8D", "#E6A23C", "#C44E52")) +
  scale_color_manual(values = c("#7F8C8D", "#E6A23C", "#C44E52")) +
  geom_label(data = challenge_annotation, aes(x = x, y = y, label = label),
             inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.15, fill = "white") +
  facet_wrap(~dataset, nrow = 1, scales = "free_y") +
  labs(title = "Normal–adenoma–CRC scores", x = NULL,
       y = "Nomogram total points", tag = "D") +
  theme_pub + theme(legend.position = "none", strip.background = element_blank(),
                    strip.text = element_text(face = "bold", size = 8))

fig <- (pA | pB) / (pC | pD)
ggsave(file.path(out, "Figure_3_Ridge_nomogram_external_validation.pdf"), fig,
       width = 12.2, height = 9.3, units = "in", device = cairo_pdf)
ggsave(file.path(out, "Figure_3_Ridge_nomogram_external_validation.png"), fig,
       width = 12.2, height = 9.3, units = "in", dpi = 400, bg = "white")

# Display-source tables are written separately from the empirical calibration
# and ROC coordinates retained by the analysis script.
source_out <- file.path(root, "manuscript/source_data")
write.table(dist_summary, file.path(source_out, "Figure_3A_nomogram_distribution_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(nom_axis_source,
            file.path(source_out, "Figure_3A_nomogram_expression_axis_source.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(nom_display_source,
            file.path(source_out, "Figure_3A_nomogram_displayed_ticks.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(example_long,
            file.path(source_out, "Figure_3A_nomogram_example_sample.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cal, file.path(source_out, "Figure_3B_external_parametric_calibration_display.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(roc, file.path(source_out, "Figure_3C_external_empirical_ROC_display.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
