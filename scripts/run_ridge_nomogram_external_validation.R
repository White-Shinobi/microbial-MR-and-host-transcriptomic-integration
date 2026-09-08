#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(pROC)
})

set.seed(20260721)
GENES <- c("GTF2IRD1", "NRXN1", "FAM135B", "KIAA1671", "WFDC2")

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork = TRUE)
OUT <- file.path(ROOT, "results/bulk/crc/ridge_nomogram_external_validation")
TABLE_MAIN <- file.path(ROOT, "manuscript/tables/main")
TABLE_SUPP <- file.path(ROOT, "manuscript/tables/supplementary")
SOURCE <- file.path(ROOT, "manuscript/source_data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_MAIN, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_SUPP, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE, row.names = FALSE, na = "")
clip <- function(p) pmin(pmax(as.numeric(p), 1e-9), 1 - 1e-9)
auc_metrics <- function(y, score) {
  roc_obj <- roc(y, score, levels = c(0, 1), direction = "<", quiet = TRUE)
  ci_obj <- as.numeric(ci.auc(roc_obj, method = "delong"))
  data.frame(AUC = as.numeric(auc(roc_obj)), AUC_CI_low = ci_obj[1],
             AUC_CI_high = ci_obj[3], stringsAsFactors = FALSE)
}

calibration_metrics <- function(y, lp) {
  p <- clip(plogis(lp))
  intercept_fit <- suppressWarnings(glm(y ~ 1, offset = lp, family = binomial(),
                                        control = glm.control(maxit = 100)))
  slope_fit <- suppressWarnings(glm(y ~ lp, family = binomial(), control = glm.control(maxit = 100)))
  data.frame(
    calibration_intercept = if (isTRUE(intercept_fit$converged)) unname(coef(intercept_fit)[1]) else NA_real_,
    calibration_slope = if (isTRUE(slope_fit$converged)) unname(coef(slope_fit)["lp"]) else NA_real_,
    Brier = mean((y - p)^2), stringsAsFactors = FALSE
  )
}

hosmer_lemeshow <- function(y, p, groups = 10L) {
  p <- clip(p)
  n <- length(y)
  groups <- min(as.integer(groups), floor(n / 5L))
  ranked <- rank(p, ties.method = "first")
  grp <- as.integer(cut(ranked, breaks = groups, labels = FALSE, include.lowest = TRUE))
  tab <- aggregate(cbind(observed_events = y, expected_events = p), list(group = grp), sum)
  tab$n <- as.integer(table(grp)[as.character(tab$group)])
  tab$observed_nonevents <- tab$n - tab$observed_events
  tab$expected_nonevents <- tab$n - tab$expected_events
  denominator <- tab$expected_events * (1 - tab$expected_events / tab$n)
  statistic <- sum((tab$observed_events - tab$expected_events)^2 / denominator)
  df <- groups - 2L
  list(statistic = statistic, df = df, p_value = pchisq(statistic, df = df, lower.tail = FALSE), table = tab)
}

smooth_calibration <- function(y, p, dataset, evaluation, span = 0.75) {
  p <- clip(p)
  fit <- tryCatch(loess(y ~ p, span = span, degree = 1,
                        control = loess.control(surface = "direct")), error = function(e) NULL)
  grid <- seq(max(min(p), 0.001), min(max(p), 0.999), length.out = 151L)
  if (is.null(fit)) {
    low <- lowess(p, y, f = span, iter = 0)
    actual <- approx(low$x, low$y, xout = grid, rule = 2, ties = mean)$y
  } else {
    actual <- as.numeric(predict(fit, newdata = data.frame(p = grid)))
  }
  keep <- is.finite(actual)
  data.frame(dataset = dataset, evaluation = evaluation,
             predicted_probability = grid[keep], actual_probability = pmin(pmax(actual[keep], 0), 1),
             stringsAsFactors = FALSE)
}

roc_coordinates <- function(y, score, dataset) {
  roc_obj <- roc(y, score, levels = c(0, 1), direction = "<", quiet = TRUE)
  xy <- coords(roc_obj, x = "all", ret = c("threshold", "specificity", "sensitivity"), transpose = FALSE)
  data.frame(dataset = dataset, model = "Ridge nomogram", threshold = xy$threshold,
             false_positive_rate = 1 - xy$specificity, true_positive_rate = xy$sensitivity,
             stringsAsFactors = FALSE)
}

patient_folds <- function(patient_ids, k = 10L) {
  ids <- unique(patient_ids)
  assignment <- sample(rep(seq_len(k), length.out = length(ids)))
  unname(setNames(assignment, ids)[patient_ids])
}

fit_ridge_cv <- function(x, y, pair_ids, seed) {
  set.seed(seed)
  foldid <- patient_folds(pair_ids, 10L)
  cv <- cv.glmnet(as.matrix(x), y, family = "binomial", alpha = 0,
                  standardize = FALSE, foldid = foldid, type.measure = "deviance", nlambda = 200)
  m <- as.matrix(coef(cv, s = "lambda.min"))
  beta <- m[c("(Intercept)", GENES), 1]
  names(beta)[1] <- "Intercept"
  list(cv = cv, beta = beta)
}

linear_predictor <- function(x, beta) {
  as.numeric(beta["Intercept"] + as.matrix(x[, GENES, drop = FALSE]) %*% beta[GENES])
}

series_table_skip <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con)); n <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE); n <- n + 1L
    if (!length(line)) stop("Missing series matrix table marker: ", path)
    if (identical(line, "!series_matrix_table_begin")) return(n)
  }
}

read_series_matrix <- function(path) {
  tab <- fread(cmd = paste("gzip -cd", shQuote(path)), skip = series_table_skip(path),
               header = TRUE, data.table = FALSE, fill = TRUE, showProgress = FALSE)
  names(tab)[1] <- "ID_REF"
  tab <- tab[tab$ID_REF != "!series_matrix_table_end" & nzchar(tab$ID_REF), , drop = FALSE]
  rownames(tab) <- tab$ID_REF; tab$ID_REF <- NULL; tab[] <- lapply(tab, as.numeric)
  as.matrix(tab)
}

read_affymetrix_annotation <- function(path) {
  lines <- readLines(gzfile(path, "rt"), warn = FALSE)
  start <- grep("^ID\\tGene title\\tGene symbol", lines)[1]
  end <- grep("^!platform_table_end", lines)[1]
  fields <- strsplit(lines[(start + 1L):(end - 1L)], "\\t")
  data.frame(ID = vapply(fields, function(x) if (length(x) >= 1L) x[1] else "", character(1)),
             `Gene symbol` = vapply(fields, function(x) if (length(x) >= 3L) x[3] else "", character(1)),
             check.names = FALSE, stringsAsFactors = FALSE)
}

probe_ids_for_gene <- function(annotation, gene) {
  symbols <- strsplit(ifelse(is.na(annotation[["Gene symbol"]]), "", annotation[["Gene symbol"]]),
                      "///", fixed = TRUE)
  keep <- vapply(symbols, function(x) gene %in% trimws(x), logical(1))
  unique(annotation$ID[keep])
}

aggregate_external_genes <- function(accession, matrix_path, audit_path, annotation,
                                     allow_unmeasured = character()) {
  expression <- read_series_matrix(matrix_path)
  transformed <- FALSE
  if (quantile(expression, 0.99, na.rm = TRUE) > 100) {
    expression <- log2(expression + 1); transformed <- TRUE
  }
  sample_audit <- read_tsv(audit_path)
  result <- data.frame(sample_id = sample_audit$sample_id, stringsAsFactors = FALSE)
  probe_audit <- data.frame()
  for (gene in GENES) {
    probes <- intersect(probe_ids_for_gene(annotation, gene), rownames(expression))
    if (!length(probes) && !gene %in% allow_unmeasured)
      stop(accession, ": no measured probe for ", gene)
    result[[gene]] <- if (length(probes)) {
      apply(expression[probes, sample_audit$sample_id, drop = FALSE], 2, median, na.rm = TRUE)
    } else {
      rep(NA_real_, nrow(sample_audit))
    }
    probe_audit <- rbind(probe_audit, data.frame(
      dataset = accession, gene_symbol = gene, n_mapped_measured_probes = length(probes),
      probe_ids = paste(probes, collapse = ";"),
      aggregation_rule = ifelse(length(probes),
                                "median across all officially annotated measured probes",
                                "not measured; neutral standardized-value imputation (z = 0)"),
      log2_transformed_before_aggregation = transformed, stringsAsFactors = FALSE))
  }
  joined <- merge(sample_audit, result, by = "sample_id", sort = FALSE)
  joined <- joined[match(sample_audit$sample_id, joined$sample_id), ]
  list(data = joined, probe_audit = probe_audit)
}

read_series_metadata_field <- function(path, field) {
  lines <- readLines(gzfile(path, "rt"), warn = FALSE)
  hit <- lines[startsWith(lines, paste0(field, "\t"))]
  if (length(hit) != 1L) stop("Expected one ", field, " line in ", path)
  values <- strsplit(hit, "\t", fixed = TRUE)[[1]][-1]
  sub('"$', "", sub('^"', "", values))
}

prepare_gse41258 <- function(matrix_path, annotation) {
  sample_ids <- read_series_metadata_field(matrix_path, "!Sample_geo_accession")
  titles <- read_series_metadata_field(matrix_path, "!Sample_title")
  group <- ifelse(startsWith(titles, "Primary Tumor"), "Tumor",
                  ifelse(startsWith(titles, "Normal Colon"), "Normal",
                         ifelse(grepl("^(Polyp|Microadenoma)", titles), "Adenoma", "Excluded")))
  technical_repeat <- grepl("(_ez|_ren|_2|rehyb)$", titles, ignore.case = TRUE)
  audit <- data.frame(
    sample_id = sample_ids, title = titles, group = group, stage = NA_character_,
    technical_repeat = technical_repeat,
    include_in_figure3 = group != "Excluded" & !technical_repeat,
    exclusion_reason = ifelse(group == "Excluded", "sample type outside primary CRC/normal/adenoma analysis",
                              ifelse(technical_repeat, "explicitly labelled repeat/rehybridized array", "")),
    stringsAsFactors = FALSE)
  write_tsv(audit, file.path(OUT, "GSE41258_sample_audit.tsv"))
  included <- audit[audit$include_in_figure3, c("sample_id", "title", "group", "stage")]
  audit_path <- file.path(OUT, "GSE41258_included_sample_manifest.tsv")
  write_tsv(included, audit_path)
  aggregate_external_genes("GSE41258", matrix_path, audit_path,
                           annotation = annotation,
                           allow_unmeasured = c("FAM135B", "KIAA1671"))
}

prepare_gse37364 <- function() {
  matrix_path <- file.path(ROOT, "data/raw/omics/crc_bulk/GSE37364/GSE37364_series_matrix.txt.gz")
  audit_path <- file.path(ROOT, "results/bulk/crc/external_validation/two_cohort_fig3",
                          "GSE37364_sample_audit.tsv")
  annotation <- read_affymetrix_annotation(
    file.path(ROOT, "data/raw/omics/crc_bulk/GSE39582/GPL570.annot.gz"))
  aggregate_external_genes("GSE37364", matrix_path, audit_path, annotation)
}

standardize_cohort <- function(frame) {
  means <- sapply(frame[, GENES, drop = FALSE], mean)
  sds <- sapply(frame[, GENES, drop = FALSE], sd)
  missing_gene <- !is.finite(means) | !is.finite(sds) | sds == 0
  safe_means <- means; safe_means[missing_gene] <- 0
  safe_sds <- sds; safe_sds[missing_gene] <- 1
  z <- sweep(sweep(as.matrix(frame[, GENES, drop = FALSE]), 2, safe_means, "-"), 2, safe_sds, "/")
  z[, missing_gene] <- 0
  list(z = as.data.frame(z, check.names = FALSE),
       parameters = data.frame(gene_symbol = GENES, mean = means, standard_deviation = sds,
                               standardization_rule = ifelse(missing_gene,
                                 "unmeasured on platform; neutral standardized-value imputation (z = 0)",
                                 "within-cohort z score without outcome labels"),
                               stringsAsFactors = FALSE))
}

nomogram_scale <- function(beta, z_min = -2.5, z_max = 2.5) {
  max_span <- max(abs(beta[GENES]) * (z_max - z_min)); rows <- data.frame()
  for (gene in GENES) {
    min_contribution <- min(beta[gene] * z_min, beta[gene] * z_max)
    for (z in seq(z_min, z_max, by = 0.5)) {
      rows <- rbind(rows, data.frame(gene_symbol = gene, z_score = z,
                                     points = (beta[gene] * z - min_contribution) / max_span * 100,
                                     coefficient = beta[gene], stringsAsFactors = FALSE))
    }
  }
  attr(rows, "max_span") <- max_span; attr(rows, "z_min") <- z_min; attr(rows, "z_max") <- z_max
  rows
}

total_points_from_z <- function(z, beta, point_scale) {
  total <- rep(0, nrow(z)); z_min <- attr(point_scale, "z_min"); z_max <- attr(point_scale, "z_max")
  max_span <- attr(point_scale, "max_span")
  for (gene in GENES) {
    zz <- pmin(pmax(z[[gene]], z_min), z_max)
    min_contribution <- min(beta[gene] * z_min, beta[gene] * z_max)
    total <- total + (beta[gene] * zz - min_contribution) / max_span * 100
  }
  total
}

# Development data.
long <- read_tsv(file.path(ROOT, "manuscript/source_data/Figures_2H_2I_expression_source.tsv"))
long <- long[long$dataset == "GSE44076" & long$gene_symbol %in% GENES, ]
training <- reshape(long[, c("sample_id", "group", "gene_symbol", "expression")],
                    idvar = c("sample_id", "group"), timevar = "gene_symbol", direction = "wide")
names(training) <- sub("^expression\\.", "", names(training))
manifest <- read_tsv(file.path(ROOT, "results/bulk/crc/discovery/GSE44076_sample_manifest.tsv"))
training$individual_id <- manifest$individual_id[match(training$sample_id, manifest$sample_id)]
training <- training[order(training$individual_id, training$group), ]
stopifnot(nrow(training) == 196L, length(unique(training$individual_id)) == 98L)

training_mean <- sapply(training[, GENES, drop = FALSE], mean)
training_sd <- sapply(training[, GENES, drop = FALSE], sd)
training_z <- as.data.frame(sweep(sweep(as.matrix(training[, GENES]), 2, training_mean, "-"),
                                       2, training_sd, "/"), check.names = FALSE)
training_y <- as.integer(training$group == "Tumor")

primary <- fit_ridge_cv(training_z, training_y, training$individual_id, 20260721)
beta <- primary$beta
training_lp <- linear_predictor(training_z, beta)
training_p <- clip(plogis(training_lp))
training_roc <- roc(training_y, training_lp, levels = c(0, 1), direction = "<", quiet = TRUE)
best_cut <- coords(training_roc, x = "best", best.method = "youden",
                   ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
training_threshold_lp <- as.numeric(best_cut$threshold[1])
training_threshold_probability <- plogis(training_threshold_lp)
point_scale <- nomogram_scale(beta)

training_predictions <- cbind(training[, c("sample_id", "individual_id", "group")], training_z)
training_predictions$observed_tumor <- training_y
training_predictions$linear_predictor <- training_lp
training_predictions$predicted_probability <- training_p
training_predictions$nomogram_total_points <- total_points_from_z(training_z, beta, point_scale)

# Frozen external evaluation. In both cohorts, the binary target is CRC versus
# the combined normal-mucosa and adenoma group. Coefficients and threshold are
# retained from GSE44076 without external refitting.
external_results <- list(); external_performance <- data.frame(); external_rocs <- data.frame()
external_curves <- data.frame(); hl_groups <- data.frame()
for (accession in c("GSE41258", "GSE37364")) {
  if (accession == "GSE41258") {
    matrix_path <- file.path(ROOT, "data/raw/omics/crc_bulk/GSE41258/GSE41258_series_matrix.txt.gz")
    annotation <- read_affymetrix_annotation(file.path(ROOT, "data/raw/omics/crc_bulk/GSE41258/GPL96.annot.gz"))
    aggregated <- prepare_gse41258(matrix_path, annotation)
  } else aggregated <- prepare_gse37364()
  scaled <- standardize_cohort(aggregated$data)
  frame <- cbind(aggregated$data[, c("sample_id", "title", "group", "stage")], scaled$z)
  frame$dataset <- accession
  frame$ridge_linear_predictor <- linear_predictor(frame, beta)
  frame$ridge_predicted_probability <- clip(plogis(frame$ridge_linear_predictor))
  frame$ridge_nomogram_total_points <- total_points_from_z(frame[, GENES], beta, point_scale)
  binary <- frame$group %in% c("Tumor", "Normal", "Adenoma")
  y <- as.integer(frame$group[binary] == "Tumor")
  lp <- frame$ridge_linear_predictor[binary]; p <- frame$ridge_predicted_probability[binary]
  auc_row <- auc_metrics(y, lp); cal_row <- calibration_metrics(y, lp)
  hl <- hosmer_lemeshow(y, p, 10L)
  pred <- lp >= training_threshold_lp
  sensitivity <- mean(pred[y == 1]); specificity <- mean(!pred[y == 0])
  external_performance <- rbind(external_performance, data.frame(
    dataset = accession,
    role = "external CRC-versus-normal-plus-adenoma validation cohort",
    outcome_definition = "CRC versus combined normal mucosa and adenoma",
    model = "Ridge nomogram at lambda.min", n_CRC = sum(y == 1),
    n_normal = sum(frame$group[binary] == "Normal"),
    n_adenoma = sum(frame$group[binary] == "Adenoma"), n_non_CRC = sum(y == 0),
    coefficient_refit = FALSE, AUC = auc_row$AUC, AUC_CI_low = auc_row$AUC_CI_low,
    AUC_CI_high = auc_row$AUC_CI_high, frozen_threshold_linear_predictor = training_threshold_lp,
    sensitivity = sensitivity, specificity = specificity,
    balanced_accuracy = mean(c(sensitivity, specificity)), Brier = cal_row$Brier,
    calibration_intercept = cal_row$calibration_intercept,
    calibration_slope = cal_row$calibration_slope,
    Hosmer_Lemeshow_groups = nrow(hl$table), Hosmer_Lemeshow_chisq = hl$statistic,
    Hosmer_Lemeshow_df = hl$df, Hosmer_Lemeshow_p = hl$p_value, stringsAsFactors = FALSE))
  hlt <- hl$table; hlt$dataset <- accession
  hl_groups <- rbind(hl_groups, hlt[, c("dataset", "group", "n", "observed_events",
                                       "expected_events", "observed_nonevents", "expected_nonevents")])
  external_rocs <- rbind(external_rocs, roc_coordinates(y, lp, accession))
  external_curves <- rbind(external_curves, smooth_calibration(y, p, accession, "Frozen external"))
  write_tsv(frame, file.path(OUT, paste0(accession, "_frozen_predictions.tsv")))
  write_tsv(aggregated$probe_audit, file.path(OUT, paste0(accession, "_probe_aggregation_audit.tsv")))
  scaled$parameters$dataset <- accession
  write_tsv(scaled$parameters, file.path(OUT, paste0(accession, "_zscore_parameters.tsv")))
  external_results[[accession]] <- frame
}

# Disease-state score distributions and comparisons in both cohorts.
group_order <- c("Normal", "Adenoma", "Tumor")
challenge_frame <- do.call(rbind, lapply(names(external_results), function(ds) {
  d <- external_results[[ds]]; d$group <- factor(d$group, levels = group_order); d
}))
challenge_summary <- data.frame(); pairwise_rows <- data.frame()
for (ds in names(external_results)) {
  d <- challenge_frame[challenge_frame$dataset == ds, ]
  kw <- kruskal.test(ridge_nomogram_total_points ~ group, data = d)
  pw <- pairwise.wilcox.test(d$ridge_nomogram_total_points, d$group,
                            p.adjust.method = "BH", exact = FALSE)
  challenge_summary <- rbind(challenge_summary, do.call(rbind, lapply(group_order, function(g) {
    v <- d$ridge_nomogram_total_points[d$group == g]
    data.frame(dataset = ds, group = g, n = length(v), median_total_points = median(v),
               IQR_total_points = IQR(v), min_total_points = min(v), max_total_points = max(v),
               Kruskal_Wallis_p = kw$p.value, stringsAsFactors = FALSE)
  })))
  for (i in seq_len(nrow(pw$p.value))) for (j in seq_len(ncol(pw$p.value))) {
    if (!is.na(pw$p.value[i, j])) pairwise_rows <- rbind(pairwise_rows, data.frame(
      dataset = ds, group_1 = rownames(pw$p.value)[i], group_2 = colnames(pw$p.value)[j],
      BH_adjusted_p = pw$p.value[i, j], stringsAsFactors = FALSE))
  }
}

coefficient_table <- data.frame(term = names(beta), coefficient = as.numeric(beta),
                                odds_ratio = exp(as.numeric(beta)), lambda_min = primary$cv$lambda.min,
                                lambda_1se = primary$cv$lambda.1se,
                                model = "Ridge logistic regression", stringsAsFactors = FALSE)
external_audit <- data.frame(
  dataset = c("GSE41258", "GSE37364"),
  role = "external CRC-versus-normal-plus-adenoma validation cohort",
  platform = c("GPL96", "GPL570"),
  measured_predictors = c("3/5; FAM135B and KIAA1671 not measured and assigned neutral z = 0", "5/5"),
  n_CRC = sapply(c("GSE41258", "GSE37364"), function(ds) sum(external_results[[ds]]$group == "Tumor")),
  n_normal = sapply(c("GSE41258", "GSE37364"), function(ds) sum(external_results[[ds]]$group == "Normal")),
  n_adenoma = sapply(c("GSE41258", "GSE37364"), function(ds) sum(external_results[[ds]]$group == "Adenoma")),
  probe_aggregation = c("median across officially annotated probes; neutral z = 0 for unmeasured FAM135B and KIAA1671",
                        "median across all officially annotated measured probes"),
  alignment = "outcome-blind within-cohort gene-wise z standardization",
  model = "Ridge logistic nomogram at lambda.min", coefficient_refit = FALSE,
  threshold_source = "frozen GSE44076 Youden threshold",
  binary_validation_samples = "CRC versus combined normal mucosa and adenoma", stringsAsFactors = FALSE)
challenge_table <- rbind(
  data.frame(record_type = "group_summary", dataset = challenge_summary$dataset,
             group = challenge_summary$group, n = challenge_summary$n,
             median_total_points = challenge_summary$median_total_points,
             IQR_total_points = challenge_summary$IQR_total_points,
             min_total_points = challenge_summary$min_total_points,
             max_total_points = challenge_summary$max_total_points,
             Kruskal_Wallis_p = challenge_summary$Kruskal_Wallis_p,
             comparison = NA_character_, BH_adjusted_p = NA_real_, stringsAsFactors = FALSE),
  data.frame(record_type = "pairwise_test", dataset = pairwise_rows$dataset,
             group = NA_character_, n = NA_integer_,
             median_total_points = NA_real_, IQR_total_points = NA_real_, min_total_points = NA_real_,
             max_total_points = NA_real_, Kruskal_Wallis_p = NA_real_,
             comparison = paste(pairwise_rows$group_1, "versus", pairwise_rows$group_2),
             BH_adjusted_p = pairwise_rows$BH_adjusted_p, stringsAsFactors = FALSE))

# Auditable results, manuscript tables and figure source data.
write_tsv(coefficient_table, file.path(OUT, "Ridge_nomogram_coefficients.tsv"))
write_tsv(point_scale, file.path(OUT, "Ridge_nomogram_point_scale.tsv"))
write_tsv(training_predictions, file.path(OUT, "GSE44076_training_predictions.tsv"))
write_tsv(data.frame(lambda = primary$cv$lambda, cv_deviance = primary$cv$cvm,
                     cv_standard_error = primary$cv$cvsd,
                     is_lambda_min = primary$cv$lambda == primary$cv$lambda.min,
                     is_lambda_1se = primary$cv$lambda == primary$cv$lambda.1se),
          file.path(OUT, "GSE44076_patient_grouped_ridge_CV.tsv"))
write_tsv(external_performance, file.path(OUT, "Ridge_nomogram_external_performance.tsv"))
write_tsv(external_curves, file.path(OUT, "external_smooth_calibration_curves.tsv"))
write_tsv(hl_groups, file.path(OUT, "external_Hosmer_Lemeshow_groups.tsv"))
write_tsv(external_rocs, file.path(OUT, "external_ROC_coordinates.tsv"))
write_tsv(challenge_frame, file.path(OUT, "two_cohort_disease_state_individual_scores.tsv"))
write_tsv(challenge_summary, file.path(OUT, "two_cohort_disease_state_summary.tsv"))
write_tsv(pairwise_rows, file.path(OUT, "two_cohort_disease_state_pairwise_tests.tsv"))
write_tsv(external_performance, file.path(TABLE_MAIN, "Table_5_Ridge_nomogram_external_validation.tsv"))
write_tsv(external_audit, file.path(TABLE_SUPP, "Table_S15_external_tissue_cohort_audit.tsv"))
write_tsv(external_performance, file.path(TABLE_SUPP, "Table_S16_Ridge_external_calibration_HL.tsv"))
write_tsv(challenge_table, file.path(TABLE_SUPP, "Table_S17_two_cohort_disease_state_challenge.tsv"))
write_tsv(point_scale, file.path(SOURCE, "Figure_3A_Ridge_nomogram_point_scale.tsv"))
write_tsv(external_curves, file.path(SOURCE, "Figure_3B_external_smooth_calibration.tsv"))
write_tsv(external_rocs, file.path(SOURCE, "Figure_3C_external_ROC_source.tsv"))
write_tsv(challenge_frame[, c("dataset", "sample_id", "title", "group", "ridge_nomogram_total_points")],
          file.path(SOURCE, "Figure_3D_two_cohort_disease_state_score_source.tsv"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

cat("\nRidge coefficients\n"); print(coefficient_table)
cat("\nExternal performance and Hosmer-Lemeshow tests\n"); print(external_performance)
cat("\nDisease-state challenge\n"); print(challenge_summary); print(pairwise_rows)
