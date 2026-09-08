#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(igraph))
args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep("^--file=", args_all)][1]
root <- normalizePath(
  file.path(dirname(sub("^--file=", "", script_arg)), ".."),
  mustWork = TRUE
)
out <- file.path(
  root, "results", "single_cell", "GSE178341",
  "label_free_epithelial"
)
edges <- read.delim(
  file.path(out, "epithelial_sketch_graph_edges.tsv.gz")
)
meta <- read.delim(
  file.path(out, "epithelial_label_free_sketch_metadata.tsv.gz")
)
n <- nrow(meta)
graph <- graph_from_data_frame(
  data.frame(
    source = as.character(edges$source + 1L),
    target = as.character(edges$target + 1L),
    weight = edges$weight
  ),
  directed = FALSE,
  vertices = data.frame(name = as.character(seq_len(n)))
)
resolutions <- c(0.05, 0.08, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50)
set.seed(20260723)
memberships <- data.frame(sketch_row = 0:(n - 1L))
audit <- list()
for (resolution in resolutions) {
  fit <- cluster_leiden(
    graph,
    objective_function = "modularity",
    weights = E(graph)$weight,
    resolution_parameter = resolution,
    n_iterations = 10
  )
  membership_i <- membership(fit)
  memberships[[sprintf("leiden_r%.2f", resolution)]] <-
    as.integer(membership_i)
  sizes <- table(membership_i)
  audit[[length(audit) + 1L]] <- data.frame(
    resolution = resolution,
    clusters = length(sizes),
    smallest_cluster = min(sizes),
    median_cluster = median(as.numeric(sizes)),
    largest_cluster = max(sizes),
    modularity = modularity(
      graph, membership_i, weights = E(graph)$weight,
      resolution = resolution
    )
  )
}
connection <- gzfile(
  file.path(out, "epithelial_leiden_memberships.tsv.gz"), "wt"
)
write.table(
  memberships, connection, sep = "\t",
  quote = FALSE, row.names = FALSE
)
close(connection)
write.table(
  do.call(rbind, audit),
  file.path(out, "epithelial_leiden_resolution_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("Completed epithelial Leiden grid for", n, "cells\n")
