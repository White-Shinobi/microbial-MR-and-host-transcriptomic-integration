#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(igraph)
})

root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[1]), ".."),
                      mustWork = FALSE)
# Rscript does not expose the script path consistently through commandArgs()[1].
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(file_arg) == 1L) {
  root <- normalizePath(file.path(dirname(file_arg), ".."), mustWork = TRUE)
}

out <- file.path(root, "results", "single_cell", "GSE178341",
                 "label_free_atlas")
edge_path <- file.path(out, "label_free_sketch_umap_graph_edges.tsv.gz")
meta_path <- file.path(out, "label_free_sketch_metadata.tsv.gz")

edges <- read.delim(edge_path, check.names = FALSE)
meta <- read.delim(meta_path, check.names = FALSE)
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

resolutions <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00)
set.seed(20260723)
memberships <- data.frame(sketch_row = 0:(n - 1L))
summary_rows <- vector("list", length(resolutions))

for (i in seq_along(resolutions)) {
  resolution <- resolutions[[i]]
  fit <- cluster_leiden(
    graph,
    objective_function = "modularity",
    weights = E(graph)$weight,
    resolution_parameter = resolution,
    n_iterations = 10
  )
  membership_i <- membership(fit)
  label <- sprintf("leiden_r%.2f", resolution)
  memberships[[label]] <- as.integer(membership_i)
  sizes <- table(membership_i)
  summary_rows[[i]] <- data.frame(
    resolution = resolution,
    clusters = length(sizes),
    smallest_cluster = min(sizes),
    median_cluster = median(as.numeric(sizes)),
    largest_cluster = max(sizes),
    modularity = modularity(
      graph,
      membership_i,
      weights = E(graph)$weight,
      resolution = resolution
    )
  )
}

membership_connection <- gzfile(
  file.path(out, "label_free_sketch_leiden_memberships.tsv.gz"),
  open = "wt"
)
write.table(
  memberships,
  membership_connection,
  sep = "\t", quote = FALSE, row.names = FALSE
)
close(membership_connection)
write.table(
  do.call(rbind, summary_rows),
  file.path(out, "label_free_sketch_leiden_resolution_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat(sprintf("Completed Leiden grid for %s sketch cells\n",
            format(n, big.mark = ",")))
