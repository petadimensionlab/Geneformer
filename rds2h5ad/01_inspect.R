#!/usr/bin/env Rscript
# Phase 1: inspect the Seurat object so we know what we are working with
# before committing to a subsample/export strategy.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

args <- commandArgs(trailingOnly = TRUE)
rds_path <- args[1]

cat("=== reading (this takes a few minutes) ===\n")
cat("file:", rds_path, "\n")
t0 <- Sys.time()
obj <- readRDS(rds_path)
cat("read time:", format(Sys.time() - t0), "\n\n")

cat("=== object ===\n")
print(obj)
cat("\nclass:", class(obj), "\n")
cat("Seurat object version:", as.character(obj@version), "\n")
cat("cells:", ncol(obj), " features:", nrow(obj), "\n\n")

cat("=== assays ===\n")
print(Assays(obj))
cat("default assay:", DefaultAssay(obj), "\n")
for (a in Assays(obj)) {
  cat("\n-- assay:", a, "class:", class(obj[[a]]), "\n")
  layers <- tryCatch(SeuratObject::Layers(obj[[a]]), error = function(e) NA)
  cat("   layers:", paste(layers, collapse = ", "), "\n")
  cat("   dim:", paste(dim(obj[[a]]), collapse = " x "), "\n")
}

cat("\n=== reductions ===\n")
print(Reductions(obj))

cat("\n=== meta.data columns ===\n")
md <- obj@meta.data
cat("n columns:", ncol(md), "\n")
print(data.frame(
  column = colnames(md),
  class  = sapply(md, function(x) class(x)[1]),
  n_uniq = sapply(md, function(x) length(unique(x))),
  example = sapply(md, function(x) paste(utils::head(unique(as.character(x)), 3), collapse = " | ")),
  row.names = NULL
))

cat("\n=== candidate donor/sample columns (n_uniq 2..60, not numeric) ===\n")
cand <- colnames(md)[sapply(md, function(x) {
  u <- length(unique(x))
  u >= 2 && u <= 60 && !is.numeric(x)
})]
for (cc in cand) {
  cat("\n--", cc, "--\n")
  print(utils::head(sort(table(md[[cc]]), decreasing = TRUE), 25))
}

cat("\n=== gene name format (first 20 rownames) ===\n")
print(utils::head(rownames(obj), 20))
cat("\nlooks like Ensembl IDs? ",
    mean(grepl("^ENS[A-Z]*G[0-9]{6,}", rownames(obj))) > 0.5, "\n")

cat("\n=== counts layer sanity ===\n")
cnt <- tryCatch(
  SeuratObject::LayerData(obj, assay = DefaultAssay(obj), layer = "counts"),
  error = function(e) GetAssayData(obj, assay = DefaultAssay(obj), slot = "counts")
)
cat("counts class:", class(cnt), " dim:", paste(dim(cnt), collapse = " x "), "\n")
sub <- cnt[1:min(500, nrow(cnt)), 1:min(500, ncol(cnt))]
vals <- if (inherits(sub, "sparseMatrix")) sub@x else as.numeric(sub)
vals <- vals[vals != 0]
cat("nonzero sampled:", length(vals), "\n")
cat("integer-like: ", isTRUE(all.equal(vals, round(vals), tolerance = 1e-8)), "\n")
cat("range:", if (length(vals)) paste(range(vals), collapse = " .. ") else "NA", "\n")
cat("first 10 values:", paste(utils::head(vals, 10), collapse = ", "), "\n")

cat("\n=== DONE ===\n")
