#!/usr/bin/env Rscript
# Phase 2: subsample a Seurat object and export raw counts in a
# language-neutral form (MatrixMarket + CSV) for assembly into .h5ad.
#
# Usage:
#   Rscript 02_export.R <rds> <outdir> <donor_col> <celltype_col> <max_per_group> <assay>
#
# Design notes:
#  - Exports the RNA assay counts layer. NOT SCT: the SCT "counts" layer holds
#    sctransform-corrected counts (depth-normalized), which is not valid
#    Geneformer tokenizer input. Geneformer needs true raw UMI counts.
#  - Subsampling happens BEFORE export so we never write a 100 GB matrix.
#  - MatrixMarket is used instead of SeuratDisk/sceasy because those break on
#    Seurat v5 multi-layer assays; writeMM never does.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(jsonlite)
})

args        <- commandArgs(trailingOnly = TRUE)
rds_path    <- args[1]
out_dir     <- args[2]
donor_col   <- args[3]
ct_col      <- args[4]
max_per_grp <- as.integer(args[5])
assay       <- if (length(args) >= 6) args[6] else "RNA"

set.seed(42)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("reading:", rds_path, "\n")
t0 <- Sys.time()
obj <- readRDS(rds_path)
cat("read in", format(Sys.time() - t0), "-", ncol(obj), "cells\n")

md <- obj@meta.data
stopifnot(donor_col %in% colnames(md), ct_col %in% colnames(md))
stopifnot(assay %in% Assays(obj))

# --- raw counts from the RNA assay (v5 layer API, v4 slot API fallback) ---
cnt <- tryCatch(
  SeuratObject::LayerData(obj, assay = assay, layer = "counts"),
  error = function(e) GetAssayData(obj, assay = assay, slot = "counts")
)
cat("assay:", assay, " counts:", paste(dim(cnt), collapse = " x "),
    " class", class(cnt)[1], "\n")
stopifnot(ncol(cnt) == nrow(md))

# --- stratified subsample: <= max_per_grp cells per donor x celltype ---
md$.grp <- paste(md[[donor_col]], md[[ct_col]], sep = "||")
keep <- unlist(lapply(split(rownames(md), md$.grp), function(ix) {
  if (length(ix) <= max_per_grp) ix else sample(ix, max_per_grp)
}), use.names = FALSE)
keep <- keep[!is.na(keep)]
cat("selected", length(keep), "of", nrow(md), "cells across",
    length(unique(md$.grp)), "donor x celltype groups\n")

cnt_sub <- cnt[, keep, drop = FALSE]
md_sub  <- md[keep, , drop = FALSE]
md_sub$.grp <- NULL

# --- drop genes that are all-zero after subsampling ---
nz <- Matrix::rowSums(cnt_sub) > 0
cat("genes retained:", sum(nz), "of", length(nz), "\n")
cnt_sub <- cnt_sub[nz, , drop = FALSE]

# --- integer check: SoupX correction can leave fractional values ---
vals <- if (inherits(cnt_sub, "sparseMatrix")) cnt_sub@x else as.numeric(cnt_sub)
is_int <- isTRUE(all.equal(vals, round(vals), tolerance = 1e-8))
cat("counts integer-like:", is_int, " range:", paste(range(vals), collapse = " .. "), "\n")
if (!is_int) {
  cat("NOTE: fractional counts (SoupX corrected); rounding for tokenizer.\n")
  if (inherits(cnt_sub, "sparseMatrix")) cnt_sub@x <- round(cnt_sub@x) else cnt_sub <- round(cnt_sub)
  cnt_sub <- Matrix::drop0(cnt_sub)
}

# --- write ---
cat("writing to", out_dir, "\n")
Matrix::writeMM(as(cnt_sub, "dgCMatrix"), file.path(out_dir, "counts.mtx"))
write.csv(data.frame(gene = rownames(cnt_sub)),
          file.path(out_dir, "genes.csv"), row.names = FALSE)
write.csv(data.frame(cell = colnames(cnt_sub)),
          file.path(out_dir, "cells.csv"), row.names = FALSE)
write.csv(md_sub, file.path(out_dir, "obs.csv"))

meta <- list(
  source_rds    = rds_path,
  assay         = assay,
  donor_col     = donor_col,
  celltype_col  = ct_col,
  max_per_group = max_per_grp,
  n_cells       = ncol(cnt_sub),
  n_genes       = nrow(cnt_sub),
  n_donors      = length(unique(md_sub[[donor_col]])),
  n_celltypes   = length(unique(md_sub[[ct_col]])),
  counts_were_integer = is_int,
  organism_guess = if (mean(grepl("^[A-Z][a-z]", rownames(cnt_sub))) > 0.7) "mouse (symbol case)" else "unknown",
  seurat_version = as.character(obj@version),
  exported_at   = as.character(Sys.time())
)
writeLines(toJSON(meta, auto_unbox = TRUE, pretty = TRUE),
           file.path(out_dir, "export_meta.json"))

cat("=== DONE ===", ncol(cnt_sub), "cells x", nrow(cnt_sub), "genes\n")
