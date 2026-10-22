#!/usr/bin/env Rscript
# ==============================================================================
# 00a_export_normalized_sparse.R  --  SHARED MAGIC input export
#
# Per major type (annot_4th; Oral+LN+MR; zero-metabolic-gene cells removed),
# write the FULL normalized (data-layer) matrix as sparse .mtx so MAGIC can be
# run ONCE on the full cell set (00b). The single imputed matrix feeds BOTH
# runs -- all_sites_bulk_after (Run A) and all_sites_bulk_first (Run B) -- so
# imputation is held identical and the only A-vs-B difference is pseudobulk
# timing.
#
# Also writes raw total-UMI per cell: Run A selects its training subsample by
# top-UMI within strata, and UMI MUST come from RAW counts (imputation destroys
# count depth). The imputed profiles of the selected cells are what train scFEA.
#
# Sparse .mtx keeps even the T combo (244K cells) compact on disk. All genes are
# written because MAGIC builds its graph on the full transcriptome (matches
# scFEA's internal MAGIC, which imputes all input genes before subsetting).
#
# Output (data/input/all_sites/shared/<MajorType>/):
#   normalized.mtx        genes x cells, log-normalized (MAGIC input)
#   genes.txt             gene names, row order of the mtx
#   barcodes.txt          cell barcodes, column order of the mtx
#   cell_metadata.csv     per-cell: annot_2nd, annot_4th, finaltype, site,
#                         site_group, hpv, study, tissue, patient, sample, umi
#                         (order matches barcodes.txt)
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

SEURAT_PATH <- "/project2/dechenli_1477/data/obj.final.2.0.rda"
MODULE_PATH <- "/project2/dechenli_1477/metabolism/scFea_analysis/scFEA/data/module_gene_m168.csv"
OUTPUT_BASE <- "/project2/dechenli_1477/ttilton/spring26/scfea/data/input/all_sites/shared"
KEEP_SITES  <- c("Oral", "LN", "MR")
MIN_CELLS   <- 1000

cat("============================================\n")
cat("00a: export full normalized matrices (sparse) for shared MAGIC\n")
cat("Started:", format(Sys.time()), "\n")
cat("============================================\n\n")

cat("Loading Seurat object...\n")
load(SEURAT_PATH)                       # -> obj
cat("Loaded:", ncol(obj), "cells,", nrow(obj), "genes\n\n")

cat("Checking data layer is log-normalized...\n")
.probe <- GetAssayData(obj, assay = "RNA", layer = "data")[1:min(500, nrow(obj)),
                                                           1:min(50, ncol(obj))]
cat(sprintf("  data-layer probe max=%.3f (expect < ~15)\n\n", max(.probe)))
if (max(.probe) > 50) stop("data layer looks like raw counts (max > 50). Aborting.")
rm(.probe)

module_genes <- read.csv(MODULE_PATH, row.names = 1)
scfea_genes  <- unique(na.omit(as.character(unlist(module_genes))))
scfea_genes  <- scfea_genes[scfea_genes != ""]
met_overlap  <- intersect(scfea_genes, rownames(obj))
cat("scFEA metabolic genes in atlas:", length(met_overlap), "/", length(scfea_genes), "\n\n")

obj$site_group <- ifelse(
  obj@meta.data$site %in% c("CA", "Distal_normal", "Tumor_adjacent_tissue", "Premalignant"), "Oral",
  ifelse(obj@meta.data$site == "LN", "LN",
  ifelse(obj@meta.data$site == "M/R", "MR", "Other")))

# pull layers once, then free the 17GB object
cat("Extracting data layer, metabolic counts, and per-cell UMI, then releasing object...\n")
data_all       <- GetAssayData(obj, assay = "RNA", layer = "data")                 # genes x cells (log-norm)
counts_all     <- GetAssayData(obj, assay = "RNA", layer = "counts")               # genes x cells (raw)
met_sums_all   <- Matrix::colSums(counts_all[met_overlap, ])                        # zero-metabolic filter basis
umi_all        <- Matrix::colSums(counts_all)                                       # total UMI per cell (raw)
meta_all       <- obj@meta.data
gene_names     <- rownames(data_all)
rm(obj, counts_all); gc(verbose = FALSE)

meta_all$met_nonzero <- met_sums_all > 0
meta_all$umi         <- as.numeric(umi_all)

in_scope  <- meta_all$site_group %in% KEEP_SITES
ct_counts <- table(meta_all$annot_4th[in_scope])
viable    <- names(ct_counts[ct_counts >= MIN_CELLS])
viable    <- viable[order(-ct_counts[viable])]
cat("Viable major types:", length(viable), "\n")
print(data.frame(major_type = viable, pooled_cells = as.integer(ct_counts[viable])), row.names = FALSE)
cat("\n")

META_COLS <- c("annot_2nd", "annot_4th", "finaltype", "site", "site_group", "hpv",
               "study", "tissue", "patient", "sample", "umi")

for (i in seq_along(viable)) {
  ct <- viable[i]
  cat("--------------------------------------------\n")
  cat(sprintf("[%d/%d] %s\n", i, length(viable), ct))

  outdir <- file.path(OUTPUT_BASE, ct)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  keep <- in_scope & meta_all$annot_4th == ct & meta_all$met_nonzero
  n_keep <- sum(keep)
  cat("  Cells (site + zero-metabolic filter):", n_keep, "\n")

  mat <- data_all[, keep, drop = FALSE]                # genes x cells, log-norm, ALL genes
  bc  <- colnames(mat)
  md  <- meta_all[keep, META_COLS, drop = FALSE]
  rownames(md) <- bc

  Matrix::writeMM(mat, file.path(outdir, "normalized.mtx"))
  writeLines(gene_names, file.path(outdir, "genes.txt"))
  writeLines(bc,         file.path(outdir, "barcodes.txt"))
  write.csv(md, file.path(outdir, "cell_metadata.csv"))

  gb <- file.info(file.path(outdir, "normalized.mtx"))$size / (1024^3)
  cat(sprintf("  Wrote normalized.mtx (%d genes x %d cells, %.2f GB) + genes/barcodes/metadata\n",
              nrow(mat), ncol(mat), gb))
  cat(sprintf("  UMI per cell: min=%.0f median=%.0f max=%.0f\n",
              min(md$umi), median(md$umi), max(md$umi)))

  rm(mat); gc(verbose = FALSE)
  cat("  Done.\n\n")
}

cat("============================================\n")
cat("00a complete:", format(Sys.time()), "\n")
cat("Next: 00b_magic_impute.py per major type (test on Myeloid/T to size the solver).\n")
cat("============================================\n")

