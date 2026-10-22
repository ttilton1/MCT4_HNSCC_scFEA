#!/usr/bin/env Rscript
# ==============================================================================
# 01_export_bulk_after.R  --  Run: all_sites/bulk_after  ( mean(flux) )
#
# Pseudobulk AFTER flux. This script only builds the TRAINING subsample; the
# predict step reads the full shared imputed matrix directly (no re-export).
#
# Per major type, from the shared imputed matrix (661 module genes x cells):
#   - if n_cells <= TRAIN_N (10000): use ALL cells (small types: Mast, Plasma)
#   - else: stratify by (tissue, hpv), allocate 10K proportionally across strata,
#     and within each stratum take the TOP-UMI cells (Chi's high-UMI training
#     selection; deterministic, no random seed). UMI is RAW (from cell_metadata),
#     but the values written are the cells' IMPUTED profiles, so train and
#     predict share one imputation and scFEA runs --sc_imputation False in both.
#
# Strata key = (tissue, hpv): tissue subsumes site (each tissue nests under one
# site) and keeps M vs R and the normal subtypes separate; site_group is dropped
# (redundant with site, hence with tissue). Decided from misc_inspect_metadata.
#
# Output (per type):
#   data/input/all_sites/bulk_after/<Type>/train_expression.csv   661 genes x <=10K cells
#   data/metadata/all_sites/bulk_after/<Type>/metadata_train.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

BASE        <- "/project2/dechenli_1477/ttilton/spring26/scfea"
SHARED_BASE <- file.path(BASE, "data/input/all_sites/shared")
INPUT_OUT   <- file.path(BASE, "data/input/all_sites/bulk_after")
META_OUT    <- file.path(BASE, "data/metadata/all_sites/bulk_after")
MODULE_PATH <- file.path(BASE, "scFEA_src/data/module_gene_m168.csv")
TYPES   <- c("T","Myeloid","Fibroblast","Epithelial","Endothelial","B","NK","Mast","Plasma")
TRAIN_N <- 10000
STRAT   <- c("tissue","hpv")

# fread mangles column 1 of the wide 244K-col imputed CSV, so reconstruct gene
# names from 00a's genes.txt (atlas order) intersected with module genes -- the
# exact row order of the imputed matrix (verified against a clean-read type).
scfea_module_genes <- unique(na.omit(as.character(unlist(
  read.csv(MODULE_PATH, row.names = 1)))))
scfea_module_genes <- scfea_module_genes[scfea_module_genes != "" & scfea_module_genes != "nan"]
clean_gene_names <- function(ct) {
  atlas <- readLines(file.path(SHARED_BASE, ct, "genes.txt"))
  atlas[atlas %in% scfea_module_genes]
}

# ---- proportional allocation + top-UMI within stratum ----
select_train <- function(meta, n_target, strat_cols, umi_col = "umi") {
  if (nrow(meta) <= n_target) { cat("    using ALL", nrow(meta), "cells (<= TRAIN_N)\n"); return(rownames(meta)) }
  meta$.strat <- apply(meta[, strat_cols, drop = FALSE], 1, paste, collapse = "_")
  gc_counts <- table(meta$.strat)
  props <- gc_counts / sum(gc_counts)
  alloc <- round(props * n_target); alloc[alloc == 0] <- 1
  while (sum(alloc) > n_target) { b <- which.max(alloc);              alloc[b] <- alloc[b] - 1 }
  while (sum(alloc) < n_target) { b <- which.max(gc_counts - alloc);  alloc[b] <- alloc[b] + 1 }
  sel <- character(0)
  for (g in names(alloc)) {
    cells <- rownames(meta[meta$.strat == g, , drop = FALSE])
    ord   <- cells[order(-meta[cells, umi_col])]            # top-UMI first
    sel   <- c(sel, head(ord, min(alloc[g], length(cells))))
  }
  sel
}

cat("============================================\n")
cat("bulk_after export: top-UMI-within-(tissue,hpv) training subsample\n")
cat("Started:", format(Sys.time()), " | TRAIN_N =", TRAIN_N, "\n")
cat("============================================\n\n")

sel_summary <- data.frame()

for (ct in TYPES) {
  cat("--------------------------------------------\n[", ct, "]\n")
  imp_path  <- file.path(SHARED_BASE, ct, "imputed_expression.csv")
  meta_path <- file.path(SHARED_BASE, ct, "cell_metadata.csv")
  if (!file.exists(imp_path) || !file.exists(meta_path)) { cat("  MISSING shared input; skipping.\n\n"); next }

  tall_path <- file.path(SHARED_BASE, ct, "imputed_tall.csv")
  if (!file.exists(tall_path)) stop(ct, ": missing ", tall_path, " (run 00c_transpose_imputed first)")
  tt    <- fread(tall_path)
  bc    <- as.character(tt[[1]])                # barcodes
  Mt    <- as.matrix(tt[, -1])                  # cells x genes
  genes <- colnames(Mt)                         # clean gene names from tall header
  M     <- t(Mt); colnames(M) <- bc; rm(tt, Mt) # -> genes x cells
  cg <- clean_gene_names(ct)
  if (!identical(rownames(M), cg))
    stop(ct, ": tall gene header does not match genes.txt-intersect-module order")
  if (max(M) > 100 || min(M) < -1e-6)
    stop(ct, ": imputed values out of expected range (max=", max(M), ") -- bad read")

  meta <- as.data.frame(fread(meta_path)); setnames(meta, 1, "barcode")
  rownames(meta) <- meta$barcode
  meta <- meta[bc, , drop = FALSE]
  if (any(is.na(meta$umi))) stop("barcode mismatch / missing umi for ", ct)

  cat(sprintf("  %d cells; strata by (%s)\n", ncol(M), paste(STRAT, collapse = ", ")))
  sel <- select_train(meta, TRAIN_N, STRAT)
  cat(sprintf("  selected %d training cells\n", length(sel)))

  idir <- file.path(INPUT_OUT, ct); mdir <- file.path(META_OUT, ct)
  dir.create(idir, recursive = TRUE, showWarnings = FALSE)
  dir.create(mdir, recursive = TRUE, showWarnings = FALSE)
  train_path <- file.path(idir, "train_expression.csv")
  data.table::fwrite(data.table(gene = rownames(M),
                                as.data.frame(M[, sel, drop = FALSE], check.names = FALSE)),
                     train_path)
  .raw <- readBin(train_path, "raw", n = file.info(train_path)$size)
  if (any(as.integer(.raw) > 127L))
    stop("non-ASCII byte in ", train_path, " -- write corruption; try a different node/path")
  rm(.raw)
  write.csv(meta[sel, , drop = FALSE], file.path(mdir, "metadata_train.csv"))
  cat(sprintf("  wrote train_expression.csv (661 x %d)\n", length(sel)))

  used_all <- ncol(M) <= TRAIN_N
  sel_summary <- rbind(sel_summary, data.frame(
    major_type = ct, n_cells = ncol(M), n_train = length(sel),
    used_all_cells = used_all,
    n_strata = if (used_all) NA else length(unique(apply(meta[, STRAT], 1, paste, collapse = "_")))))

  rm(M, meta); gc(verbose = FALSE)
  cat("  done.\n\n")
}

cat("============================================\nTRAINING SELECTION SUMMARY:\n")
print(sel_summary, row.names = FALSE)
dir.create(INPUT_OUT, recursive = TRUE, showWarnings = FALSE)
write.csv(sel_summary, file.path(INPUT_OUT, "_train_selection_summary.csv"), row.names = FALSE)
cat("\nNext: 02_train_bulk_after.slurm -> 03_predict_bulk_after*.slurm\n")
cat("Predict reads the FULL shared imputed matrix directly (not re-exported here).\n")
cat("============================================\n")

