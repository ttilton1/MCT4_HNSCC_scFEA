#!/usr/bin/env Rscript
# ==============================================================================
# 01_export_bulk_before.R  --  Run: all_sites/bulk_before  ( flux(mean) )
#
# Pseudobulk BEFORE flux. For each major type, read the SHARED imputed matrix
# (data/input/all_sites/shared/<Type>/imputed_expression.csv, 661 module genes x
# cells) and collapse to (sample x finaltype) MEANS. Each pseudobulk unit is one
# "cell" for scFEA. Cell set + imputation are identical to bulk_after (both read
# the same shared matrix), so the ONLY difference between the two runs is when
# pseudobulk happens relative to flux.
#
# Input  (per type): data/input/all_sites/shared/<Type>/
#          imputed_expression.csv   661 genes x cells (from 00b MAGIC)
#          cell_metadata.csv        per-cell covariates incl umi (from 00a)
# Output (per type):
#          data/input/all_sites/bulk_before/<Type>/train_expression.csv
#              661 genes x units  (scFEA convention: genes=rows, units=cols)
#          data/metadata/all_sites/bulk_before/<Type>/metadata_units.csv
#              one row per unit: sample, finaltype, patient, study, tissue, site,
#              site_group, hpv, annot_2nd, annot_4th, n_cells, median_umi
#
# scFEA then trains == predicts on this matrix (02_scfea_bulk_before.slurm,
# --mode train --sc_imputation False): the training-byproduct flux IS the result,
# one row per (sample x finaltype) unit.
#
# Light + fast (imputed matrices are small; T's is 2.6 GB, read with fread).
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)   # fread: fast read of the imputed CSVs (esp. T, 2.6 GB)
  library(Matrix)
})

BASE        <- "/project2/dechenli_1477/ttilton/spring26/scfea"
SHARED_BASE <- file.path(BASE, "data/input/all_sites/shared")
INPUT_OUT   <- file.path(BASE, "data/input/all_sites/bulk_before")
META_OUT    <- file.path(BASE, "data/metadata/all_sites/bulk_before")
MODULE_PATH <- file.path(BASE, "scFEA_src/data/module_gene_m168.csv")
TYPES <- c("T","Myeloid","Fibroblast","Epithelial","Endothelial","B","NK","Mast","Plasma")

COVARS <- c("sample","finaltype","patient","study","tissue","site","site_group",
            "hpv","annot_2nd","annot_4th")

# ---- authoritative gene names ----------------------------------------------
# fread mangles column 1 of the wide 244K-col imputed CSV (multithread chunk
# reassembly + a couple of decode faults), so we do NOT trust the CSV's gene
# column. The imputed matrix rows are atlas-order intersected with module genes;
# 00a's per-type genes.txt gives atlas order, so this reconstructs the exact row
# names (verified: identical to a clean-read type's names).
scfea_module_genes <- unique(na.omit(as.character(unlist(
  read.csv(MODULE_PATH, row.names = 1)))))
scfea_module_genes <- scfea_module_genes[scfea_module_genes != "" & scfea_module_genes != "nan"]
clean_gene_names <- function(ct) {
  atlas <- readLines(file.path(SHARED_BASE, ct, "genes.txt"))
  atlas[atlas %in% scfea_module_genes]
}

cat("============================================\n")
cat("bulk_before export: pseudobulk (sample x finaltype) of the shared imputed matrix\n")
cat("Started:", format(Sys.time()), "\n")
cat("============================================\n\n")

unit_summary <- data.frame()

for (ct in TYPES) {
  cat("--------------------------------------------\n")
  cat("[", ct, "]\n")

  imp_path  <- file.path(SHARED_BASE, ct, "imputed_expression.csv")
  meta_path <- file.path(SHARED_BASE, ct, "cell_metadata.csv")
  if (!file.exists(imp_path) || !file.exists(meta_path)) {
    cat("  MISSING shared input; skipping.\n\n"); next
  }

  # read the TALL imputed (cells x genes) written by 00c -- fread mis-parses the
  # 244K-wide genes x cells CSV, but reads the tall shape reliably.
  tall_path <- file.path(SHARED_BASE, ct, "imputed_tall.csv")
  if (!file.exists(tall_path)) stop(ct, ": missing ", tall_path, " (run 00c_transpose_imputed first)")
  tt    <- fread(tall_path)
  bc    <- as.character(tt[[1]])                # barcodes (tall row index)
  Mt    <- as.matrix(tt[, -1])                  # cells x genes
  genes <- colnames(Mt)                         # clean gene names from tall header
  M     <- t(Mt); colnames(M) <- bc; rm(tt, Mt) # -> genes x cells
  # cross-check against the authoritative gene list (order + count)
  cg <- clean_gene_names(ct)
  if (!identical(rownames(M), cg))
    stop(ct, ": tall gene header does not match genes.txt-intersect-module order")
  if (max(M) > 100 || min(M) < -1e-6)
    stop(ct, ": imputed values out of expected range [0,~5] (max=", max(M), ") -- bad read")

  # metadata: rows = barcodes
  meta <- fread(meta_path)
  setnames(meta, 1, "barcode")
  meta <- as.data.frame(meta)
  rownames(meta) <- meta$barcode
  meta <- meta[bc, , drop = FALSE]      # align to matrix column order
  if (any(is.na(meta$sample))) stop("barcode mismatch between imputed matrix and cell_metadata for ", ct)

  cat(sprintf("  %d genes x %d cells\n", nrow(M), ncol(M)))

  # ---- group key = sample x finaltype ----
  grp <- factor(paste(meta$sample, meta$finaltype, sep = "__"))
  n_units <- nlevels(grp)
  cat("  pseudobulk units:", n_units, "\n")

  # column means per unit via sparse indicator: (genes x cells) %*% (cells x units)
  ind <- sparseMatrix(i = seq_along(grp), j = as.integer(grp), x = 1,
                      dims = c(length(grp), n_units), dimnames = list(NULL, levels(grp)))
  sizes  <- Matrix::colSums(ind)
  pb     <- (M %*% ind)                 # genes x units (sum of imputed values)
  pb     <- sweep(as.matrix(pb), 2, sizes, "/")   # -> means
  rownames(pb) <- genes
  colnames(pb) <- levels(grp)

  # ---- per-unit metadata (covariates constant within sample x finaltype) ----
  first_idx <- match(levels(grp), grp)
  mu <- meta[first_idx, COVARS, drop = FALSE]
  mu$n_cells    <- as.integer(sizes)
  mu$median_umi <- as.numeric(tapply(meta$umi, grp, median))
  rownames(mu)  <- levels(grp)

  for (cv in c("site","hpv","patient","study","site_group")) {
    nu <- tapply(meta[[cv]], grp, function(z) length(unique(z)))
    if (any(nu > 1, na.rm = TRUE))
      cat(sprintf("  WARNING: %d units have >1 %s value (unexpected)\n", sum(nu > 1), cv))
  }
  cat(sprintf("  n_cells/unit: min=%d median=%d max=%d  (downstream gates at 8/15)\n",
              min(sizes), as.integer(median(sizes)), max(sizes)))

  # ---- write ----
  idir <- file.path(INPUT_OUT, ct); mdir <- file.path(META_OUT, ct)
  dir.create(idir, recursive = TRUE, showWarnings = FALSE)
  dir.create(mdir, recursive = TRUE, showWarnings = FALSE)
  out_expr <- file.path(idir, "train_expression.csv")
  # fwrite (data.table C writer) -- write.csv corrupted T's 30MB matrix (byte-stream
  # buffering bug on large output). fwrite is faster and clean. Guard: read back and
  # fail loudly if any non-ASCII byte slipped in, rather than shipping it downstream.
  data.table::fwrite(data.table(gene = rownames(pb),
                                as.data.frame(pb, check.names = FALSE)),
                     out_expr)
  .raw <- readBin(out_expr, "raw", n = file.info(out_expr)$size)
  if (any(as.integer(.raw) > 127L))
    stop("non-ASCII byte in ", out_expr, " -- write corruption; try a different node/path")
  rm(.raw)
  write.csv(mu, file.path(mdir, "metadata_units.csv"))
  cat(sprintf("  wrote %s (%d genes x %d units, %.1f MB)\n",
              out_expr, nrow(pb), ncol(pb), file.info(out_expr)$size / 1e6))

  unit_summary <- rbind(unit_summary, data.frame(
    major_type = ct, n_cells = ncol(M), n_units = n_units,
    min_cells_per_unit = min(sizes), med_cells_per_unit = as.integer(median(sizes)),
    max_cells_per_unit = max(sizes)))

  rm(M, ind, pb, meta); gc(verbose = FALSE)
  cat("  done.\n\n")
}

cat("============================================\n")
cat("UNIT SUMMARY:\n"); print(unit_summary, row.names = FALSE)
dir.create(INPUT_OUT, recursive = TRUE, showWarnings = FALSE)
write.csv(unit_summary, file.path(INPUT_OUT, "_unit_summary.csv"), row.names = FALSE)
cat("\nbulk_before export complete:", format(Sys.time()), "\n")
cat("Next: 02_scfea_bulk_before.slurm (scFEA --mode train --sc_imputation False)\n")
cat("============================================\n")

