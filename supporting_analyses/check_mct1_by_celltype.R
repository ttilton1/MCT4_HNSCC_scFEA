#!/usr/bin/env Rscript
# ==============================================================================
# check_mct1_by_celltype.R  --  is MCT1 concentrated on epithelial cells?
#
#   Rscript check_mct1_by_celltype.R --arm bulk_before
#
# The glioma study of Zhu et al. reports SLC16A1 on malignant cells and SLC16A3
# on immune cells. We report the MCT4 half; this asks whether our data also shows
# the MCT1 half, so the discussion can either claim the full parallel or state
# which part we cannot confirm.
#
# Expression is comparable between cell types in a way flux is not: the
# underlying matrix is CP10K normalized across all cells before imputation,
# whereas scFEA normalizes flux to the expression mean within each run.
#
# Carcinoma units only, since that is where the comparison is meaningful.
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(data.table); library(dplyr) })

CANCER_T <- c("CA","CA-E","CA-A","CA-c","CA-p")
CTS <- c("Myeloid","Epithelial","Fibroblast","Endothelial","B","T","NK","Mast","Plasma")

res <- rbindlist(lapply(CTS, function(ct) {
  ep <- file.path(ROOT, "data/input/all_sites", ARM, ct, "train_expression.csv")
  mp <- file.path(ROOT, "metadata/all_sites", ARM, ct, "metadata_units.csv")
  if (!file.exists(ep) || !file.exists(mp)) return(NULL)
  e <- fread(ep); gn <- as.character(e[[1]])
  if (!all(c("SLC16A3","SLC16A1") %in% gn)) return(NULL)
  E <- as.matrix(e[, -1]); rownames(E) <- gn
  mu <- fread(mp); setnames(mu, 1, "id"); mu <- as.data.frame(mu)
  rownames(mu) <- as.character(mu$id); mu <- mu[colnames(E), ]
  k <- mu$tissue %in% CANCER_T & mu$n_cells >= GATES$min_cells
  if (!sum(k)) return(NULL)
  data.table(celltype = ct, units = sum(k),
             MCT1 = median(E["SLC16A1", k]),
             MCT4 = median(E["SLC16A3", k]))
}))
res[, `:=`(ratio = MCT4 / MCT1,
           mct1_share = 100 * MCT1 / (MCT1 + MCT4))]

say("=== transporter expression by cell type, carcinoma units ===")
print(as.data.frame(res[order(-MCT1)] %>%
  mutate(across(c(MCT1, MCT4, ratio, mct1_share), ~round(.x, 3)))), row.names = FALSE)

say("\n  MCT1 is highest in: %s", res[which.max(MCT1)]$celltype)
say("  MCT4 is highest in: %s", res[which.max(MCT4)]$celltype)
say("  epithelial MCT1 rank: %d of %d cell types",
    which(res[order(-MCT1)]$celltype == "Epithelial"), nrow(res))
say("\n  If epithelial ranks first or near first on MCT1, our data shows the same")
say("  split the glioma study reports. If not, we report only the MCT4 half.")
