#!/usr/bin/env Rscript
# ==============================================================================
# 32a_build_cillo.R  --  assemble GSE139324 (Cillo 2020) into a Seurat object
#
#   Rscript 32a_build_cillo.R
#
# GEO gives flat, GSM-prefixed 10x triplets:
#   GSM4138110_HNSCC_1_PBMC_{genes,barcodes}.tsv.gz + _matrix.mtx.gz
# Read10X cannot read that layout (prefixes + flat dir + old-style 'genes'), so
# we pair the three files per sample by hand, read each matrix, tag every cell
# with sample / patient / tissue parsed from the filename, merge, normalize, and
# assign broad lineages from canonical markers.
#
# TISSUE for the tumor-vs-normal test (C4): the right normal control for a tumor
# microenvironment claim is HD_Tonsil (tissue-resident, oropharyngeal), NOT
# HD_PBMC (blood myeloid are a different population). So C4 = TIL vs HD_Tonsil.
#
# Output -> data/independent/cillo_GSE139324.rds
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(stringr)
})
say <- function(...) message(sprintf(...))

DIR <- paste0(Sys.getenv("SCFEA_ROOT", unset = "/Users/Tommy/dechenli_lab_local/scfea26"), "/data/independent/GSE139324_RAW")
OUT <- paste0(Sys.getenv("SCFEA_ROOT", unset = "/Users/Tommy/dechenli_lab_local/scfea26"), "/data/independent/cillo_GSE139324.rds")

# ---- discover samples from the matrix files ----------------------------------
mtx <- list.files(DIR, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE)
say("found %d matrix files", length(mtx))
stopifnot(length(mtx) > 0)

sample_of <- function(p) sub("_matrix\\.mtx\\.gz$", "", basename(p))  # GSM..._HNSCC_1_PBMC
samples <- sample_of(mtx)

read_one <- function(stub) {
  mf <- file.path(DIR, paste0(stub, "_matrix.mtx.gz"))
  bf <- file.path(DIR, paste0(stub, "_barcodes.tsv.gz"))
  gf <- file.path(DIR, paste0(stub, "_genes.tsv.gz"))
  if (!all(file.exists(mf, bf, gf))) { say("  [skip] incomplete triplet: %s", stub); return(NULL) }
  m  <- readMM(mf)
  bc <- readLines(bf)
  gn <- read.delim(gf, header = FALSE, stringsAsFactors = FALSE)
  # genes.tsv is usually 2 cols (ensembl, symbol); use symbol (col 2) if present
  sym <- if (ncol(gn) >= 2) gn[[2]] else gn[[1]]
  rownames(m) <- make.unique(sym)
  colnames(m) <- paste0(stub, "|", bc)
  m
}

# parse metadata from the stub: GSM4138110_HNSCC_1_PBMC  /  GSM4138168_HD_1_Tonsil
parse_meta <- function(stub) {
  s <- sub("^GSM[0-9]+_", "", stub)            # HNSCC_1_PBMC | HD_1_Tonsil | HD_1_PBMC
  parts <- strsplit(s, "_")[[1]]
  cohort <- parts[1]                            # HNSCC or HD
  pid    <- paste(parts[1], parts[2], sep = "_")# HNSCC_1 / HD_1
  tis    <- tail(parts, 1)                      # PBMC / TIL / Tonsil
  tissue <- dplyr::case_when(
    cohort == "HNSCC" & tis == "TIL"  ~ "tumor_TIL",
    cohort == "HNSCC" & tis == "PBMC" ~ "tumor_PBMC",
    cohort == "HD"    & tis == "Tonsil" ~ "normal_tonsil",
    cohort == "HD"    & tis == "PBMC"  ~ "normal_PBMC",
    TRUE ~ "other")
  list(patient = pid, sample = s, tissue = tissue, cohort = cohort)
}

# ---- read + merge (memory-conscious: build list, then merge once) -----------
objs <- list()
for (i in seq_along(samples)) {
  st <- samples[i]
  m <- read_one(st); if (is.null(m)) next
  md <- parse_meta(st)
  o <- CreateSeuratObject(counts = m, project = md$sample, min.cells = 3, min.features = 200)
  o$patient <- md$patient; o$sample <- md$sample
  o$tissue  <- md$tissue;  o$cohort <- md$cohort
  objs[[st]] <- o
  if (i %% 10 == 0) say("  read %d / %d samples", i, length(samples))
}
say("merging %d samples...", length(objs))
obj <- if (length(objs) == 1) objs[[1]] else
       merge(objs[[1]], y = objs[-1], add.cell.ids = NULL)
rm(objs); gc()
say("merged: %d genes x %d cells", nrow(obj), ncol(obj))

# ---- basic QC + normalize ----------------------------------------------------
obj[["pct_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & pct_mt < 20)
say("after QC: %d cells", ncol(obj))
obj <- NormalizeData(obj, verbose = FALSE)
# Seurat v5 keeps one layer per merged sample; join so GetAssayData spans all cells
obj <- tryCatch(JoinLayers(obj), error = function(e) { say("  (JoinLayers skipped: %s)", conditionMessage(e)); obj })

# ---- lineage from canonical markers (marker-score, no clustering needed) -----
# CD45-sorted immune dataset, so mostly T/NK/B/myeloid; still score broadly.
sig <- list(
  Myeloid    = c("LYZ","CD68","CD14","FCGR3A","ITGAX","C1QA","C1QB","CLEC10A","S100A8"),
  T          = c("CD3D","CD3E","CD8A","CD4","IL7R","TRAC"),
  NK         = c("NKG7","GNLY","KLRD1","NCAM1"),
  B          = c("CD79A","CD79B","MS4A1","IGHM"),
  Plasma     = c("MZB1","IGHG1","XBP1","JCHAIN"),
  Epithelial = c("KRT5","KRT14","EPCAM","SFN"))
sig <- lapply(sig, function(g) intersect(g, rownames(obj)))
E <- GetAssayData(obj, layer = "data")
score <- sapply(sig, function(g) if (length(g)) Matrix::colMeans(E[g,,drop=FALSE]) else rep(0, ncol(obj)))
obj$lineage <- colnames(score)[max.col(score, ties.method = "first")]
# collapse NK into T/NK for the myeloid-vs-lymphoid contrast
obj$lineage_broad <- ifelse(obj$lineage %in% c("T","NK"), "T/NK", obj$lineage)

say("\nlineage assignment:")
print(table(obj$lineage_broad))
say("\ntissue x lineage (myeloid coverage is what C1/C2 need):")
print(table(obj$tissue, obj$lineage_broad))

# quick sanity: MCT genes present?
for (g in c("SLC16A3","SLC16A1","SLC16A4","SPP1","MMP11","CD68"))
  say("  %-9s %s", g, ifelse(g %in% rownames(obj), "present", "MISSING"))

saveRDS(obj, OUT)
say("\n-> %s", OUT)
say("\nNEXT:")
say("  Rscript 32_independent_cohort.R --obj %s \\", OUT)
say("      --name cillo_GSE139324 --celltype lineage_broad --sample sample --tissue tissue")
