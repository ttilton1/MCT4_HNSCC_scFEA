#!/usr/bin/env Rscript
# ==============================================================================
# 32b_build_kurten.R  --  assemble GSE164690 (Kurten 2021) into a Seurat object
#
#   Rscript 32b_build_kurten.R
#
# WHY THIS COHORT. Whole-tumor HNSCC with THREE fractions per patient:
#   PBL    = blood         -> DROP (not the tumor microenvironment)
#   CD45p  = immune         -> myeloid + T (re-tests the Cillo result)
#   CD45n  = non-immune     -> FIBROBLASTS + epithelial + endothelial
# The CD45n fraction is the point: it lets us test FIBROBLAST MCT4 specificity,
# which the CD45+-sorted Cillo cohort could not.
#
# HPV. Filenames carry only HN01..HN18, no HPV label. Kurten reported HPV status
# per patient. FILL IN THE MAP BELOW from the paper's Table S1 before trusting the
# HPV split; until then every patient is tagged NA and the HPV test is skipped.
#
# Output -> data/independent/kurten_GSE164690.rds
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(stringr)
})
say <- function(...) message(sprintf(...))

DIR <- "/Users/Tommy/dechenli_lab_local/scfea26/data/independent/GSE164690_RAW"
OUT <- "/Users/Tommy/dechenli_lab_local/scfea26/data/independent/kurten_GSE164690.rds"

KEEP_FRACTIONS <- c("CD45p", "CD45n")   # drop PBL (blood)

# -------- HPV MAP: fill from Kurten 2021 Table S1 (TRUE = HPV+, FALSE = HPV-) --
# Leave as NA to skip HPV; edit these once you have the paper's assignments.
HPV_MAP <- c(
  HN01 = NA, HN02 = NA, HN03 = NA, HN04 = NA, HN05 = NA, HN06 = NA,
  HN07 = NA, HN08 = NA, HN09 = NA, HN10 = NA, HN11 = NA, HN12 = NA,
  HN13 = NA, HN14 = NA, HN15 = NA, HN16 = NA, HN17 = NA, HN18 = NA)

# ---- discover samples --------------------------------------------------------
mtx <- list.files(DIR, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE)
stub_of <- function(p) sub("_matrix\\.mtx\\.gz$", "", basename(p))
stubs <- stub_of(mtx)
# keep only wanted fractions
frac_of <- function(s) tail(strsplit(s, "_")[[1]], 1)
stubs <- stubs[sapply(stubs, frac_of) %in% KEEP_FRACTIONS]
say("samples to load (CD45p + CD45n): %d", length(stubs))

read_one <- function(stub) {
  mf <- file.path(DIR, paste0(stub, "_matrix.mtx.gz"))
  bf <- file.path(DIR, paste0(stub, "_barcodes.tsv.gz"))
  ff <- file.path(DIR, paste0(stub, "_features.tsv.gz"))
  if (!all(file.exists(mf, bf, ff))) { say("  [skip] incomplete: %s", stub); return(NULL) }
  m  <- readMM(mf)
  bc <- readLines(bf)
  ft <- read.delim(ff, header = FALSE, stringsAsFactors = FALSE)
  sym <- if (ncol(ft) >= 2) ft[[2]] else ft[[1]]      # features: ensembl, symbol, type
  rownames(m) <- make.unique(sym)
  colnames(m) <- paste0(stub, "|", bc)
  m
}

parse_meta <- function(stub) {
  s <- sub("^GSM[0-9]+_", "", stub)          # HN01_CD45p
  pid  <- strsplit(s, "_")[[1]][1]            # HN01
  frac <- tail(strsplit(s, "_")[[1]], 1)      # CD45p / CD45n
  list(patient = pid, sample = s, fraction = frac,
       hpv = if (pid %in% names(HPV_MAP)) HPV_MAP[[pid]] else NA)
}

objs <- list()
for (i in seq_along(stubs)) {
  st <- stubs[i]; m <- read_one(st); if (is.null(m)) next
  md <- parse_meta(st)
  o <- CreateSeuratObject(m, project = md$sample, min.cells = 3, min.features = 200)
  o$patient <- md$patient; o$sample <- md$sample
  o$fraction <- md$fraction
  o$hpv <- ifelse(is.na(md$hpv), NA_character_, ifelse(md$hpv, "Positive", "Negative"))
  o$tissue <- "tumor"
  objs[[st]] <- o
  if (i %% 10 == 0) say("  read %d / %d", i, length(stubs))
}
say("merging %d samples...", length(objs))
obj <- if (length(objs) == 1) objs[[1]] else merge(objs[[1]], objs[-1])
rm(objs); gc()
say("merged: %d genes x %d cells", nrow(obj), ncol(obj))

obj[["pct_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & pct_mt < 20)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- tryCatch(JoinLayers(obj), error = function(e) obj)   # Seurat v5
say("after QC: %d cells", ncol(obj))

# ---- lineage: broad, marker-score. CD45n gives us FIBROBLASTS ---------------
sig <- list(
  Myeloid     = c("LYZ","CD68","CD14","FCGR3A","ITGAX","C1QA","C1QB","CLEC10A"),
  T           = c("CD3D","CD3E","CD8A","CD4","IL7R","TRAC"),
  NK          = c("NKG7","GNLY","KLRD1"),
  B           = c("CD79A","MS4A1","CD79B"),
  Plasma      = c("MZB1","IGHG1","JCHAIN"),
  Fibroblast  = c("COL1A1","COL1A2","DCN","LUM","PDGFRB","ACTA2","MMP11","FAP"),
  Endothelial = c("PECAM1","VWF","CLDN5","CDH5"),
  Epithelial  = c("KRT5","KRT14","EPCAM","SFN","KRT6A"))
sig <- lapply(sig, function(g) intersect(g, rownames(obj)))
E <- GetAssayData(obj, layer = "data")
score <- sapply(sig, function(g) if (length(g)) Matrix::colMeans(E[g,,drop=FALSE]) else rep(0, ncol(obj)))
obj$lineage <- colnames(score)[max.col(score, ties.method = "first")]
obj$lineage_broad <- dplyr::recode(obj$lineage, NK = "T/NK", T = "T/NK")

say("\nlineage assignment:"); print(table(obj$lineage_broad))
say("\nfraction x lineage (CD45n should carry the fibroblasts):")
print(table(obj$fraction, obj$lineage_broad))
say("\nHPV labelling: %s",
    if (all(is.na(obj$hpv))) "NONE set -- edit HPV_MAP from Kurten Table S1" else "present")
if (!all(is.na(obj$hpv))) print(table(obj$hpv, obj$lineage_broad))

for (g in c("SLC16A3","SLC16A1","SLC16A4","SPP1","MMP11","CD68","COL1A1"))
  say("  %-9s %s", g, ifelse(g %in% rownames(obj), "present", "MISSING"))

saveRDS(obj, OUT)
say("\n-> %s", OUT)
say("\nNEXT (fibroblast MCT4 specificity is the new test here):")
say("  Rscript 32_independent_cohort.R --obj %s \\", OUT)
say("      --name kurten_GSE164690 --celltype lineage_broad --sample sample")
