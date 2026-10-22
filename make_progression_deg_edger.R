#!/usr/bin/env Rscript
# ==============================================================================
# make_progression_deg_edger.R  --  paired pseudobulk DE with edgeR
#
#   R --vanilla
#   source("make_progression_deg_edger.R")
#
# Replaces the CP10K-then-average approach with the standard pseudobulk route:
# sum raw counts within each patient x stage x cell type, then model with edgeR
# using a paired design (~ patient + stage).
#
# WHY. Averaging normalized values discards the count structure that negative
# binomial models use, weights a 15-cell profile the same as a 500-cell one, and
# leaves a paired Wilcoxon on ten patients that floors at p = 0.002. Summing raw
# counts and fitting edgeR is what Squair et al. do and what a reader expects.
#
# Patients are restricted to the flux analysis's set, as before, so every paired
# comparison in this work rests on the same patients.
#
# Output -> results/deg/  (edgeR files are suffixed _edger)
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(readr); library(edgeR)
})

ROOT     <- Sys.getenv("SCFEA_ROOT", unset = "/Users/Tommy/dechenli_lab_local/scfea26")
OBJ      <- file.path(ROOT, "data/obj.final.2.0.rda")
OUT      <- file.path(ROOT, "results/deg"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
FLUXSET  <- file.path(ROOT, "results/all_sites/bulk_before/fig4/fig4D_expression_ratio_paired.csv")
LINEAGES <- c("Myeloid","Fibroblast","T")
NORMAL_T <- c("NL","Normal"); CANCER_T <- c("CA","CA-E","CA-A","CA-c","CA-p")
PROG     <- c("He","Kim")
MIN_CELLS <- 10

if (!exists(".obj_deg")) {
  cat("loading Seurat object (slow, once)...\n")
  load(OBJ); .obj_deg <<- get(ls()[sapply(ls(), function(x) inherits(get(x), "Seurat"))][1])
}
obj <- .obj_deg
md  <- obj@meta.data
md$stage <- ifelse(md$tissue %in% NORMAL_T, "normal",
            ifelse(md$tissue %in% CANCER_T, "carcinoma", NA))

flux_pt <- if (file.exists(FLUXSET)) read.csv(FLUXSET)$patient else NULL
if (is.null(flux_pt)) cat("[!] flux patient list not found; using all paired patients\n")

for (lin in LINEAGES) {
  cat(sprintf("\n=== %s ===\n", lin))
  keep <- !is.na(md$stage) & md$annot_4th == lin & md$study %in% PROG
  if (!sum(keep)) { cat("  no cells\n"); next }
  m <- md[keep, ]
  m$key <- paste(m$patient, m$stage, sep = "__")

  # ---- sum raw counts within each patient x stage ---------------------------
  C <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")[, rownames(m), drop = FALSE]
  ind <- Matrix::sparse.model.matrix(~ 0 + factor(m$key))
  colnames(ind) <- levels(factor(m$key))
  pb <- as.matrix(C %*% ind)                       # genes x profiles, summed counts
  ncell <- as.integer(table(factor(m$key))[colnames(pb)])

  info <- tibble(key = colnames(pb), n_cells = ncell) %>%
    mutate(patient = sub("__.*", "", key), stage = sub(".*__", "", key)) %>%
    filter(n_cells >= MIN_CELLS)
  paired <- info %>% count(patient) %>% filter(n == 2) %>% pull(patient)
  if (!is.null(flux_pt)) {
    dropped <- setdiff(paired, flux_pt)
    if (length(dropped)) cat(sprintf("  restricting to the flux patient set: dropping %s\n",
                                     paste(dropped, collapse = ", ")))
    paired <- intersect(paired, flux_pt)
  }
  cat(sprintf("  paired patients: %d\n", length(paired)))
  if (length(paired) < 4) { cat("  too few, skipping\n"); next }

  info <- info %>% filter(patient %in% paired) %>% arrange(patient, stage)
  Y <- pb[, info$key, drop = FALSE]

  # ---- edgeR, paired design -------------------------------------------------
  d <- DGEList(counts = Y)
  grp <- factor(info$stage, levels = c("normal","carcinoma"))
  pat <- factor(info$patient)
  design <- model.matrix(~ pat + grp)
  keepg <- filterByExpr(d, design = design)
  cat(sprintf("  genes after filterByExpr: %d of %d\n", sum(keepg), nrow(d)))
  d <- d[keepg, , keep.lib.sizes = FALSE]
  d <- calcNormFactors(d)
  d <- estimateDisp(d, design)
  fit <- glmQLFit(d, design)
  qlf <- glmQLFTest(fit, coef = ncol(design))     # carcinoma vs normal
  tt  <- topTags(qlf, n = Inf)$table %>% tibble::rownames_to_column("gene")

  # per-patient direction, for the "n of n rising" column the text uses
  cpm_m <- edgeR::cpm(d, log = FALSE)
  rise <- sapply(seq_len(nrow(cpm_m)), function(i) {
    v <- tapply(cpm_m[i, ], list(info$patient, info$stage), mean)
    sum(v[, "carcinoma"] > v[, "normal"], na.rm = TRUE)
  })
  res <- tt %>% mutate(n_rising = rise[match(gene, rownames(cpm_m))],
                       n_patients = length(paired)) %>%
    select(gene, logFC, logCPM, F, PValue, FDR, n_rising, n_patients) %>%
    arrange(PValue)

  f <- file.path(OUT, sprintf("DEG_%s_normal_vs_carcinoma_edger.csv", tolower(sub(" .*","",lin))))
  write_csv(res, f)
  cat(sprintf("  -> %s\n  %d genes | FDR < 0.05: %d | FDR < 0.10: %d\n",
              f, nrow(res), sum(res$FDR < .05), sum(res$FDR < .10)))

  KEY <- c("SLC16A3","SLC16A1","LDHA","HK2")
  cat("  key genes:\n")
  print(as.data.frame(res %>% filter(gene %in% KEY) %>%
    transmute(gene, logFC = round(logFC, 3),
              rising = sprintf("%d/%d", n_rising, n_patients),
              PValue = signif(PValue, 3), FDR = signif(FDR, 3))), row.names = FALSE)
}
cat("\nCompare these against the CP10K-average versions before changing the text.\n")
