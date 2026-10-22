#!/usr/bin/env Rscript
# ==============================================================================
# 32_independent_cohort.R  --  does the MCT4 finding hold in data we never touched?
#
#   Rscript 32_independent_cohort.R --obj /path/to/seurat.rds [--tumor_only]
#
# WHY THIS IS THE STRONGEST REMAINING VALIDATION.
# Every result so far comes from the five cohorts in our atlas (He, Kim, Koh,
# Pittet, Tirosh) -- the same data that GENERATED the hypothesis. A dataset that
# played no role in generating it is the only way to test whether the finding
# travels. CPTAC gave protein-level support but is bulk; this gives cell-type
# resolution, which is what the claim is actually about.
#
# WHY NO scFEA HERE.
# The grounding (07c) showed M_78 flux tracks SLC16A3 at rho = 0.84 (bulk_after)
# and 0.89 (bulk_before) in myeloid, ~3x above SLC16A1. For this cell type the
# flux IS an MCT4 readout. So the claims can be tested directly on SLC16A3
# expression without re-running the pipeline -- which also removes scFEA as a
# shared dependency between discovery and validation, making the test cleaner
# rather than weaker.
#
# PRE-SPECIFIED CLAIMS (state them before looking):
#   C1. MCT4 is MYELOID-BIASED: SLC16A3 higher in myeloid than in T cells.
#   C2. ISOFORM SPECIFICITY: within myeloid, SLC16A3 > SLC16A1; the ordering does
#       NOT hold in T cells. (The single strongest result in the study.)
#   C3. NICHE ASSOCIATION: SLC16A3 co-expresses with SPP1 / MMP11 across samples.
#   C4. TUMOR ENRICHMENT: myeloid SLC16A3 higher in tumor than normal tissue.
#       (Only testable if the object has a normal/control arm.)
#
# A claim that fails here is not fatal -- annotation and platform differ -- but a
# claim that PASSES in independent data is the strongest evidence available
# short of a bench experiment.
#
# Output -> results/independent/<dataset>/
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(stringr)
})
say <- function(...) message(sprintf(...))

# ---- args -------------------------------------------------------------------
a <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, d = NULL) if (k %in% a) a[which(a == k) + 1] else d
OBJ   <- getarg("--obj")
NAME  <- getarg("--name", "independent")
CTCOL <- getarg("--celltype")      # metadata column with cell type
SMPL  <- getarg("--sample")        # metadata column with sample/patient
TISS  <- getarg("--tissue")        # metadata column with tumor/normal (optional)
TUMOR_ONLY <- "--tumor_only" %in% a || "--tumour_only" %in% a  # restrict to tumor tissue
if (is.null(OBJ)) stop("usage: --obj <seurat.rds> [--celltype col] [--sample col] [--tissue col]")

ROOT <- "/Users/Tommy/dechenli_lab_local/scfea26"
OUT  <- file.path(ROOT, "results", "independent", NAME)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MCT4 <- "SLC16A3"; MCT1 <- "SLC16A1"; MCT_ORPHAN <- "SLC16A4"
NICHE <- c("SPP1", "MMP11", "CD68", "ITGAX", "CLEC10A")

say("loading %s ...", OBJ)
obj <- if (grepl("\\.rda$|\\.RData$", OBJ)) { nm <- load(OBJ); get(nm[1]) } else readRDS(OBJ)
say("  %d genes x %d cells", nrow(obj), ncol(obj))
md <- obj@meta.data

# ---- resolve metadata columns ------------------------------------------------
guess <- function(pats, given) {
  if (!is.null(given) && given %in% names(md)) return(given)
  hit <- names(md)[grepl(paste(pats, collapse="|"), names(md), ignore.case = TRUE)]
  if (!length(hit)) return(NA_character_)
  hit[1]
}
CTCOL <- guess(c("celltype","cell_type","annot","cluster","ident","label"), CTCOL)
SMPL  <- guess(c("sample","patient","donor","orig.ident"), SMPL)
TISS  <- guess(c("tissue","type","condition","group","site","status"), TISS)
say("  celltype col: %s | sample col: %s | tissue col: %s", CTCOL, SMPL, TISS)
if (is.na(CTCOL) || is.na(SMPL)) {
  say("\n  available metadata columns:"); print(names(md)); stop("specify --celltype / --sample")
}
say("\n  cell types present:")
print(sort(table(md[[CTCOL]]), decreasing = TRUE))

# ---- restrict to tumor tissue ------------------------------------------------
# Transporter ratios must be compared on like tissue across cohorts. Cillo pools
# tumor TIL with blood PBMC and healthy tonsil, and blood myeloid cells are a
# different population from tumor macrophages. Kurten is tumor-only by
# construction (PBL dropped at build) and our atlas is restricted to carcinoma,
# so this flag brings Cillo onto the same footing.
if (TUMOR_ONLY) {
  if (is.na(TISS)) stop("--tumor_only requires --tissue")
  tv <- as.character(md[[TISS]])
  keep <- grepl("tumou?r|TIL|cancer|carcinoma|malig", tv, ignore.case = TRUE) &
         !grepl("PBMC|PBL|blood|normal|healthy|tonsil|adjacent", tv, ignore.case = TRUE)
  say("\n  --tumor_only: keeping %d of %d cells", sum(keep), length(keep))
  say("    kept   : %s", paste(sort(unique(tv[keep])),  collapse = ", "))
  say("    dropped: %s", paste(sort(unique(tv[!keep])), collapse = ", "))
  if (!sum(keep)) stop("no cells matched the tumor filter -- check --tissue values")
  obj <- obj[, keep]; md <- obj@meta.data
  say("  cell types after filtering:")
  print(sort(table(md[[CTCOL]]), decreasing = TRUE))
}


# ---- gene availability -------------------------------------------------------
G <- c(MCT4, MCT1, MCT_ORPHAN, NICHE)
have <- intersect(G, rownames(obj))
miss <- setdiff(G, rownames(obj))
say("\ngenes present: %s", paste(have, collapse=", "))
if (length(miss)) say("  missing: %s", paste(miss, collapse=", "))
if (!MCT4 %in% have) stop("SLC16A3 not in this object -- cannot test the claim")

DefaultAssay(obj) <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
E <- GetAssayData(obj, layer = "data")[have, , drop = FALSE]

D <- data.frame(t(as.matrix(E)), check.names = FALSE)
D$celltype <- as.character(md[[CTCOL]])
D$sample   <- as.character(md[[SMPL]])
if (!is.na(TISS)) D$tissue <- as.character(md[[TISS]])

# broad lineage buckets from whatever labels the authors used
lineage <- function(x) {
  x <- tolower(x)
  dplyr::case_when(
    grepl("mono|macro|myelo|dc$|dendritic|cdc|pdc|mast|neutro|tam", x) ~ "Myeloid",
    grepl("^t |t.cell|cd4|cd8|treg|tcell|nk", x)                      ~ "T/NK",
    grepl("fibro|caf|stromal|myofib", x)                              ~ "Fibroblast",
    grepl("^b |b.cell|plasma", x)                                     ~ "B/Plasma",
    grepl("epithe|tumor|malign|carcin|keratin", x)                    ~ "Epithelial",
    grepl("endo|vascular|lymphatic", x)                               ~ "Endothelial",
    TRUE ~ "Other")
}
D$lineage <- lineage(D$celltype)
say("\nlineage assignment:"); print(table(D$lineage))
if (sum(D$lineage == "Myeloid") < 200 || sum(D$lineage == "T/NK") < 200)
  say("  !! thin myeloid or T compartment -- C1/C2 will be underpowered")

# pseudobulk to sample x lineage (avoids treating cells as independent)
PB <- D %>% group_by(sample, lineage) %>%
  summarise(across(all_of(have), mean), n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 10)
say("pseudobulk units (sample x lineage, >=10 cells): %d", nrow(PB))

# ================= C1: is MCT4 myeloid-biased? ================================
say("\n============================================================")
say("C1 -- is SLC16A3 (MCT4) higher in MYELOID than in T cells?")
say("============================================================")
c1 <- PB %>% filter(lineage %in% c("Myeloid","T/NK")) %>%
  group_by(lineage) %>%
  summarise(n = n(), median_MCT4 = median(.data[[MCT4]]),
            mean_MCT4 = mean(.data[[MCT4]]), .groups="drop")
print(as.data.frame(c1), row.names = FALSE)
if (nrow(c1) == 2) {
  w <- wilcox.test(PB[[MCT4]][PB$lineage=="Myeloid"], PB[[MCT4]][PB$lineage=="T/NK"])
  say("  Wilcoxon myeloid vs T: p = %.3g", w$p.value)
  say("  VERDICT C1: %s",
      ifelse(c1$median_MCT4[c1$lineage=="Myeloid"] >
             c1$median_MCT4[c1$lineage=="T/NK"] & w$p.value < .05,
             "SUPPORTED", "not supported"))
}

# ================= C2: isoform specificity (the key test) =====================
say("\n============================================================")
say("C2 -- ISOFORM SPECIFICITY: MCT4 > MCT1 in myeloid, not in T?")
say("  (the strongest result in the study -- does it travel?)")
say("============================================================")
if (MCT1 %in% have) {
  c2 <- PB %>% group_by(lineage) %>%
    summarise(n = n(),
              MCT4 = median(.data[[MCT4]]), MCT1 = median(.data[[MCT1]]),
              ratio_MCT4_MCT1 = MCT4 / pmax(MCT1, 1e-9),
              ratio_median_of_ratios = median(.data[[MCT4]] / pmax(.data[[MCT1]], 1e-9)),
              p_paired = tryCatch(
                wilcox.test(.data[[MCT4]], .data[[MCT1]], paired = TRUE)$p.value,
                error = function(e) NA_real_), .groups="drop") %>%
    arrange(desc(ratio_MCT4_MCT1))
  print(as.data.frame(c2 %>% mutate(across(where(is.numeric), ~signif(.x,3)))),
        row.names = FALSE)
  write_csv(c2, file.path(OUT, "C2_isoform_specificity.csv"))
  m_r <- c2$ratio_MCT4_MCT1[c2$lineage=="Myeloid"]
  t_r <- c2$ratio_MCT4_MCT1[c2$lineage=="T/NK"]
  if (length(m_r) && length(t_r))
    say("  VERDICT C2: %s  (myeloid ratio %.2f vs T ratio %.2f)",
        ifelse(m_r > 1 & m_r > t_r, "SUPPORTED", "not supported"), m_r, t_r)
} else say("  SLC16A1 not measured -- C2 untestable")

# ================= C3: niche association ======================================
say("\n============================================================")
say("C3 -- does MCT4 co-express with the niche markers?")
say("============================================================")
myl <- PB %>% filter(lineage == "Myeloid")
if (nrow(myl) >= 8) {
  c3 <- lapply(intersect(NICHE, have), function(g) {
    r <- suppressWarnings(cor(myl[[MCT4]], myl[[g]], method="spearman",
                              use="complete.obs"))
    tibble(gene = g, rho = r, n = nrow(myl))
  }) %>% bind_rows() %>% arrange(desc(rho))
  print(as.data.frame(c3 %>% mutate(rho = round(rho,3))), row.names = FALSE)
  write_csv(c3, file.path(OUT, "C3_niche_coexpression.csv"))
} else say("  too few myeloid pseudobulk units (n=%d)", nrow(myl))

# ================= C4: tumor enrichment =======================================
say("\n============================================================")
say("C4 -- is myeloid MCT4 higher in TUMOR than normal?")
say("============================================================")
if (!is.na(TISS) && "tissue" %in% names(PB)) {
  PBt <- PB %>% mutate(grp = case_when(
    grepl("tumou?r|cancer|hnscc|malig", tolower(tissue)) ~ "tumor",
    grepl("normal|healthy|tonsil|adjacent|control|pbmc") ~ "normal",
    TRUE ~ NA_character_)) %>% filter(!is.na(grp), lineage == "Myeloid")
  if (n_distinct(PBt$grp) == 2 && nrow(PBt) >= 10) {
    w <- wilcox.test(PBt[[MCT4]] ~ PBt$grp)
    s <- PBt %>% group_by(grp) %>% summarise(n=n(), med=median(.data[[MCT4]]), .groups="drop")
    print(as.data.frame(s), row.names = FALSE)
    say("  Wilcoxon tumor vs normal: p = %.3g", w$p.value)
  } else say("  no usable tumor/normal split in this object")
} else say("  no tissue column -- C4 untestable (tumor-only dataset)")

# ---- figure ------------------------------------------------------------------
pl <- PB %>% select(lineage, all_of(intersect(c(MCT4, MCT1, MCT_ORPHAN), have))) %>%
  pivot_longer(-lineage, names_to = "gene", values_to = "expr") %>%
  mutate(gene = recode(gene,
    SLC16A3 = "SLC16A3 (MCT4, exporter)",
    SLC16A1 = "SLC16A1 (MCT1, importer)",
    SLC16A4 = "SLC16A4 (orphan)"))
p <- ggplot(pl, aes(reorder(lineage, expr, median), expr, fill = gene)) +
  geom_boxplot(outlier.size = .5, position = position_dodge(.8)) +
  coord_flip() +
  scale_fill_manual(values = c("#b2182b","#4393c3","grey70")) +
  labs(title = sprintf("MCT isoform expression by lineage (%s)", NAME),
       subtitle = "independent cohort, pseudobulked to sample x lineage.\nMCT4 above MCT1 in myeloid but not T = the specificity claim replicating.",
       x = NULL, y = "mean log-normalized expression", fill = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "fig_mct_isoforms.png"), p, width = 9, height = 5.5,
       dpi = 150, bg = "white")

write_csv(PB, file.path(OUT, "pseudobulk_expression.csv"))
say("\n-> %s", OUT)
say("\nREAD: C2 is the one that matters. If MCT4 exceeds MCT1 in myeloid and the")
say("ordering does not hold in T cells, the isoform specificity has replicated in")
say("data that had no role in generating the hypothesis.")
