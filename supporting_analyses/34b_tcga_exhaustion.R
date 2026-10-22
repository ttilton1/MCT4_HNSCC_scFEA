#!/usr/bin/env Rscript
# ==============================================================================
# 34b_tcga_exhaustion.R  --  bulk cross-check: export signature vs exhaustion
#
#   Rscript 34b_tcga_exhaustion.R
#
# The single-cell version (34) is the real test -- it separates myeloid export
# from T-cell exhaustion by measuring them in different cells. This bulk version
# is a SUPPORTING cross-check only: across TCGA tumors, does the export signature
# correlate with a T-cell exhaustion signature?
#
# BULK CAVEAT (why this is support, not primary): bulk cannot separate "T cells
# are more exhausted" from "there are more T cells" or "there is more myeloid".
# A positive correlation is consistent with the single-cell finding; it cannot
# establish the cell-level relationship on its own. Reported as corroboration.
#
# Output -> results/tcell_dysfunction/tcga/
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(ggplot2); library(readr) })
say <- function(...) message(sprintf(...))
ROOT <- "/Users/Tommy/dechenli_lab_local/scfea26"
OUT  <- file.path(ROOT,"results/tcell_dysfunction/tcga"); dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
EXPR <- file.path(ROOT,"tcga/hnsc_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt")

EXPORT_SIG <- c("SLC16A3","LDHA","HK2")
EXHAUSTION <- c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","ENTPD1","TOX")

ex <- fread(EXPR)[Hugo_Symbol!="" & !is.na(Hugo_Symbol)][!duplicated(Hugo_Symbol)]
M <- as.matrix(ex[,-(1:2)]); rownames(M) <- ex$Hugo_Symbol; M[!is.finite(M)] <- 0
Z <- t(scale(t(log2(M+1))))
sc_ <- function(g) colMeans(Z[intersect(g,rownames(Z)),,drop=FALSE],na.rm=TRUE)

df <- data.frame(export = sc_(EXPORT_SIG), exhaustion = sc_(EXHAUSTION),
                 mct4 = Z["SLC16A3",])
r_sig <- cor.test(df$export, df$exhaustion, method="spearman")
r_m   <- cor.test(df$mct4,   df$exhaustion, method="spearman")
say("=== TCGA bulk (SUPPORTING cross-check) ===")
say("  export signature vs exhaustion: rho = %.3f, p = %.3g, n = %d",
    r_sig$estimate, r_sig$p.value, nrow(df))
say("  SLC16A3 alone   vs exhaustion: rho = %.3f, p = %.3g", r_m$estimate, r_m$p.value)

write_csv(df, file.path(OUT,"tcga_export_exhaustion.csv"))
p <- ggplot(df, aes(export, exhaustion)) +
  geom_point(alpha=.4, size=1) + geom_smooth(method="lm", colour="#b2182b") +
  labs(title="Export signature vs T-cell exhaustion, TCGA-HNSC (bulk cross-check)",
       subtitle=sprintf("rho = %.2f, p = %.3g, n = %d. BULK -- cannot separate exhaustion from T-cell abundance.",
                        r_sig$estimate, r_sig$p.value, nrow(df)),
       x="export signature (SLC16A3+LDHA+HK2)", y="exhaustion signature") +
  theme_bw(base_size=11)
ggsave(file.path(OUT,"fig_tcga_exhaustion.png"), p, width=7, height=5, dpi=150, bg="white")
say("-> %s", OUT)
