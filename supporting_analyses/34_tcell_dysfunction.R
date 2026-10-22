#!/usr/bin/env Rscript
# ==============================================================================
# 34_tcell_dysfunction.R  --  does myeloid lactate export track T-cell exhaustion?
#
#   Rscript 34_tcell_dysfunction.R --cohort atlas   --arm bulk_before
#   Rscript 34_tcell_dysfunction.R --cohort cillo   --obj .../cillo_GSE139324.rds
#   Rscript 34_tcell_dysfunction.R --cohort kurten  --obj .../kurten_GSE164690.rds
#
# THE CLAIM THIS ADDRESSES.
# The niche is called "immunosuppressive", but so far we have only shown the
# export EXISTS and sits near SPP1+/MMP11+ cells. We have NOT shown it is
# associated with actual T-cell dysfunction. Lactate is a known immunosuppressant
# (acidifies the TME, impairs T-cell proliferation/cytokine output, linked to
# exhaustion). Prediction: samples with high MYELOID export have MORE EXHAUSTED
# T cells.
#
# THE DESIGN THAT MAKES IT A MICROENVIRONMENT CLAIM (not a within-cell artifact).
# One point per SAMPLE. x = export measured in the MYELOID cells; y = exhaustion
# measured in the T cells. x and y come from DIFFERENT cells in the same sample,
# so a correlation means "samples where myeloid export is high are samples where
# T cells are exhausted" -- a tissue-level association, not a single-cell one.
#
# WHAT "EXPORT" IS, PER COHORT.
#   atlas  : real scFEA M_78 flux SHARE in myeloid (the actual quantity).
#   cillo  : SLC16A3 (MCT4) expression in myeloid -- the grounded proxy
#            (07c: M_78 flux tracks SLC16A3 at rho 0.84-0.89 in myeloid).
#   kurten : same MCT4 proxy.
# Labelled in output so the comparison is honest.
#
# EXHAUSTION SCORE = mean z of PDCD1, HAVCR2, LAG3, TIGIT, CTLA4, ENTPD1, TOX
# in T/NK cells, per sample.
#
# HONEST LIMIT (stated in output): correlation, not causation. High export with
# high exhaustion is CONSISTENT with lactate-driven suppression but does not prove
# it -- reverse causation and a shared driver (hypoxia) are not excluded.
#
# Output -> results/tcell_dysfunction/<cohort>/
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(stringr)
})
say <- function(...) message(sprintf(...))
a <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, d=NULL) if (k %in% a) a[which(a==k)+1] else d
COHORT <- getarg("--cohort"); OBJ <- getarg("--obj"); ARM <- getarg("--arm","bulk_before")
if (is.null(COHORT)) stop("--cohort atlas|cillo|kurten required")

ROOT <- "/Users/Tommy/dechenli_lab_local/scfea26"
OUT  <- file.path(ROOT, "results", "tcell_dysfunction", COHORT)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

EXHAUSTION <- c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","ENTPD1","TOX")

# ============================================================================
# per-cohort: build a per-sample table of (myeloid_export, tcell_exhaustion)
# ============================================================================
build_sc <- function(obj_path) {
  suppressPackageStartupMessages(library(Seurat))
  o <- readRDS(obj_path)
  o <- tryCatch(SeuratObject::JoinLayers(o), error=function(e) o)
  md <- o@meta.data
  ct <- if ("lineage_broad" %in% names(md)) "lineage_broad" else "lineage"
  smp <- "sample"
  ex_genes <- intersect(EXHAUSTION, rownames(o))
  need <- c("SLC16A3", ex_genes)
  E <- GetAssayData(o, layer="data")[intersect(need, rownames(o)), , drop=FALSE]
  D <- data.frame(t(as.matrix(E)), check.names=FALSE)
  D$lin <- as.character(md[[ct]]); D$sample <- as.character(md[[smp]])

  # myeloid export = mean SLC16A3 in myeloid cells, per sample
  myl <- D %>% filter(lin=="Myeloid") %>% group_by(sample) %>%
    summarise(export = mean(SLC16A3), n_myeloid = n(), .groups="drop") %>%
    filter(n_myeloid >= 10)
  # T-cell exhaustion = mean z of markers within T cells, per sample
  tz <- D %>% filter(lin %in% c("T/NK","T")) %>%
    mutate(across(all_of(ex_genes), ~ as.numeric(scale(.x)))) %>%
    group_by(sample) %>%
    summarise(exhaustion = mean(c_across(all_of(ex_genes)), na.rm=TRUE),
              n_tcell = n(), .groups="drop") %>%
    filter(n_tcell >= 10)
  inner_join(myl, tz, by="sample") %>% mutate(export_type="SLC16A3 (MCT4) expr")
}

build_atlas <- function() {
  # myeloid export = M_78 flux SHARE in myeloid pseudobulk units, averaged/sample
  source("config_all_sites.R")
  cw <- readRDS(file.path(PATHS$cw_dir,"cw_finaltype.rds"))
  mc <- mcols_of(cw)
  myl <- cw %>% filter(annot_4th=="Myeloid", n_cells>=GATES$min_cells)
  X <- as.matrix(myl[,mc]); myl$export <- (X/rowSums(X))[,"M_78"]
  myl_s <- myl %>% group_by(sample=patient) %>%     # atlas 'sample'~patient
    summarise(export = weighted.mean(export, n_cells), n_myeloid=sum(n_cells), .groups="drop")

  # T-cell exhaustion from the imputed expression used for grounding
  ep <- file.path(ROOT,"data/input/all_sites/bulk_before","T","train_expression.csv")
  mp <- file.path(ROOT,"metadata/all_sites/bulk_before","T","metadata_units.csv")
  e <- data.table::fread(ep); g <- as.character(e[[1]])
  ex_g <- intersect(EXHAUSTION, g)
  Em <- as.matrix(e[,-1]); rownames(Em) <- g
  mu <- data.table::fread(mp); data.table::setnames(mu,1,"id")
  mu <- as.data.frame(mu); rownames(mu) <- as.character(mu$id)
  ids <- colnames(Em); mu <- mu[ids,]
  Z <- t(scale(t(Em[ex_g,,drop=FALSE])))
  tz <- tibble(sample=mu$patient, exh=colMeans(Z,na.rm=TRUE), n=mu$n_cells) %>%
    group_by(sample) %>% summarise(exhaustion=weighted.mean(exh,n),
                                   n_tcell=sum(n), .groups="drop")
  inner_join(myl_s, tz, by="sample") %>% mutate(export_type="M_78 flux share (real)")
}

DF <- if (COHORT=="atlas") build_atlas() else build_sc(OBJ)
say("cohort %s | samples with both myeloid & T: %d", COHORT, nrow(DF))
if (nrow(DF) < 6) stop("too few samples for a correlation")

# ---- the correlation --------------------------------------------------------
ct <- cor.test(DF$export, DF$exhaustion, method="spearman")
say("\n=== myeloid export vs T-cell exhaustion (%s) ===", COHORT)
say("  export measure: %s", DF$export_type[1])
say("  Spearman rho = %.3f | p = %.4g | n = %d samples",
    ct$estimate, ct$p.value, nrow(DF))
# bootstrap CI (small n -> report uncertainty honestly)
set.seed(1)
bs <- replicate(2000, {
  i <- sample(nrow(DF), replace=TRUE)
  suppressWarnings(cor(DF$export[i], DF$exhaustion[i], method="spearman"))
})
say("  bootstrap 95%% CI: [%.2f, %.2f]", quantile(bs,.025,na.rm=T), quantile(bs,.975,na.rm=T))

write_csv(DF, file.path(OUT,"per_sample.csv"))

p <- ggplot(DF, aes(export, exhaustion)) +
  geom_smooth(method="lm", se=TRUE, colour="#b2182b", fill="#f4c9c4") +
  geom_point(size=2.4, alpha=.8) +
  labs(title=sprintf("Myeloid lactate export vs T-cell exhaustion (%s)", COHORT),
       subtitle=sprintf("one point = one sample. x from myeloid cells, y from T cells.\nexport = %s | Spearman rho = %.2f, p = %.3g, n = %d",
                        DF$export_type[1], ct$estimate, ct$p.value, nrow(DF)),
       x="myeloid export (per sample)", y="T-cell exhaustion score") +
  theme_bw(base_size=11)
ggsave(file.path(OUT,"fig_export_vs_exhaustion.png"), p, width=7, height=5.5, dpi=150, bg="white")
say("-> %s", OUT)
say("\nLIMIT: correlational. Consistent with lactate-driven suppression; does not")
say("exclude reverse causation or a shared driver (e.g. hypoxia). x and y are from")
say("DIFFERENT cells in the same sample, so this is a tissue-level association.")
