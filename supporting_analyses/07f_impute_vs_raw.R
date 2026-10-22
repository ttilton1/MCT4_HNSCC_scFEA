#!/usr/bin/env Rscript
# ==============================================================================
# 07f_impute_vs_raw.R  --  why do 07d and 07e disagree?
#
#   Rscript 07f_impute_vs_raw.R --arm bulk_before
#
# 07d (MAGIC-imputed expression, scFEA's own 661-gene input): HK2 co-moves with
#     export flux in myeloid at rho = 0.64, clears the null.
# 07e (raw log-normalized expression from the Seurat object): the same HK2 sits
#     inside the null; nothing external clears it.
#
# Two possible explanations, with very different consequences:
#   (A) IMPUTATION. MAGIC fills dropouts, so imputed single-gene expression
#       correlates cleanly; raw expression is mostly zeros and the correlation is
#       crushed by dropout noise. -> the co-movement is REAL but only visible
#       after imputation; a defensible one-line caveat.
#   (B) SOMETHING ELSE (different units, different cells, a join bug). -> we do
#       not trust 07d and have to re-examine.
#
# This isolates it: take the SAME units, correlate export flux against HK2
# computed THREE ways -- (1) MAGIC-imputed (from train_expression.csv, exactly
# what 07d used), (2) raw log-norm pseudobulk (what 07e used), (3) raw but only
# using cells with nonzero HK2 (removes the dropout mass). If (1) clears the null
# and (2) does not, and (3) moves toward (1), then it is DROPOUT -- explanation A.
#
# Run for a few genes so it is not a single-gene fluke: HK2 (the disputed one),
# SLC16A3 (MCT4, the headline enzyme), HIF1A (the failed hypoxia test), ACTB (neg).
#
# Output -> results/all_sites/<ARM>/07f_impute_vs_raw/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(Matrix); library(stringr) })

OUT <- file.path(PATHS$results, "07f_impute_vs_raw")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
MIN_N_COHORT <- 8; N_NEG <- 200; set.seed(1)

GENES  <- c("HK2","SLC16A3","HIF1A","GLUT1_SLC2A1"="SLC2A1","ACTB")
CT     <- "Myeloid"                 # the cell type where 07d/07e disagree most
MODULE <- "M_78"

# ---- (1) imputed expression: exactly the matrix 07d/07c used -----------------
say("loading MAGIC-imputed expression (scFEA input) for %s...", CT)
# NOTE (arm-safety): both arms consume the SAME MAGIC-imputed expression -- the
# imputation ran once in shared/ and the arms differ only in whether pseudobulking
# happens before or after scFEA. Unit-level expression is therefore identical, so
# grounding reads the single bulk_before copy regardless of --arm. This also
# isolates the comparison: any arm difference in enzyme correlation reflects the
# FLUX method, not different input expression.
ep <- file.path(ROOT, "data/input/all_sites", "bulk_before", CT, "train_expression.csv")
mp <- file.path(file.path(ROOT, "metadata", "all_sites", "bulk_before"), CT, "metadata_units.csv")
e  <- fread(ep); g <- as.character(e[[1]]); Eimp <- as.matrix(e[,-1]); rownames(Eimp) <- g
mu <- fread(mp); setnames(mu,1,"id"); mu[, id:=as.character(id)]
mu <- as.data.frame(mu); rownames(mu) <- mu$id; mu <- mu[colnames(Eimp),,drop=FALSE]
imp <- as_tibble(t(Eimp[intersect(unname(GENES), g),,drop=FALSE])) %>%
  mutate(patient=mu$patient, tissue=mu$tissue, finaltype=mu$finaltype,
         study=mu$study, n_cells=mu$n_cells) %>%
  group_by(patient, tissue, finaltype, study) %>%
  summarise(across(any_of(unname(GENES)), ~ weighted.mean(.x, n_cells)),
            n_cells=sum(n_cells), .groups="drop")
say("  imputed genes available: %s", paste(intersect(unname(GENES), g), collapse=", "))

# ---- (2,3) raw expression from the object, pseudobulked two ways --------------
say("loading full object for RAW expression (this is the slow step)...")
suppressPackageStartupMessages(library(Seurat))
nm <- load(PATHS$seurat); obj <- get(nm[1])
DefaultAssay(obj) <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
md <- obj@meta.data; md$barcode <- rownames(md)
md <- md[md$annot_4th == CT, ]
raw <- GetAssayData(obj, layer = "data")[intersect(unname(GENES), rownames(obj)), md$barcode, drop=FALSE]

md$unit <- paste(md$patient, md$tissue, md$finaltype, sep="||")
units <- split(md$barcode, md$unit)
# (2) mean over ALL cells (includes dropout zeros) ; (3) mean over NONZERO cells only
raw_all  <- sapply(units, function(bc) Matrix::rowMeans(raw[, bc, drop=FALSE]))
raw_nz   <- sapply(units, function(bc){ m<-raw[,bc,drop=FALSE]
  apply(m,1,function(r){ nz<-r[r>0]; if(length(nz)) mean(nz) else 0 }) })
mk_df <- function(M, suf){ uk<-do.call(rbind, strsplit(colnames(M),"\\|\\|"))
  as_tibble(t(M)) %>% rename_with(~paste0(.x,suf)) %>%
    mutate(patient=uk[,1], tissue=uk[,2], finaltype=uk[,3]) }
R2 <- mk_df(raw_all, "_rawAll"); R3 <- mk_df(raw_nz, "_rawNZ")

# ---- flux, same units --------------------------------------------------------
FL <- readRDS(file.path(PATHS$cw_dir,"cw_finaltype.rds")) %>%
  filter(n_cells >= GATES$min_cells, annot_4th==CT)
mcols <- mcols_of(FL); Fm <- as.matrix(FL[,mcols])
FL$export <- (Fm/rowSums(Fm))[,MODULE]
key <- c("patient","tissue","finaltype")
D <- FL %>% select(all_of(key), study, export) %>%
  left_join(imp, by=c(key,"study")) %>%
  left_join(R2, by=key) %>% left_join(R3, by=key)

# ---- empirical null (raw all-cell, matched to the raw pipeline) ---------------
say("building null from raw expression...")
poolg <- sample(setdiff(rownames(obj), unname(GENES)), N_NEG)
rawp <- GetAssayData(obj, layer="data")[poolg, md$barcode, drop=FALSE]
np <- sapply(units, function(bc) Matrix::rowMeans(rawp[,bc,drop=FALSE]))
NP <- mk_df(np, "")
Dn <- FL %>% select(all_of(key), study, export) %>% left_join(NP, by=key)
cor_meta <- function(df, a, b){
  pc <- df %>% group_by(study) %>% filter(n()>=MIN_N_COHORT) %>%
    summarise(r=suppressWarnings(cor(.data[[a]],.data[[b]],method="spearman",
              use="complete.obs")), n=n(), .groups="drop") %>% filter(is.finite(r))
  if(!nrow(pc)) return(NA_real_)
  z<-atanh(pmin(pmax(pc$r,-.999),.999)); w<-pc$n-3; tanh(sum(w*z)/sum(w)) }
null_rho <- sapply(poolg, function(g) cor_meta(Dn, "export", g))
bar <- quantile(null_rho, .999, na.rm=TRUE)
say("  raw null 99.9th pctl bar = %.3f", bar)

# ---- the comparison ----------------------------------------------------------
say("\n============================================================")
say("HK2 (and friends) vs export flux, THREE expression pipelines, %s", CT)
say("bar (raw null 99.9th) = %.2f", bar)
say("============================================================")
tab <- lapply(intersect(unname(GENES), colnames(imp)), function(gn){
  tibble(gene=gn,
    imputed   = cor_meta(D, "export", gn),
    raw_all   = cor_meta(D, "export", paste0(gn,"_rawAll")),
    raw_nz    = cor_meta(D, "export", paste0(gn,"_rawNZ")))
}) %>% bind_rows() %>%
  mutate(across(where(is.numeric), ~round(.x,3)),
         imputed_beats = imputed > bar, rawall_beats = raw_all > bar)
print(as.data.frame(tab), row.names=FALSE)
write_csv(tab, file.path(OUT,"impute_vs_raw.csv"))

say("\nREAD:")
say("  If imputed clears the bar and raw_all does not, but raw_NZ (dropping the")
say("  dropout zeros) moves back UP toward imputed -> the discrepancy is DROPOUT.")
say("  The co-movement is real; it is just invisible in raw, dropout-limited data.")
say("  Then: cite 07d (imputed), and state the raw-expression caveat honestly.")
say("  If raw_NZ stays low too -> the imputed signal may be a MAGIC smoothing")
say("  artifact, and the HK2 co-movement should be treated with caution.")

# dropout rate context
say("\ndropout context (fraction of %s cells with zero counts):", CT)
for (gn in intersect(unname(GENES), rownames(raw)))
  say("  %-8s zero in %.0f%% of cells", gn, 100*mean(raw[gn,]==0))

# ---- figure ------------------------------------------------------------------
pl <- tab %>% select(gene, imputed, raw_all, raw_nz) %>%
  pivot_longer(-gene, names_to="pipeline", values_to="rho") %>%
  mutate(pipeline=factor(pipeline, levels=c("imputed","raw_nz","raw_all"),
                         labels=c("MAGIC imputed\n(07d)","raw, nonzero cells","raw, all cells\n(07e)")))
p <- ggplot(pl, aes(pipeline, rho, group=gene, colour=gene)) +
  geom_hline(yintercept=bar, linetype=2, colour="grey40") +
  annotate("text", x=0.7, y=bar, label="null bar", vjust=-0.4, size=3, colour="grey40") +
  geom_line(linewidth=1) + geom_point(size=3) +
  labs(title=sprintf("Why 07d and 07e disagree (%s, export flux vs gene)", CT),
       subtitle="Same units, three expression pipelines. If imputed clears the bar and raw does not,\nbut raw-nonzero recovers, the gap is dropout -- the co-movement is real, just dropout-limited.",
       x=NULL, y="Spearman rho (export share ~ gene)", colour="gene") +
  theme_bw(base_size=11)
ggsave(file.path(OUT,"fig_impute_vs_raw.png"), p, width=9, height=6, dpi=150, bg="white")
say("\n-> %s", OUT)
