#!/usr/bin/env Rscript
# ==============================================================================
# 01g_dominance_vs_agreement.R
#   For which modules does flux carry information beyond gene expression?
#
#   Rscript 01g_dominance_vs_agreement.R --arm bulk_before
#
# MOTIVATION. Module flux SHARE tracks module gene-expression SHARE closely for
# some modules (M_78 rho = 0.89, M_6 0.91, M_2 0.87, M_51 0.91) but not others
# (M_71 glucose import, rho = 0.25). The modules that agree closely tend to be
# dominated by one highly expressed gene; M_71 divides across four GLUT genes.
#
# HYPOTHESIS: scFEA flux approximates a rescaled expression measurement when one
# gene dominates a module's inputs, and diverges from it only when several genes
# contribute comparably and the network structure has something to arbitrate.
#
# If that holds, it is a practical statement about where the method adds
# information, and it identifies which modules a user can interpret as more than
# an expression readout.
#
# METRIC. For each module with >=2 genes present:
#   dominance = (mean expression of the single largest gene) / (sum over module genes)
#               1.0 = one gene supplies everything; 1/n = perfectly even
#   agreement = Spearman(flux share, module expression share) across units
#
# Run per major cell type, since dominance is cell-type specific (SLC16A3 is 92%
# of the lactate module in myeloid but 53% in T cells).
#
# Output -> results/all_sites/<ARM>/01g_dominance/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(stringr) })

OUT <- file.path(PATHS$results, "01g_dominance")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MGF <- file.path(ROOT, "reference", "module_gene_m168.csv")
if (!file.exists(MGF)) MGF <- file.path(ROOT, "module_gene_m168.csv")
mg <- read.csv(MGF, row.names = 1, stringsAsFactors = FALSE)
MR <- tryCatch(load_modref(), error = function(e) NULL)

CTS <- c("Myeloid","T","Fibroblast","Epithelial","B","Endothelial")

cw_all <- as.data.frame(readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))) %>%
  filter(n_cells >= GATES$min_cells)
mc_all <- mcols_of(cw_all)

per_ct <- function(ct) {
  ep <- file.path(ROOT, "data/input/all_sites/bulk_before", ct, "train_expression.csv")
  mp <- file.path(ROOT, "metadata/all_sites/bulk_before", ct, "metadata_units.csv")
  if (!file.exists(ep) || !file.exists(mp)) return(NULL)
  cw <- cw_all %>% filter(annot_4th == ct)
  if (nrow(cw) < 20) return(NULL)
  S <- as.matrix(cw[, mc_all]); S <- S / rowSums(S)

  e <- fread(ep); g <- as.character(e[[1]])
  E <- as.matrix(e[, -1]); rownames(E) <- g
  mu <- fread(mp); setnames(mu, 1, "id"); mu <- as.data.frame(mu)
  rownames(mu) <- as.character(mu$id); mu <- mu[colnames(E), ]
  key <- paste(mu$patient, mu$tissue, mu$finaltype)
  cwk <- paste(cw$patient, cw$tissue, cw$finaltype)
  tot <- colSums(E)

  do.call(rbind, lapply(mc_all, function(m) {
    gs <- as.character(unlist(mg[m, ]))
    gs <- intersect(gs[!is.na(gs) & gs != "" & gs != "nan"], g)
    if (length(gs) < 2) return(NULL)
    gm <- rowMeans(E[gs, , drop = FALSE])
    if (sum(gm) <= 0) return(NULL)
    es <- colSums(E[gs, , drop = FALSE]) / tot
    eu <- tapply(es, key, mean)
    r <- suppressWarnings(cor(S[, m], eu[match(cwk, names(eu))],
                              method = "spearman", use = "complete.obs"))
    data.frame(celltype = ct, module = m, n_genes = length(gs),
               dominance = max(gm) / sum(gm), top_gene = names(which.max(gm)),
               agreement = r, stringsAsFactors = FALSE)
  }))
}

say("computing dominance and flux-expression agreement per cell type...")
D <- do.call(rbind, lapply(CTS, per_ct))
D <- D[is.finite(D$agreement), ]
say("modules x cell types evaluated: %s", format(nrow(D), big.mark = ","))

# ---- the central relationship ------------------------------------------------
say("\n=== does single-gene dominance predict flux-expression agreement? ===")
for (ct in unique(D$celltype)) {
  d <- D[D$celltype == ct, ]
  say("  %-12s n=%3d  rho(dominance, agreement) = %+.2f",
      ct, nrow(d), cor(d$dominance, d$agreement, method = "spearman"))
}
say("  pooled: rho = %+.2f", cor(D$dominance, D$agreement, method = "spearman"))
say("\n  POSITIVE => flux approximates an expression readout for single-gene-")
say("  dominated modules, and diverges where several genes contribute. That")
say("  identifies which modules the model adds information to.")

# ---- binned view -------------------------------------------------------------
D$bin <- cut(D$dominance, breaks = c(0, .4, .6, .8, .95, 1),
             labels = c("<0.4","0.4-0.6","0.6-0.8","0.8-0.95",">0.95"))
say("\n=== agreement by dominance bin (all cell types) ===")
b <- aggregate(agreement ~ bin, D, function(x) c(n = length(x), med = median(x)))
print(data.frame(bin = b$bin, n = b$agreement[, "n"],
                 median_agreement = round(b$agreement[, "med"], 3)), row.names = FALSE)

# ---- where does the model add most? -----------------------------------------
say("\n=== modules where flux DIVERGES most from expression (model contributes) ===")
lo <- D[order(D$agreement), ][1:12, ]
print(data.frame(celltype = lo$celltype, module = lo$module, genes = lo$n_genes,
                 dominance = round(lo$dominance, 2), top = lo$top_gene,
                 agreement = round(lo$agreement, 2)), row.names = FALSE)

say("\n=== key modules in context ===")
KEY <- c("M_78","M_6","M_2","M_51","M_71","M_113","M_116")
k <- D[D$module %in% KEY & D$celltype %in% c("Myeloid","T"), ]
print(data.frame(celltype = k$celltype, module = k$module, genes = k$n_genes,
                 dominance = round(k$dominance, 2), top = k$top_gene,
                 agreement = round(k$agreement, 2)), row.names = FALSE)

write_csv(D, file.path(OUT, "dominance_vs_agreement.csv"))

p <- ggplot(D, aes(dominance, agreement)) +
  geom_point(alpha = .45, size = 1.3) +
  geom_smooth(method = "loess", se = FALSE, colour = "#b2182b") +
  facet_wrap(~ celltype) +
  labs(title = sprintf("Where does flux carry information beyond expression? (%s)", ARM),
       subtitle = "Each point is one module in one cell type. Dominance = share of module expression from its\nlargest gene. Agreement = correlation between flux share and module expression share.",
       x = "single-gene dominance of module expression",
       y = "agreement between flux share and expression share") +
  theme_bw(base_size = 10)
ggsave(file.path(OUT, "fig_dominance.png"), p, width = 10, height = 6.5, dpi = 150, bg = "white")
say("\n-> %s", OUT)
