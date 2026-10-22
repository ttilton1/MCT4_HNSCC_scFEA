#!/usr/bin/env Rscript
# ==============================================================================
# 01f_flux_expression_coupling.R
#   how tightly does total flux track total metabolic-gene expression, in OUR data?
#
#   Rscript 01f_flux_expression_coupling.R --arm bulk_before
#
# WHY THIS IS NEEDED. The draft currently asserts, from reading scFEA's loss
# function, that per-cell total flux is "anchored" to total metabolic-gene
# expression and that magnitude comparisons therefore "largely restate" an
# expression difference. That argument is vulnerable in three ways, and the
# method's authors would raise all three:
#
#   1. The expression-coherence term is a SOFT penalty, one of four, weighted
#      equally with flux balance. The optimum trades off between them. "Anchored"
#      overstates what a soft constraint guarantees.
#   2. We inferred a data property from source code rather than measuring it.
#   3. Flux tracking enzyme abundance is the biological premise of the method,
#      not a defect. Calling it redundant misreads a deliberate design choice.
#
# The defensible version of the claim is empirical and specific to our data:
# report the observed correlation between per-unit total flux and per-unit total
# metabolic-gene expression, and state what follows for the comparison we make.
# If that correlation is very high, "a magnitude comparison closely tracks an
# expression comparison" is an observation, not an interpretation of the loss.
# If it is moderate, the claim must be weakened accordingly.
#
# OUTPUTS
#   A. rho(total flux, total metabolic-gene expression) per cell type
#   B. the same for the carcinoma-vs-normal CONTRAST: does the fold-change in
#      total flux track the fold-change in expression sum?
#   C. variance in total flux explained by expression sum (R^2)
#
# Output -> results/all_sites/<ARM>/01f_coupling/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })

OUT <- file.path(PATHS$results, "01f_coupling")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

NORMAL_T <- c("NL","NL-contact","Normal")
CANCER_T <- c("CA","CA-E","CA-A","CA-c","CA-p")

# ---- flux totals per unit ----------------------------------------------------
cw <- as.data.frame(readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))) %>%
  filter(n_cells >= GATES$min_cells)
mc <- mcols_of(cw)
cw$total_flux <- rowSums(as.matrix(cw[, mc]))

# ---- expression totals per unit, from the matrix scFEA consumed --------------
say("reading the expression matrices supplied to scFEA...")
expr_tot <- list()
for (ct in TYPES) {
  ep <- file.path(ROOT, "data/input/all_sites/bulk_before", ct, "train_expression.csv")
  mp <- file.path(ROOT, "metadata/all_sites/bulk_before", ct, "metadata_units.csv")
  if (!file.exists(ep)) next
  e <- fread(ep); g <- as.character(e[[1]])
  E <- as.matrix(e[, -1]); rownames(E) <- g
  mu <- fread(mp); setnames(mu, 1, "id"); mu <- as.data.frame(mu)
  rownames(mu) <- as.character(mu$id); mu <- mu[colnames(E), ]
  expr_tot[[ct]] <- tibble(patient = mu$patient, tissue = mu$tissue,
                           finaltype = mu$finaltype, n = mu$n_cells,
                           expr_sum = colSums(E))
}
EX <- bind_rows(expr_tot) %>% group_by(patient, tissue, finaltype) %>%
  summarise(expr_sum = weighted.mean(expr_sum, n), .groups = "drop")

D <- cw %>% select(patient, tissue, finaltype, annot_4th, study, total_flux) %>%
  inner_join(EX, by = c("patient","tissue","finaltype"))
say("units with both flux and expression totals: %s", format(nrow(D), big.mark=","))

# ============ A: per-cell-type coupling ======================================
say("\n=== A: rho(total flux, total metabolic-gene expression) ===")
a <- D %>% group_by(annot_4th) %>%
  summarise(n = n(),
            rho = cor(total_flux, expr_sum, method = "spearman"),
            pearson = cor(total_flux, expr_sum),
            r2 = summary(lm(total_flux ~ expr_sum))$r.squared, .groups = "drop") %>%
  arrange(desc(rho))
print(as.data.frame(a %>% mutate(across(where(is.numeric), ~round(.x, 3)))),
      row.names = FALSE)
overall <- cor(D$total_flux, D$expr_sum, method = "spearman")
say("  pooled across cell types: rho = %.3f", overall)
say("\n  READ: this is the number the claim must rest on. If rho > 0.9, a")
say("  magnitude comparison closely tracks an expression comparison and can be")
say("  stated as an observation. If it is moderate, weaken the claim.")

# ============ B: does the CONTRAST track too? ================================
say("\n=== B: carcinoma vs normal, fold change in flux vs in expression ===")
fc <- D %>% mutate(grp = case_when(tissue %in% NORMAL_T ~ "normal",
                                   tissue %in% CANCER_T ~ "cancer", TRUE ~ NA_character_)) %>%
  filter(!is.na(grp)) %>%
  group_by(annot_4th, grp) %>%
  summarise(flux = median(total_flux), expr = median(expr_sum), .groups="drop") %>%
  pivot_wider(names_from = grp, values_from = c(flux, expr)) %>%
  mutate(fc_flux = flux_cancer / flux_normal,
         fc_expr = expr_cancer / expr_normal)
print(as.data.frame(fc %>% transmute(celltype = annot_4th,
      fc_flux = round(fc_flux, 3), fc_expr = round(fc_expr, 3),
      ratio = round(fc_flux / fc_expr, 3))), row.names = FALSE)
if (sum(!is.na(fc$fc_flux)) >= 4) {
  r_fc <- cor(fc$fc_flux, fc$fc_expr, use = "complete.obs")
  say("  correlation of the two fold-changes across cell types: r = %.3f", r_fc)
}
say("  READ: ratio near 1 means the flux fold-change is the expression")
say("  fold-change. That is the specific claim being made in Section 1.")

write_csv(a, file.path(OUT, "coupling_by_celltype.csv"))
write_csv(fc, file.path(OUT, "foldchange_comparison.csv"))

# ---- figure ------------------------------------------------------------------
p <- ggplot(D, aes(expr_sum, total_flux)) +
  geom_point(alpha = .25, size = .7) +
  geom_smooth(method = "lm", colour = "#b2182b", se = FALSE) +
  facet_wrap(~ annot_4th, scales = "free") +
  labs(title = "Total flux against total metabolic-gene expression, per unit",
       subtitle = sprintf("pooled Spearman rho = %.2f. The expression-coherence term in scFEA's objective\npenalises deviation between these two quantities.", overall),
       x = "sum of metabolic-gene expression", y = "sum of flux across 168 modules") +
  theme_bw(base_size = 10)
ggsave(file.path(OUT, "fig_coupling.png"), p, width = 11, height = 7, dpi = 150, bg = "white")
say("\n-> %s", OUT)
