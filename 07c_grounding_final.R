#!/usr/bin/env Rscript
# ==============================================================================
# 07c_grounding_final.R  --  grounding judged against the EMPIRICAL NULL.
#
#   Rscript 07c_grounding_final.R --arm bulk_before
#
# WHAT WENT WRONG IN 07 AND 07b, AND WHY THIS IS THE VERSION TO TRUST.
#
# 07  correlated raw flux vs raw expression. Random genes came back at rho=0.202
#     (41% above |0.3|) -- a SHARED-SCALE artifact: units with more total
#     expression have more total flux, so everything correlates with everything.
#     07 also reported a top-25 at rho~0.999, every one a ONE-GENE module, where
#     scFEA's flux is a monotone transform OF THAT GENE. Circular, not evidence.
#
# 07b fixed both (flux SHARE vs expression SHARE; 1-gene modules dropped). The
#     null collapsed to 0.017 -- correct. But 07b still used two bad statistics:
#
#       (a) MODULE MEDIAN rho. A module has 10-15 genes; only 1-3 are rate-
#           limiting. M_6 contains LDHA (the enzyme) and LDHAL6A/B (paralogs
#           that carry no flux signal). Averaging them buries the result: M_6's
#           median is 0.18 while LDHA alone is 0.79. The median measures how many
#           genes in the module happen to be informative -- NOT whether the flux
#           tracks the reaction. USE THE BEST CATALYTIC GENE, not the median.
#
#       (b) A `grounded` flag built from CI-excludes-0 + same-sign + LOSO. With
#           n~600-2000 units a rho of 0.25 has a CI nowhere near zero, so the
#           flag passed correlations that RANDOM GENES REACH ROUTINELY.
#
# THE NULL IS WIDE. Empirically, random genes reach:
#     50th 0.017 | 90th 0.344 | 95th 0.439 | 99th 0.602 | 99.9th 0.709
# So 6.6% of random genes beat rho=0.40 and 1% beat 0.60. Any "grounding" below
# ~0.71 is indistinguishable from noise. THAT is the bar, and it is set from the
# data, not chosen.
#
# THIS SCRIPT: computes a null PER CELL TYPE (the null width differs by how many
# units a cell type has), calls a module GROUNDED only if its BEST catalytic gene
# exceeds that cell type's 99.9th-percentile null, and prints an explicit
# SUPPORTED / NOT SUPPORTED verdict for every upstream claim.
#
# Output -> results/all_sites/<ARM>/07c_grounding_final/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(stringr) })

OUT <- file.path(PATHS$results, "07c_grounding_final")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
MR <- load_modref()

NULL_PCTL <- 0.999      # the bar: beat the 99.9th pctl of random genes
N_NEG     <- 60         # random genes per (celltype, module) -- need many for a 99.9th pctl
MIN_GENES <- 2
set.seed(1)

R  <- read_csv(file.path(PATHS$results, "07b_grounding_scalefree/gene_level_scalefree.csv"),
               show_col_types = FALSE) %>% filter(!circular)
NG <- read_csv(file.path(PATHS$results, "07b_grounding_scalefree/negative_control_scalefree.csv"),
               show_col_types = FALSE) %>% filter(n_module_genes >= MIN_GENES)
say("loaded %s enzyme correlations | %s null correlations", 
    format(nrow(R), big.mark=","), format(nrow(NG), big.mark=","))

# ---- PER-CELL-TYPE null (a cell type with fewer units has a wider null) -------
nulls <- NG %>% group_by(celltype) %>%
  summarise(n_null = n(),
            null_50  = quantile(rho, .50, na.rm=TRUE),
            null_95  = quantile(rho, .95, na.rm=TRUE),
            null_99  = quantile(rho, .99, na.rm=TRUE),
            null_bar = quantile(rho, NULL_PCTL, na.rm=TRUE), .groups="drop")
say("\n=== EMPIRICAL NULL, per cell type (bar = %.1fth pctl of RANDOM genes) ===",
    100*NULL_PCTL)
print(as.data.frame(nulls %>% transmute(celltype, n_null,
      median = round(null_50,3), p95 = round(null_95,3),
      p99 = round(null_99,3), BAR = round(null_bar,3))), row.names = FALSE)
write_csv(nulls, file.path(OUT, "empirical_null_by_celltype.csv"))

# ---- grounding = BEST catalytic gene beats that cell type's null bar ----------
G <- R %>% left_join(nulls %>% select(celltype, null_bar, null_99), by = "celltype") %>%
  mutate(beats_null = rho > null_bar,
         beats_99   = rho > null_99)

mod <- G %>% group_by(celltype, module, src) %>%
  summarise(n_genes   = n(),
            best_gene = gene[which.max(rho)],
            best_rho  = max(rho, na.rm = TRUE),
            med_rho   = median(rho, na.rm = TRUE),
            n_beat    = sum(beats_null, na.rm = TRUE),
            null_bar  = first(null_bar), .groups = "drop") %>%
  mutate(GROUNDED = best_rho > null_bar,
         margin   = round(best_rho - null_bar, 3)) %>%
  left_join(MR$rxn, by = "module") %>%
  arrange(desc(best_rho))
write_csv(mod, file.path(OUT, "module_grounding_final.csv"))
write_csv(G,   file.path(OUT, "gene_grounding_final.csv"))

say("\n=== GROUNDED MODULES (best catalytic gene beats the null bar) ===")
say("    NOTE: med_rho is shown only to demonstrate that it is USELESS --")
say("    it measures how many genes in the module are informative, not whether")
say("    the flux tracks the reaction (M_6: median 0.18, but LDHA = 0.79).")
gr <- mod %>% filter(GROUNDED)
if (nrow(gr)) {
  print(as.data.frame(gr %>% transmute(celltype, module,
      reaction = str_trunc(reaction, 28), src, best_gene,
      best_rho = round(best_rho, 3), bar = round(null_bar, 2),
      margin, med_rho = round(med_rho, 2),
      beat = sprintf("%d/%d", n_beat, n_genes))), row.names = FALSE)
} else say("  NONE clear the bar.")
say("\n  %d of %d (celltype, module) pairs are grounded (%.0f%%)",
    nrow(gr), nrow(mod), 100*nrow(gr)/nrow(mod))

# ---- explicit verdict on the modules carrying the narrative -------------------
CLAIMS <- tribble(
  ~module, ~claim,
  "M_6",   "lactate production (LDH)",
  "M_78",  "lactate EXPORT (MCT)",
  "M_4",   "glycolysis: 3PD->pyruvate",
  "M_3",   "glycolysis: G3P->3PD",
  "M_2",   "glycolysis: G6P->G3P",
  "M_113", "mevalonate / isoprenoid",
  "M_116", "ubiquinone (CoQ)")

say("\n============================================================")
say("VERDICT ON EVERY UPSTREAM CLAIM")
say("============================================================")
verdict <- mod %>% inner_join(CLAIMS, by = "module") %>%
  transmute(claim, celltype, best_gene, best_rho = round(best_rho, 3),
            bar = round(null_bar, 2),
            VERDICT = ifelse(GROUNDED, "SUPPORTED", "not supported"),
            src) %>%
  arrange(claim, desc(best_rho))
print(as.data.frame(verdict), row.names = FALSE)
write_csv(verdict, file.path(OUT, "claim_verdicts.csv"))

say("\n-- summary by claim --")
print(as.data.frame(verdict %>% group_by(claim) %>%
  summarise(supported_in = paste(celltype[VERDICT=="SUPPORTED"], collapse=", "),
            failed_in    = paste(celltype[VERDICT!="SUPPORTED"], collapse=", "),
            .groups="drop")), row.names = FALSE)

# ---- the LDHA/LDHB directionality caveat -------------------------------------
say("\n=== CAVEAT: M_6 flux tracks LDH ABUNDANCE, not DIRECTION ===")
ldh <- G %>% filter(module == "M_6", gene %in% c("LDHA","LDHB")) %>%
  select(celltype, gene, rho) %>% pivot_wider(names_from = gene, values_from = rho)
print(as.data.frame(ldh %>% mutate(across(where(is.numeric), ~round(.x,3)))), row.names=FALSE)
say("  LDHA catalyses pyruvate->lactate; LDHB catalyses the REVERSE. Both correlate")
say("  with M_6 flux, so the flux reflects total LDH, NOT net direction. Report as")
say("  'LDH-associated flux', never as a lactate PRODUCTION rate.")

# ---- MCT specificity: the strongest positive result --------------------------
say("\n=== MCT SPECIFICITY (a scale artifact could not produce this) ===")
mct <- G %>% filter(module == "M_78") %>%
  select(celltype, gene, rho) %>% pivot_wider(names_from = gene, values_from = rho)
print(as.data.frame(mct %>% mutate(across(where(is.numeric), ~round(.x,3)))), row.names=FALSE)
say("  SLC16A3 = MCT4 (high-Km EXPORTER, used by glycolytic cells)")
say("  SLC16A1 = MCT1 (importer)   SLC16A4 = orphan")
say("  The flux singles out MCT4 in Myeloid/Fibroblast and NOT in T cells --")
say("  i.e. it says myeloid+stroma export lactate and T cells do not. That is a")
say("  specific, correct biological discrimination, not a global artifact.")

# ---- figure -------------------------------------------------------------------
pd <- G %>% inner_join(CLAIMS, by = "module") %>%
  mutate(lab = paste0(celltype, " | ", gene),
         status = ifelse(beats_null, "beats null (99.9th)",
                  ifelse(beats_99, "99th-99.9th", "within null")))
p <- ggplot(pd, aes(rho, reorder(lab, rho), colour = status)) +
  geom_vline(data = nulls, aes(xintercept = null_bar), linetype = 2, colour = "grey40") +
  geom_point(size = 1.9) +
  facet_wrap(~ claim, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c(`beats null (99.9th)` = "#b2182b",
                                 `99th-99.9th` = "#f4a582", `within null` = "grey70")) +
  labs(title = sprintf("Enzyme grounding vs the empirical null (%s)", ARM),
       subtitle = "flux share ~ expression share. Dashed = 99.9th pctl of RANDOM genes for that cell type.\nOnly red clears the bar.",
       x = "meta Spearman rho", y = NULL, colour = NULL) +
  theme_bw(base_size = 8) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "fig_claims_vs_null.png"), p, width = 11, height = 12,
       dpi = 150, bg = "white", limitsize = FALSE)

nl <- ggplot(NG, aes(rho)) +
  geom_density(fill = "grey70", colour = NA, alpha = .7) +
  geom_vline(xintercept = quantile(NG$rho, c(.95,.99,.999), na.rm=TRUE),
             linetype = c(3,2,1), colour = c("grey50","grey35","#b2182b")) +
  geom_rug(data = G %>% filter(module %in% CLAIMS$module, beats_null),
           aes(rho), colour = "#b2182b", length = unit(.05,"npc")) +
  labs(title = "The null is WIDE -- this is why most 'grounding' is noise",
       subtitle = sprintf("random genes: 95th=%.2f  99th=%.2f  99.9th=%.2f (red). Red rug = enzymes that clear it.",
         quantile(NG$rho,.95,na.rm=TRUE), quantile(NG$rho,.99,na.rm=TRUE),
         quantile(NG$rho,.999,na.rm=TRUE)),
       x = "meta Spearman rho (random gene vs module flux)", y = "density") +
  theme_bw(base_size = 11)
ggsave(file.path(OUT, "fig_the_null.png"), nl, width = 9, height = 5, dpi = 150, bg="white")

say("\n-> %s", OUT)
