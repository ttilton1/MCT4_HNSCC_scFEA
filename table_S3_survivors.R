#!/usr/bin/env Rscript
# ==============================================================================
# table_S3_survivors.R  --  Supplementary Table S3
#
#   Rscript table_S3_survivors.R --arm bulk_before
#
# All nineteen comparisons surviving the progression filters, across every cell
# type, with per-cohort effect sizes and flux-enzyme criterion status. The paper
# follows two of them (M_78 and M_6 in cDC2); this table exists so the others are
# visible rather than absent, and so a reader can see which were set aside and on
# what basis.
#
# Per-cohort deltas are shown because cross-cohort concordance was required to
# survive, which constrains direction but not magnitude. Several T-cell results
# are much larger in Kim than in He, and that is visible here rather than hidden
# inside a meta-analysed effect.
#
# Output -> results/all_sites/<ARM>/tables/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages(library(dplyr))

OUT <- file.path(PATHS$results, "tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

meta <- read_csv(file.path(PATHS$results, "05c_progression_meta/finaltype/meta_all.csv"),
                 show_col_types = FALSE)
gr   <- read_csv(file.path(PATHS$results, "07c_grounding_final/module_grounding_final.csv"),
                 show_col_types = FALSE) %>%
        select(celltype, module, best_gene, best_rho, null_bar, GROUNDED)

tab <- meta %>% filter(survivor) %>%
  left_join(gr, by = c("celltype","module")) %>%
  mutate(criterion = case_when(is.na(GROUNDED) ~ "not tested",
                               GROUNDED        ~ "meets",
                               TRUE            ~ "does not meet"),
         in_paper = ifelse(unit == "DC_C2_cDC2" & module %in% c("M_78","M_6"),
                           "yes", "")) %>%
  arrange(factor(criterion, levels = c("meets","does not meet","not tested")),
          desc(abs(delta_meta))) %>%
  transmute(`cell type` = celltype, subtype = unit, module, reaction,
            `Δ meta` = round(delta_meta, 4),
            `Δ He` = round(delta_He, 4), `Δ Kim` = round(delta_Kim, 4),
            `n He` = n_pt_He, `n Kim` = n_pt_Kim,
            FDR = signif(fdr, 2),
            `best gene` = best_gene, `ρ` = round(best_rho, 3),
            threshold = round(null_bar, 2), criterion,
            `in paper` = in_paper)

say("=== Supplementary Table S3: %d progression survivors ===", nrow(tab))
print(as.data.frame(tab), row.names = FALSE)
write_csv(tab, file.path(OUT, "table_S3_progression_survivors.csv"))

say("\nby criterion status:")
print(as.data.frame(tab %>% count(criterion)), row.names = FALSE)
say("\nNK results are present in this table for completeness but are not")
say("interpreted anywhere: the NK empirical null rests on 71 draws, too few to")
say("estimate a 99.9th percentile.")
say("\n-> %s", file.path(OUT, "table_S3_progression_survivors.csv"))
