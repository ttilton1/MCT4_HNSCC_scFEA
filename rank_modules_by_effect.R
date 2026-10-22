#!/usr/bin/env Rscript
# ==============================================================================
# rank_modules_by_effect.R  --  where lactate export sits among all modules
#
#   Rscript rank_modules_by_effect.R --arm bulk_before
#   Rscript rank_modules_by_effect.R --arm bulk_after
#
# WHY THIS EXISTS. With ten paired patients the Wilcoxon signed-rank test floors
# at p = 0.002, which every module returns if all ten patients move the same way
# regardless of how large the move is. In the myeloid panel 93 of 162 modules sit
# at that floor and 120 clear FDR < 0.05. So a p-value says only that the
# direction was consistent, and quoting one for a single module overstates what
# it establishes.
#
# The informative quantity is effect size among modules that rise in every
# patient, and whether that ranking survives the second pseudobulking arm.
#
# THE RESULT THIS PRODUCES. Lactate export (M_78) ranks second in both arms:
# 2 of 33 in bulk_before, 2 of 18 in bulk_after. It is the highest-ranked such
# module meeting the flux-enzyme criterion. The module ranked above it differs
# between arms -- M_4 (3PD to pyruvate) in bulk_before, M_116 (farnesyl-PP to
# ubiquinone) in bulk_after -- and neither meets the criterion in myeloid
# (PGAM1 0.534 and COQ6 0.548, against a 0.70 threshold).
#
# Writes a CSV for plotting later.
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(dplyr) })

OUT <- file.path(PATHS$results, "module_ranks")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
NL <- "NL"; CA <- c("CA","CA-E","CA-A","CA-c","CA-p")
u  <- readRDS(file.path(PATHS$results, "03_qc", "usable_modules.rds"))
mt <- read_csv(file.path(PATHS$results, "05c_progression_meta/finaltype/meta_all.csv"),
               show_col_types = FALSE) %>% distinct(module, reaction, supermodule)
gr <- read_csv(file.path(PATHS$results, "07c_grounding_final/module_grounding_final.csv"),
               show_col_types = FALSE) %>% filter(celltype == "Myeloid") %>%
      select(module, best_gene, best_rho, null_bar, GROUNDED)

cw <- as.data.frame(readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))) %>%
  filter(annot_4th == "Myeloid", n_cells >= GATES$min_cells,
         study %in% PROG_COHORTS, tissue %in% c(NL, CA)) %>%
  mutate(stage = ifelse(tissue %in% NL, "normal", "carcinoma"))
mc <- intersect(mcols_of(cw), u)
S  <- as.matrix(cw[, mc]); S <- S / rowSums(S)

res <- lapply(mc, function(m) {
  d <- cw %>% mutate(v = 100 * S[, m]) %>%
    group_by(patient, stage) %>%
    summarise(v = weighted.mean(v, n_cells), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = stage, values_from = v) %>%
    filter(!is.na(normal), !is.na(carcinoma))
  if (nrow(d) < 4) return(NULL)
  tibble(module = m, n = nrow(d),
         normal = median(d$normal), carcinoma = median(d$carcinoma),
         change = median(d$carcinoma) - median(d$normal),
         n_rising = sum(d$carcinoma > d$normal),
         p = suppressWarnings(wilcox.test(d$carcinoma, d$normal, paired = TRUE))$p.value)
}) %>% bind_rows() %>%
  mutate(FDR = p.adjust(p, "BH")) %>%
  left_join(mt, by = "module") %>% left_join(gr, by = "module") %>%
  mutate(meets_criterion = ifelse(is.na(GROUNDED), NA, GROUNDED),
         consistent = n_rising == n)

say("=== ARM %s: %d modules ===", ARM, nrow(res))
say("  at the Wilcoxon floor p = %.4f: %d | FDR < 0.05: %d",
    min(res$p), sum(res$p == min(res$p)), sum(res$FDR < 0.05))
say("  -> a p-value here says the direction was consistent, nothing more")

up <- res %>% filter(consistent, change > 0) %>% arrange(desc(change)) %>%
  mutate(rank = row_number())
say("\n  modules rising in ALL %d patients: %d", max(res$n), nrow(up))
print(as.data.frame(up %>% slice_head(n = 12) %>%
  transmute(rank, module, reaction = substr(reaction, 1, 32),
            change = round(change, 3), gene = best_gene,
            rho = round(best_rho, 3), meets = meets_criterion)), row.names = FALSE)

r78 <- up$rank[up$module == "M_78"]
say("\n  M_78 rank %d of %d", r78, nrow(up))
first_ok <- up %>% filter(meets_criterion) %>% slice_head(n = 1)
say("  highest-ranked module meeting the criterion: %s (rank %d)",
    first_ok$module, first_ok$rank)
above <- up %>% filter(rank < r78)
if (nrow(above)) {
  say("  ranked above M_78:")
  for (i in seq_len(nrow(above)))
    say("    %s  %s  change %.3f  %s (rho %.3f vs %.2f)  meets: %s",
        above$module[i], substr(above$reaction[i], 1, 26), above$change[i],
        above$best_gene[i], above$best_rho[i], above$null_bar[i],
        ifelse(is.na(above$meets_criterion[i]), "not tested", above$meets_criterion[i]))
}

write_csv(up,  file.path(OUT, sprintf("module_ranks_%s.csv", ARM)))
write_csv(res, file.path(OUT, sprintf("module_all_%s.csv", ARM)))
say("\n-> %s", OUT)
say("   module_ranks_%s.csv  consistent risers, ranked, for plotting", ARM)
say("   module_all_%s.csv    every module, for the p-value floor argument", ARM)
