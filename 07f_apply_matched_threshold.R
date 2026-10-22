#!/usr/bin/env Rscript
# ==============================================================================
# 07f_apply_matched_threshold.R  --  what changes under the matched null
#
#   Rscript 07f_apply_matched_threshold.R --arm bulk_before
#
# 07e computes a threshold matched to the statistic actually tested: the maximum
# correlation across a module's genes, so the null takes the maximum of k random
# draws rather than a single draw. This script reports every number in the paper
# that the change touches, without overwriting anything.
#
# A SCOPE LIMIT WORTH STATING. A matched threshold at percentile P needs the null
# pool to support the quantile P^(1/k). With 300 random genes per module the pool
# supports the median module size in Myeloid and T cells but not in B, Epithelial,
# Fibroblast, Endothelial or NK, where the resampled maximum saturates at the pool
# maximum. Thresholds for those cell types are lower bounds, and this script marks
# them so they are not quoted as estimates.
#
# Output -> results/all_sites/<ARM>/07e_null_max/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(dplyr) })

OUT <- file.path(PATHS$results, "07e_null_max")
NULL_PCTL <- 0.999

th <- read_csv(file.path(OUT, "null_max_matched_thresholds.csv"), show_col_types = FALSE)
NG <- read_csv(file.path(PATHS$results, "07b_grounding_scalefree/negative_control_scalefree.csv"),
               show_col_types = FALSE)

# which cell types have a pool large enough for their own module sizes
support <- th %>% group_by(celltype) %>%
  summarise(n_null = sum(NG$celltype == first(celltype)),
            median_k = median(n_genes_tested),
            max_k = max(n_genes_tested), .groups = "drop") %>%
  mutate(needed_at_median = NULL_PCTL^(1/median_k),
         supported = needed_at_median < 1 - 1/n_null)
say("=== is the matched threshold estimable at each cell type's median module size? ===")
print(as.data.frame(support %>% transmute(celltype, n_null, median_k, max_k,
  `quantile needed` = round(needed_at_median, 5),
  `pool supports` = ifelse(supported, "yes", "NO - threshold is a lower bound"))),
  row.names = FALSE)

th <- th %>% left_join(support %>% select(celltype, supported), by = "celltype") %>%
  mutate(pass_matched = best_rho > bar_ct)

say("\n=== criterion counts ===")
say("  old single-draw threshold : %d of %d", sum(th$pass_old), nrow(th))
say("  matched threshold         : %d of %d", sum(th$pass_matched), nrow(th))
say("  matched, restricted to cell types where the pool supports it: %d of %d",
    sum(th$pass_matched & th$supported), sum(th$supported))

# ---- the 19 progression survivors ------------------------------------------
meta <- read_csv(file.path(PATHS$results, "05c_progression_meta/finaltype/meta_all.csv"),
                 show_col_types = FALSE) %>% filter(survivor)
sv <- meta %>%
  left_join(th %>% select(celltype, module, best_gene, best_rho, old_bar, bar_ct,
                          pass_old, pass_matched, supported),
            by = c("celltype","module")) %>%
  mutate(old_status = ifelse(is.na(pass_old), "not tested",
                      ifelse(pass_old, "meets", "does not meet")),
         new_status = ifelse(is.na(pass_matched), "not tested",
                      ifelse(pass_matched, "meets", "does not meet")),
         changed = old_status != new_status) %>%
  arrange(desc(changed), desc(abs(delta_meta)))

say("\n=== the 19 progression survivors under each threshold ===")
print(as.data.frame(sv %>% transmute(celltype, unit, module,
  delta = round(delta_meta, 4), gene = best_gene, rho = round(best_rho, 3),
  old = round(old_bar, 3), matched = round(bar_ct, 3),
  old_status, new_status, pool_ok = supported)), row.names = FALSE)

say("\n  meeting the criterion, old threshold     : %d of %d",
    sum(sv$old_status == "meets"), nrow(sv))
say("  meeting the criterion, matched threshold : %d of %d",
    sum(sv$new_status == "meets"), nrow(sv))
if (any(sv$changed)) {
  say("\n  survivors changing status:")
  for (i in which(sv$changed))
    say("    %s %s %s: %s -> %s (rho %.3f, old %.3f, matched %.3f)",
        sv$celltype[i], sv$unit[i], sv$module[i], sv$old_status[i], sv$new_status[i],
        sv$best_rho[i], sv$old_bar[i], sv$bar_ct[i])
}
write_csv(sv, file.path(OUT, "survivors_under_matched_threshold.csv"))

# ---- M_78 by cell type, which Figure 5A draws ------------------------------
say("\n=== M_78 across cell types, which Figure 5A shows ===")
print(as.data.frame(th %>% filter(module == "M_78") %>%
  transmute(celltype, k = n_genes_tested, gene = best_gene,
            rho = round(best_rho, 3), old = round(old_bar, 3),
            matched = round(bar_ct, 3),
            old_status = ifelse(pass_old, "meets", "fails"),
            new_status = ifelse(pass_matched, "meets", "fails"),
            pool_ok = supported) %>% arrange(desc(rho))), row.names = FALSE)

say("\n-> %s/survivors_under_matched_threshold.csv", OUT)
