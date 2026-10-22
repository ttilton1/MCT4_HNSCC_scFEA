#!/usr/bin/env Rscript
# ==============================================================================
# 07g_patch_matched_threshold.R  --  rewrite module_grounding_final.csv with the matched
#                          threshold, so everything downstream picks it up
#
#   Rscript 07g_patch_matched_threshold.R --arm bulk_before
#
# 07c computes null_bar as the 99.9th percentile of INDIVIDUAL random-gene
# correlations. The statistic actually tested is the maximum across a module's
# genes, so that threshold does not model the right quantity (see 07e). This
# script replaces null_bar with the matched threshold and recomputes GROUNDED.
#
# It preserves every column and the file path, so fig2, fig4, fig5,
# table_S3_survivors and make_supplementary_tables need no changes -- rerun them
# and they use the new criterion.
#
# The original is backed up alongside, since the change is not reversible from
# the output.
#
# Requires 07e_null_max_matched.R to have run.
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(dplyr) })

GR <- file.path(PATHS$results, "07c_grounding_final", "module_grounding_final.csv")
NB <- file.path(PATHS$results, "07c_grounding_final", "empirical_null_by_celltype.csv")
TH <- file.path(PATHS$results, "07e_null_max", "null_max_matched_thresholds.csv")
stopifnot(file.exists(GR), file.exists(TH))

g   <- read_csv(GR, show_col_types = FALSE)
th  <- read_csv(TH, show_col_types = FALSE) %>%
       select(celltype, module, n_genes_tested, bar_matched = bar_ct)

bak <- sub("\\.csv$", "_singledraw_backup.csv", GR)
if (!file.exists(bak)) { write_csv(g, bak); say("backed up original -> %s", basename(bak)) }

before <- sum(g$GROUNDED, na.rm = TRUE)
out <- g %>%
  left_join(th, by = c("celltype","module")) %>%
  mutate(null_bar_singledraw = null_bar,
         null_bar = ifelse(is.na(bar_matched), null_bar, bar_matched),
         k_tested = n_genes_tested,
         GROUNDED = best_rho > null_bar) %>%
  select(-bar_matched, -n_genes_tested)
after <- sum(out$GROUNDED, na.rm = TRUE)

say("=== module_grounding_final.csv rewritten with the matched threshold ===")
say("  rows: %d | meeting the criterion: %d -> %d", nrow(out), before, after)
say("  modules with no matched threshold available: %d",
    sum(is.na(out$k_tested)))

say("\n  by cell type:")
print(as.data.frame(out %>% group_by(celltype) %>%
  summarise(n = n(),
            single_draw = round(first(null_bar_singledraw), 3),
            matched_median = round(median(null_bar), 3),
            met_before = sum(best_rho > null_bar_singledraw, na.rm = TRUE),
            met_now = sum(GROUNDED, na.rm = TRUE), .groups = "drop")),
  row.names = FALSE)

write_csv(out, GR)

# the per-cell-type null file is also read by fig5's threshold segments; add the
# matched median so that panel can draw the right line
if (file.exists(NB)) {
  nb <- read_csv(NB, show_col_types = FALSE) %>%
    left_join(out %>% group_by(celltype) %>%
                summarise(null_bar_matched_median = median(null_bar), .groups="drop"),
              by = "celltype")
  write_csv(nb, NB)
  say("\n  empirical_null_by_celltype.csv gains null_bar_matched_median")
}

say("\nNow rerun, in any order:")
say("  fig2_allocation_panels.R, fig4_progression_all_panels.R,")
say("  fig5_transporter_all_panels.R, table_S3_survivors.R,")
say("  make_supplementary_tables.R")
say("\nTo revert: copy %s back over %s", basename(bak), basename(GR))
