#!/usr/bin/env Rscript
# ==============================================================================
# check_threshold_sensitivity.R
#   does the choice of 99.9th percentile change what is retained?
#
#   Rscript check_threshold_sensitivity.R --arm bulk_before
#
# The threshold is a judgement call. With 254 comparisons, the 99.9th percentile
# is roughly Bonferroni-like; the 99th would admit around 2.5 expected false
# positives and might retain results currently set aside. Rather than defending
# one cut, report how much turns on it.
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages(library(readr))

n <- read_csv(file.path(PATHS$results, "07c_grounding_final/empirical_null_by_celltype.csv"),
              show_col_types = FALSE)
g <- read_csv(file.path(PATHS$results, "07c_grounding_final/module_grounding_final.csv"),
              show_col_types = FALSE)
# module_grounding_final already carries null_bar; only the lower percentiles
# need joining, or the join silently creates null_bar.x / null_bar.y
g <- g %>% left_join(n %>% select(celltype, null_95, null_99), by = "celltype")
say("rows in module_grounding_final: %d | GROUNDED = TRUE: %d",
    nrow(g), sum(g$GROUNDED, na.rm = TRUE))
say("  (the paper quotes 59 of 254 -- if these differ, trace which subset that came from)")
say("  src values present: %s", paste(unique(g$src), collapse = " | "))
say("  by src:")
print(as.data.frame(g %>% count(src, GROUNDED)), row.names = FALSE)

say("=== retained at each threshold, of %d cell-type x module pairs ===", nrow(g))
for (col in c("null_95","null_99","null_bar")) {
  lab <- c(null_95 = "95th", null_99 = "99th", null_bar = "99.9th")[col]
  say("  %-6s : %3d retained (%.0f%%)", lab, sum(g$best_rho > g[[col]], na.rm = TRUE),
      100*mean(g$best_rho > g[[col]], na.rm = TRUE))
}

say("\n=== the key modules at each threshold ===")
KEY <- c("M_78","M_6","M_2","M_113","M_115","M_116")
k <- g %>% filter(module %in% KEY, celltype %in% c("Myeloid","T","Fibroblast")) %>%
  transmute(celltype, module, best_gene, rho = round(best_rho,3),
            p95 = best_rho > null_95, p99 = best_rho > null_99, p999 = best_rho > null_bar)
print(as.data.frame(k), row.names = FALSE)

say("\n=== would any of the 19 progression survivors change status? ===")
s <- read_csv(file.path(PATHS$results, "05c_progression_meta/finaltype/meta_survivors.csv"),
              show_col_types = FALSE) %>%
     left_join(g %>% select(celltype, module, best_rho, null_95, null_99, null_bar),
               by = c("celltype","module"))
tab <- s %>% transmute(celltype, unit, module,
                       rho = round(best_rho,3),
                       at_99  = best_rho > null_99,
                       at_999 = best_rho > null_bar) %>%
  filter(!is.na(rho))
print(as.data.frame(tab[order(-tab$rho), ]), row.names = FALSE)
say("\n  retained at 99.9th: %d | at 99th: %d",
    sum(tab$at_999, na.rm=TRUE), sum(tab$at_99, na.rm=TRUE))
say("  changed by relaxing to the 99th percentile: %s",
    paste(with(tab, paste0(unit, " ", module))[which(tab$at_99 & !tab$at_999)], collapse = ", "))
