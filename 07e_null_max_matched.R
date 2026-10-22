#!/usr/bin/env Rscript
# ==============================================================================
# 07e_null_max_matched.R  --  a null that models the statistic we actually test
#
#   Rscript 07e_null_max_matched.R --arm bulk_before
#
# THE PROBLEM. We call a module grounded if its BEST-correlating catalytic gene
# beats a threshold. That statistic is a maximum over the module's 10 to 15
# genes. The null in 07c is built from correlations of INDIVIDUAL random genes,
# so it models a single draw, not a maximum over k draws.
#
# A maximum over k draws is stochastically larger than one draw. Concretely, if
# the threshold is the 99.9th percentile of single draws, a null module of 15
# genes exceeds it with probability 1 - 0.999^15 = 1.5 percent, not 0.1 percent.
# The 99.9th percentile makes the test stringent but it does not make it a test
# of the right quantity.
#
# THE FIX. For each module, draw k random-gene correlations where k is that
# module's gene count, take the maximum, and repeat. The 99.9th percentile of
# those maxima is the matched threshold. The statistic in the null now has the
# same form as the statistic being tested.
#
# TWO POOLS ARE COMPARED, because the choice is not obvious:
#   celltype  draw from all random-gene correlations for that cell type. More
#             distinct values, but assumes the null is exchangeable across
#             modules within a cell type, which is what 07c already assumes.
#   module    draw from that module's own 60 random-gene correlations. Correct
#             in principle, but 60 distinct values limits the tail estimate.
#
# This script does not overwrite 07c. It writes its own thresholds and a
# comparison, so the change in status of every module is visible before anything
# downstream is rerun.
#
# Output -> results/all_sites/<ARM>/07e_null_max/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

OUT <- file.path(PATHS$results, "07e_null_max")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
B         <- 20000   # resamples per module
NULL_PCTL <- 0.999
MIN_GENES <- 2
set.seed(1)

R  <- read_csv(file.path(PATHS$results, "07b_grounding_scalefree/gene_level_scalefree.csv"),
               show_col_types = FALSE) %>% filter(!circular)
NG <- read_csv(file.path(PATHS$results, "07b_grounding_scalefree/negative_control_scalefree.csv"),
               show_col_types = FALSE) %>% filter(n_module_genes >= MIN_GENES)
say("null correlations available: %s", format(nrow(NG), big.mark = ","))

# a max-of-k threshold at percentile P needs the pool to support the quantile
# P^(1/k). Report where that fails rather than letting the resampling saturate
# silently at the pool maximum, which is what a smaller pool did.
say("\n=== does the pool support the required quantile? ===")
for (ct in sort(unique(NG$celltype))) {
  n <- sum(NG$celltype == ct)
  ks <- c(2, 5, 10, 15, 25, 40)
  ok <- ks[NULL_PCTL^(1/ks) < 1 - 1/n]
  say("  %-11s n = %6s | supported k: %s", ct, format(n, big.mark = ","),
      if (length(ok)) paste(ok, collapse = ", ") else "none")
}

# the tested statistic: best catalytic gene per (celltype, module)
best <- R %>% group_by(celltype, module) %>%
  summarise(k = first(n_module_genes), n_genes_tested = n(),
            best_gene = gene[which.max(rho)], best_rho = max(rho, na.rm = TRUE),
            .groups = "drop")
say("modules with a tested best gene: %d across %d cell types",
    nrow(best), n_distinct(best$celltype))

# the old threshold, for comparison
old <- NG %>% group_by(celltype) %>%
  summarise(old_bar = quantile(rho, NULL_PCTL, na.rm = TRUE), n_null = n(), .groups = "drop")

# ---- matched null: max of k draws -------------------------------------------
pool_ct <- split(NG$rho, NG$celltype)
pool_md <- split(NG$rho, paste(NG$celltype, NG$module, sep = "|"))

max_thresh <- function(pool, k, B) {
  if (length(pool) < 2 || is.na(k) || k < 2) return(NA_real_)
  m <- matrix(sample(pool, k * B, replace = TRUE), nrow = k)
  unname(quantile(apply(m, 2, max), NULL_PCTL, na.rm = TRUE))
}

say("\nresampling %s maxima per module...", format(B, big.mark = ","))
res <- best %>%
  rowwise() %>%
  # k is the module's gene count, but the maximum is taken over the genes that
  # actually had a correlation computed. Where those differ (36 of 304 modules,
  # almost always by one gene) the null must use the tested count, otherwise it
  # draws more values than the statistic maximises over and is conservative.
  mutate(bar_ct = max_thresh(pool_ct[[celltype]], n_genes_tested, B),
         bar_md = max_thresh(pool_md[[paste(celltype, module, sep = "|")]],
                             n_genes_tested, B)) %>%
  ungroup() %>%
  left_join(old, by = "celltype") %>%
  mutate(pass_old = best_rho > old_bar,
         pass_ct  = best_rho > bar_ct,
         pass_md  = best_rho > bar_md)

say("\nmodules where the module gene count differs from the tested count: %d of %d",
    sum(res$k != res$n_genes_tested), nrow(res))
say("  the null uses the tested count, which is what the maximum is taken over")

say("\n=== thresholds by cell type ===")
print(as.data.frame(res %>% group_by(celltype) %>%
  summarise(modules = n(), median_k = median(n_genes_tested),
            `old (single draw)` = round(first(old_bar), 3),
            `matched, celltype pool` = round(median(bar_ct), 3),
            `matched, module pool`   = round(median(bar_md), 3), .groups = "drop")),
  row.names = FALSE)

say("\n=== how many modules pass ===")
say("  old threshold           : %d of %d", sum(res$pass_old), nrow(res))
say("  matched, celltype pool  : %d of %d", sum(res$pass_ct,  na.rm = TRUE), nrow(res))
say("  matched, module pool    : %d of %d", sum(res$pass_md,  na.rm = TRUE), nrow(res))

say("\n=== modules that change status (celltype pool) ===")
chg <- res %>% filter(pass_old != pass_ct) %>% arrange(celltype, desc(best_rho))
if (nrow(chg)) {
  print(as.data.frame(chg %>% transmute(celltype, module, k, best_gene,
        best_rho = round(best_rho, 3), old = round(old_bar, 3),
        matched = round(bar_ct, 3),
        change = ifelse(pass_old & !pass_ct, "now fails", "now passes"))), row.names = FALSE)
} else say("  none")

say("\n=== the modules the paper depends on ===")
KEY <- tibble(celltype = c("Myeloid","Myeloid","Myeloid","T","Fibroblast"),
              module   = c("M_78","M_6","M_8","M_3","M_78"))
print(as.data.frame(res %>% inner_join(KEY, by = c("celltype","module")) %>%
  transmute(celltype, module, k = n_genes_tested, best_gene, best_rho = round(best_rho, 3),
            old = round(old_bar, 3), matched_ct = round(bar_ct, 3),
            matched_md = round(bar_md, 3),
            verdict = ifelse(pass_ct, "still meets", "NO LONGER MEETS"))), row.names = FALSE)

write_csv(res, file.path(OUT, "null_max_matched_thresholds.csv"))
say("\n-> %s", OUT)
say("\nIf the counts differ materially, 07c's null_bar should be replaced by the")
say("matched threshold and everything downstream of the criterion rerun.")
