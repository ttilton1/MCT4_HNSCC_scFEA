#!/usr/bin/env Rscript
# ==============================================================================
# check_m78_all_lineages.R  --  is the export rise myeloid-specific?
#
#   Rscript check_m78_all_lineages.R --arm bulk_before
#   Rscript check_m78_all_lineages.R --arm bulk_after
#
# THE QUESTION AND WHY IT MATTERS NOW. The paper says myeloid cells reprogram
# toward lactate export during progression. Two observations threaten the
# "specifically": FB_C1_MMP11 M_78 is among the nineteen progression survivors,
# and the fibroblast progression DEG has SLC16A3 at +2.92 log2FC with HIF1A up.
# If export allocation rises in every lineage, the finding is a tissue-wide
# reorganisation in which myeloid cells are highest, not a myeloid-specific one.
# That is still a claim, and arguably a larger one, but it is a different paper.
#
# WHAT THIS IS NOT. Figure 3A already shows myeloid highest in carcinoma against
# every lineage, paired within patient. That is a comparison of levels. This asks
# whether the CHANGE from normal to carcinoma is myeloid-specific, which nothing
# has tested.
#
# TWO CONSTRAINTS ON READING THE OUTPUT.
#
#   The criterion is per lineage. M_78 meets it in myeloid (SLC16A3 rho = 0.891),
#   fibroblast (0.766) and B cells (0.706), and fails in T cells, where the best
#   correlate is SLC16A1 at 0.480. It was never tested in epithelium. A lineage
#   failing the criterion cannot contribute an interpretable answer, and that is
#   part of the result rather than a gap in it.
#
#   Both arms. M_57 rose 10 of 10 in myeloid and fell in every other lineage in
#   bulk_before, then lost that specificity in bulk_after. A specificity claim
#   tested in one arm is worth little.
#
# THE THREE OUTCOMES.
#   Only myeloid rises            -> "myeloid specifically" holds, paper unchanged.
#   Myeloid rises most, others too -> reframe: the tissue reorganises lactate
#                                    handling across compartments, myeloid most.
#   All rise similarly            -> the myeloid emphasis is arbitrary and the
#                                    paper should rest on the level comparison
#                                    in Fig 3A, not on progression.
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(dplyr) })

NL <- "NL"; CA <- c("CA","CA-E","CA-A","CA-c","CA-p")
M  <- "M_78"
u  <- readRDS(file.path(PATHS$results, "03_qc", "usable_modules.rds"))
gr <- read_csv(file.path(PATHS$results, "07c_grounding_final/module_grounding_final.csv"),
               show_col_types = FALSE) %>% filter(module == M) %>%
      select(lineage = celltype, best_gene, best_rho, null_bar, GROUNDED)

cwa <- as.data.frame(readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))) %>%
  filter(n_cells >= GATES$min_cells, study %in% PROG_COHORTS, tissue %in% c(NL, CA)) %>%
  mutate(stage = ifelse(tissue %in% NL, "normal", "carcinoma"))
say("=== ARM %s | gate >= %d cells ===", ARM, GATES$min_cells)

out <- lapply(sort(unique(cwa$annot_4th)), function(lin) {
  cw <- cwa %>% filter(annot_4th == lin)
  mc <- intersect(mcols_of(cw), u)
  if (!M %in% mc) return(NULL)
  S  <- as.matrix(cw[, mc]); S <- S / rowSums(S)
  d  <- cw %>% mutate(v = 100 * S[, M]) %>%
    group_by(study, patient, stage) %>%
    summarise(v = weighted.mean(v, n_cells), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = stage, values_from = v)
  if (!all(c("normal","carcinoma") %in% names(d))) return(NULL)
  d <- d %>% filter(!is.na(normal), !is.na(carcinoma))
  if (nrow(d) < 4) return(NULL)
  per <- sapply(split(d, d$study), function(z)
    if (nrow(z) >= 2) median(z$carcinoma) - median(z$normal) else NA)
  tibble(lineage = lin, n = nrow(d),
         normal = median(d$normal), carcinoma = median(d$carcinoma),
         change = median(d$carcinoma) - median(d$normal),
         fold = median(d$carcinoma) / median(d$normal),
         rising = sum(d$carcinoma > d$normal),
         He = unname(per["He"]), Kim = unname(per["Kim"]),
         p = suppressWarnings(wilcox.test(d$carcinoma, d$normal, paired = TRUE))$p.value)
}) %>% bind_rows() %>%
  left_join(gr, by = "lineage") %>%
  mutate(criterion = ifelse(is.na(GROUNDED), "not tested",
                     ifelse(GROUNDED, "meets", "fails"))) %>%
  arrange(desc(change))

say("\n=== M_78 allocation, normal mucosa to carcinoma, by lineage ===")
print(as.data.frame(out %>% transmute(lineage, n,
  normal = round(normal, 2), carcinoma = round(carcinoma, 2),
  change = round(change, 2), fold = round(fold, 2),
  rising = sprintf("%d/%d", rising, n),
  He = round(He, 2), Kim = round(Kim, 2), p = signif(p, 2),
  gene = best_gene, rho = round(best_rho, 3), criterion)), row.names = FALSE)

say("\n  lineages rising: %d of %d | rising in ALL patients: %d",
    sum(out$change > 0), nrow(out), sum(out$rising == out$n))
say("  of those meeting the criterion: %d rise, %d rise in all patients",
    sum(out$change > 0 & out$criterion == "meets"),
    sum(out$rising == out$n & out$criterion == "meets"))

my <- out %>% filter(lineage == "Myeloid")
if (nrow(my)) {
  say("\n  myeloid change %.2f points (%.2f-fold), %d of %d patients",
      my$change, my$fold, my$rising, my$n)
  ot <- out %>% filter(lineage != "Myeloid", criterion == "meets")
  if (nrow(ot)) {
    say("  other criterion-meeting lineages:")
    for (i in seq_len(nrow(ot)))
      say("    %-11s %+.2f points (%.2f-fold), %d of %d  -- %.0f%% of the myeloid change",
          ot$lineage[i], ot$change[i], ot$fold[i], ot$rising[i], ot$n[i],
          100 * ot$change[i] / my$change)
  }
}

say("\n=== is the myeloid change larger, tested within patient? ===")
say("  paired difference in the change itself, myeloid minus other lineage")
delta <- function(lin) {
  cw <- cwa %>% filter(annot_4th == lin)
  mc <- intersect(mcols_of(cw), u); if (!M %in% mc) return(NULL)
  S <- as.matrix(cw[, mc]); S <- S / rowSums(S)
  cw %>% mutate(v = 100 * S[, M]) %>%
    group_by(patient, stage) %>%
    summarise(v = weighted.mean(v, n_cells), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = stage, values_from = v) %>%
    filter(!is.na(normal), !is.na(carcinoma)) %>%
    transmute(patient, d = carcinoma - normal)
}
dm <- delta("Myeloid")
for (lin in setdiff(out$lineage, "Myeloid")) {
  do <- delta(lin); if (is.null(do)) next
  j <- inner_join(dm, do, by = "patient", suffix = c("_mye","_oth"))
  if (nrow(j) < 4) { say("  vs %-11s only %d shared patients", lin, nrow(j)); next }
  say("  vs %-11s n=%2d | myeloid change larger in %d of %d | median gap %+.2f | p = %.3f",
      lin, nrow(j), sum(j$d_mye > j$d_oth), nrow(j),
      median(j$d_mye - j$d_oth),
      suppressWarnings(wilcox.test(j$d_mye, j$d_oth, paired = TRUE))$p.value)
}
say("\n  Run under both arms before concluding anything about specificity.")
