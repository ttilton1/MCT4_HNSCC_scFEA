#!/usr/bin/env Rscript
# ==============================================================================
# check_export_production_ratio.R
#
#   Rscript check_export_production_ratio.R --arm bulk_before
#
# Myeloid leads every lineage on lactate export (M_78) but not on lactate
# production (M_6), where fibroblasts and NK cells allocate more. Two things
# need checking before that becomes a claim.
#
# 1. THE SCALE OBJECTION. Allocation is a share of each cell's own total, so a
#    lineage with higher total flux could exceed another in absolute terms while
#    allocating less. Each cross-lineage comparison is vulnerable to this.
#
#    The ratio M_78 / M_6 is not. Both are shares of the same denominator, so the
#    total cancels: (M_78/T) / (M_6/T) = M_78/M_6. Whatever a lineage's flux
#    scale, its export-to-production ratio is unaffected. If myeloid's ratio
#    exceeds fibroblasts', that lineage tips toward export relative to
#    production, and no assumption about scale can reverse it.
#
# 2. DOES THE SPLIT CHANGE WITH STAGE? If it is only a carcinoma property, it is
#    a description of the tumour. If it holds in normal mucosa too, it is a
#    property of the lineages.
#
# What this still cannot say: which lineage releases more lactate into tissue.
# That needs absolute rates and cell abundance, neither of which this design has.
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages(library(dplyr))

u  <- readRDS(file.path(PATHS$results, "03_qc", "usable_modules.rds"))
cw <- as.data.frame(readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))) %>%
  filter(n_cells >= GATES$min_cells)
mc <- intersect(mcols_of(cw), u)
S  <- as.matrix(cw[, mc]); S <- S / rowSums(S)
cw$M_6 <- 100 * S[, "M_6"]; cw$M_78 <- 100 * S[, "M_78"]
cw$ratio <- cw$M_78 / cw$M_6           # scale-free: the total cancels

CA <- c("CA","CA-E","CA-A","CA-c","CA-p"); NL <- c("NL","Normal")

paired <- function(d, col, label) {
  w <- d %>% group_by(patient, annot_4th) %>%
    summarise(v = mean(.data[[col]]), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = annot_4th, values_from = v)
  others <- setdiff(names(w), c("patient","Myeloid"))
  out <- lapply(others, function(ct) {
    x <- w %>% filter(!is.na(Myeloid), !is.na(.data[[ct]]))
    if (nrow(x) < 8) return(NULL)
    dd <- x$Myeloid - x[[ct]]
    tibble(other = ct, n = nrow(x), pct_higher = round(100*mean(dd > 0)),
           med_diff = round(median(dd), 3),
           p = signif(suppressWarnings(wilcox.test(dd)$p.value), 2))
  }) %>% bind_rows() %>% arrange(desc(med_diff))
  say("\n  -- %s --", label); print(as.data.frame(out), row.names = FALSE)
  out
}

say("=== 1. export-to-production ratio, carcinoma, paired within patient ===")
say("    (scale-free: total flux cancels, so this is immune to the objection")
say("     that lineages are not on a common magnitude scale)")
r_ca <- paired(cw %>% filter(tissue %in% CA), "ratio",
               "M_78 / M_6, myeloid minus other lineage")

say("\n=== 2. the same ratio in NORMAL mucosa ===")
say("    if the split holds here too it is a lineage property, not a tumour one")
r_nl <- paired(cw %>% filter(tissue %in% NL), "ratio",
               "M_78 / M_6, normal mucosa")

say("\n=== 3. ratio levels by lineage and tissue ===")
lv <- cw %>% filter(tissue %in% c(NL, CA)) %>%
  mutate(tc = ifelse(tissue %in% NL, "normal", "carcinoma")) %>%
  group_by(annot_4th, tc) %>%
  summarise(n = n(), median_ratio = round(median(ratio), 2), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = tc, values_from = c(n, median_ratio)) %>%
  arrange(desc(median_ratio_carcinoma))
print(as.data.frame(lv), row.names = FALSE)
say("\n  ratio > 1 means a lineage allocates more to exporting lactate than to")
say("  making it; < 1 the reverse. This is a within-cell-type quantity.")

say("\n=== 4. does the myeloid-fibroblast contrast hold in both tissues? ===")
for (tis in list(list(NL,"normal"), list(CA,"carcinoma"))) {
  d <- cw %>% filter(tissue %in% tis[[1]], annot_4th %in% c("Myeloid","Fibroblast"))
  w <- d %>% group_by(patient, annot_4th) %>% summarise(v = mean(ratio), .groups="drop") %>%
    tidyr::pivot_wider(names_from = annot_4th, values_from = v) %>%
    filter(!is.na(Myeloid), !is.na(Fibroblast))
  if (nrow(w) >= 5) {
    say("  %s: n = %d paired patients, myeloid ratio %.2f vs fibroblast %.2f, p = %.3g",
        tis[[2]], nrow(w), median(w$Myeloid), median(w$Fibroblast),
        suppressWarnings(wilcox.test(w$Myeloid, w$Fibroblast, paired = TRUE))$p.value)
  } else say("  %s: only %d paired patients, not tested", tis[[2]], nrow(w))
}
