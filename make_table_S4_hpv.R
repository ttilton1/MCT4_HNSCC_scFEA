#!/usr/bin/env Rscript
# ==============================================================================
# make_table_S4_hpv.R  --  Supplementary Table S4, HPV-stratified comparisons
#
#   Rscript make_table_S4_hpv.R --arm bulk_before
#
# The 39 myeloid lactate-handling comparisons (13 subtypes x 3 modules), plus
# the two comparisons that met the retention criteria anywhere in the HPV screen,
# recorded below a divider.
#
# TWO THINGS THE CODE SETTLES, both stated in the caption:
#   1. The test is the Wilcoxon rank-sum test within each cohort (06_hpv.R
#      line 51), requiring at least five patients per HPV arm in that cohort.
#   2. Single-cohort comparisons were TESTED AND REPORTED, not excluded, but
#      could never be retained: concordance is defined as
#      (all(z>0)|all(z<0)) & n()>=2, so a comparison present in one cohort is
#      assigned concordant = FALSE and fails the retention rule whatever its
#      p-value. 34 of the 39 are in that position.
#
# Output -> results/all_sites/<ARM>/supplementary/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(dplyr) })

OUT <- file.path(PATHS$results, "supplementary")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
HPV <- file.path(PATHS$results, "06_hpv", "finaltype")
LACT <- c(M_4 = "M_4 glycolysis (3PD to pyruvate)",
          M_6 = "M_6 lactate production (pyruvate to lactate)",
          M_78 = "M_78 lactate export")

all_ <- read_csv(file.path(HPV, "hpv_all.csv"), show_col_types = FALSE)
say("hpv_all.csv: %s rows | columns: %s", format(nrow(all_), big.mark=","),
    paste(names(all_), collapse = ", "))

# the celltype column is logical in some outputs (a column collision upstream);
# fall back to the subtype prefix when that happens
if (is.logical(all_$celltype)) {
  say("[!] celltype is logical in this file; deriving it from the subtype prefix")
  all_$celltype <- ifelse(grepl("^(MQ_|DC_|Mono_|Neutrophil)", all_$unit), "Myeloid",
                   ifelse(grepl("^(CD4_|CD8_|Treg)", all_$unit), "T", NA_character_))
}

# ---- patient counts per cohort x subtype x arm, from the crosswalk ----------
# matches 06_hpv.R, which now restricts to carcinoma before averaging each
# patient (line 69), so every patient value here is carcinoma-only
cw <- as.data.frame(readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))) %>%
  filter(n_cells >= GATES$min_cells,
         tissue %in% names(STAGE_MAP)[STAGE_MAP == "CA"],
         hpv %in% c("Negative","Positive"))
pt <- cw %>% distinct(study, finaltype, patient, hpv) %>%
  count(study, finaltype, hpv) %>%
  tidyr::pivot_wider(names_from = hpv, values_from = n, values_fill = 0) %>%
  rename(unit = finaltype, n_pos = Positive, n_neg = Negative) %>%
  filter(n_pos >= 5, n_neg >= 5) %>%          # the arms that could be tested
  mutate(txt = sprintf("%s (%d+/%d-)", study, n_pos, n_neg))
cohort_n <- pt %>% group_by(unit) %>%
  summarise(arms = paste(txt, collapse = "; "), .groups = "drop")

# ---- the 39 -----------------------------------------------------------------
tbl <- all_ %>%
  filter(celltype == "Myeloid", module %in% names(LACT)) %>%
  left_join(cohort_n, by = "unit") %>%
  mutate(`cell subtype`            = unit,
         `module`                  = LACT[module],
         `cohorts testable`        = cohorts,
         `patients per arm`        = ifelse(is.na(arms), "-", arms),
         `effect (pct points)`     = round(100 * delta_meta, 4),
         `higher in`               = higher_in,
         `relative shift`          = round(rel_effect, 3),
         `same direction`          = ifelse(n_cohorts < 2, "NA (one cohort)",
                                     ifelse(concordant, "yes", "no")),
         `combined P`              = signif(p_meta, 3),
         `FDR`                     = signif(fdr, 3),
         `>= 10% shift`            = ifelse(pass_rel, "yes", "no"),
         `retained`                = ifelse(survivor, "yes", "no")) %>%
  arrange(fdr) %>%
  select(`cell subtype`, module, `cohorts testable`, `patients per arm`,
         `effect (pct points)`, `higher in`, `relative shift`, `same direction`,
         `combined P`, FDR, `>= 10% shift`, retained)

say("\n=== the %d myeloid lactate-handling comparisons ===", nrow(tbl))
say("  minimum FDR: %.3f | retained: %d", min(tbl$FDR), sum(tbl$retained == "yes"))
say("  testable in one cohort only: %d of %d",
    sum(tbl$`same direction` == "NA (one cohort)"), nrow(tbl))
print(as.data.frame(head(tbl, 8)), row.names = FALSE)

# ---- the comparisons that met the criteria anywhere in the screen -----------
surv <- all_ %>% filter(survivor) %>%
  left_join(cohort_n, by = "unit") %>%
  mutate(`cell subtype`        = unit,
         module                = sprintf("%s %s", module, reaction),
         `cohorts testable`    = cohorts,
         `patients per arm`    = ifelse(is.na(arms), "-", arms),
         `effect (pct points)` = round(100 * delta_meta, 6),
         `higher in`           = higher_in,
         `relative shift`      = round(rel_effect, 3),
         `same direction`      = ifelse(concordant, "yes", "no"),
         `combined P`          = signif(p_meta, 3),
         FDR                   = signif(fdr, 3),
         `>= 10% shift`        = ifelse(pass_rel, "yes", "no"),
         retained              = "yes") %>%
  arrange(fdr) %>% select(names(tbl))
say("\n=== comparisons meeting the criteria anywhere in the screen: %d ===", nrow(surv))
if (nrow(surv)) {
  say("  NOTE: a relative shift far above 1 means the module's own flux is near zero,")
  say("  so the ratio is unstable. Check the absolute allocation before interpreting:")
  for (i in seq_len(nrow(surv)))
    say("    %s %s: effect %.2e pct points, relative shift %.1f",
        surv$`cell subtype`[i], surv$module[i],
        as.numeric(surv$`effect (pct points)`[i]), as.numeric(surv$`relative shift`[i]))
}
print(as.data.frame(surv), row.names = FALSE)

divider <- tbl[1, ]; divider[] <- ""
divider$`cell subtype` <- sprintf(
  "--- met the retention criteria elsewhere in the screen (%d of %s comparisons) ---",
  nrow(surv), format(nrow(all_), big.mark = ","))
# the surviving row(s) clear the relative-effect floor only because the floor is
# a ratio and the module's own flux is near zero. That is the artifact 06_hpv.R
# warns about at the end of its run, so the note travels with the table rather
# than living only in the console output.
note <- tbl[1, ]; note[] <- ""
note$`cell subtype` <- paste(
  "Note: the row(s) above clear the 10% relative floor because that floor is a",
  "ratio and the module's flux is near zero. M_19 has a median allocation of",
  "-5.6e-08% in MQ_C2_APOE, six orders of magnitude below M_78 (0.032%) and",
  "negative, which is optimization residue rather than a measurable quantity.",
  "This is recorded for completeness and is not interpreted as a finding.")

out <- bind_rows(mutate(tbl, across(everything(), as.character)),
                 mutate(divider, across(everything(), as.character)),
                 mutate(surv, across(everything(), as.character)),
                 mutate(note, across(everything(), as.character)))

CAP <- sprintf(paste0(
  "Wilcoxon rank-sum test within each cohort, at least five patients per HPV arm; ",
  "cohort z-scores combined by Stouffer's method weighted by the square root of ",
  "the smaller arm; Benjamini-Hochberg across all %s cell subtype by module ",
  "comparisons. Retention required FDR < %.2f, a shift of at least %d%% of the ",
  "module's flux, and the same direction in both cohorts. Comparisons testable in ",
  "only one cohort were reported but could not be retained, since concordance ",
  "requires two cohorts. Analysis is restricted to carcinoma units; each patient ",
  "contributes one value. The minimum FDR among the %d lactate-handling ",
  "comparisons was %.2f."),
  format(nrow(all_), big.mark = ","), FDR_THRESH, round(100 * REL_FLOOR),
  nrow(tbl), min(tbl$FDR))

write_csv(out, file.path(OUT, "TableS4_hpv_myeloid.csv"))
writeLines(CAP, file.path(OUT, "TableS4_hpv_caption.txt"))

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "S4 HPV myeloid")
  openxlsx::writeData(wb, 1, out, startRow = 1)
  openxlsx::writeData(wb, 1, CAP, startRow = nrow(out) + 3)
  openxlsx::addStyle(wb, 1, openxlsx::createStyle(textDecoration = "bold"),
                     rows = 1, cols = seq_len(ncol(out)))
  openxlsx::setColWidths(wb, 1, cols = seq_len(ncol(out)), widths = "auto")
  openxlsx::saveWorkbook(wb, file.path(OUT, "TableS4_hpv_myeloid.xlsx"), overwrite = TRUE)
  say("\n-> TableS4_hpv_myeloid.xlsx and .csv")
} else {
  say("\n[!] openxlsx not installed; CSV written only.")
  say("    install.packages(\"openxlsx\")  then re-run for the xlsx")
}
say("\nCaption:\n%s", CAP)
say("\nFor Methods:")
say("  (1) test: Wilcoxon rank-sum within cohort (06_hpv.R line 51)")
say("  (2) single-cohort comparisons were reported, not excluded, but could not")
say("      be retained: concordance requires n_cohorts >= 2 by definition.")
say("  (3) the gate is five PATIENTS per arm (06_hpv.R groups by patient before")
say("      testing), restricted to carcinoma units, so the comparison matches the")
say("      carcinoma framing used in every other cross-cohort analysis.")
