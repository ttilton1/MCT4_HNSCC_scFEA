#!/usr/bin/env Rscript
# ==============================================================================
# 03b_ln_offset_check.R  --  is the LN dominance (QC Q2) GLOBAL or MODULE-SPECIFIC?
#
#   Rscript 03b_ln_offset_check.R --arm bulk_before
#
# QC Q2 found LN has the highest mean flux for 73-86% of modules in the immune
# types (Myeloid .78, NK .73, B .82, Plasma .86, T .61). Two very different
# explanations, with opposite consequences:
#
#   GLOBAL OFFSET  -- LN units simply have higher TOTAL metabolic-gene expression,
#     so scFEA's geneExprScale lifts EVERY module together. Near-uniform, not
#     pathway-specific. CONSEQUENCE: double-centering (specialization) removes it
#     entirely by construction, so specialization is SAFE. Only ABSOLUTE
#     cross-site flux comparisons are affected.
#
#   MODULE-SPECIFIC -- LN differentially inflates particular pathways. CONSEQUENCE:
#     a real batch/biology confound that survives centering, and site would have to
#     be modeled explicitly.
#
# Three tests:
#   T1. Per-unit TOTAL flux (sum across modules) by site. A global offset shows up
#       here as a clean LN > Oral shift.
#   T2. After removing each unit's total (row-centering, i.e. flux share), does LN
#       still dominate? If dominance collapses to ~chance, it was global.
#   T3. Correlation of per-module LN-vs-Oral log-ratio with the module's mean flux.
#       A global scaling produces a flat, near-constant log-ratio across modules.
#
# Output -> results/all_sites/<ARM>/03b_ln_offset/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2) })

OUT <- file.path(PATHS$results, "03b_ln_offset"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
usable <- readRDS(file.path(PATHS$results, "03_qc", "usable_modules.rds"))

cw <- readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds")) %>% filter(n_cells >= GATES$min_cells)
mcols <- intersect(mcols_of(cw), usable)
say("\n=== LN offset check | %s | %s units | %d modules ===",
    ARM, format(nrow(cw), big.mark=","), length(mcols))

X <- as.matrix(cw[, mcols])

# ---- T1: total flux per unit, by site ----------------------------------------
say("\n-- T1: per-unit TOTAL flux by site (global offset shows up here) --")
t1 <- cw %>% mutate(total_flux = rowSums(X)) %>%
  group_by(annot_4th, site_group) %>%
  summarise(n = n(), med_total = median(total_flux), .groups = "drop") %>%
  pivot_wider(id_cols = annot_4th, names_from = site_group,
              values_from = med_total, names_prefix = "tot_") %>%
  mutate(LN_over_Oral = round(tot_LN / tot_Oral, 2))
print(as.data.frame(t1), row.names = FALSE)
write_csv(t1, file.path(OUT, "T1_total_flux_by_site.csv"))
say("  READ: LN_over_Oral >> 1 across the board = LN units carry a GLOBAL flux lift.")

# ---- T2: does LN dominance survive row-centering (flux SHARE)? ----------------
say("\n-- T2: site dominance BEFORE vs AFTER removing each unit's total --")
share <- X / rowSums(X)                       # each unit's flux SHARE per module
cw_sh <- bind_cols(cw %>% select(annot_4th, site_group), as_tibble(share))

dom <- function(d) {
  d %>% pivot_longer(all_of(mcols), names_to = "module", values_to = "v") %>%
    group_by(annot_4th, module, site_group) %>%
    summarise(m = mean(v), .groups = "drop") %>%
    group_by(annot_4th, module) %>% slice_max(m, n = 1, with_ties = FALSE) %>% ungroup() %>%
    count(annot_4th, site_group) %>% group_by(annot_4th) %>%
    mutate(frac = n / sum(n)) %>% ungroup()
}
before <- dom(bind_cols(cw %>% select(annot_4th, site_group), as_tibble(X))) %>%
  filter(site_group == "LN") %>% transmute(annot_4th, LN_frac_raw = round(frac, 2))
after  <- dom(cw_sh) %>% filter(site_group == "LN") %>%
  transmute(annot_4th, LN_frac_share = round(frac, 2))
t2 <- full_join(before, after, by = "annot_4th") %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0)),
         collapsed = LN_frac_share < 0.5) %>%
  arrange(desc(LN_frac_raw))
print(as.data.frame(t2), row.names = FALSE)
write_csv(t2, file.path(OUT, "T2_dominance_before_after.csv"))
say("  READ: if LN_frac_share collapses toward ~1/3 (chance), the dominance was a")
say("        GLOBAL offset -> double-centering removes it -> SPECIALIZATION IS SAFE.")
say("        If it stays high, LN inflates specific pathways -> must model site.")

# ---- T3: is the LN/Oral log-ratio flat across modules? ------------------------
say("\n-- T3: per-module LN/Oral log-ratio (flat = global scaling) --")
t3 <- cw %>% select(annot_4th, site_group, all_of(mcols)) %>%
  pivot_longer(all_of(mcols), names_to = "module", values_to = "flux") %>%
  filter(site_group %in% c("LN","Oral")) %>%
  group_by(annot_4th, module, site_group) %>%
  summarise(m = mean(flux) + 1e-12, .groups = "drop") %>%
  pivot_wider(names_from = site_group, values_from = m) %>%
  mutate(log2_LN_Oral = log2(LN / Oral))
t3_sum <- t3 %>% group_by(annot_4th) %>%
  summarise(med_log2 = round(median(log2_LN_Oral, na.rm = TRUE), 2),
            iqr_log2 = round(IQR(log2_LN_Oral, na.rm = TRUE), 2),
            frac_pos = round(mean(log2_LN_Oral > 0, na.rm = TRUE), 2), .groups = "drop") %>%
  arrange(desc(med_log2))
print(as.data.frame(t3_sum), row.names = FALSE)
write_csv(t3, file.path(OUT, "T3_module_log_ratio.csv"))
say("  READ: med_log2 >0 with a SMALL iqr = uniform lift (global). A large iqr means")
say("        modules are shifted by different amounts = module-specific.")

p <- ggplot(t3, aes(log2_LN_Oral, annot_4th)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_jitter(height = .2, size = .6, alpha = .4) +
  labs(title = sprintf("LN vs Oral, per module (%s)", ARM),
       subtitle = "tight cluster off zero = GLOBAL offset (centering removes it). Wide spread = module-specific.",
       x = "log2(mean LN flux / mean Oral flux)", y = NULL) +
  theme_bw(base_size = 11)
ggsave(file.path(OUT, "fig_ln_vs_oral.png"), p, width = 8, height = 5, dpi = 150, bg = "white")

say("\nVERDICT depends on T2: LN_frac_share collapsing = specialization unaffected.")
say("-> %s", OUT)
