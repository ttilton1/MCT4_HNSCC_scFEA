#!/usr/bin/env Rscript
# ==============================================================================
# 03_qc_flux.R  --  Tier 1. Do NOT trust any downstream result until this passes.
#
#   Rscript 03_qc_flux.R --arm bulk_before
#
# Four checks, each of which can invalidate a finding:
#
#  Q1 DEPTH / SIZE CONFOUND. Does flux track n_cells (bulk_before) or median_umi?
#     If a module's flux correlates with how many cells the unit was built from,
#     then any contrast where group sizes differ (CA units are bigger than Normal
#     units!) is partly a size artifact. This is the new-arm analogue of the old
#     raw-counts depth confound. Reports the fraction of modules with |rho| > 0.3
#     per major type -- the number to watch.
#
#  Q2 BATCH / SITE STRUCTURE. Which site is "highest" for each module, per major
#     type (your 08 diagnostic). LN dominance was the raw-counts batch signature;
#     with normalized+MAGIC input it should be diffuse. A type where one site wins
#     >50% of modules is suspect.
#
#  Q3 MODULE USABILITY. Modules that are ~constant or ~zero everywhere cannot
#     support a contrast; they inflate the multiple-testing burden and produce
#     rank-test artifacts. Emits a keep-list used by downstream scripts.
#
#  Q4 STUDY SEPARATION. Per major type, how much of flux variance is explained by
#     study alone (one-way, per module)? High = cross-cohort claims are unsafe
#     (they already are -- this quantifies it).
#
# Output -> results/all_sites/<ARM>/03_qc/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(purrr) })

OUT <- file.path(PATHS$results, "03_qc"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
MR  <- load_modref()

cw <- readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds")) %>%
  filter(n_cells >= GATES$min_cells)
mcols <- mcols_of(cw)
say("\n=== QC | ARM = %s | gate n_cells >= %d | %s units ===",
    ARM, GATES$min_cells, format(nrow(cw), big.mark = ","))

long <- cw %>%
  select(all_of(c(UNIT_KEYS, "site","site_group","hpv","study","annot_4th","n_cells")), all_of(mcols)) %>%
  pivot_longer(all_of(mcols), names_to = "module", values_to = "flux")

# ---- Q1: size / depth confound ------------------------------------------------
say("\n-- Q1: does flux track unit size (n_cells)? --")
q1 <- long %>% group_by(annot_4th, module) %>%
  filter(n_distinct(n_cells) > 2) %>%
  summarise(rho = suppressWarnings(cor(flux, n_cells, method = "spearman", use = "complete.obs")),
            .groups = "drop")
q1_sum <- q1 %>% group_by(annot_4th) %>%
  summarise(n_modules = n(),
            med_rho = round(median(rho, na.rm = TRUE), 3),
            frac_abs_gt_0.3 = round(mean(abs(rho) > 0.3, na.rm = TRUE), 3),
            frac_abs_gt_0.5 = round(mean(abs(rho) > 0.5, na.rm = TRUE), 3),
            .groups = "drop") %>% arrange(desc(frac_abs_gt_0.3))
print(as.data.frame(q1_sum), row.names = FALSE)
write_csv(q1, file.path(OUT, "Q1_flux_vs_ncells_rho.csv"))
write_csv(q1_sum, file.path(OUT, "Q1_summary.csv"))
say("  READ: frac_abs_gt_0.3 is the alarm. If a major type is high, its contrasts")
say("        between groups of different unit SIZE (CA units are bigger than Normal)")
say("        are partly a size artifact -> report n_cells as a covariate / check.")

p1 <- ggplot(q1, aes(rho, annot_4th)) +
  geom_vline(xintercept = c(-0.3, 0, 0.3), linetype = c(2,1,2), colour = c("grey60","grey30","grey60")) +
  geom_jitter(height = 0.2, size = 0.6, alpha = 0.4) +
  labs(title = sprintf("Q1. Flux vs unit size (%s)", ARM),
       subtitle = "each point = one module. Mass outside +/-0.3 = size confound in that major type.",
       x = "Spearman rho (flux ~ n_cells)", y = NULL) +
  theme_bw(base_size = 11)
ggsave(file.path(OUT, "Q1_flux_vs_ncells.png"), p1, width = 8, height = 5, dpi = 150, bg = "white")

# ---- Q2: which site dominates each module ------------------------------------
say("\n-- Q2: site dominance per major type (batch signature) --")
q2 <- long %>% group_by(annot_4th, module, site_group) %>%
  summarise(m = mean(flux, na.rm = TRUE), .groups = "drop") %>%
  group_by(annot_4th, module) %>% slice_max(m, n = 1, with_ties = FALSE) %>% ungroup() %>%
  count(annot_4th, site_group, name = "n_modules_highest") %>%
  group_by(annot_4th) %>% mutate(frac = round(n_modules_highest / sum(n_modules_highest), 2)) %>%
  ungroup() %>% arrange(annot_4th, desc(frac))
print(as.data.frame(q2 %>% pivot_wider(id_cols = annot_4th, names_from = site_group,
                                       values_from = frac, values_fill = 0)), row.names = FALSE)
write_csv(q2, file.path(OUT, "Q2_site_dominance.csv"))
say("  READ: in the RAW-counts run LN won 50-68%% of modules (batch signature).")
say("        Diffuse (~0.2-0.4 each) = normalization worked. One site > 0.5 = suspect.")

# ---- Q3: module usability ------------------------------------------------------
say("\n-- Q3: module usability --")
q3 <- long %>% group_by(module) %>%
  summarise(mean_flux = mean(flux, na.rm = TRUE),
            sd_flux   = sd(flux, na.rm = TRUE),
            frac_zero = mean(abs(flux) < 1e-9, na.rm = TRUE),
            cv        = ifelse(mean_flux > 0, sd_flux / mean_flux, NA_real_), .groups = "drop") %>%
  mutate(usable = frac_zero < 0.5 & sd_flux > 1e-8) %>%
  arrange(desc(frac_zero))
n_drop <- sum(!q3$usable)
say("  %d / %d modules unusable (>=50%% zeros or no variance) -> excluded downstream",
    n_drop, nrow(q3))
if (n_drop) print(as.data.frame(q3 %>% filter(!usable) %>%
  left_join(MR$sm, by = "module") %>% select(module, supermodule, mean_flux, frac_zero)), row.names = FALSE)
write_csv(q3 %>% left_join(MR$sm, by = "module"), file.path(OUT, "Q3_module_usability.csv"))
saveRDS(q3$module[q3$usable], file.path(OUT, "usable_modules.rds"))

# ---- Q4: study separation ------------------------------------------------------
say("\n-- Q4: variance explained by study (cross-cohort safety) --")
q4 <- long %>% group_by(annot_4th, module) %>%
  filter(n_distinct(study) >= 2) %>%
  summarise(r2_study = {
    fit <- try(summary(lm(flux ~ study))$r.squared, silent = TRUE)
    if (inherits(fit, "try-error")) NA_real_ else fit
  }, .groups = "drop")
q4_sum <- q4 %>% group_by(annot_4th) %>%
  summarise(med_r2_study = round(median(r2_study, na.rm = TRUE), 3),
            frac_gt_0.25 = round(mean(r2_study > 0.25, na.rm = TRUE), 3), .groups = "drop") %>%
  arrange(desc(med_r2_study))
print(as.data.frame(q4_sum), row.names = FALSE)
write_csv(q4, file.path(OUT, "Q4_study_r2.csv"))
say("  READ: this is WHY progression is run WITHIN cohort. High r2_study just")
say("        confirms the design; it does not invalidate within-cohort trends.")

say("\nQC done -> %s", OUT)
say("Gate downstream on: usable_modules.rds (Q3) + n_cells >= %d.", GATES$min_cells)
