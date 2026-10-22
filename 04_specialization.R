#!/usr/bin/env Rscript
# ==============================================================================
# 04_specialization.R  --  Tier 2. WHO over-allocates to which pathway.
#
#   Rscript 04_specialization.R --arm bulk_before
#
# Specialization = double-centered residual of flux within each (patient x major
# type) block: rows = the patient's finaltype (or subtype) units, cols = modules.
# Removing the row mean (this unit's overall metabolic activity) and the column
# mean (this module's typical level in that patient) leaves "does THIS unit
# over-allocate to THIS module, relative to its peers in the same patient".
# Patient-internal => cohort/batch offsets cancel by construction.
#
# WITHIN-MODEL: every block is one major type = one trained scFEA network. Flux is
# never compared across major types. (Cross-SITE is now fine -- the all_sites arms
# train one model spanning Oral/LN/MR -- so unlike the old pipeline we do NOT
# restrict to Oral.)
#
# Per (major type x unit x module): one-sample Wilcoxon of the residual vs 0
# across patients, + FDR, + DATA-DERIVED effect floor, + LOSO across studies,
# + cross-cohort direction consistency. Fully general: no hardcoded modules,
# no hardcoded focal cell types. Ranked hit table out.
#
# Output -> results/all_sites/<ARM>/04_specialization/<level>/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(purrr); library(stringr) })

LEVELS <- c("finaltype", "subtype")
MR     <- load_modref()
usable <- readRDS(file.path(PATHS$results, "03_qc", "usable_modules.rds"))
say("usable modules: %d / 168 (from QC Q3)", length(usable))

# double-center one (patient x major type) block; NA if too few units
double_center <- function(w, min_units) {
  if (nrow(w) < min_units) return(NULL)
  u <- w$unit; X <- as.matrix(w[, setdiff(names(w), "unit"), drop = FALSE])
  R <- sweep(X, 1, rowMeans(X)); R <- sweep(R, 2, colMeans(X)) + mean(X)
  as_tibble(R) %>% mutate(unit = u) %>%
    pivot_longer(-unit, names_to = "module", values_to = "spec")
}

run_level <- function(level) {
  OUT <- file.path(PATHS$results, "04_specialization", level)
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  say("\n================ LEVEL: %s ================", level)

  cw <- readRDS(file.path(PATHS$cw_dir,
                          sprintf("cw_%s.rds", if (level == "finaltype") "finaltype" else "subtype"))) %>%
    filter(n_cells >= GATES$min_cells)
  unit_col <- if (level == "finaltype") "finaltype" else "annot_2nd"
  mcols <- intersect(mcols_of(cw), usable)

  d <- cw %>% rename(celltype = annot_4th, unit = !!unit_col) %>% filter(!is.na(unit))
  say("units: %s | patients: %d | modules: %d",
      format(nrow(d), big.mark=","), n_distinct(d$patient), length(mcols))

  # ---- residuals, per (patient x celltype) block ----
  resid <- d %>%
    select(patient, celltype, unit, study, hpv, site_group, all_of(mcols)) %>%
    group_by(patient, celltype) %>%
    group_modify(function(g, k) {
      w <- g %>% select(unit, all_of(mcols))
      r <- double_center(w, GATES$min_units_dc)
      if (is.null(r)) return(tibble())
      r %>% left_join(distinct(g, unit, study, hpv, site_group), by = "unit")
    }) %>% ungroup()

  n_blocks <- resid %>% distinct(patient, celltype) %>% nrow()
  say("double-centered blocks (>= %d units): %d | residual rows: %s",
      GATES$min_units_dc, n_blocks, format(nrow(resid), big.mark=","))
  if (!nrow(resid)) { say("  no testable blocks at this level -- skipping"); return(NULL) }

  # ---- core test: is this unit's residual consistently != 0 across patients? ----
  core <- resid %>%
    group_by(celltype, unit, module) %>%
    summarise(n_pt = n_distinct(patient),
              med_spec = median(spec, na.rm = TRUE),
              wilcox_p = if (n_distinct(patient) >= GATES$min_patients)
                suppressWarnings(wilcox.test(spec, mu = 0)$p.value) else NA_real_,
              .groups = "drop") %>%
    filter(!is.na(wilcox_p)) %>%
    mutate(fdr = p.adjust(wilcox_p, "BH"))

  # ---- DATA-DERIVED effect floor (arm-specific scale; see config) ----
  floor_abs <- derive_floor(core$med_spec)
  say("data-derived effect floor (%.0fth pctl of |median residual|): %.5f",
      100*FLOOR_PCTL, floor_abs)

  # ---- LOSO: does the sign survive dropping each study? ----
  studies <- sort(unique(resid$study))
  loso <- map_dfr(studies, function(s) {
    resid %>% filter(study != s) %>%
      group_by(celltype, unit, module) %>%
      summarise(med = median(spec, na.rm = TRUE), n_pt = n_distinct(patient), .groups = "drop") %>%
      filter(n_pt >= GATES$min_patients) %>% mutate(dropped = s)
  })
  loso_sum <- loso %>%
    inner_join(core %>% select(celltype, unit, module, med_full = med_spec),
               by = c("celltype","unit","module")) %>%
    group_by(celltype, unit, module) %>%
    summarise(n_holdouts = n(),
              loso_stable = mean(sign(med) == sign(med_full), na.rm = TRUE),
              .groups = "drop")

  # ---- cross-cohort consistency: same direction inside every cohort? ----
  cons <- resid %>%
    group_by(celltype, unit, module, study) %>%
    summarise(cmed = median(spec, na.rm = TRUE), n_pt = n_distinct(patient), .groups = "drop") %>%
    filter(n_pt >= 2) %>%
    group_by(celltype, unit, module) %>%
    summarise(n_cohorts = n(),
              consistent = (all(cmed > 0) | all(cmed < 0)) & n() >= 2,
              .groups = "drop")

  res <- core %>%
    left_join(loso_sum, by = c("celltype","unit","module")) %>%
    left_join(cons,     by = c("celltype","unit","module")) %>%
    left_join(MR$sm,  by = "module") %>% left_join(MR$rxn, by = "module") %>%
    mutate(
      n_cohorts   = coalesce(n_cohorts, 0L),
      consistent  = coalesce(consistent, FALSE),
      loso_stable = coalesce(loso_stable, NA_real_),
      direction   = ifelse(med_spec > 0, "over", "under"),
      pass_fdr    = fdr < FDR_THRESH,
      pass_effect = abs(med_spec) >= floor_abs,
      pass_loso   = !is.na(loso_stable) & loso_stable == 1,
      pass_multi  = n_cohorts >= 2,
      notable     = pass_fdr & pass_effect & pass_multi,
      survivor    = notable & consistent & pass_loso) %>%
    arrange(desc(survivor), desc(notable), fdr, desc(abs(med_spec)))

  write_csv(res, file.path(OUT, "specialization_all.csv"))
  write_csv(res %>% filter(survivor), file.path(OUT, "specialization_survivors.csv"))

  funnel <- tibble(
    step = factor(c("all tests","FDR<0.10","+ effect floor","+ >=2 cohorts",
                    "+ consistent","+ LOSO stable (survivors)"),
                  levels = c("all tests","FDR<0.10","+ effect floor","+ >=2 cohorts",
                             "+ consistent","+ LOSO stable (survivors)")),
    n = c(nrow(res), sum(res$pass_fdr),
          sum(res$pass_fdr & res$pass_effect),
          sum(res$notable),
          sum(res$notable & res$consistent),
          sum(res$survivor)))
  write_csv(funnel, file.path(OUT, "funnel.csv"))
  say("\nfunnel: %s", paste(sprintf("%s=%d", funnel$step, funnel$n), collapse = " -> "))

  if (sum(res$survivor)) {
    say("\n-- TOP SURVIVORS (%s) --", level)
    print(as.data.frame(res %>% filter(survivor) %>%
      arrange(desc(abs(med_spec))) %>% slice_head(n = 20) %>%
      transmute(celltype, unit, module, reaction, supermodule, direction,
                med_spec = round(med_spec, 4), fdr = signif(fdr, 3),
                n_pt, n_cohorts)), row.names = FALSE)
  } else say("\n-- no survivors at %s level --", level)

  # ---- figures ----
  pF <- ggplot(funnel, aes(step, n)) +
    geom_col(fill = "#4393c3", width = .7) +
    geom_text(aes(label = n), vjust = -0.4, fontface = "bold", size = 4) +
    labs(title = sprintf("Specialization funnel (%s, %s)", ARM, level),
         subtitle = sprintf("effect floor = %.4f (data-derived, %.0fth pctl)", floor_abs, 100*FLOOR_PCTL),
         x = NULL, y = "# tests") +
    theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(file.path(OUT, "fig_funnel.png"), pF, width = 8, height = 5, dpi = 150, bg = "white")

  top <- res %>% filter(notable) %>% arrange(desc(abs(med_spec))) %>% slice_head(n = 25)
  if (nrow(top)) {
    top <- top %>% mutate(lab = str_trunc(paste0(unit, " | ", reaction), 42))
    pT <- ggplot(top, aes(reorder(lab, abs(med_spec)), med_spec, fill = direction)) +
      geom_col(width = .7) + coord_flip() + geom_hline(yintercept = 0, colour = "grey50") +
      geom_text(data = filter(top, survivor), aes(label = "*"), hjust = -0.3, size = 5) +
      scale_fill_manual(values = c(over = "#b2182b", under = "#2166ac"), name = NULL) +
      labs(title = sprintf("Top specialization (%s, %s)", ARM, level),
           subtitle = "median double-centered residual | * = survivor (FDR + effect + >=2 cohorts + consistent + LOSO)",
           x = NULL, y = "specialization residual") +
      theme_bw(base_size = 10)
    ggsave(file.path(OUT, "fig_top_hits.png"), pT, width = 10, height = 7, dpi = 150, bg = "white")
  }

  say("-> %s", OUT)
  res
}

out <- lapply(LEVELS, run_level)
names(out) <- LEVELS
say("\n=== SUMMARY (%s) ===", ARM)
for (lv in LEVELS) if (!is.null(out[[lv]]))
  say("  %-10s: %d notable | %d survivors", lv,
      sum(out[[lv]]$notable), sum(out[[lv]]$survivor))
say("\nNEXT: 05_progression.R, 06_hpv.R. Then 07_marker_grounding.R will take the")
say("TOP SURVIVORS from here and ground them in enzyme expression (generalized 23/23b/23c).")
