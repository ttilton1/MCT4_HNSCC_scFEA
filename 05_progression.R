#!/usr/bin/env Rscript
# ==============================================================================
# 05_progression.R  --  Tier 2. Does flux trend along NL -> Pre -> CA?
#
#   Rscript 05_progression.R --arm bulk_before
#
# AXIS: ORAL MUCOSA ONLY (NL -> Pre/LP -> CA*). LN-normal and NL-contact are
# EXCLUDED (see config): LN-normal is a normal lymph node that inherits
# site_group="Oral" from site="Distal_normal", and LN units carry a large global
# flux offset -- folding it into "Normal" would manufacture a trend. NL-contact is
# field-cancerized, not distal normal.
#
# WITHIN COHORT ALWAYS. He is primary (spans all 3 stages); Kim is REPLICATION by
# direction agreement, not by pooling. Studies are never pooled for this axis --
# QC Q4 shows study explains 6-45% of flux variance depending on major type.
#
# TEST: Jonckheere-Terpstra (ordered-alternative trend), computed across PATIENTS
# per (major type x unit x module). One value per patient per stage first, so no
# pseudoreplication. Reports the paired-patient overlap as a diagnostic -- JT
# assumes independent groups, and He's design may put a patient in >1 stage.
#
# Also: He-only CA-E -> CA-A extension (staged carcinoma), where available.
#
# Output -> results/all_sites/<ARM>/05_progression/<level>/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(purrr); library(stringr) })

LEVELS <- c("finaltype", "subtype")
MR     <- load_modref()
usable <- readRDS(file.path(PATHS$results, "03_qc", "usable_modules.rds"))

# ---- Jonckheere-Terpstra, normal approximation with tie correction ------------
jt_test <- function(x, g) {                       # x numeric, g ORDERED factor
  ok <- is.finite(x) & !is.na(g); x <- x[ok]; g <- droplevels(factor(g[ok], levels = STAGE_LEVELS))
  ns <- table(g); k <- length(ns); N <- sum(ns)
  if (k < 3 || any(ns < 2)) return(c(z = NA_real_, p = NA_real_))
  U <- 0
  lv <- levels(g)
  for (i in 1:(k-1)) for (j in (i+1):k) {
    a <- x[g == lv[i]]; b <- x[g == lv[j]]
    U <- U + sum(outer(a, b, "<")) + 0.5 * sum(outer(a, b, "=="))
  }
  EU <- (N^2 - sum(ns^2)) / 4
  VU <- (N^2 * (2*N + 3) - sum(ns^2 * (2*ns + 3))) / 72
  if (VU <= 0) return(c(z = NA_real_, p = NA_real_))
  z <- (U - EU) / sqrt(VU)
  c(z = z, p = 2 * pnorm(-abs(z)))                # two-sided; sign of z = direction
}

run_level <- function(level) {
  OUT <- file.path(PATHS$results, "05_progression", level)
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  say("\n================ LEVEL: %s ================", level)

  cw <- readRDS(file.path(PATHS$cw_dir,
                          sprintf("cw_%s.rds", if (level == "finaltype") "finaltype" else "subtype"))) %>%
    filter(n_cells >= GATES$min_cells)
  unit_col <- if (level == "finaltype") "finaltype" else "annot_2nd"
  mcols <- intersect(mcols_of(cw), usable)

  d <- cw %>%
    filter(!tissue %in% STAGE_EXCLUDED) %>%
    mutate(stage = unname(STAGE_MAP[tissue])) %>%
    filter(!is.na(stage), study %in% PROG_COHORTS) %>%
    rename(celltype = annot_4th, unit = !!unit_col) %>% filter(!is.na(unit))

  say("oral-axis units: %s | %s", format(nrow(d), big.mark=","),
      paste(sprintf("%s=%d", names(table(d$study)), table(d$study)), collapse = " "))
  say("stages: %s", paste(sprintf("%s=%d", names(table(d$stage)), table(d$stage)), collapse = " "))

  # one value per PATIENT per stage (n_cells-weighted), then JT across patients
  pt <- d %>% group_by(study, celltype, unit, patient, stage) %>%
    summarise(across(all_of(mcols), ~ weighted.mean(.x, n_cells)),
              n_cells = sum(n_cells), .groups = "drop") %>%
    mutate(stage = factor(stage, levels = STAGE_LEVELS, ordered = TRUE))

  # pseudoreplication / pairing diagnostic
  pair <- pt %>% distinct(study, patient, stage) %>% count(study, patient) %>%
    group_by(study) %>% summarise(patients = n(), in_multi_stage = sum(n > 1),
                                  pct = round(100*mean(n > 1)), .groups = "drop")
  say("\n-- pairing diagnostic (JT assumes INDEPENDENT groups) --")
  print(as.data.frame(pair), row.names = FALSE)
  if (any(pair$in_multi_stage > 0))
    say("  CAUTION: %d patients appear in >1 stage -> JT p-values are anticonservative.",
        sum(pair$in_multi_stage))
  write_csv(pair, file.path(OUT, "pairing_diagnostic.csv"))

  long <- pt %>% pivot_longer(all_of(mcols), names_to = "module", values_to = "flux")

  # ---- JT within each cohort ----
  jt <- long %>% group_by(study, celltype, unit, module) %>%
    filter(n_distinct(stage) == 3) %>%
    summarise(n_pt = n_distinct(patient),
              n_min_stage = min(table(droplevels(stage))),
              tt = list(jt_test(flux, stage)),
              med_Normal = median(flux[stage == "Normal"]),
              med_Pre    = median(flux[stage == "Pre"]),
              med_CA     = median(flux[stage == "CA"]),
              .groups = "drop") %>%
    mutate(z = map_dbl(tt, "z"), p = map_dbl(tt, "p"), tt = NULL,
           delta = med_CA - med_Normal,
           direction = ifelse(z > 0, "up", "down")) %>%
    filter(!is.na(p), n_min_stage >= 2)

  # He primary + FDR within He
  he <- jt %>% filter(study == "He") %>% mutate(fdr = p.adjust(p, "BH"))
  km <- jt %>% filter(study == "Kim") %>%
    select(celltype, unit, module, z_kim = z, p_kim = p, delta_kim = delta,
           dir_kim = direction, n_pt_kim = n_pt)

  floor_abs <- derive_floor(he$delta)
  say("\ndata-derived effect floor on |CA - Normal|: %.5f", floor_abs)

  res <- he %>% left_join(km, by = c("celltype","unit","module")) %>%
    left_join(MR$sm, by = "module") %>% left_join(MR$rxn, by = "module") %>%
    mutate(
      replicated = !is.na(dir_kim) & dir_kim == direction,      # DIRECTION agreement
      rel_effect = abs(delta) / pmax(abs(med_Normal), 1e-9),
      pass_fdr    = fdr < FDR_THRESH,
      pass_effect = abs(delta) >= floor_abs,
      pass_rel    = rel_effect >= REL_FLOOR,
      notable     = pass_fdr & pass_effect,
      survivor    = notable & replicated) %>%
    arrange(desc(survivor), fdr, desc(abs(delta)))

  write_csv(res, file.path(OUT, "progression_all.csv"))
  write_csv(res %>% filter(survivor), file.path(OUT, "progression_survivors.csv"))

  funnel <- tibble(
    step = factor(c("He tests","FDR<0.10","+ effect floor","+ Kim direction (survivors)"),
                  levels = c("He tests","FDR<0.10","+ effect floor","+ Kim direction (survivors)")),
    n = c(nrow(res), sum(res$pass_fdr), sum(res$notable), sum(res$survivor)))
  write_csv(funnel, file.path(OUT, "funnel.csv"))
  say("funnel: %s", paste(sprintf("%s=%d", funnel$step, funnel$n), collapse = " -> "))
  say("Kim-testable fraction: %.0f%% of He tests had a Kim counterpart",
      100*mean(!is.na(res$dir_kim)))

  if (sum(res$survivor)) {
    say("\n-- TOP PROGRESSION SURVIVORS (%s) --", level)
    print(as.data.frame(res %>% filter(survivor) %>% arrange(desc(abs(delta))) %>%
      slice_head(n = 20) %>%
      transmute(celltype, unit, module, reaction, supermodule, direction,
                delta = signif(delta, 3), fdr = signif(fdr, 3),
                n_pt, kim_dir = dir_kim)), row.names = FALSE)
  } else say("\n-- no survivors at %s level --", level)

  pF <- ggplot(funnel, aes(step, n)) + geom_col(fill = "#5aae61", width = .7) +
    geom_text(aes(label = n), vjust = -0.4, fontface = "bold") +
    labs(title = sprintf("Progression funnel (%s, %s)", ARM, level),
         subtitle = "He primary (JT trend, oral axis NL->Pre->CA); Kim = direction replication",
         x = NULL, y = "# tests") +
    theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
  ggsave(file.path(OUT, "fig_funnel.png"), pF, width = 8, height = 5, dpi = 150, bg = "white")

  top <- res %>% filter(notable) %>% arrange(desc(abs(delta))) %>% slice_head(n = 25)
  if (nrow(top)) {
    top <- top %>% mutate(lab = str_trunc(paste0(celltype, "/", unit, " | ", reaction), 46))
    pT <- ggplot(top, aes(reorder(lab, abs(delta)), delta, fill = direction)) +
      geom_col(width = .7) + coord_flip() + geom_hline(yintercept = 0, colour = "grey50") +
      geom_text(data = filter(top, survivor), aes(label = "*"), hjust = -0.3, size = 5) +
      scale_fill_manual(values = c(up = "#b2182b", down = "#2166ac"), name = NULL) +
      labs(title = sprintf("Progression NL->CA (%s, %s)", ARM, level),
           subtitle = "* = replicated in Kim by direction",
           x = NULL, y = "median CA - median Normal (He)") +
      theme_bw(base_size = 10)
    ggsave(file.path(OUT, "fig_top_hits.png"), pT, width = 10, height = 7, dpi = 150, bg = "white")
  }

  say("-> %s", OUT)
  res
}

out <- lapply(LEVELS, run_level); names(out) <- LEVELS
say("\n=== SUMMARY (%s) ===", ARM)
for (lv in LEVELS) if (!is.null(out[[lv]]))
  say("  %-10s: %d notable | %d replicated survivors", lv,
      sum(out[[lv]]$notable), sum(out[[lv]]$survivor))
