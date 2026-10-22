#!/usr/bin/env Rscript
# ==============================================================================
# 05c_progression_meta.R  --  two-cohort meta-analysis + He fine-stage extension
#
#   Rscript 05c_progression_meta.R --arm bulk_before
#
# WHY. 05 (absolute flux) was ~95% "up" -- a GLOBAL scale shift, not biology.
# 05b showed total flux rises with stage in every cell type (CA/Normal 1.08-1.55;
# that is ONE finding), and that on scale-free flux SHARE the direction balance
# is restored and only a handful of real allocation shifts survive. This script
# is the proper inferential version of 05b:
#
#   * ALWAYS on flux SHARE (module flux / that unit's total). Absolute flux is
#     not used for hypothesis testing -- it only carries the global shift.
#
#   * BOTH cohorts are DISCOVERY. 05/05b used He (n=8) as primary and Kim
#     (n=14-22) merely as a direction check, which wastes the larger cohort and
#     rests every FDR on 8 patients. Here each cohort gets its own JT trend test,
#     then STOUFFER meta-analysis weighted by sqrt(n_patients). Concordant
#     direction is REQUIRED -- a hit driven by one cohort alone is not a hit.
#
#   * HE FINE-STAGE EXTENSION. He is the only cohort with staged carcinoma
#     (CA-E early -> CA-A advanced), so He alone runs the 4-level axis
#     NL -> Pre -> CA-E -> CA-A. That is what He is FOR; Kim cannot do it.
#
# Output -> results/all_sites/<ARM>/05c_progression_meta/<level>/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(purrr); library(stringr) })

LEVELS <- c("finaltype", "subtype")
MR     <- load_modref()
usable <- readRDS(file.path(PATHS$results, "03_qc", "usable_modules.rds"))

# Jonckheere-Terpstra for an arbitrary ordered factor (3 or 4 levels)
jt_test <- function(x, g, lv) {
  ok <- is.finite(x) & !is.na(g); x <- x[ok]
  g <- droplevels(factor(g[ok], levels = lv))
  ns <- table(g); k <- length(ns); N <- sum(ns)
  if (k < 3 || any(ns < 2)) return(c(z = NA_real_, p = NA_real_))
  U <- 0; L <- levels(g)
  for (i in 1:(k-1)) for (j in (i+1):k) {
    a <- x[g == L[i]]; b <- x[g == L[j]]
    U <- U + sum(outer(a, b, "<")) + 0.5 * sum(outer(a, b, "=="))
  }
  EU <- (N^2 - sum(ns^2)) / 4
  VU <- (N^2*(2*N+3) - sum(ns^2*(2*ns+3))) / 72
  if (VU <= 0) return(c(z = NA_real_, p = NA_real_))
  z <- (U - EU)/sqrt(VU); c(z = z, p = 2*pnorm(-abs(z)))
}

# ---- Confirm what each cohort has via printing them out. ---------------------
cw0 <- readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds")) %>% filter(n_cells >= GATES$min_cells)
say("\n=== tissue x study (patients) -- which cohort can stage carcinoma? ===")
xt <- cw0 %>% distinct(study, tissue, patient) %>% count(study, tissue) %>%
  pivot_wider(names_from = study, values_from = n, values_fill = 0)
print(as.data.frame(xt), row.names = FALSE)
write_csv(xt, file.path(PATHS$results, "05c_tissue_by_study.csv"))

HE_FINE <- c("NL", "Pre", "CA-E", "CA-A")     # He's 4-level axis (verify above)
FINE_LV <- c("Normal", "Pre", "CA-E", "CA-A")

run_level <- function(level) {
  OUT <- file.path(PATHS$results, "05c_progression_meta", level)
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  say("\n================ LEVEL: %s ================", level)

  cw <- readRDS(file.path(PATHS$cw_dir,
         sprintf("cw_%s.rds", if (level=="finaltype") "finaltype" else "subtype"))) %>%
    filter(n_cells >= GATES$min_cells)
  unit_col <- if (level=="finaltype") "finaltype" else "annot_2nd"
  mcols <- intersect(mcols_of(cw), usable)

  base <- cw %>% filter(!tissue %in% STAGE_EXCLUDED) %>%
    rename(celltype = annot_4th, unit = !!unit_col) %>% filter(!is.na(unit))

  # ---------- A. 3-stage axis, BOTH cohorts, on flux SHARE ----------
  d3 <- base %>% mutate(stage = unname(STAGE_MAP[tissue])) %>%
    filter(!is.na(stage), study %in% PROG_COHORTS)

  pt <- d3 %>% group_by(study, celltype, unit, patient, stage) %>%
    summarise(across(all_of(mcols), ~ weighted.mean(.x, n_cells)),
              n_cells = sum(n_cells), .groups = "drop")
  X <- as.matrix(pt[, mcols]); tot <- rowSums(X)
  ps <- bind_cols(pt %>% select(study, celltype, unit, patient, stage),
                  as_tibble(X / tot))                      # SHARE

  long <- ps %>% pivot_longer(all_of(mcols), names_to = "module", values_to = "share") %>%
    mutate(stage = factor(stage, levels = STAGE_LEVELS, ordered = TRUE))

  per_cohort <- long %>% group_by(study, celltype, unit, module) %>%
    filter(n_distinct(stage) == 3, min(table(droplevels(stage))) >= 2) %>%
    summarise(n_pt = n_distinct(patient),
              tt = list(jt_test(share, stage, STAGE_LEVELS)),
              med_N = median(share[stage=="Normal"]),
              med_CA = median(share[stage=="CA"]), .groups = "drop") %>%
    mutate(z = map_dbl(tt,"z"), p = map_dbl(tt,"p"), tt = NULL,
           delta = med_CA - med_N) %>% filter(!is.na(z))

  wide <- per_cohort %>%
    select(study, celltype, unit, module, z, p, delta, n_pt) %>%
    pivot_wider(names_from = study, values_from = c(z, p, delta, n_pt))

  need <- c("z_He","z_Kim","n_pt_He","n_pt_Kim","delta_He","delta_Kim")
  for (nm in need) if (!nm %in% names(wide)) wide[[nm]] <- NA_real_

  # ---- STOUFFER meta (weights = sqrt(n_patients)); concordance REQUIRED ----
  meta <- wide %>%
    filter(!is.na(z_He), !is.na(z_Kim)) %>%
    mutate(
      w_He  = sqrt(n_pt_He), w_Kim = sqrt(n_pt_Kim),
      z_meta = (w_He*z_He + w_Kim*z_Kim) / sqrt(w_He^2 + w_Kim^2),
      p_meta = 2*pnorm(-abs(z_meta)),
      fdr    = p.adjust(p_meta, "BH"),
      concordant  = sign(z_He) == sign(z_Kim),
      delta_meta  = (w_He*delta_He + w_Kim*delta_Kim) / (w_He + w_Kim),
      direction   = ifelse(z_meta > 0, "up", "down")) %>%
    left_join(MR$sm, by="module") %>% left_join(MR$rxn, by="module")

  floor_abs <- derive_floor(meta$delta_meta, 0.99)
  meta <- meta %>%
    mutate(pass_fdr = fdr < FDR_THRESH,
           pass_effect = abs(delta_meta) >= floor_abs,
           survivor = pass_fdr & pass_effect & concordant) %>%
    arrange(desc(survivor), fdr, desc(abs(delta_meta)))

  say("\n-- meta-analysable tests (present in BOTH cohorts): %s of %s --",
      format(nrow(meta), big.mark=","), format(nrow(wide), big.mark=","))
  say("effect floor |delta share| (99th pctl): %.5f", floor_abs)
  say("direction balance: up=%d down=%d  (a global shift cannot survive in SHARE)",
      sum(meta$direction=="up"), sum(meta$direction=="down"))
  say("concordant across cohorts: %.0f%%", 100*mean(meta$concordant))

  funnel <- tibble(
    step = factor(c("meta tests","FDR<0.10","+ effect floor","+ concordant (survivors)"),
      levels = c("meta tests","FDR<0.10","+ effect floor","+ concordant (survivors)")),
    n = c(nrow(meta), sum(meta$pass_fdr),
          sum(meta$pass_fdr & meta$pass_effect), sum(meta$survivor)))
  say("funnel: %s", paste(sprintf("%s=%d", funnel$step, funnel$n), collapse=" -> "))
  write_csv(funnel, file.path(OUT, "funnel.csv"))
  write_csv(meta, file.path(OUT, "meta_all.csv"))
  write_csv(meta %>% filter(survivor), file.path(OUT, "meta_survivors.csv"))

  if (sum(meta$survivor)) {
    say("\n-- META SURVIVORS (%s): both cohorts, same direction, scale-free --", level)
    print(as.data.frame(meta %>% filter(survivor) %>%
      arrange(desc(abs(delta_meta))) %>% slice_head(n=25) %>%
      transmute(celltype, unit, module, reaction, direction,
                d = signif(delta_meta,3), fdr = signif(fdr,3),
                z_He = round(z_He,1), z_Kim = round(z_Kim,1),
                n_He = n_pt_He, n_Kim = n_pt_Kim)), row.names=FALSE)
  } else say("\n-- no meta survivors at %s --", level)

  # forest-ish: per-cohort z for the survivors
  sv <- meta %>% filter(survivor) %>% arrange(desc(abs(delta_meta))) %>% slice_head(n=20)
  if (nrow(sv)) {
    fp <- sv %>% mutate(lab = str_trunc(paste0(celltype,"/",unit," | ",reaction), 44)) %>%
      select(lab, delta_meta, He = z_He, Kim = z_Kim) %>%
      pivot_longer(c(He, Kim), names_to = "cohort", values_to = "z")
    p1 <- ggplot(fp, aes(z, reorder(lab, abs(delta_meta)), colour = cohort)) +
      geom_vline(xintercept = 0, colour = "grey40") +
      geom_point(size = 2.6, alpha = .85) +
      scale_colour_manual(values = c(He = "#d6604d", Kim = "#4393c3")) +
      labs(title = sprintf("Progression meta-survivors (%s, %s)", ARM, level),
           subtitle = "JT trend z on flux SHARE, per cohort. Both sides of zero = concordant.",
           x = "JT z (Normal -> Pre -> CA)", y = NULL) +
      theme_bw(base_size = 10)
    ggsave(file.path(OUT, "fig_forest.png"), p1, width = 10, height = 7, dpi = 150, bg = "white")
  }

  # ---------- B. HE FINE-STAGE EXTENSION: descriptive, NOT discovery ----------
  # He has only ~3 patients per stage (NL=3, Pre=3, CA-E=3, CA-A=5), so JT across
  # 4 levels has almost no resolution: p-values become discrete and tied, and any
  # near-zero module that happens to order monotonically looks "significant"
  # (an earlier version of this block returned 2053 "hits" with effects ~1e-5).
  # So He is NOT used for discovery here. Instead: take the modules that ALREADY
  # survived the two-cohort meta, and describe what they do across
  # NL -> Pre -> CA-E -> CA-A. No new tests, no new multiple-testing burden --
  # just "does the established trend continue into ADVANCED carcinoma?"
  say("\n-- He fine-stage EXTENSION of the meta survivors (descriptive) --")
  sv_keys <- meta %>% filter(survivor) %>% select(celltype, unit, module)
  dHe <- base %>% filter(study == "He", tissue %in% HE_FINE) %>%
    mutate(fstage = factor(recode(tissue, NL = "Normal", Pre = "Pre"),
                           levels = FINE_LV, ordered = TRUE))

  if (!nrow(sv_keys) || n_distinct(dHe$fstage) < 4) {
    say("  no survivors or He lacks all 4 levels -- skipping")
  } else {
    ptH <- dHe %>% group_by(celltype, unit, patient, fstage) %>%
      summarise(across(all_of(mcols), ~ weighted.mean(.x, n_cells)),
                n_cells = sum(n_cells), .groups = "drop")
    XH <- as.matrix(ptH[, mcols])
    psH <- bind_cols(ptH %>% select(celltype, unit, patient, fstage),
                     as_tibble(XH / rowSums(XH)))

    ext <- psH %>% pivot_longer(all_of(mcols), names_to="module", values_to="share") %>%
      inner_join(sv_keys, by = c("celltype","unit","module")) %>%
      group_by(celltype, unit, module, fstage) %>%
      summarise(med_share = median(share), n_pt = n_distinct(patient), .groups="drop")

    tab <- ext %>% select(-n_pt) %>%
      pivot_wider(names_from = fstage, values_from = med_share) %>%
      left_join(MR$rxn, by="module") %>%
      mutate(across(any_of(FINE_LV), ~ signif(.x, 3)),
             step_CAE_to_CAA = signif(`CA-A` - `CA-E`, 3),
             continues = sign(`CA-A` - `CA-E`) == sign(`CA-E` - Normal))
    say("  survivors traceable across all 4 He stages: %d of %d",
        sum(complete.cases(tab %>% select(any_of(FINE_LV)))), nrow(sv_keys))
    say("  'continues' = the CA-E -> CA-A step keeps the same direction as NL -> CA-E")
    print(as.data.frame(tab %>% filter(complete.cases(select(., any_of(FINE_LV)))) %>%
      arrange(desc(abs(step_CAE_to_CAA))) %>% slice_head(n = 20) %>%
      select(celltype, unit, reaction, any_of(FINE_LV), step_CAE_to_CAA, continues)),
      row.names = FALSE)
    write_csv(tab, file.path(OUT, "he_fine_stage_extension.csv"))

    pE <- ext %>% inner_join(MR$rxn, by = "module") %>%
      mutate(lab = str_trunc(paste0(celltype,"/",unit," | ",reaction), 40))
    p2 <- ggplot(pE, aes(fstage, med_share, group = lab, colour = lab)) +
      geom_line(alpha = .8) + geom_point(size = 1.8) +
      guides(colour = "none") + facet_wrap(~lab, scales = "free_y", ncol = 4) +
      labs(title = sprintf("Meta-survivors across He's 4 stages (%s, %s)", ARM, level),
           subtitle = "descriptive; He n~3/stage. Does the trend continue into ADVANCED carcinoma?",
           x = NULL, y = "median flux share") +
      theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
    ggsave(file.path(OUT, "fig_he_fine_stage.png"), p2,
           width = 13, height = 2 + 2.1*ceiling(n_distinct(pE$lab)/4), dpi = 150,
           bg = "white", limitsize = FALSE)
  }

  say("-> %s", OUT)
  meta
}

out <- lapply(LEVELS, run_level); names(out) <- LEVELS
say("\n=== SUMMARY (%s) ===", ARM)
for (lv in LEVELS) if (!is.null(out[[lv]]))
  say("  %-10s: %d meta survivors (both cohorts, concordant, scale-free)",
      lv, sum(out[[lv]]$survivor))
