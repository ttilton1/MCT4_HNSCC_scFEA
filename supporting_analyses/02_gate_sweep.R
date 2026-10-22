#!/usr/bin/env Rscript
# ==============================================================================
# 02_gate_sweep.R  --  choose the n_cells gate from the data, not by guessing.
#
#   Rscript 02_gate_sweep.R --arm bulk_before
#
# The gate trades unit QUALITY against POWER. In bulk_before a low-n_cells unit
# is not merely noisy -- it is a different object (n_cells=1 is a single cell's
# imputed profile handed to scFEA, i.e. a per-cell flux wearing a pseudobulk
# label). So the floor must be high enough that every unit is a real mean, but
# low enough that the progression stages keep >= min_patients.
#
# For each candidate gate this reports, per major type:
#   - units kept / dropped, and the fraction of units that were degenerate (<3 cells)
#   - patients kept
#   - PROGRESSION power: patients per stage (Normal/Pre/CA) within He and Kim
#   - HPV power: patients per HPV arm
#   - finaltypes still testable (>= min_patients)
#
# Pick the largest gate that keeps progression and HPV testable, then set
# GATES$min_cells_* in config_all_sites.R.
# Output -> results/all_sites/<ARM>/02_gate_sweep/
# ==============================================================================
source("config_all_sites.R")

OUT <- file.path(PATHS$results, "02_gate_sweep"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
GATES_TRY <- c(0, 3, 5, 8, 10, 15, 20, 30)

cw <- readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))
cw <- cw %>% mutate(stage = unname(STAGE_MAP[tissue]))

say("\n=== gate sweep | ARM = %s | %s finaltype units ===", ARM, format(nrow(cw), big.mark=","))

# ---- 0. how degenerate is the ungated table? ---------------------------------
degen <- cw %>% group_by(annot_4th) %>%
  summarise(units = n(),
            n1  = sum(n_cells == 1), n_lt3 = sum(n_cells < 3), n_lt8 = sum(n_cells < 8),
            pct_lt3 = round(100*mean(n_cells < 3), 1),
            pct_lt8 = round(100*mean(n_cells < 8), 1),
            med_cells = as.integer(median(n_cells)), .groups = "drop") %>%
  arrange(desc(pct_lt3))
say("\n-- degenerate units (n_cells=1 is a SINGLE CELL, not a pseudobulk) --")
print(as.data.frame(degen), row.names = FALSE)
write_csv(degen, file.path(OUT, "degenerate_units.csv"))

# ---- 1. sweep -----------------------------------------------------------------
sweep_one <- function(g) {
  d <- cw %>% filter(n_cells >= g)

  base <- d %>% group_by(annot_4th) %>%
    summarise(units = n(), patients = n_distinct(patient),
              finaltypes = n_distinct(finaltype), .groups = "drop")

  # progression power: patients per stage, per cohort (He / Kim), per major type
  prog <- d %>% filter(study %in% PROG_COHORTS, !is.na(stage)) %>%
    distinct(study, annot_4th, stage, patient) %>%
    count(study, annot_4th, stage, name = "n_pt") %>%
    pivot_wider(names_from = stage, values_from = n_pt, values_fill = 0)
  for (s in STAGE_LEVELS) if (!s %in% names(prog)) prog[[s]] <- 0L
  prog <- prog %>%
    mutate(testable = Normal >= GATES$min_patients &
                      Pre    >= GATES$min_patients &
                      CA     >= GATES$min_patients)

  # HPV power: patients per arm
  hpv <- d %>% filter(hpv %in% c("Positive","Negative")) %>%
    distinct(annot_4th, hpv, patient) %>% count(annot_4th, hpv, name = "n_pt") %>%
    pivot_wider(names_from = hpv, values_from = n_pt, values_fill = 0)
  for (h in c("Positive","Negative")) if (!h %in% names(hpv)) hpv[[h]] <- 0L
  hpv <- hpv %>% mutate(hpv_testable = Positive >= GATES$min_per_group_hpv &
                                       Negative >= GATES$min_per_group_hpv)

  list(
    summary = base %>% mutate(gate = g),
    prog    = prog %>% mutate(gate = g),
    hpv     = hpv  %>% mutate(gate = g),
    top     = tibble(gate = g,
                     units = nrow(d),
                     patients = n_distinct(d$patient),
                     finaltypes = n_distinct(d$finaltype),
                     prog_testable = sum(prog$testable),
                     hpv_testable  = sum(hpv$hpv_testable)))
}

res <- lapply(GATES_TRY, sweep_one)
top  <- bind_rows(lapply(res, `[[`, "top"))
prog <- bind_rows(lapply(res, `[[`, "prog"))
hpvt <- bind_rows(lapply(res, `[[`, "hpv"))
summ <- bind_rows(lapply(res, `[[`, "summary"))

say("\n-- HEADLINE: what each gate costs --")
say("   prog_testable = # (cohort x major type) with >=%d patients in ALL 3 stages",
    GATES$min_patients)
say("   hpv_testable  = # major types with >=%d patients in BOTH HPV arms",
    GATES$min_per_group_hpv)
print(as.data.frame(top), row.names = FALSE)

say("\n-- progression power per cohort x major type (patients per stage) --")
print(as.data.frame(prog %>% filter(gate %in% c(0, 8, 15, 20)) %>%
                    arrange(gate, study, annot_4th)), row.names = FALSE)

say("\n-- HPV power per major type --")
print(as.data.frame(hpvt %>% filter(gate %in% c(0, 8, 15, 20)) %>%
                    arrange(gate, annot_4th)), row.names = FALSE)

write_csv(top,  file.path(OUT, "gate_sweep_headline.csv"))
write_csv(prog, file.path(OUT, "gate_sweep_progression_power.csv"))
write_csv(hpvt, file.path(OUT, "gate_sweep_hpv_power.csv"))
write_csv(summ, file.path(OUT, "gate_sweep_by_majortype.csv"))

say("\n=== HOW TO READ ===")
say("Find the LARGEST gate where prog_testable and hpv_testable have not yet")
say("collapsed. That is the floor: every surviving unit is a real pseudobulk,")
say("and you still have the patients to test. Set GATES$min_cells_* in config.")
say("If a gate of even 3-5 already guts progression, the honest read is that the")
say("thin units ARE your power -- report the gate and the sensitivity, don't hide it.")
say("\nOutputs -> %s", OUT)
