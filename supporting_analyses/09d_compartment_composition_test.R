#!/usr/bin/env Rscript
# ==============================================================================
# 09d_compartment_composition_test.R
#   WHY does myeloid export allocation rise from normal mucosa to carcinoma?
#
#   Rscript 09d_compartment_composition_test.R --arm bulk_before
#   Rscript 09d_compartment_composition_test.R --arm bulk_after
#
# WHAT WE KNOW. In the 12 patients with both normal and carcinoma tissue, myeloid
# export allocation roughly doubled (0.0170 -> 0.0331, p = 5e-4). We also know it
# is not attributable to any single state: cDC2 rose in step with the compartment
# rather than beyond it (ratio to patient myeloid mean unchanged, p = 0.63).
#
# WHAT WE DO NOT KNOW is why the compartment rises. Two explanations:
#
#   (A) COMPOSITION. High-export states appear during progression. SPP1+ TAM goes
#       from 0-2% of myeloid cells in normal tissue to 7-9% in carcinoma, and has
#       the highest export allocation of any myeloid state. Adding such cells
#       raises the compartment mean with no individual cell changing.
#
#   (B) COORDINATED REPROGRAMMING. The states present in normal tissue each shift
#       upward, plausibly driven by a shared environmental cue.
#
# These are distinguishable. If the rise is compositional, removing the emergent
# states should abolish it. If resident states are also reprogramming, the rise
# should persist among states present in BOTH normal and carcinoma.
#
# TESTS
#   A. recompute the paired normal-vs-carcinoma rise EXCLUDING SPP1+ TAM
#   B. recompute it using ONLY states present in both stages ("resident" states)
#   C. decompose the change into a composition term and a within-state term
#      (a standard shift-share decomposition):
#         delta_total = SUM_s [ w_s(CA) - w_s(N) ] * x_s(N)      <- composition
#                     + SUM_s   w_s(N) * [ x_s(CA) - x_s(N) ]    <- within-state
#                     + interaction
#      where w_s is state s's share of myeloid cells and x_s its export allocation.
#
# DECISION: if the within-state term dominates, report coordinated reprogramming.
# If composition dominates, report emergence of high-export states. If mixed,
# report both with their relative contributions.
#
# Output -> results/all_sites/<ARM>/09d_composition/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(purrr) })

OUT <- file.path(PATHS$results, "09d_composition")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
usable <- readRDS(file.path(PATHS$results, "03_qc", "usable_modules.rds"))

NORMAL_T <- c("NL","Normal"); CANCER_T <- c("CA","CA-E","CA-A","CA-c","CA-p")
EMERGENT <- "MQ_C4_SPP1"

cw <- as.data.frame(readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))) %>%
  filter(annot_4th == "Myeloid", study %in% PROG_COHORTS, n_cells >= GATES$min_cells) %>%
  mutate(stage = case_when(tissue %in% NORMAL_T ~ "normal",
                           tissue %in% CANCER_T ~ "carcinoma", TRUE ~ NA_character_)) %>%
  filter(!is.na(stage))
mc <- intersect(mcols_of(cw), usable)
S <- as.matrix(cw[, mc]); S <- S / rowSums(S)
cw$export <- S[, "M_78"]

paired_test <- function(d, label) {
  if (!nrow(d)) { say("  %-42s no units", label); return(NULL) }
  pw <- d %>% group_by(study, patient, stage) %>%
    summarise(e = weighted.mean(export, n_cells), .groups = "drop") %>%
    pivot_wider(names_from = stage, values_from = e)
  # a state set may be absent from one stage entirely (emergent states have no
  # normal-tissue units by definition), in which case that column never appears
  if (!all(c("normal", "carcinoma") %in% names(pw))) {
    missing <- setdiff(c("normal", "carcinoma"), names(pw))
    say("  %-42s not present in %s tissue -- no paired test possible", label, missing)
    return(NULL)
  }
  pw <- pw %>% filter(!is.na(normal), !is.na(carcinoma))
  if (nrow(pw) < 5) { say("  %-42s too few paired patients (%d)", label, nrow(pw)); return(NULL) }
  w <- suppressWarnings(wilcox.test(pw$carcinoma, pw$normal, paired = TRUE))
  say("  %-42s n=%2d  %.4f -> %.4f  (%.2fx)  p = %.3g",
      label, nrow(pw), median(pw$normal), median(pw$carcinoma),
      median(pw$carcinoma)/median(pw$normal), w$p.value)
  tibble(set = label, n = nrow(pw), normal = median(pw$normal),
         carcinoma = median(pw$carcinoma), p = w$p.value)
}

# ---- which states are present in both stages? -------------------------------
present <- cw %>% group_by(finaltype, stage) %>% summarise(n = n(), .groups="drop") %>%
  pivot_wider(names_from = stage, values_from = n, values_fill = 0)
resident <- present %>% filter(normal > 0, carcinoma > 0) %>% pull(finaltype)
emergent_states <- present %>% filter(normal == 0, carcinoma > 0) %>% pull(finaltype)
say("=== myeloid states present in NORMAL tissue (resident): %d ===", length(resident))
say("  %s", paste(resident, collapse = ", "))
say("=== states appearing only in CARCINOMA (emergent): %d ===", length(emergent_states))
say("  %s", paste(emergent_states, collapse = ", "))

# ---- A and B: paired rise under different state sets ------------------------
say("\n=== paired normal vs carcinoma, myeloid export allocation ===")
res <- bind_rows(
  paired_test(cw, "all myeloid states"),
  paired_test(cw %>% filter(finaltype != EMERGENT), "excluding SPP1+ TAM"),
  paired_test(cw %>% filter(!finaltype %in% emergent_states), "resident states only"),
  paired_test(cw %>% filter(finaltype %in% emergent_states), "emergent states only"))
if (!is.null(res)) write_csv(res, file.path(OUT, "paired_by_state_set.csv"))
say("\n  READ: if the rise persists when emergent states are removed, resident")
say("  myeloid cells are also reprogramming. If it collapses, the compartment")
say("  change is driven by the appearance of high-export states.")

# ---- C: shift-share decomposition -------------------------------------------
say("\n=== decomposition of the compartment change ===")
comp <- cw %>% group_by(study, stage, finaltype) %>%
  summarise(cells = sum(n_cells), x = weighted.mean(export, n_cells), .groups="drop") %>%
  group_by(study, stage) %>% mutate(w = cells / sum(cells)) %>% ungroup()

dec <- map_dfr(PROG_COHORTS, function(s) {
  n <- comp %>% filter(study == s, stage == "normal") %>% select(finaltype, wN = w, xN = x)
  c_ <- comp %>% filter(study == s, stage == "carcinoma") %>% select(finaltype, wC = w, xC = x)
  j <- full_join(n, c_, by = "finaltype") %>%
    mutate(across(c(wN, wC), ~replace_na(.x, 0)),
           xN = ifelse(is.na(xN), 0, xN), xC = ifelse(is.na(xC), 0, xC))
  tibble(study = s,
         total_change  = sum(j$wC * j$xC) - sum(j$wN * j$xN),
         composition   = sum((j$wC - j$wN) * j$xN),
         within_state  = sum(j$wN * (j$xC - j$xN)),
         interaction   = sum((j$wC - j$wN) * (j$xC - j$xN)))
})
print(as.data.frame(dec %>% mutate(across(where(is.numeric), ~signif(.x, 3)),
  pct_composition = round(100*composition/total_change),
  pct_within      = round(100*within_state/total_change),
  pct_interaction = round(100*interaction/total_change))), row.names = FALSE)
write_csv(dec, file.path(OUT, "shift_share_decomposition.csv"))
say("\n  composition  = change from cell states shifting in relative abundance")
say("  within_state = change from states altering their own allocation")
say("  interaction  = joint term (states that both grew and changed)")
say("  Percentages can exceed 100 or be negative when terms oppose.")

# ---- figure ------------------------------------------------------------------
pd <- dec %>% select(study, composition, within_state, interaction) %>%
  pivot_longer(-study, names_to = "term", values_to = "value")
p <- ggplot(pd, aes(term, value, fill = term)) +
  geom_col() + facet_wrap(~ study) +
  geom_hline(yintercept = 0, colour = "grey40") +
  scale_fill_manual(values = c(composition = "#b2182b", within_state = "#2166ac",
                               interaction = "grey70"), guide = "none") +
  labs(title = sprintf("What drives the rise in myeloid export allocation? (%s)", ARM),
       subtitle = "Decomposition of the normal-to-carcinoma change into shifts in cell-state abundance\n(composition) and changes within states (within-state).",
       x = NULL, y = "contribution to change in export allocation") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(OUT, "fig_decomposition.png"), p, width = 8, height = 5, dpi = 150, bg = "white")
say("\n-> %s", OUT)
