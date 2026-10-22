#!/usr/bin/env Rscript
# ==============================================================================
# fig1_atlas_overview.R  --  Figure 1: the atlas
#
#   Rscript fig1_atlas_overview.R --arm bulk_before
#
# Replaces the coupling figure as the opening display item. A reader meeting the
# methods argument first has no picture of what data exists; this establishes the
# resource and the vocabulary (unit, cell subtype, tissue class) before either is
# used. The flux-expression coupling is a single number in the text and a
# supplementary panel, not two main-text panels.
#
#   1A  what each cohort contributes, by tissue class. Shows immediately that
#       staged tissue exists in only two cohorts, which is why progression
#       analyses are restricted to them.
#   1B  cells by major cell type, and the myeloid subtypes resolved within the compartment
#       the paper concerns.
#   1C  units per cohort and tissue class. This is the data-density table: where
#       a claim rests on many units and where it rests on few.
#
# Output -> results/all_sites/<ARM>/fig1/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

OUT <- file.path(PATHS$results, "fig1")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

NL_T <- c("NL","Normal"); ADJ_T <- "NL-contact"; PRE_T <- c("Pre","LP")
CA_T <- c("CA","CA-E","CA-A","CA-c","CA-p")
TCLASS <- c("normal","tumor-adjacent","premalignant","carcinoma")   # "other" (lymph node,
# recurrent, metastatic) is excluded from every analysis and is not shown

cw_all <- as.data.frame(readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))) %>%
  filter(n_cells >= GATES$min_cells)

# Atlas totals are reported before restricting to the oral axis, since the figure
# describes the resource. Lymph node, recurrent and metastatic tissue are excluded
# from every analysis in this work and are not shown in the panels.
say("=== atlas, all tissue ===")
say("  units %s | patients %d | cohorts %d | cell subtypes %d",
    format(nrow(cw_all), big.mark = ","), n_distinct(cw_all$patient),
    n_distinct(cw_all$study), n_distinct(cw_all$finaltype))
say("  cells in retained units %s (645,351 across all 8,639 possible units)",
    format(sum(cw_all$n_cells), big.mark = ","))

cw <- cw_all %>%
  mutate(tclass = case_when(tissue %in% NL_T  ~ "normal",
                            tissue == ADJ_T   ~ "tumor-adjacent",
                            tissue %in% PRE_T ~ "premalignant",
                            tissue %in% CA_T  ~ "carcinoma",
                            TRUE ~ NA_character_)) %>%
  filter(!is.na(tclass)) %>%
  mutate(tclass = factor(tclass, levels = TCLASS))

say("\n=== oral axis only (shown in panels) ===")
say("  units %s | patients %d | cells %s",
    format(nrow(cw), big.mark = ","), n_distinct(cw$patient),
    format(sum(cw$n_cells), big.mark = ","))
say("  excluded: %s units of lymph node, recurrent or metastatic tissue",
    format(nrow(cw_all) - nrow(cw), big.mark = ","))

# ============================ PANEL 1A ======================================
# cells contributed per cohort per tissue class
a <- cw %>% group_by(study, tclass) %>%
  summarise(cells = sum(n_cells), .groups = "drop")
staged <- cw %>% group_by(study) %>%
  summarise(has_norm = any(tclass == "normal"), has_ca = any(tclass == "carcinoma"),
            .groups = "drop") %>% filter(has_norm & has_ca) %>% pull(study)
say("\n=== 1A: cohorts with both normal and carcinoma: %s ===", paste(staged, collapse = ", "))
print(as.data.frame(a %>% pivot_wider(names_from = tclass, values_from = cells,
                                      values_fill = 0)), row.names = FALSE)

a <- a %>% mutate(lab = ifelse(study %in% staged, paste0(study, " *"), study))
A <- ggplot(a, aes(cells/1000, reorder(lab, cells, sum), fill = tclass)) +
  geom_col(width = .68) +
  scale_fill_brewer(palette = "RdYlBu", direction = -1, name = NULL, drop = FALSE) +
  labs(title = "A   Cell contribution by cohort and tissue class",
       x = "cells (thousands)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())

# ============================ PANEL 1B ======================================
b_lin <- cw %>% group_by(annot_4th) %>%
  summarise(cells = sum(n_cells), subtypes = n_distinct(finaltype), .groups = "drop") %>%
  arrange(desc(cells))
say("\n=== 1B: cells by major cell type ===")
print(as.data.frame(b_lin), row.names = FALSE)

B <- ggplot(b_lin, aes(cells/1000, reorder(annot_4th, cells))) +
  geom_col(fill = "grey60", width = .68) +
  geom_text(aes(label = sprintf("%d subtype%s", subtypes, ifelse(subtypes == 1, "", "s"))), hjust = -0.12,
            size = 2.8, colour = "grey35") +
  scale_x_continuous(expand = expansion(mult = c(0, .18))) +
  labs(title = "B   Cells and subtypes by major cell type",
       x = "cells (thousands)", y = NULL) +
  theme_bw(base_size = 10) + theme(panel.grid.major.y = element_blank())

# ============================ PANEL 1C ======================================
cc <- cw %>% count(study, tclass, name = "units")
say("\n=== 1C: units per cohort x tissue class ===")
print(as.data.frame(cc %>% pivot_wider(names_from = tclass, values_from = units,
                                       values_fill = 0)), row.names = FALSE)

C <- ggplot(cc, aes(tclass, study, fill = units)) +
  geom_tile(colour = "white", linewidth = .8) +
  geom_text(aes(label = units, colour = units > 250), size = 3, show.legend = FALSE) +
  scale_fill_gradient(low = "grey93", high = "#2166ac", name = "units") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey25")) +
  labs(title = "C   Units per cohort and tissue class",
       x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1))

fig <- (A | B) / C + plot_layout(heights = c(1, 0.85))
ggsave(file.path(OUT, "fig1_atlas.png"), fig, width = 12.5, height = 8.5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig1_atlas.pdf"), fig, width = 12.5, height = 8.5, bg = "white")
write_csv(cc, file.path(OUT, "fig1C_units.csv"))
say("\n-> %s", OUT)

say("\nCaption draft:")
say("  Our cell atlas. 645,351 cells from 114 patients across 5 published cohorts,")
say("  provided as 167 samples; patients contributing more than one tissue are")
say("  represented separately for each. The unit of analysis for the study, and")
say("  depicted in panel C, is one value for each patient, tissue and cell subtype")
say("  combination (patient x tissue x cell subtype), requiring at least ten cells,")
say("  which gives 5,459 units across all tissues. Panels show the oral axis (3,818")
say("  units, 453,411 cells from 82 patients); lymph node, recurrent and metastatic")
say("  tissue are excluded from every analysis in this work.")
say("  (A) Cells contributed by each cohort, by tissue class, counting only cells in")
say("  units that met the ten-cell threshold; asterisks mark the two cohorts containing")
say("  both normal mucosa and carcinoma, to which progression analyses are restricted.")
say("  (B) Cells and cell subtypes by major cell type. (C) Units available for each")
say("  cohort and tissue class, showing where analyses rest on many units and where on")
say("  few. Empty cells mean the combination was not sampled or did not reach the")
say("  ten-cell threshold. Cells were never pooled across cohorts: statistics were")
say("  computed within cohort and combined afterwards.")
