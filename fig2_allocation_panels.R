#!/usr/bin/env Rscript
# ==============================================================================
# fig2_allocation_panels.R  --  Figure 2: "the move that makes it visible"
#
#   Rscript fig2_allocation_panels.R --arm bulk_before
#
# Two panels making the case for analysing allocation rather than magnitude.
#
#   2A  direction balance. On magnitude nearly every change is an increase, which
#       is what an analysis anchored to total metabolic-gene expression produces.
#       On allocation the balance restores. This is the clearest single argument
#       for the approach: a magnitude analysis cannot show reductions when the
#       quantity it measures rises almost everywhere.
#
#   2B  the allocation result. Volcano of the two-cohort meta-analysis, with the
#       19 surviving comparisons labelled.
#
# NOTE ON THE TWO DIRECTION FIGURES. The balance across all tested comparisons
# (38% up / 62% down) and the direction of the survivors (18 of 19 up) describe
# different things and both belong in the figure. Shares sum to one, so a few
# modules gaining substantial share requires many others to lose a little; the
# gains concentrate and clear an effect threshold, the losses spread and do not.
# Panel A shows the first, panel B the second.
#
# Sources: 05_progression (magnitude), 05c_progression_meta (allocation).
#
# Output -> results/all_sites/<ARM>/fig2/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(stringr) })

OUT <- file.path(PATHS$results, "fig2")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MAG  <- read_csv(file.path(PATHS$results, "05_progression/finaltype/progression_all.csv"),
                 show_col_types = FALSE)
ALLO <- read_csv(file.path(PATHS$results, "05c_progression_meta/finaltype/meta_all.csv"),
                 show_col_types = FALSE)

# ============================ PANEL 2A ======================================
# The magnitude table is per cohort; use He (the discovery cohort) so the two
# panels describe comparable numbers of tests.
mag_dir <- MAG %>% filter(study == "He", !is.na(direction)) %>% count(direction) %>%
  mutate(pct = 100 * n / sum(n), analysis = "magnitude\n(raw flux)")
allo_dir <- ALLO %>% filter(!is.na(direction)) %>% count(direction) %>%
  mutate(pct = 100 * n / sum(n), analysis = "allocation\n(flux share)")
A_dat <- bind_rows(mag_dir, allo_dir) %>%
  mutate(analysis = factor(analysis, levels = c("magnitude\n(raw flux)","allocation\n(flux share)")),
         direction = factor(direction, levels = c("down","up")))

say("=== 2A direction balance ===")
print(as.data.frame(A_dat %>% transmute(analysis = gsub("\n"," ",analysis),
                                        direction, n, pct = round(pct))), row.names = FALSE)

A <- ggplot(A_dat, aes(analysis, pct, fill = direction)) +
  geom_col(width = .6) +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            position = position_stack(vjust = .5), colour = "white", size = 4.2) +
  scale_fill_manual(values = c(up = "#b2182b", down = "#4393c3"),
                    labels = c(up = "increase", down = "decrease"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .02))) +
  labs(title = "A   Direction of change, magnitude vs allocation",
       x = NULL, y = "% of comparisons") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

# ============================ PANEL 2B ======================================
# Survivors are distinguished by whether their module meets the flux-enzyme
# criterion in that cell type (07c). Ten of the nineteen do not. Because that
# criterion is applied as an exclusion rule before interpretation, the figure
# should not give those ten equal weight. The two lactate-export survivors are annotated, since they carry the
# analyses that follow. SPP1+ macrophages do not appear here and cannot: they have
# no normal-tissue units, so no progression comparison exists for them.
GR <- read_csv(file.path(PATHS$results, "07c_grounding_final/module_grounding_final.csv"),
               show_col_types = FALSE) %>%
  select(celltype, module, GROUNDED, best_gene, best_rho, null_bar)

B_dat <- ALLO %>% filter(is.finite(delta_meta), is.finite(fdr)) %>%
  left_join(GR, by = c("celltype","module")) %>%
  mutate(GROUNDED = ifelse(is.na(GROUNDED), FALSE, GROUNDED),
         class = case_when(!survivor ~ "not significant",
                           survivor & GROUNDED ~ "meets flux-enzyme criterion",
                           TRUE ~ "does not meet criterion"),
         lab = ifelse(survivor, paste0(str_trunc(unit, 14), " · ", module), NA),
         is_export = survivor & module == "M_78")

say("\n=== 2B: %d tested | %d survivors (%d meet criterion, %d do not) ===",
    nrow(B_dat), sum(B_dat$survivor),
    sum(B_dat$survivor & B_dat$GROUNDED), sum(B_dat$survivor & !B_dat$GROUNDED))
say("  lactate-export survivors: %s",
    paste(B_dat$unit[B_dat$is_export], collapse = ", "))

B <- ggplot(B_dat, aes(delta_meta, -log10(fdr))) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_hline(yintercept = -log10(FDR_THRESH), linetype = 2, colour = "grey60") +
  geom_point(data = filter(B_dat, !survivor), colour = "grey85", size = 1, alpha = .5) +
  geom_point(data = filter(B_dat, survivor),
             aes(fill = class), shape = 21, colour = "grey30", size = 2.6, stroke = .4) +
  geom_point(data = filter(B_dat, is_export), shape = 21, fill = NA,
             colour = "#b2182b", size = 5, stroke = 1) +
  ggrepel::geom_text_repel(aes(label = lab, fontface = ifelse(is_export, "bold", "plain")),
                           size = 2.5, max.overlaps = 25,
                           segment.colour = "grey70", min.segment.length = 0) +
  scale_fill_manual(values = c(`meets flux-enzyme criterion` = "#b2182b",
                               `does not meet criterion` = "grey65"), name = NULL) +
  labs(title = "B   Changes in allocation from normal mucosa to carcinoma",
       x = "change in allocation (carcinoma minus normal)", y = expression(-log[10]~FDR)) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

# ============================ PANEL 2C ======================================
# The 19 survivors, each shown against the flux-enzyme criterion. Panel B encodes
# this as fill colour; here it is the subject. A survivor is interpreted only if
# its module's flux tracks its own catalytic genes above that cell type's null
# threshold, and ten of the nineteen do not.
C_dat <- B_dat %>% filter(survivor) %>%
  mutate(lab2 = paste0(str_trunc(unit, 16), " \u00b7 ", module),
         gene_lab = ifelse(is.na(best_gene), "not tested",
                    sprintf("%s %.2f / %.2f", best_gene, best_rho, null_bar))) %>%
  arrange(GROUNDED, desc(abs(delta_meta))) %>%
  mutate(lab2 = factor(lab2, levels = lab2))

say("\n=== 2C: %d survivors | %d meet the criterion, %d do not ===",
    nrow(C_dat), sum(C_dat$GROUNDED), sum(!C_dat$GROUNDED))

Cp <- ggplot(C_dat, aes(delta_meta, lab2, fill = class)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_col(width = .68) +
  geom_text(aes(label = gene_lab,
                hjust = ifelse(delta_meta > 0, -0.08, 1.08)),
            size = 2.4, colour = "grey30") +
  scale_fill_manual(values = c(`meets flux-enzyme criterion` = "#b2182b",
                               `does not meet criterion` = "grey65"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(.35, .35))) +
  labs(title = "C   Criterion status of the 19 surviving comparisons",
       x = "change in allocation (carcinoma minus normal)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 7))

fig <- (A + B + plot_layout(widths = c(1, 1.9))) / Cp + plot_layout(heights = c(1, 1.15))
ggsave(file.path(OUT, "fig2_allocation.png"), fig, width = 12.5, height = 10.5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig2_allocation.pdf"), fig, width = 12.5, height = 10.5, bg = "white")
write_csv(C_dat %>% select(celltype, unit, module, reaction, delta_meta, fdr,
                           best_gene, best_rho, null_bar, GROUNDED),
          file.path(OUT, "fig2C_survivor_criterion.csv"))
write_csv(A_dat, file.path(OUT, "fig2A_direction_balance.csv"))
say("\n-> %s", OUT)

say("\nCaption draft:")
say("  Analysing allocation rather than magnitude changes what is visible.")
say("  (A) On raw flux, %.0f%% of tested changes are increases, consistent with an", 
    mag_dir$pct[mag_dir$direction=="up"])
say("  estimate anchored to total metabolic-gene expression. On allocation the balance")
say("  restores to %.0f%% increases and %.0f%% decreases.",
    allo_dir$pct[allo_dir$direction=="up"], allo_dir$pct[allo_dir$direction=="down"])
say("  (B) Meta-analysis of allocation change across the two staged cohorts. Nineteen")
say("  comparisons survive the effect floor, FDR threshold and concordance requirement,")
say("  of which %d are increases: shares sum to one, so gains concentrate in a few",
    sum(B_dat$survivor & B_dat$delta_meta > 0))
say("  modules while losses spread thinly and do not clear an effect threshold. %d of",
    sum(B_dat$survivor & B_dat$GROUNDED))
say("  the nineteen meet the flux-enzyme criterion (Methods); the two lactate-export")
say("  comparisons, circled, are among them and carry the analyses that follow.")
say("  (C) The nineteen survivors against the flux-enzyme criterion. Each is annotated")
say("  with its best-correlating catalytic gene, that correlation, and the null")
say("  threshold for its cell type. %d do not clear the threshold and are not",
    sum(B_dat$survivor & !B_dat$GROUNDED))
say("  interpreted. Units are patient x tissue x cell subtype pseudobulks of at least")
say("  10 cells.")
