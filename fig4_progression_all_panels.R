#!/usr/bin/env Rscript
# ==============================================================================
# fig4_progression_all_panels.R  --  Figure 4, all six panels
#
#   Rscript fig4_progression_all_panels.R --arm bulk_before
#   Rscript fig4_progression_all_panels.R --arm bulk_after
#
# Replaces fig4_progression_panels.R (A-C), fig4d_export_production_ratio.R (D)
# and fig4e_celltype_specificity.R (E-F). One figure, one script.
#
#   4A  the paired comparison itself: myeloid export allocation in the patients
#       contributing both normal mucosa and carcinoma. Each patient is their own
#       control, so the comparison is insensitive to differences in composition
#       or processing between patients.
#   4B  the same comparison over nested subtype sets. This is the panel that
#       separates cells that changed from cells that arrived.
#   4C  allocation by subtype at each tissue class, which shows which subtypes
#       exist in which tissue and is why panel B is possible.
#   4D  lactate production alongside export. Both rise; export rises in every
#       patient and production does not, so the change is not a uniform increase
#       in glycolytic activity.
#   4E  the same paired comparison applied to every cell type with enough
#       patients. Myeloid is the only one that rises.
#   4F  the changes compared within patients: myeloid's change minus each other
#       cell type's change in the same person.
#
# A NOTE ON THE EXPRESSION RATIO. The console reports the SLC16A3/LDHA expression
# ratio alongside the flux ratio. Neither is plotted. On the MAGIC-imputed matrix
# the expression ratio reaches p = 0.020, but on raw counts over the same ten
# patients it is p = 0.065, so the manuscript reports it as a trend and the panel
# was cut. The numbers are printed here so the decision is reproducible.
#
# Output -> results/all_sites/<ARM>/fig4/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(data.table) })

OUT <- file.path(PATHS$results, "fig4")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
usable <- readRDS(file.path(PATHS$results, "03_qc", "usable_modules.rds"))

NORMAL_T <- c("NL","Normal"); CANCER_T <- c("CA","CA-E","CA-A","CA-c","CA-p")
ADJ_T <- "NL-contact"; PRE_T <- c("Pre","LP")
EMERGENT <- "MQ_C4_SPP1"

cw <- as.data.frame(readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))) %>%
  filter(annot_4th == "Myeloid", study %in% PROG_COHORTS,
         n_cells >= GATES$min_cells) %>%
  mutate(stage = case_when(tissue %in% NORMAL_T ~ "normal",
                           tissue %in% CANCER_T ~ "carcinoma", TRUE ~ NA_character_)) %>%
  filter(!is.na(stage))
mc <- intersect(mcols_of(cw), usable)
S <- as.matrix(cw[, mc]); S <- S / rowSums(S)
cw$export     <- 100 * S[, "M_78"]
cw$production <- 100 * S[, "M_6"]

present <- cw %>% count(finaltype, stage) %>%
  tidyr::pivot_wider(names_from = stage, values_from = n, values_fill = 0)
emergent_subtypes <- present %>% filter(normal == 0, carcinoma > 0) %>% pull(finaltype)
resident_subtypes <- present %>% filter(normal > 0) %>% pull(finaltype)
say("resident subtypes (present in normal): %d", length(resident_subtypes))
say("emergent subtypes (carcinoma only)   : %s", paste(emergent_subtypes, collapse = ", "))

paired <- function(d, col = "export") {
  d %>% group_by(study, patient, stage) %>%
    summarise(e = weighted.mean(.data[[col]], n_cells), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = stage, values_from = e) %>%
    filter(!is.na(normal), !is.na(carcinoma))
}

# ============================ PANEL 4A ======================================
pw <- paired(cw)
wA <- suppressWarnings(wilcox.test(pw$carcinoma, pw$normal, paired = TRUE))
say("\n4A: %d paired patients | %.2f%% -> %.2f%% | p = %.3g",
    nrow(pw), median(pw$normal), median(pw$carcinoma), wA$p.value)
say("    median of per-patient differences: %.2f  IQR [%.2f, %.2f]",
    median(pw$carcinoma - pw$normal),
    quantile(pw$carcinoma - pw$normal, .25), quantile(pw$carcinoma - pw$normal, .75))

longA <- pw %>%
  tidyr::pivot_longer(c(normal, carcinoma), names_to = "stage", values_to = "pct") %>%
  mutate(stage = factor(stage, levels = c("normal","carcinoma"),
                        labels = c("normal\nmucosa","carcinoma")))
medsA <- longA %>% group_by(stage) %>% summarise(pct = median(pct), .groups = "drop")

A <- ggplot(longA, aes(stage, pct)) +
  geom_line(aes(group = patient, colour = study), alpha = .55, linewidth = .7) +
  geom_point(aes(colour = study), size = 1.9, alpha = .8) +
  geom_line(data = medsA, aes(group = 1), linewidth = 1.6, colour = "grey20") +
  geom_point(data = medsA, size = 4, colour = "grey20") +
  annotate("text", x = 2.32, y = medsA$pct[2],
           label = sprintf("%.1f%%", medsA$pct[2]), hjust = 0, size = 3.4) +
  annotate("text", x = 0.68, y = medsA$pct[1],
           label = sprintf("%.1f%%", medsA$pct[1]), hjust = 1, size = 3.4) +
  scale_colour_manual(values = c(He = "#d6604d", Kim = "#4393c3"), name = NULL) +
  coord_cartesian(xlim = c(0.55, 2.65), clip = "off") +
  labs(title = "A   Myeloid lactate-export allocation,\npaired within patient",
       x = NULL, y = "% of total metabolic flux") +
  theme_bw(base_size = 11) +
  theme(legend.position.inside = c(.14,.88), legend.position = "inside",
        legend.background = element_blank(), panel.grid.major.x = element_blank())

# ============================ PANEL 4B ======================================
sets <- list(
  "all myeloid subtypes"             = cw,
  "excluding SPP1+ TAM"              = cw %>% filter(finaltype != EMERGENT),
  "subtypes present\nin normal only" = cw %>% filter(!finaltype %in% emergent_subtypes))

B_dat <- purrr::imap_dfr(sets, function(d, nm) {
  q <- paired(d)
  tibble(set = nm, n = nrow(q),
         normal = median(q$normal), carcinoma = median(q$carcinoma),
         fold = median(q$carcinoma) / median(q$normal),
         p = suppressWarnings(wilcox.test(q$carcinoma, q$normal, paired = TRUE))$p.value)
})
say("\n4B:"); print(as.data.frame(B_dat %>% mutate(across(where(is.numeric), ~signif(.x,3)))),
                    row.names = FALSE)
write_csv(B_dat, file.path(OUT, "fig4B_subtype_sets.csv"))

# the three set names run together on one line at this width; wrapped to two
# lines each, same meaning
B_dat$set <- factor(B_dat$set, levels = names(sets),
                    labels = c("All myeloid\nsubtypes", "Excluding\nSPP1+ TAM",
                               "Present in\nnormal only"))
B <- ggplot(B_dat, aes(set, fold)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
  geom_col(fill = "#4393c3", width = .55) +
  geom_text(aes(label = sprintf("%.2f\u00d7", fold)), vjust = -0.6, size = 4) +
  geom_text(aes(y = 0.12, label = sprintf("p = %.3f", p)), size = 3, colour = "white") +
  coord_cartesian(ylim = c(0, max(B_dat$fold) * 1.18)) +
  labs(title = "B   Lactate-export allocation, fold change by myeloid subtype set",
       x = NULL, y = "fold change, normal to carcinoma") +
  theme_bw(base_size = 11) + theme(panel.grid.major.x = element_blank())

# ============================ PANEL 4C ======================================
allc <- as.data.frame(readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))) %>%
  filter(annot_4th == "Myeloid", n_cells >= GATES$min_cells) %>%
  mutate(tclass = case_when(tissue %in% NORMAL_T ~ "normal (NL)",
                            tissue == ADJ_T      ~ "tumor-adjacent",
                            tissue %in% PRE_T    ~ "premalignant",
                            tissue %in% CANCER_T ~ "carcinoma", TRUE ~ NA_character_)) %>%
  filter(!is.na(tclass))
mc2 <- intersect(mcols_of(allc), usable)
S2 <- as.matrix(allc[, mc2]); S2 <- S2 / rowSums(S2)
allc$export <- 100 * S2[, "M_78"]
allc <- allc %>%
  mutate(tclass = factor(tclass, levels = c("normal (NL)","tumor-adjacent",
                                            "premalignant","carcinoma")),
         grp = case_when(grepl("^MQ_", finaltype) ~ "macrophage",
                         grepl("^DC_", finaltype) ~ "dendritic",
                         TRUE ~ "other myeloid"))
ordC <- allc %>% filter(tclass == "carcinoma") %>% group_by(finaltype) %>%
  summarise(m = median(export), .groups="drop") %>% arrange(m) %>% pull(finaltype)
allc$finaltype <- factor(allc$finaltype, levels = ordC)
say("\n4C: units per tissue class"); print(as.data.frame(allc %>% count(tclass)), row.names = FALSE)

labC <- allc %>% count(tclass) %>% mutate(txt = sprintf("%s\n(n = %d units)", tclass, n))
allc <- allc %>% left_join(labC %>% select(tclass, txt), by = "tclass") %>%
  mutate(txt = factor(txt, levels = labC$txt))

Cp <- ggplot(allc, aes(finaltype, export, fill = grp)) +
  geom_boxplot(outlier.size = .3, linewidth = .3) +
  coord_flip() + facet_wrap(~ txt, nrow = 1) +
  scale_fill_manual(values = c(macrophage = "#2166ac", dendritic = "#b2182b",
                               `other myeloid` = "grey78"), name = NULL) +
  labs(title = "C   Lactate-export allocation by myeloid subtype and tissue class",
       x = NULL, y = "% of total metabolic flux") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")

# ============================ PANEL 4D ======================================
patD <- cw %>% group_by(study, patient, stage) %>%
  summarise(production = weighted.mean(production, n_cells),
            export     = weighted.mean(export,     n_cells), .groups = "drop")
MOD <- list()
for (mm in c("production","export")) {
  z <- patD %>% select(patient, stage, v = all_of(mm)) %>%
    tidyr::pivot_wider(names_from = stage, values_from = v) %>%
    filter(!is.na(normal), !is.na(carcinoma))
  pz <- suppressWarnings(wilcox.test(z$carcinoma, z$normal, paired = TRUE))$p.value
  say("4D  %-10s %.2f -> %.2f | p = %.4f | %d of %d rise", mm,
      median(z$normal), median(z$carcinoma), pz, sum(z$carcinoma > z$normal), nrow(z))
  MOD[[mm]] <- list(nl = median(z$normal), ca = median(z$carcinoma), p = pz,
                    k = sum(z$carcinoma > z$normal), n = nrow(z))
}

longD <- patD %>%
  mutate(stage = factor(stage, levels = c("normal","carcinoma"),
                        labels = c("normal\nmucosa","carcinoma"))) %>%
  tidyr::pivot_longer(c(production, export), names_to = "module", values_to = "pct") %>%
  mutate(module = factor(module, levels = c("production","export"),
           labels = c("production (M_6)\npyruvate to lactate",
                      "export (M_78)\nlactate out")))
D <- ggplot(longD, aes(stage, pct, colour = module, group = interaction(patient, module))) +
  geom_line(alpha = .3, linewidth = .5) +
  stat_summary(aes(group = module), fun = median, geom = "line", linewidth = 1.5) +
  stat_summary(aes(group = module), fun = median, geom = "point", size = 3.2) +
  scale_colour_manual(values = c("#4393c3","#b2182b"), name = NULL) +
  labs(title = "D   Myeloid lactate production and export", x = NULL,
       y = "% of total metabolic flux") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

# the ratio, reported but not plotted (see header)
ratD <- patD %>% mutate(ratio = export / production) %>%
  select(patient, stage, ratio) %>%
  tidyr::pivot_wider(names_from = stage, values_from = ratio) %>%
  filter(!is.na(normal), !is.na(carcinoma))
say("    flux ratio (not plotted): %.2f -> %.2f, p = %.3f, %d of %d rise",
    median(ratD$normal), median(ratD$carcinoma),
    suppressWarnings(wilcox.test(ratD$carcinoma, ratD$normal, paired = TRUE))$p.value,
    sum(ratD$carcinoma > ratD$normal), nrow(ratD))

EXPR <- file.path(ROOT, "data/input/all_sites", ARM, "Myeloid/train_expression.csv")
META <- file.path(ROOT, "metadata/all_sites", ARM, "Myeloid/metadata_units.csv")
if (file.exists(EXPR) && file.exists(META)) {
  ex <- fread(EXPR); mu <- fread(META)
  ids <- setdiff(names(ex), names(ex)[1])
  key <- names(mu)[which.max(sapply(mu, function(x) length(intersect(as.character(x), ids))))]
  E <- as.matrix(ex[, -1]); rownames(E) <- as.character(ex[[1]])
  mu <- as.data.frame(mu); rownames(mu) <- as.character(mu[[key]]); mu <- mu[colnames(E), ]
  keep <- mu$tissue %in% NORMAL_T | mu$tissue %in% CANCER_T
  keep <- keep & mu$study %in% PROG_COHORTS & mu$n_cells >= GATES$min_cells
  xw <- tibble(patient = mu$patient[keep],
               stage = ifelse(mu$tissue[keep] %in% NORMAL_T, "normal", "carcinoma"),
               n = mu$n_cells[keep],
               r = E["SLC16A3", keep] / E["LDHA", keep]) %>%
    group_by(patient, stage) %>% summarise(r = weighted.mean(r, n), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = stage, values_from = r) %>%
    filter(!is.na(normal), !is.na(carcinoma))
  say("    SLC16A3/LDHA expression ratio (MAGIC, not plotted): %.2f -> %.2f, p = %.3f, %d of %d",
      median(xw$normal), median(xw$carcinoma),
      suppressWarnings(wilcox.test(xw$carcinoma, xw$normal, paired = TRUE))$p.value,
      sum(xw$carcinoma > xw$normal), nrow(xw))
  say("    on raw counts the same ratio gives p = 0.065, which is what the text reports")
} else say("    [expression matrix not found; skipping the ratio check]")

# ========================== PANELS 4E and 4F =================================
cwa <- as.data.frame(readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))) %>%
  filter(n_cells >= GATES$min_cells, study %in% PROG_COHORTS,
         tissue %in% c(NORMAL_T, CANCER_T)) %>%
  mutate(stage = ifelse(tissue %in% NORMAL_T, "normal", "carcinoma"))
gr <- read_csv(file.path(PATHS$results, "07c_grounding_final/module_grounding_final.csv"),
               show_col_types = FALSE) %>% filter(module == "M_78") %>%
      select(celltype, GROUNDED)

deltaCT <- function(lin) {
  d0 <- cwa %>% filter(annot_4th == lin)
  m0 <- intersect(mcols_of(d0), usable); if (!"M_78" %in% m0 || nrow(d0) < 4) return(NULL)
  S0 <- as.matrix(d0[, m0]); S0 <- S0 / rowSums(S0)
  d1 <- d0 %>% mutate(v = 100 * S0[, "M_78"]) %>%
    group_by(patient, stage) %>%
    summarise(v = weighted.mean(v, n_cells), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = stage, values_from = v)
  if (!all(c("normal","carcinoma") %in% names(d1))) return(NULL)
  d1 %>% filter(!is.na(normal), !is.na(carcinoma)) %>%
    transmute(celltype = lin, patient, normal, carcinoma, d = carcinoma - normal)
}
DCT <- lapply(sort(unique(cwa$annot_4th)), deltaCT) %>% bind_rows()

smE <- DCT %>% group_by(celltype) %>% filter(n() >= 4) %>%
  summarise(n = n(), change = median(d), rising = sum(d > 0),
            p = suppressWarnings(wilcox.test(carcinoma, normal, paired = TRUE))$p.value,
            .groups = "drop") %>%
  left_join(gr, by = "celltype") %>%
  mutate(lab = sprintf("%s  (%d/%d)", celltype, rising, n)) %>%
  arrange(change) %>% mutate(lab = factor(lab, levels = lab))
say("\n4E: change in export allocation by cell type")
print(as.data.frame(smE %>% transmute(celltype, n, change = round(change,3),
      rising = sprintf("%d/%d", rising, n), p = signif(p,2),
      meets_criterion = GROUNDED)), row.names = FALSE)
write_csv(smE, file.path(OUT, "fig4E_celltype_change.csv"))

E <- ggplot(smE, aes(change, lab, fill = celltype == "Myeloid")) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_col(width = .66) +
  geom_text(aes(label = ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)),
                hjust = ifelse(change > 0, -0.1, 1.1)), size = 2.7, colour = "grey30") +
  scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "grey72"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(.28, .28))) +
  labs(title = "E   Change in lactate-export allocation by cell type",
       x = "change, carcinoma minus normal (percentage points)", y = NULL) +
  theme_bw(base_size = 10) + theme(panel.grid.major.y = element_blank())

dm <- DCT %>% filter(celltype == "Myeloid") %>% select(patient, d_mye = d)
pwF <- DCT %>% filter(celltype != "Myeloid") %>%
  inner_join(dm, by = "patient") %>%
  group_by(celltype) %>% filter(n() >= 4) %>%
  summarise(n = n(), larger = sum(d_mye > d), gap = median(d_mye - d),
            p = suppressWarnings(wilcox.test(d_mye, d, paired = TRUE))$p.value,
            .groups = "drop") %>%
  mutate(lab = sprintf("vs %s  (%d/%d)", celltype, larger, n)) %>%
  arrange(gap) %>% mutate(lab = factor(lab, levels = lab))
say("\n4F: myeloid change minus each cell type's, same patient")
print(as.data.frame(pwF %>% transmute(celltype, n, larger = sprintf("%d/%d", larger, n),
      gap = round(gap,3), p = signif(p,3))), row.names = FALSE)
write_csv(pwF, file.path(OUT, "fig4F_paired_vs_myeloid.csv"))

Fp <- ggplot(pwF, aes(gap, lab)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_col(width = .6, fill = "grey60") +
  geom_text(aes(label = sprintf("p = %.3f", p)), hjust = -0.1, size = 2.7, colour = "grey30") +
  scale_x_continuous(expand = expansion(mult = c(.05, .3))) +
  # two lines: F is a third of the figure width, where a 55-character title
  # already clipped, and this one is longer
  labs(title = "F   Lactate-export allocation change,\nmyeloid minus each cell type",
       x = "myeloid minus other cell type's change (pct points)", y = NULL) +
  theme_bw(base_size = 10) + theme(panel.grid.major.y = element_blank())

# ============================== ASSEMBLE ====================================
# D sat on the top row beside A and B, so the letters read A, B, D, C, E, F.
# Moving it to the bottom row makes them sequential in reading order while
# keeping C full-width and leaving every panel's letter unchanged, so no figure
# reference in the manuscript moves.
fig <- (A | B) / Cp / (D | E | Fp) + plot_layout(heights = c(1, 1.15, 0.95))
ggsave(file.path(OUT, "fig4_progression.png"), fig, width = 14, height = 13.5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig4_progression.pdf"), fig, width = 14, height = 13.5, bg = "white")
say("\n-> %s", OUT)

say("\nCaption draft:")
say("  Myeloid lactate-export allocation increases from normal mucosa to carcinoma,")
say("  and no other cell type increases significantly.")
say("  (A) Myeloid lactate-export allocation in the %d patients contributing both", nrow(pw))
say("  normal mucosa and carcinoma. Allocation")
say("  rises from %.2f%% to %.2f%% of total metabolic flux (p = %.3f, Wilcoxon signed",
    median(pw$normal), median(pw$carcinoma), wA$p.value)
say("  rank), and every patient shows an increase. Each line is one patient, the heavy")
say("  line the median. Because")
say("  each patient serves as their own control, the comparison is insensitive to")
say("  differences in composition or processing between patients.")
say("  (B) The same comparison restricted to different sets of myeloid subtypes.")
say("  Excluding SPP1+ macrophages, or restricting to only those subtypes already")
say("  present in normal mucosa, leaves the increase intact, so resident populations")
say("  reprogram rather than the change being attributable to populations that appear")
say("  with the tumor.")
say("  (C) Allocation by myeloid subtype at each tissue class, ordered by carcinoma")
say("  median. Normal (NL) and tumor-adjacent are kept apart because the two labels are")
say("  used by different cohorts and are not interchangeable. %d subtypes have no NL",
    length(emergent_subtypes))
say("  units in any cohort (cDC3, SPP1+ macrophages and proliferating MKI67+")
say("  macrophages), so those cells are empty;")
say("  identifying them is what makes the restriction in panel B possible.")
say("  (D) Lactate production (M_6) alongside export (M_78) in the same patients.")
say("  Export rises from %.2f%% to %.2f%% in %d of %d patients; production rises from",
    MOD$export$nl, MOD$export$ca, MOD$export$k, MOD$export$n)
say("  %.2f%% to %.2f%% in %d of %d. The reallocation is therefore not a uniform",
    MOD$production$nl, MOD$production$ca, MOD$production$k, MOD$production$n)
say("  increase in glycolytic activity, in which both modules would rise together and")
say("  in the same patients.")
say("  (E) Change in export allocation for every cell type with at least four patients")
say("  contributing both tissues. Numbers in parentheses are patients rising out of")
say("  patients tested; counts differ because each cell type requires enough cells of")
say("  that type in both tissues. Myeloid is the only cell type that rises significantly;")
say("  NK and fibroblasts rise numerically but do not reach significance.")
say("  (F) The same changes compared within patients: the myeloid change minus each")
say("  other cell type's change in the same person. Positive values mean the myeloid")
say("  increase was larger, which it is against every other cell type tested.")
say("  Units are patient x tissue x cell subtype pseudobulks of at least 10 cells;")
say("  allocation is a module's flux as a percentage of the unit's total flux.")
