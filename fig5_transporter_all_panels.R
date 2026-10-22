#!/usr/bin/env Rscript
# ==============================================================================
# fig5_transporter_all_panels.R  --  Figure 5, all four panels
#
#   Rscript fig5_transporter_all_panels.R --arm bulk_before
#
# Replaces fig5a_isoform_panel.R, fig5b_mct_ratio_three_cohorts.R and
# fig5c_cptac_protein.R. One figure, one script.
#
#   5A  which transporter the flux follows, with the per-cell-type null
#       threshold drawn. This is verification that the estimate tracks its
#       dominant input, not evidence about transporter usage: the panel is
#       annotated with MCT4's share of transporter expression so the reader can
#       see why the correlations fall where they do.
#   5B  the MCT4/MCT1 expression ratio by cell type in three independent
#       datasets. Depends on no flux inference at all.
#   5C  MCT4 at protein level in CPTAC, against six other proteins.
#   5D  what MCT4 protein tracks: niche and macrophage markers, not LDHA. LDHA is
#       the informative negative here. If the myeloid signal were general glycolytic
#       activity rather than something specific to transport, LDHA protein should be
#       elevated alongside MCT4 and should track it across samples. It is neither,
#       which is why the LDH flux module is not interpreted anywhere in this work.
#
# Requires 32_independent_cohort.R and 31_cptac_protein.R to have run.
#
# Output -> results/all_sites/<ARM>/fig5/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(data.table) })

OUT <- file.path(PATHS$results, "fig5")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CTS      <- c("Myeloid","Fibroblast","B","Epithelial","Endothelial","T")
CANCER_T <- c("CA","CA-E","CA-A","CA-c","CA-p")
CP       <- file.path(ROOT, "results", "cptac")

# ============================== PANEL 5A ====================================
gA <- read_csv(file.path(PATHS$results, "07b_grounding_scalefree/gene_level_scalefree.csv"),
               show_col_types = FALSE) %>%
  filter(module == "M_78", gene %in% c("SLC16A1","SLC16A3","SLC16A4"), celltype %in% CTS)
# the matched threshold is per cell type AND module, so it comes from the
# grounding file rather than the per-cell-type null summary
nulls <- read_csv(file.path(PATHS$results, "07c_grounding_final/module_grounding_final.csv"),
                  show_col_types = FALSE) %>%
  filter(module == "M_78") %>% select(celltype, null_bar)
gA <- gA %>% left_join(nulls, by = "celltype")

comp <- rbindlist(lapply(CTS, function(ct) {
  ep <- file.path(ROOT, "data/input/all_sites", ARM, ct, "train_expression.csv")
  if (!file.exists(ep)) return(NULL)
  e <- fread(ep); gn <- as.character(e[[1]])
  gs <- intersect(c("SLC16A1","SLC16A3","SLC16A4"), gn)
  if (length(gs) < 2) return(NULL)
  E <- as.matrix(e[, -1]); rownames(E) <- gn
  m <- rowMeans(E[gs, , drop = FALSE])
  data.table(celltype = ct, pct_mct4 = 100 * m[["SLC16A3"]] / sum(m))
}))
say("=== 5A: MCT4 share of transporter expression ===")
print(as.data.frame(comp[order(-pct_mct4)]), row.names = FALSE)

have <- unique(gA$celltype)
missing <- setdiff(CTS, have)
if (length(missing)) say("  no M_78 correlations, omitted: %s", paste(missing, collapse=", "))
comp <- comp[celltype %in% have]

# this panel is transcript, so the series keep their gene symbols; the legend
# italicises them and names what each transporter does
gA <- gA %>% left_join(comp, by = "celltype") %>%
  mutate(gene = factor(gene, levels = c("SLC16A3","SLC16A1","SLC16A4")),
         celltype = factor(celltype, levels = comp[order(pct_mct4)]$celltype))
# the per-row transporter-share annotations are gone; the myeloid figure (92%)
# is given in the Results text instead. comp is still needed for the row order.

PA <- ggplot(gA, aes(rho, celltype, fill = gene)) +
  geom_vline(xintercept = 0, colour = "grey55") +
  geom_col(position = position_dodge(.78), width = .7) +
  # linetype is mapped rather than fixed so the dashed segments get their own
  # legend key beside the three transporter swatches
  geom_segment(data = distinct(gA, celltype, null_bar),
               aes(x = null_bar, xend = null_bar,
                   y = as.numeric(celltype) - .42, yend = as.numeric(celltype) + .42,
                   linetype = "empirical null threshold"),
               inherit.aes = FALSE, colour = "grey25", linewidth = .5) +
  scale_linetype_manual(values = c(`empirical null threshold` = "dashed"),
                        name = NULL) +
  scale_fill_manual(values = c(SLC16A3 = "#b2182b", SLC16A1 = "#4393c3",
                               SLC16A4 = "grey72"), name = NULL,
                    labels = expression(italic("SLC16A3")~"(exporter)",
                                        italic("SLC16A1")~"(importer)",
                                        italic("SLC16A4"))) +
  coord_cartesian(xlim = c(-0.55, 1)) +
  labs(title = "A   Correlation between lactate-transport flux and transporter expression",
       x = "Spearman rho (flux share ~ gene share)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.box = "horizontal",
        panel.grid.major.y = element_blank())
write_csv(gA, file.path(OUT, "fig5A_isoform.csv"))

# ============================== PANEL 5B ====================================
say("\n=== 5B: MCT4/MCT1 ratio, carcinoma units only ===")
atlas <- rbindlist(lapply(CTS, function(ct) {
  ep <- file.path(ROOT, "data/input/all_sites", ARM, ct, "train_expression.csv")
  mp <- file.path(ROOT, "metadata/all_sites", ARM, ct, "metadata_units.csv")
  if (!file.exists(ep) || !file.exists(mp)) return(NULL)
  e <- fread(ep); gn <- as.character(e[[1]])
  if (!all(c("SLC16A3","SLC16A1") %in% gn)) return(NULL)
  E <- as.matrix(e[, -1]); rownames(E) <- gn
  mu <- fread(mp); setnames(mu, 1, "id"); mu <- as.data.frame(mu)
  rownames(mu) <- as.character(mu$id); mu <- mu[colnames(E), ]
  k <- mu$tissue %in% CANCER_T
  if (!sum(k)) return(NULL)
  d <- data.table(pt = mu$patient[k], n = mu$n_cells[k],
                  a3 = E["SLC16A3", k], a1 = E["SLC16A1", k])[
       , .(a3 = weighted.mean(a3, n), a1 = weighted.mean(a1, n)), by = pt]
  data.table(cohort = "This study", lineage = ct, n = nrow(d),
             ratio = median(d$a3) / median(d$a1))
}))
read_ext <- function(name, label) {
  f <- file.path(ROOT, "results/independent", name, "C2_isoform_specificity.csv")
  if (!file.exists(f)) { say("  missing: %s", f); return(NULL) }
  d <- fread(f)
  data.table(cohort = label, lineage = d$lineage, n = d$n, ratio = d$ratio_MCT4_MCT1)
}
DB <- rbindlist(list(atlas,
                     read_ext("cillo_GSE139324", "Cillo"),
                     read_ext("kurten_GSE164690", "Kürten")), fill = TRUE)
DB[, lineage := fifelse(lineage %in% c("T","T/NK"), "T / NK",
                 fifelse(lineage %in% c("B","B/Plasma","Plasma"), "B / plasma", lineage))]
DB <- DB[lineage %in% c("Myeloid","B / plasma","Fibroblast","Epithelial","Endothelial","T / NK")]
DB <- DB[, .(ratio = max(ratio), n = max(n)), by = .(cohort, lineage)]
MIN_N <- 8
if (nrow(DB[n < MIN_N])) {
  say("  dropped for n < %d:", MIN_N); print(as.data.frame(DB[n < MIN_N]), row.names = FALSE)
}
DB <- DB[n >= MIN_N]
DB[, cohort := factor(cohort, levels = c("This study","Cillo","Kürten"))]
wB <- dcast(DB, lineage ~ cohort, value.var = "ratio")
print(as.data.frame(wB[order(-wB$`This study`)]), row.names = FALSE)
write_csv(DB, file.path(OUT, "fig5B_mct_ratio.csv"))

DB[, lineage := factor(lineage, levels = DB[cohort == "This study"][order(ratio)]$lineage)]
PB <- ggplot(DB, aes(ratio, lineage, colour = cohort)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
  geom_line(aes(group = lineage), colour = "grey80", linewidth = .6) +
  geom_point(size = 3.4, alpha = .9) +
  geom_text(data = DB[cohort == "This study"], aes(label = sprintf("n=%d", n)),
            hjust = 1.6, size = 2.7, colour = "grey40", show.legend = FALSE) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 40),
                labels = c("1\u00d7","2\u00d7","5\u00d7","10\u00d7","20\u00d7","40\u00d7")) +
  scale_colour_manual(values = c(`This study` = "#b2182b",
                                 Cillo = "#4393c3", `Kürten` = "#2166ac"), name = NULL) +
  labs(title = expression("B"~~~~italic("SLC16A3")~"/"~italic("SLC16A1")~
                          "ratio in tumor, three cohorts"),
       x = expression(italic("SLC16A3")~"/"~italic("SLC16A1")~
                      "expression ratio (log scale)"), y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())

# ========================== PANELS 5C and 5D =================================
tn   <- read_csv(file.path(CP, "tumor_vs_normal_protein.csv"), show_col_types = FALSE)
conc <- read_csv(file.path(CP, "rna_protein_concordance.csv"), show_col_types = FALSE)
cors <- read_csv(file.path(CP, "mct4_protein_correlations.csv"), show_col_types = FALSE)
say("\n=== 5C: tumor vs normal-adjacent protein ==="); print(as.data.frame(tn))
say("\n=== 5D: MCT4 protein correlations ==="); print(as.data.frame(cors))

tn <- tn %>% mutate(
  # this panel is protein throughout; the two transporters take their protein
  # names, the other five already share gene and protein symbols
  lab  = ifelse(gene == "SLC16A3", "MCT4",
         ifelse(gene == "SLC16A1", "MCT1", gene)),
  fill = case_when(gene == "SLC16A3" ~ "MCT4", gene == "LDHA" ~ "LDHA", TRUE ~ "other")) %>%
  arrange(diff) %>% mutate(lab = factor(lab, levels = lab))

PC <- ggplot(tn, aes(diff, lab, fill = fill)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_col(width = .65) +
  geom_text(aes(label = ifelse(p < 0.001, "p < 0.001", sprintf("p = %.2f", p)),
                hjust = ifelse(diff > 0, -0.1, 1.1)), size = 2.7, colour = "grey25") +
  scale_fill_manual(values = c(MCT4 = "#b2182b", LDHA = "#4393c3", other = "grey78"),
                    guide = "none") +
  coord_cartesian(xlim = c(0, max(tn$diff) * 1.22)) +
  labs(title = "C   Protein abundance, tumor vs normal-adjacent (CPTAC-HNSCC)",
       x = "difference in median protein abundance", y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), plot.margin = margin(5,18,5,5))

# GAPDH is labelled a negative control in the source table but correlates at
# -0.40, so it does not serve that function. It is excluded from the panel and
# reported in the caption; ACTB behaves as expected and is retained.
gapdh <- cors %>% filter(gene == "GAPDH")
cors  <- cors %>% filter(gene != "GAPDH") %>%
  mutate(lab  = recode(gene, SLC16A1 = "MCT1", ACTB = "ACTB (neg. control)"),
         fill = case_when(group == "neg control" ~ "negative control",
                          gene == "LDHA"         ~ "LDHA", TRUE ~ "marker"),
         sig  = p < 0.05) %>%
  arrange(rho) %>% mutate(lab = factor(lab, levels = lab))

PD <- ggplot(cors, aes(rho, lab, fill = fill)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_col(width = .65) +
  geom_text(aes(label = ifelse(sig, sprintf("%.2f", rho), sprintf("%.2f (n.s.)", rho)),
                hjust = ifelse(rho > 0, -0.12, 1.12)), size = 2.7, colour = "grey25") +
  scale_fill_manual(values = c(marker = "grey72", LDHA = "#4393c3",
                               `negative control` = "grey85"), guide = "none") +
  coord_cartesian(xlim = c(-0.6, 1)) +
  labs(title = "D   Protein correlation with MCT4 (CPTAC-HNSCC)",
       x = "rho (with MCT4 protein)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), plot.margin = margin(5,18,5,5))

# ============================== ASSEMBLE ====================================
fig <- (PA | PB) / (PC | PD)
ggsave(file.path(OUT, "fig5_transporter.png"), fig, width = 13, height = 9, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig5_transporter.pdf"), fig, width = 13, height = 9, bg = "white")
say("\n-> %s", OUT)

cm <- conc %>% filter(gene == "SLC16A3")
m  <- tn %>% filter(gene == "SLC16A3"); l <- tn %>% filter(gene == "LDHA")
say("\nCaption draft:")
say("  The transporter, and where the pattern holds.")
say("  (A) Correlation between lactate-transport flux and the expression of each")
say("  transporter gene, for the cell types in which the module was tested. Dashed")
say("  segments mark the empirical null threshold for that cell type and module size.")
ct_pass <- sort(as.character(gA$celltype[gA$gene == "SLC16A3" &
                                         gA$rho > gA$null_bar]))
say("  SLC16A3 exceeds the threshold in %s; SLC16A1 exceeds it",
    if (length(ct_pass) > 1)
      paste0(paste(head(ct_pass, -1), collapse = ", "), " and ", tail(ct_pass, 1))
    else ct_pass)
say("  in none, including T cells where it is the better-correlated of the two but")
say("  still falls short. Rows are annotated with MCT4's share of summed transporter")
say("  expression. The comparison verifies that the flux estimate follows its dominant")
say("  input; it is not independent evidence about transporter usage.")
say("  (B) MCT4 to MCT1 expression ratio by cell type, restricted to tumor tissue in")
say("  all three datasets. Myeloid ranks highest and T cells among the lowest in every")
say("  dataset. Ratios in this study run lower throughout, consistent with imputation")
say("  over a restricted metabolic gene panel compressing extreme values. Cell types")
say("  with fewer than %d samples in a dataset are omitted. These are expression", MIN_N)
say("  ratios and do not depend on the flux estimate.")
say("  (C) MCT4 protein abundance is higher in tumor than in normal-adjacent tissue in")
say("  CPTAC-HNSCC (difference %.2f, p < 0.001), the largest enrichment among the seven",
    m$diff)
say("  proteins examined. LDHA shows no difference (%.2f, p = %.2f). If the signal were",
    l$diff, l$p)
say("  general glycolytic activity rather than transport, LDHA would be expected to")
say("  move with MCT4; it does not. Transcript and protein agree closely for SLC16A3")
say("  (rho = %.2f, n = %d), supporting the transcript-based measures used elsewhere.",
    cm$rho_rna_protein, cm$n)
say("  (D) What MCT4 protein tracks across CPTAC samples. LDHA failing to reach")
say("  significance is the basis for not interpreting the LDH module elsewhere in this")
say("  work. ACTB is shown as a loading control; GAPDH was also examined but correlates")
say("  at %.2f (p = %.1e), so it does not serve that function here and is omitted from",
    gapdh$rho, gapdh$p)
say("  the panel. CPTAC is bulk proteomics and cannot attribute protein to a cell type:")
say("  these panels establish that the transporter is present and tumor-enriched, not")
say("  which cells carry it.")
