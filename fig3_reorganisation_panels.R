#!/usr/bin/env Rscript
# ==============================================================================
# fig3_reorganisation_panels.R  --  Figure 3: "what the reorganisation is"
#
#   Rscript fig3_reorganisation_panels.R --arm bulk_before
#
# Self-contained: builds all three panels rather than borrowing from 04f/04h, so
# the figure has one source and one set of conventions. (04f's third panel, the
# transporter comparison, now lives in fig5a and is not part of this figure.)
#
#   3A  myeloid versus each other cell type, paired within patient, carcinoma only.
#       Patients differ substantially in overall export allocation, so cell types
#       are compared inside each patient rather than pooled. The annotation is
#       the percentage of patients in which myeloid is higher, and the full range
#       is shown: NK at 72% is the weakest comparison and belongs on the plot.
#
#   3B  export allocation by myeloid subtype in CARCINOMA. This figure is about the
#       structure of the reorganisation, not its timing, so only carcinoma is
#       shown here; the same panel across tissue classes belongs to Figure 4,
#       where change over tissue is the subject.
#
#   3C  where lactate export sits among all modules ranked by myeloid bias. This
#       keeps 3A honest: export is elevated in myeloid cells but is not the most
#       myeloid-biased process, and the distribution is centered near zero, so the
#       elevation is not a global inflation of myeloid shares.
#
# Output -> results/all_sites/<ARM>/fig3/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

OUT <- file.path(PATHS$results, "fig3")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
usable <- readRDS(file.path(PATHS$results, "03_qc", "usable_modules.rds"))
MODULE <- "M_78"

NL_T  <- c("NL","Normal"); ADJ_T <- "NL-contact"
PRE_T <- c("Pre","LP");    CA_T  <- c("CA","CA-E","CA-A","CA-c","CA-p")

cw <- as.data.frame(readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))) %>%
  filter(n_cells >= GATES$min_cells)
mc <- intersect(mcols_of(cw), usable)
S  <- as.matrix(cw[, mc]); S <- S / rowSums(S)
cw$export <- S[, MODULE]

# ============================ PANEL 3A ======================================
ca <- cw %>% filter(tissue %in% CA_T)
wide <- ca %>% group_by(patient, annot_4th) %>%
  summarise(e = mean(export), .groups = "drop") %>%
  pivot_wider(names_from = annot_4th, values_from = e)
others <- setdiff(names(wide), c("patient","Myeloid"))

pa <- lapply(others, function(ct) {
  d <- wide %>% filter(!is.na(Myeloid), !is.na(.data[[ct]]))
  if (nrow(d) < 8) return(NULL)
  tibble(other = ct, patient = d$patient, diff = 100 * (d$Myeloid - d[[ct]]))
}) %>% bind_rows()

stats <- pa %>% group_by(other) %>%
  summarise(n = n(), pct_pos = 100*mean(diff > 0), med = median(diff),
            p = suppressWarnings(wilcox.test(diff)$p.value), .groups = "drop") %>%
  arrange(desc(med))
say("=== 3A: myeloid minus each cell type, paired, carcinoma ===")
print(as.data.frame(stats %>% mutate(med = signif(med,3), pct_pos = round(pct_pos),
                                     p = signif(p,2))), row.names = FALSE)
say("  range of %% of patients myeloid-higher: %.0f%% to %.0f%%",
    min(stats$pct_pos), max(stats$pct_pos))

pa$other  <- factor(pa$other, levels = stats$other)
pa_sum <- pa %>% group_by(other) %>%
  summarise(med = median(diff), lo = quantile(diff,.25), hi = quantile(diff,.75), .groups="drop") %>%
  left_join(stats %>% select(other, pct_pos), by = "other")

A <- ggplot(pa, aes(diff, other)) +
  geom_vline(xintercept = 0, colour = "grey45", linetype = 2) +
  geom_jitter(height = .16, alpha = .3, size = .9, colour = "grey45") +
  geom_errorbarh(data = pa_sum, inherit.aes = FALSE,
                 aes(y = other, xmin = lo, xmax = hi), height = 0,
                 colour = "#b2182b", linewidth = .9) +
  geom_point(data = pa_sum, inherit.aes = FALSE, aes(x = med, y = other),
             size = 3, colour = "#b2182b") +
  geom_text(data = pa_sum, inherit.aes = FALSE,
            aes(x = max(pa$diff, na.rm = TRUE), y = other,
                label = sprintf("%.0f%%", pct_pos)),
            hjust = 1, size = 3, colour = "grey30") +
  labs(title = "A   Lactate-export allocation, myeloid versus each cell type",
       x = "difference in export allocation, percentage points (myeloid minus other)", y = NULL) +
  theme_bw(base_size = 10)

# ============================ PANEL 3C ======================================
pw  <- cw %>% filter(tissue %in% CA_T, annot_4th %in% c("Myeloid","T"))
Sp  <- S[cw$tissue %in% CA_T & cw$annot_4th %in% c("Myeloid","T"), , drop = FALSE]
grp <- pw$annot_4th
adv <- vapply(colnames(Sp), function(m) {
  a <- median(Sp[grp=="Myeloid", m]); b <- median(Sp[grp=="T", m])
  if (b <= 0) return(NA_real_); log2(a/b) }, numeric(1))
adv <- adv[is.finite(adv)]
pct <- 100 * mean(adv < adv[[MODULE]])
say("\n=== 3C: M_78 at the %.0fth percentile of %d modules by myeloid bias ===",
    pct, length(adv))

C <- ggplot(tibble(l2 = adv), aes(l2)) +
  geom_histogram(bins = 45, fill = "grey82", colour = "white") +
  geom_vline(xintercept = 0, colour = "grey45", linetype = 2) +
  geom_vline(xintercept = adv[[MODULE]], colour = "#b2182b", linewidth = 1) +
  annotate("text", x = adv[[MODULE]], y = Inf, label = "  lactate export",
           colour = "#b2182b", hjust = 0, vjust = 1.8, size = 3) +
  labs(title = "C   Myeloid-to-T-cell allocation ratio across all modules",
       x = "log2 (myeloid allocation / T-cell allocation)", y = "modules") +
  theme_bw(base_size = 10)

# ============================ PANEL 3B ======================================
mye <- cw %>% filter(annot_4th == "Myeloid", tissue %in% CA_T) %>%
  mutate(grp = case_when(grepl("^MQ_", finaltype) ~ "macrophage",
                         grepl("^DC_", finaltype) ~ "dendritic",
                         TRUE ~ "other myeloid"))
ord <- mye %>% group_by(finaltype) %>% summarise(m = median(export), .groups="drop") %>%
  arrange(m) %>% pull(finaltype)
mye$finaltype <- factor(mye$finaltype, levels = ord)
say("\n=== 3B: %d carcinoma myeloid units across %d subtypes ===",
    nrow(mye), n_distinct(mye$finaltype))
print(as.data.frame(mye %>% group_by(finaltype) %>%
  summarise(n = n(), median_pct = round(100*median(export),2), .groups="drop") %>%
  arrange(desc(median_pct))), row.names = FALSE)

B <- ggplot(mye, aes(finaltype, 100*export, fill = grp)) +
  geom_boxplot(outlier.size = .3, linewidth = .3) +
  coord_flip() +
  scale_fill_manual(values = c(macrophage = "#2166ac", dendritic = "#b2182b",
                               `other myeloid` = "grey78"), name = NULL) +
  labs(title = "B   Lactate-export allocation by myeloid subtype in carcinoma",
       x = NULL, y = "% of total metabolic flux") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")

fig <- (A | B) / C + plot_layout(heights = c(1.25, 1))
ggsave(file.path(OUT, "fig3_reorganisation.png"), fig, width = 12, height = 9, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig3_reorganisation.pdf"), fig, width = 12, height = 9, bg = "white")
write_csv(stats, file.path(OUT, "fig3A_paired_stats.csv"))
say("\n-> %s", OUT)

say("\nCaption draft:")
say("  Lactate-export allocation across cell types and myeloid subtypes.")
say("  (A) Myeloid cells allocate more flux to lactate export than every other cell")
say("  type, compared within each patient in carcinoma tissue. Each point is one")
say("  patient; the annotation gives the percentage of patients in which myeloid is")
say("  higher (%.0f-%.0f%%). The difference is small in size but holds in nearly every",
    min(stats$pct_pos), max(stats$pct_pos))
say("  patient.")
say("  (B) Allocation by myeloid subtype in carcinoma, across %d units. SPP1+", nrow(mye))
say("  macrophages are the most export-committed subtype; cDC1 and pDC the least.")
say("  (C) Ratio of myeloid to T-cell allocation, shown for the %d of 162 usable", length(adv))
say("  modules for which the ratio is defined. The remainder are excluded because")
say("  median T-cell allocation is zero or, at values near zero, slightly negative.")
say("  Lactate export sits at the %.0fth percentile, so it is elevated without being", pct)
say("  the most myeloid-biased process. The distribution is centered near zero, so")
say("  myeloid allocation is not inflated across the board.")
say("  Units are patient x tissue x cell subtype pseudobulks of at least 10 cells.")
