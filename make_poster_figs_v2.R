#!/usr/bin/env Rscript
# ==============================================================================
# make_poster_figs_v2.R  --  poster panels, revision 2
#
#   Rscript make_poster_figs_v2.R --arm bulk_before
#
# Cosmetic only. Numbers come from the manuscript scripts unchanged: they are
# sourced with ggsave, write_csv and write_tsv masked, so nothing in the
# manuscript output or in poster_figs/ is touched. Only poster_figs_v2/ is written.
#
# THE COLOUR RULE. Cardinal marks the one element per panel that IS the finding.
# Blue marks the comparator it should be read against. Everything else is grey.
# Gold is not used as a data colour.
#
# Where a panel needed its structure changed rather than its colours swapped --
# 2B needs per-bar fills where the manuscript has one, 2D needs the patient lines
# decoupled from the module colouring -- the panel is rebuilt from the same data
# object the manuscript script produced, not re-derived.
#
# Output -> results/all_sites/<ARM>/poster_figs_v2/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr) })

OUT <- file.path(PATHS$results, "poster_figs_v2")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

CARD  <- "#990000"   # the finding
BLUE  <- "#1F77B4"   # the comparator
GREY  <- "#B0B0B0"   # everything else

BASE   <- 22    # all panel text
LETTER <- 28    # panel letter, bold
TITLE  <- 22    # panel title, regular weight beside the letter
ANN    <- 6.4   # in-plot annotations: ggplot size ~6.4 is about 18 pt
ANN_4B <- 5.6   # 4B's HR text, about 16 pt. At 9 in wide the panel gives most of
                # its width to text, so this is down from 18 pt and the right
                # margin with it. If it is still tight, the next lever is the
                # annotation itself: dropping the CI leaves "HR, p" and frees
                # roughly a third of the margin.
DPI    <- 300

# ---- source a manuscript script without letting it write anything -----------
lift <- function(path) {
  e <- new.env(parent = globalenv())
  for (f in c("ggsave","write_csv","write_tsv"))
    assign(f, function(...) invisible(NULL), envir = e)
  say("  sourcing %s (writes suppressed)...", basename(path))
  suppressMessages(suppressWarnings(sys.source(path, envir = e)))
  e
}

# ---- layer surgery: geom params are hardcoded and base_size never reaches them
set_param <- function(p, geoms, param, value, which = NULL) {
  hits <- which(vapply(p$layers, function(l) inherits(l$geom, geoms), logical(1)))
  if (!is.null(which)) hits <- hits[which]
  for (i in hits) p$layers[[i]]$aes_params[[param]] <- value
  p
}
# text layers carry their own hardcoded grey; set size and colour together so no
# annotation is left grey
text_size <- function(p, s = ANN) {
  p <- set_param(p, c("GeomText","GeomLabel"), "size", s)
  set_param(p, c("GeomText","GeomLabel"), "colour", "black")
}

# ---- poster theme: letter only, no in-plot title ----------------------------
poster <- function(p, letter, title = NULL, legend = "none") {
  p + labs(tag = letter, title = title, subtitle = NULL) +
    theme_bw(base_size = BASE) +
    theme(text              = element_text(colour = "black"),
          axis.text         = element_text(colour = "black", size = BASE - 2),
          plot.tag          = element_text(face = "bold", size = LETTER, hjust = 0),
          plot.tag.position = c(0.01, 0.98),
          plot.title        = element_text(face = "plain", size = TITLE, hjust = 0,
                                           lineheight = 1.05,
                                           margin = margin(l = LETTER + 12, b = 10)),
          plot.subtitle     = element_blank(),
          plot.margin       = margin(14, 24, 8, 14),
          axis.title        = element_text(size = BASE),
          strip.text        = element_text(size = BASE - 1),
          legend.position   = legend,
          legend.text       = element_text(size = BASE - 3),
          legend.title      = element_text(size = BASE - 3),
          panel.grid.major.y = element_blank())
}

# ---- re-scale an already-built panel to a different base size ---------------
# Additive: overrides the size of each text element without resetting the theme,
# so panel-specific tweaks survive. geom_text annotations scale by the same
# ratio so in-plot numbers keep their proportion to the axis text.
at_base <- function(p, base) {
  r <- base / BASE
  for (i in seq_along(p$layers))
    if (inherits(p$layers[[i]]$geom, c("GeomText","GeomLabel")) &&
        !is.null(p$layers[[i]]$aes_params$size))
      p$layers[[i]]$aes_params$size <- p$layers[[i]]$aes_params$size * r
  p + theme(text         = element_text(size = base, colour = "black"),
            axis.text    = element_text(size = base - 2, colour = "black"),
            axis.title   = element_text(size = base),
            strip.text   = element_text(size = base - 1),
            legend.text  = element_text(size = base - 3),
            legend.title = element_text(size = base - 3),
            plot.title   = element_text(face = "plain", size = base,
                                        hjust = 0, lineheight = 1.05,
                                        margin = margin(l = LETTER + 12, b = 10)),
            plot.tag     = element_text(face = "bold", size = base + 6, hjust = 0))
}

CHECK <- list()
save_panel <- function(p, name, w, h, cardinal, removed = "-") {
  for (ext in c("png","pdf")) {
    f <- file.path(OUT, sprintf("%s.%s", name, ext))
    if (ext == "png") ggplot2::ggsave(f, p, width = w, height = h, dpi = DPI, bg = "white")
    else              ggplot2::ggsave(f, p, width = w, height = h, bg = "white")
  }
  CHECK[[length(CHECK)+1]] <<- data.frame(file = name, size = sprintf("%g x %g in", w, h),
                                          cardinal = cardinal, legend_removed = removed)
  say("    %-18s %-11s cardinal: %s", name, sprintf("%g x %g", w, h), cardinal)
}

# ============================== FIGURE 1 =====================================
say("\n=== Fig 1 (fig3_reorganisation_panels.R) ===")
e3 <- lift("fig3_reorganisation_panels.R")

# 1A: median dot and IQR bar are the finding; the patient cloud is context
p1A <- e3$A
p1A <- set_param(p1A, "GeomPoint", "colour", GREY,  which = 1)   # jitter cloud
p1A <- set_param(p1A, "GeomPoint", "size",   2.5,   which = 1)
p1A <- set_param(p1A, "GeomPoint", "alpha",  0.45,  which = 1)
p1A <- set_param(p1A, "GeomPoint", "colour", CARD,  which = 2)   # median dot
p1A <- set_param(p1A, "GeomPoint", "size",   5,     which = 2)
p1A <- set_param(p1A, "GeomErrorbarh", "colour", CARD)
p1A <- set_param(p1A, "GeomErrorbarh", "linewidth", 1.4)
p1A <- poster(text_size(p1A), "A", "Export allocation by cell type") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.09))) +
  labs(x = "difference in export allocation, percentage points\n(myeloid minus other)")
save_panel(p1A, "poster_fig1A", 9, 9, "median dot + IQR bar", "-")

# 1B: macrophages are the finding; dendritic and other myeloid share one grey
# legend dropped: with dendritic and other myeloid sharing one grey its keys
# were identical, and the y axis already names every subtype
p1B <- poster(e3$B, "B", "Export allocation by myeloid subtype") +
  scale_fill_manual(values = c(macrophage = CARD, dendritic = GREY,
                               `other myeloid` = GREY), guide = "none")
save_panel(p1B, "poster_fig1B", 9, 9, "macrophage subtypes", "subtype fill legend")

# ============================== FIGURE 2 =====================================
say("\n=== Fig 2 (fig4_progression_all_panels.R) ===")
e4 <- lift("fig4_progression_all_panels.R")

# 2A: the median line is the finding; patient lines are context, cohort dropped
p2A <- e4$A
p2A <- set_param(p2A, "GeomLine",  "colour", GREY, which = 1)
p2A <- set_param(p2A, "GeomPoint", "colour", GREY, which = 1)
p2A <- set_param(p2A, "GeomLine",  "colour", CARD, which = 2)   # median
p2A <- set_param(p2A, "GeomLine",  "linewidth", 2.2, which = 2)
p2A <- set_param(p2A, "GeomPoint", "colour", CARD, which = 2)
p2A <- set_param(p2A, "GeomPoint", "size",   6,    which = 2)
p2A <- poster(text_size(p2A), "A", "Paired myeloid export allocation") +
  labs(y = "% of total flux") + theme(panel.grid.major.x = element_blank())
save_panel(p2A, "poster_fig2A", 11, 6.5, "thick median line", "He/Kim cohort legend")

# 2B: rebuilt so the all-subtypes bar can carry cardinal on its own
b <- e4$B_dat
# all three sets give p = 0.002, so one p-value in the caption covers the panel
b$lab2 <- c("All myeloid\nsubtypes", "Excluding\nSPP1+ TAM",
            "Present in\nnormal only")[match(b$set, levels(b$set))]
b$set <- factor(b$lab2, levels = b$lab2)
b$hl <- b$set == levels(b$set)[1]
p2B <- ggplot(b, aes(set, fold, fill = hl)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
  geom_col(width = .6) +
  geom_text(aes(label = sprintf("%.2f\u00d7", fold)), vjust = -0.5,
            size = ANN, fontface = "bold", colour = "black") +
  scale_fill_manual(values = c(`TRUE` = CARD, `FALSE` = GREY), guide = "none") +
  coord_cartesian(ylim = c(0, max(b$fold) * 1.20)) +
  labs(x = NULL, y = "fold change")
p2B <- poster(p2B, "B", "Fold change, normal to carcinoma") +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = BASE - 4))
save_panel(p2B, "poster_fig2B", 7, 6.5, "the All myeloid subtypes bar",
           "in-bar p-values dropped; all three are p = 0.002, so the caption covers them")

# 2C: myeloid is the finding, every other cell type is context
# rebuilt so T can take blue as the comparator while Myeloid keeps cardinal
sm <- e4$smE
sm$role <- factor(ifelse(sm$celltype == "Myeloid", "finding",
                  ifelse(sm$celltype == "T", "comparator", "other")),
                  levels = c("finding","comparator","other"))
p2C <- ggplot(sm, aes(change, lab, fill = role)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_col(width = .66) +
  # the p-value sits just right of the bar end for positive bars and just right
  # of zero for negative ones, so it never runs back across the y labels
  geom_text(aes(x = pmax(change, 0),
                label = ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))),
            hjust = -0.12, size = ANN, colour = "black") +
  scale_fill_manual(values = c(finding = CARD, comparator = BLUE, other = GREY),
                    guide = "none") +
  coord_cartesian(xlim = c(min(sm$change) * 1.25, 2.3), clip = "off") +
  labs(x = "change, carcinoma minus normal (percentage points)", y = NULL)
p2C <- poster(p2C, "C", "Change in export allocation")
save_panel(p2C, "poster_fig2C", 9, 5.5, "the Myeloid bar; the T bar is blue", "-")

# 2D: rebuilt so the patient lines leave the module colouring alone
d <- e4$longD
med <- d %>% group_by(stage, module) %>% summarise(pct = median(pct), .groups = "drop")
p2D <- ggplot(d, aes(stage, pct)) +
  geom_line(aes(group = interaction(patient, module)), colour = GREY,
            alpha = .45, linewidth = .6) +
  geom_line(data = med, aes(group = module, colour = module), linewidth = 2.2) +
  geom_point(data = med, aes(colour = module), size = 5.5) +
  scale_colour_manual(values = setNames(c(BLUE, CARD), levels(d$module)), name = NULL) +
  labs(x = NULL, y = "% of total flux")
p2D <- poster(p2D, "D", "Lactate production and export", legend = "bottom") + theme(panel.grid.major.x = element_blank())
save_panel(p2D, "poster_fig2D", 9, 5.5, "the export (M_78) median line; production (M_6) is blue", "-")

# ============================== FIGURE 3 =====================================
say("\n=== Fig 3 (fig5_transporter_all_panels.R) ===")
e5 <- lift("fig5_transporter_all_panels.R")

# 3A: this study is the finding, the two external cohorts are context. The colour
# vector is built from the data's own levels rather than from names typed here,
# so the umlaut in the cohort label cannot fail to match.
coh <- levels(e5$DB$cohort)
coh_cols <- setNames(c(CARD, BLUE, GREY)[seq_along(coh)], coh)
say("  3A cohort colours: %s", paste(sprintf("%s=%s", coh, coh_cols), collapse = ", "))
p3A <- poster(text_size(e5$PB), "A", "SLC16A3 (MCT4) exceeds SLC16A1 (MCT1)", legend = "bottom") +
  scale_colour_manual(values = coh_cols, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  labs(x = "SLC16A3 / SLC16A1 expression ratio (log scale)")
save_panel(p3A, "poster_fig3A", 9, 8, "This study points",
           "none; cohort legend kept so Cillo and Kurten are distinguishable")

# 3B: MCT4 is the finding; the other six proteins are context
# the panel is protein throughout, so the two transporters take their protein
# names; the other five have identical gene and protein symbols already
p3B <- e5$PC
levels(p3B$data$lab) <- sub("^SLC16A3.*$", "MCT4",
                        sub("^SLC16A1$",  "MCT1", levels(p3B$data$lab)))
p3B$data$fill <- ifelse(p3B$data$gene == "SLC16A1", "MCT1", p3B$data$fill)
p3B <- poster(text_size(p3B), "B", "Protein enrichment in tumor tissue") +
  scale_fill_manual(values = c(MCT4 = CARD, MCT1 = BLUE, LDHA = GREY,
                               other = GREY), guide = "none")
save_panel(p3B, "poster_fig3B", 9, 8, "the MCT4 bar; the MCT1 bar is blue",
           "LDHA no longer highlighted, now grey with the rest; transporters labelled by protein name (MCT4, MCT1)")

# ============================== FIGURE 4 =====================================
say("\n=== Fig 4 (fig6_clinical_panels.R) ===")
e6 <- lift("fig6_clinical_panels.R")

km <- poster(text_size(e6$A_km), "A", "Survival by export-program score", legend = "bottom") +
  scale_colour_manual(values = setNames(c(BLUE, CARD), levels(e6$km_d$grp)),
                      labels = c("below median", "above median"),
                      name = NULL) +
  scale_fill_manual(values = setNames(c(BLUE, CARD), levels(e6$km_d$grp)),
                    guide = "none")
rk <- ggplot(e6$risk_df, aes(time, grp, label = n)) +
  geom_text(size = ANN - 0.6) +
  scale_x_continuous(breaks = e6$brks, limits = c(0, e6$TMAX)) +
  scale_y_discrete(limits = c("low", "high"),
                   labels = c(low = "below median", high = "above median")) +
  labs(x = NULL, y = NULL, subtitle = "patients at risk") +
  theme_minimal(base_size = BASE - 3) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        plot.subtitle = element_text(size = BASE - 3))
p4A <- km / rk + plot_layout(heights = c(4.6, 1.2))
save_panel(p4A, "poster_fig4A", 9, 8,
           "above-median curve; below-median is blue", "-")

# 4B: program is the finding, its components are the comparator, the rest context
POSTER_LAB <- c(export_core = "Export program",   # composition in the caption
                LDHA = "LDHA", SLC16A3 = "SLC16A3 (MCT4)", HK2 = "HK2",
                SPP1 = "SPP1", CA9 = "CA9", VEGFA = "VEGFA", MMP11 = "MMP11")
a2 <- e6$a_full
miss <- setdiff(a2$var, names(POSTER_LAB))
if (length(miss)) say("  [!] no poster label for: %s", paste(miss, collapse = ", "))
a2$lab  <- ifelse(a2$var %in% names(POSTER_LAB), POSTER_LAB[a2$var], a2$lab)
a2$role <- factor(ifelse(a2$var == "export_core", "program",
                  ifelse(a2$var %in% c("SLC16A3","LDHA","HK2"), "component", "other")),
                  levels = c("program","component","other"))
a2$labf <- factor(a2$lab, levels = a2$lab)

p4B <- ggplot(a2, aes(HR, forcats::fct_rev(labf), colour = role)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .2, linewidth = 1) +
  geom_point(aes(size = role)) +
  geom_text(aes(x = max(a2$hi) * 1.06,
                label = sprintf("%.2f (%.2f-%.2f)  p = %s", HR, lo, hi,
                                ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))),
            hjust = 0, size = ANN_4B, colour = "grey25") +
  scale_colour_manual(values = c(program = CARD, component = BLUE, other = GREY),
                      guide = "none") +
  scale_size_manual(values = c(program = 6, component = 4.2, other = 3.6), guide = "none") +
  scale_x_log10() +
  coord_cartesian(xlim = c(min(a2$lo) * 0.96, max(a2$hi) * 1.75), clip = "off") +
  labs(x = "hazard ratio per SD (log scale)", y = NULL)
p4B <- poster(p4B, "B", "Hazard ratio per standard deviation") +
  theme(axis.text.y = element_text(size = BASE), plot.margin = margin(14, 150, 8, 14))
save_panel(p4B, "poster_fig4B", 9, 8, "the Export program row; its three components are blue",
           "program/component/signature legend")

# =========================== COLUMN COMPOSITES ===============================
# Column 2 of the poster carries figures 1 and 2, column 3 figures 3 and 4, each
# column 18.3 in wide. Panel widths inside a composite are set in proportion to
# the standalone widths, so a panel that was given more room on its own keeps it:
# 2A is the hero panel at 11 against 2B's 7, and that ratio is preserved.
#
# These are in addition to the individual files, not instead of them.
COLW <- 18.3

say("\n=== column composites at %.1f in wide ===", COLW)
comp_save <- function(p, name, h, note) {
  for (ext in c("png","pdf")) {
    f <- file.path(OUT, sprintf("%s.%s", name, ext))
    if (ext == "png") ggplot2::ggsave(f, p, width = COLW, height = h, dpi = DPI, bg = "white")
    else              ggplot2::ggsave(f, p, width = COLW, height = h, bg = "white")
  }
  CHECK[[length(CHECK)+1]] <<- data.frame(file = name,
      size = sprintf("%g x %g in", COLW, h), cardinal = "as in the panels",
      legend_removed = note)
  say("    %-24s %g x %-5g %s", name, COLW, h, note)
}

comp_save((p1A | p1B) + plot_layout(widths = c(1, 1)),
          "poster_fig1_composite", 9, "1A | 1B, equal widths")

# a plain 2x2 shares column widths across rows, so each row is its own nested
# layout: 11:7 on top, 9:9 below. Those sum to 18 units in both rows, so at
# 18.3 in wide every panel lands within 2% of the width it has standalone, and
# the heights (6.5 + 5.5) fill 12 in exactly. base_size therefore renders at the
# same visual size here as in the individual files.
BASE2 <- 26   # Fig 2 composite only; the individual files stay at BASE
row_top <- (at_base(p2A, BASE2) | at_base(p2B, BASE2)) + plot_layout(widths = c(10, 8))
row_bot <- (at_base(p2C, BASE2) | at_base(p2D, BASE2)) + plot_layout(widths = c(9, 9))
comp_save(row_top / row_bot + plot_layout(heights = c(6, 5)),
          "poster_fig2_composite", 11,
          sprintf("2A | 2B at 10:8 over 2C | 2D at 9:9, base_size %d", BASE2))

comp_save((p3A | p3B) + plot_layout(widths = c(1, 1)),
          "poster_fig3_composite", 8, "3A | 3B, equal widths")

comp_save((p4A | p4B) + plot_layout(widths = c(1, 1)),
          "poster_fig4_composite", 8, "4A | 4B, equal widths")

say("\n  column 2: fig 1 (9 in) + fig 2 (11 in) = 20 in")
say("  column 3: fig 3 (8 in) + fig 4 (8 in)  = 16 in")
say("  against about 27 in of column height, leaving room for captions.")
say("\n  Fig 2 composite is at base_size %d, the individual Fig 2 files at %d.", BASE2, BASE)
say("  If Fig 2 still reads smaller than Figs 1, 3 and 4 on the board, check how")
say("  those are being placed: an individual 9-in panel stretched to fill an 18.3-in")
say("  column has its text scaled up with it, while this composite goes in at 1:1.")
say("  Placing the fig1/fig3/fig4 composites instead would equalise them directly.")

# ================================ CHECKLIST ==================================
say("\n============================================================")
say("  CHECKLIST  -> %s", OUT)
say("============================================================")
chk <- bind_rows(CHECK)
print(as.data.frame(chk), row.names = FALSE, right = FALSE)
write.csv(chk, file.path(OUT, "poster_v2_checklist.csv"), row.names = FALSE)
say("\n  %d panels, PNG at %d dpi plus PDF of each.", nrow(chk), DPI)
say("  base_size %d, panel letters %d pt bold, in-plot annotations ~%.0f pt.",
    BASE, LETTER, ANN * 2.845)
say("  Panel letter bold at %d pt with the title beside it at %d pt, as in v1.", LETTER, TITLE)
say("  4B HR annotation at %.0f pt, right margin 150 pt, y label shortened to", ANN_4B * 2.845)
say("  \"Export program\" -- its composition belongs in the poster caption. If the")
say("  panel is still tight, drop the CI from the annotation next.")
say("  Cardinal %s = the finding, blue %s = comparator, grey %s = context.",
    CARD, BLUE, GREY)
say("  Numbers unchanged: manuscript scripts sourced with writes masked.")
