#!/usr/bin/env Rscript
# ==============================================================================
# fig6_clinical_panels.R  --  Figure 6: clinical association in TCGA-HNSC
#
#   Rscript fig6_clinical_panels.R
#
# Standalone: reads the saved Cox outputs, does not source config_all_sites.R.
#
# Two versions are written:
#   fig6_full.png        panel A shows the program, then all three of its
#                        components, then the other genes tested.
#   fig6_restricted.png  panel A shows the program alone.
# Both share panel B, the myeloid adjustment.
#
# On the components. The program is SLC16A3 + LDHA + HK2 (EXPORT_CORE in
# 30_tcga_signature.R). HK2 was originally absent from the Cox loop, so an
# earlier version of this panel showed two of three components alongside two
# genes that are not in the program. 30_tcga_signature.R was amended to fit HK2
# and re-run. SPP1 and MMP11 are kept in a separate block: they are niche markers
# from the fuller signature, not part of the export program.
#
# On the SLC16A3 null. It does not reach significance alone (HR 1.14, p = 0.10),
# is unchanged by myeloid adjustment, and does not improve a model containing
# SPP1 (LRT p = 0.17). Two things bear on that. Bulk expression averages over all
# compartments, so a myeloid-confined program is diluted in proportion to the
# myeloid fraction. And the study has limited power at that effect size: the
# power at each observed HR is computed below and reported on the panel, rather
# than quoting the 80% bound, which would sit above LDHA's significant 1.20 and
# read as a contradiction.
#
# On LDHA carrying more of the signal. This does not conflict with excluding the
# LDH flux module elsewhere. That exclusion concerns whether scFEA's estimate of
# LDH flux is trustworthy; it is not a claim about LDHA expression and outcome.
#
# Output -> results/tcga/fig6/
# ==============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr)
  library(patchwork); library(forcats); library(survival)
})

TC  <- "/Users/Tommy/dechenli_lab_local/scfea26/results/tcga"
OUT <- file.path(TC, "fig6"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

uni <- read_csv(file.path(TC, "cox_univariate.csv"), show_col_types = FALSE)
adj <- read_csv(file.path(TC, "deconvolution/cox_adjusted.csv"), show_col_types = FALSE)
pat <- read_csv(file.path(TC, "per_patient_scores.csv"), show_col_types = FALSE)
cat("=== univariate Cox ===\n"); print(as.data.frame(uni))
cat("\n=== myeloid-adjusted ===\n"); print(as.data.frame(adj))

N        <- unique(uni$n)
N_EVENTS <- 184   # events among the 417 complete cases, from incremental_cox.txt

# Power for a continuous predictor scaled to SD, alpha 0.05 two-sided:
#   z = |log(HR)| * sqrt(events);  power = P(reject | that shift)
pow_at <- function(hr) pnorm(abs(log(hr)) * sqrt(N_EVENTS) - qnorm(0.975))
MIN_DETECTABLE <- exp((qnorm(0.975) + qnorm(0.80)) / sqrt(N_EVENTS))
POWER_AT_SLC   <- pow_at(uni$HR[uni$var == "SLC16A3"])
POWER_AT_LDHA  <- pow_at(uni$HR[uni$var == "LDHA"])
cat(sprintf("\npower with %d events: %.0f%% at SLC16A3 (HR %.2f), %.0f%% at LDHA (HR %.2f);\n  80%% power would require HR >= %.2f\n",
    N_EVENTS, 100*POWER_AT_SLC, uni$HR[uni$var=="SLC16A3"],
    100*POWER_AT_LDHA, uni$HR[uni$var=="LDHA"], MIN_DETECTABLE))

PROGRAM <- c("SLC16A3", "LDHA", "HK2")
missing <- setdiff(PROGRAM, uni$var)
if (length(missing))
  cat(sprintf("\n[!] program component(s) missing from cox_univariate.csv: %s\n",
              paste(missing, collapse = ", ")))

# sig_full = SPP1, MMP11, SLC16A3, LDHA, HK2, CA9. Panel B shows the three
# program components plus whichever of the remaining sig_full genes were fitted.
# CA9 was not in the original Cox loop; if it is still absent the panel shows two
# of three, which invites the same question HK2 did.
# sig_full = EXPORT_CORE + NICHE_ID + HYPOXIA, i.e. seven genes
# (30_tcga_signature.R line 42). The header comment on line 23 of that script
# lists only six and omits VEGFA; the code is correct, the comment is stale.
SIGFULL_EXTRA <- c("SPP1", "MMP11", "CA9", "VEGFA")
missing_sf <- setdiff(SIGFULL_EXTRA, uni$var)
if (length(missing_sf))
  cat(sprintf("[!] sig_full gene(s) missing from cox_univariate.csv: %s\n    Add to the lapply() in 30_tcga_signature.R and re-run, or state the\n    omission in the legend.\n",
              paste(missing_sf, collapse = ", ")))

# Gene symbols are italicised with ggtext where it is installed; without it the
# same labels render in plain text, so the figure still builds.
MD <- requireNamespace("ggtext", quietly = TRUE)
if (!MD) cat("[i] ggtext not installed; row labels will render without italics\n")
em <- function(g) if (MD) paste0("*", g, "*") else g
LAB <- c(sig_full    = "full signature",
         export_core = paste0("Lactate metabolic program (", em("SLC16A3"), " + ",
                              em("LDHA"), " + ", em("HK2"), ")"),
         SLC16A3     = paste0(em("SLC16A3"), " (MCT4)"),
         LDHA = em("LDHA"), HK2 = em("HK2"), SPP1 = em("SPP1"),
         MMP11 = em("MMP11"), CA9 = em("CA9"), VEGFA = em("VEGFA"))

forest <- function(d, title, sub = NULL) {
  # positions computed here rather than referenced inside aes(), which ggplot
  # discourages and which produced a warning per row
  lab_x  <- max(d$hi) * 1.05
  x_lo   <- min(d$lo) * 0.95
  x_hi   <- max(d$hi) * 2
  d$labf <- factor(d$lab, levels = d$lab)

  ggplot(d, aes(HR, fct_rev(labf))) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
    geom_errorbarh(aes(xmin = lo, xmax = hi, colour = sig), height = .18, linewidth = .7) +
    geom_point(aes(colour = sig, size = is_prog)) +
    geom_text(aes(x = lab_x,
                  label = sprintf("%.2f (%.2f-%.2f)  p = %s", HR, lo, hi,
                                  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))),
              hjust = 0, size = 2.9, colour = "grey25") +
    scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "#b2182b"),
                        drop = FALSE, guide = "none") +
    scale_size_manual(values = c(`FALSE` = 2.2, `TRUE` = 3.4),
                      drop = FALSE, guide = "none") +
    scale_x_log10() +
    coord_cartesian(xlim = c(x_lo, x_hi), clip = "off") +
    labs(title = title, x = "hazard ratio per SD (log scale)", y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.y = element_blank(), plot.margin = margin(5, 95, 5, 5)) +
    (if (MD) theme(axis.text.y = ggtext::element_markdown()) else NULL)
}

mk <- function(d) d %>% mutate(sig     = factor(p < 0.05, levels = c(FALSE, TRUE)),
                               is_prog = factor(var == "export_core", levels = c(FALSE, TRUE)))

# ---- panel A: Kaplan-Meier ---------------------------------------------------
# Median split on the export program. A dichotomised split is shown because a
# survival curve is the form a clinical reader expects; the Cox models in panels
# B and C treat the score as continuous, which is the analysis the paper reports.
# OS.time units vary by TCGA release: days in most, months in some. Infer from
# the maximum rather than assuming, since assuming months when the data is in
# days silently clips the curve.
TMAX <- 5   # x-axis limit in years; follow-up runs longer and is noted on the panel
MAXT <- max(pat$OS.time, na.rm = TRUE)
TIME_DIVISOR <- if (MAXT > 1000) 365.25 else 12
cat(sprintf("\nOS.time max = %.1f -> treating as %s (divisor %.2f)\n",
            MAXT, ifelse(TIME_DIVISOR > 100, "days", "months"), TIME_DIVISOR))

km_d <- pat %>% filter(!is.na(export_core), !is.na(OS), !is.na(OS.time), OS.time > 0) %>%
  mutate(grp = factor(ifelse(export_core > median(export_core), "high", "low"),
                      levels = c("low","high")),
         yrs = OS.time / TIME_DIVISOR)
fit <- survfit(Surv(yrs, OS) ~ grp, data = km_d)
lr  <- survdiff(Surv(yrs, OS) ~ grp, data = km_d)
LR_P <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
cat(sprintf("\nKM: n = %d (%d high, %d low), %d events, log-rank p = %.4f\n",
            nrow(km_d), sum(km_d$grp=="high"), sum(km_d$grp=="low"), sum(km_d$OS), LR_P))
cat(sprintf("  follow-up to %.1f years; axis truncated at %d\n", max(km_d$yrs), TMAX))
cat("  NOTE: 30_tcga_signature.R reports log-rank p = 0.0171 for its own KM, which\n        splits on sig_full (line 283). This panel splits on export_core, the\n        program the paper reports. Only one should appear in the text.\n")

sf_df <- broom::tidy(fit) %>%
  mutate(grp = factor(sub("grp=", "", strata), levels = c("low","high")))
brks  <- 0:TMAX
risk  <- summary(fit, times = brks, extend = TRUE)
risk_df <- tibble(time = risk$time, n = risk$n.risk,
                  grp = factor(sub("grp=", "", as.character(risk$strata)),
                               levels = c("low","high")))

# risk_df is built above so the caption can quote the number still at risk at TMAX
A_km <- ggplot(sf_df, aes(time, estimate, colour = grp)) +
  geom_step(linewidth = .8) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = grp),
              alpha = .12, colour = NA) +
  annotate("text", x = 0.15, y = 0.12, hjust = 0, size = 3.2,
           label = sprintf("log-rank p = %.3f", LR_P)) +
  scale_colour_manual(values = c(low = "#4393c3", high = "#b2182b"),
                      labels = c("below median", "above median"),
                      name = "lactate metabolic program") +
  scale_fill_manual(values = c(low = "#4393c3", high = "#b2182b"), guide = "none") +
  scale_x_continuous(breaks = brks, limits = c(0, TMAX)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "A   Overall survival by lactate metabolic program score",
       x = "years from diagnosis", y = "overall survival") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom",
                                   legend.margin = margin(t = -4))

A_risk <- ggplot(risk_df, aes(time, fct_rev(grp), label = n)) +
  geom_text(size = 3) +
  scale_x_continuous(breaks = brks, limits = c(0, TMAX)) +
  scale_y_discrete(labels = c(high = "above median", low = "below median")) +
  labs(x = NULL, y = NULL, subtitle = "patients at risk") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank())

KM <- wrap_elements(A_km / A_risk + plot_layout(heights = c(5, 1)))

# ---- panel B, full ----------------------------------------------------------
a_full <- uni %>% filter(var != "sig_full") %>%
  mutate(block = case_when(var == "export_core" ~ 1L,
                           var %in% PROGRAM     ~ 2L,
                           TRUE                 ~ 3L)) %>%
  arrange(block, desc(HR)) %>%
  # the pipeline suffixes (component / sig_full only) are dropped; block still
  # orders the rows, and the caption says which genes belong to which set
  mutate(lab = unname(LAB[var])) %>%
  mk()

A_full <- forest(a_full,
  "B   Hazard ratio per standard deviation, program and components",
  sprintf("univariate Cox, TCGA-HNSC, n = %d, adjusted for age and stratified by stage. SLC16A3 alone does not reach\nsignificance, but with %d events a hazard ratio of %.2f would be detected in only about %.0f%% of studies this\nsize, so that result is uninformative rather than negative. A composite of correlated genes carries less\nmeasurement noise than any one of them, so part of its advantage is arithmetic; the program is reported as\nthe unit of association, not as evidence that the genes act jointly.",
          N, N_EVENTS, uni$HR[uni$var=="SLC16A3"], 100*POWER_AT_SLC))

# ---- panel A, restricted ----------------------------------------------------
a_rest <- uni %>% filter(var == "export_core") %>% mutate(lab = LAB[var]) %>% mk()

A_rest <- forest(a_rest,
  "B   Hazard ratio per standard deviation, lactate metabolic program",
  sprintf("univariate Cox, TCGA-HNSC, n = %d, adjusted for age and stratified by stage.\nComponents are reported in the supplement.", N))

# ---- panel B, shared --------------------------------------------------------
b <- adj %>%
  mutate(var = model,
         # italicised with the same em() helper as 6B, so the two panels agree.
         # ggtext renders a line break as <br>, not \n.
         lab = recode(model,
           "SLC16A3 alone"      = paste0(em("SLC16A3"), " alone"),
           "SLC16A3 | +myeloid" = paste0(em("SLC16A3"), ", adjusted for",
                                         if (MD) "<br>" else "\n", "myeloid content"),
           "myeloid | +SLC16A3" = paste0("myeloid content, adjusted",
                                         if (MD) "<br>" else "\n", "for ",
                                         em("SLC16A3")))) %>%
  arrange(desc(HR)) %>% mk()

B <- forest(b,
  "C   Hazard ratio before and after adjustment for myeloid content",
  "SLC16A3's hazard ratio is unchanged when myeloid content is added to the model, and myeloid content carries\nno signal of its own once SLC16A3 is included.")

for (nm in c("full","restricted")) {
  h   <- if (nm == "full") 12 else 10
  fig <- KM / (if (nm == "full") A_full else A_rest) / B +
         plot_layout(heights = c(1.5, if (nm == "full") 1.25 else .6, 1))
  ggsave(file.path(OUT, sprintf("fig6_%s.png", nm)), fig, width = 9.5, height = h,
         dpi = 300, bg = "white")
  ggsave(file.path(OUT, sprintf("fig6_%s.pdf", nm)), fig, width = 9.5, height = h, bg = "white")
}
cat("\n->", OUT, "\n")

ec <- uni %>% filter(var == "export_core"); s3 <- uni %>% filter(var == "SLC16A3")
ld <- uni %>% filter(var == "LDHA");        hk <- uni %>% filter(var == "HK2")
sf <- uni %>% filter(var == "sig_full")

cat(sprintf(
"\nCaption draft:\n  Clinical association in TCGA-HNSC.\n  (A) Overall survival split at the median of the lactate metabolic program score (log-rank\n  p = %.3f), with numbers at risk beneath. The curve is unadjusted and includes all\n  %d patients with survival data (%d events); the Cox models in B and C are\n  complete-case after age and stage (n = %d, %d events), which is why the counts\n  differ. The split is for display only; hazard ratios come from the continuous\n  models. Follow-up extends to %.1f years and the axis is truncated at %d.\n\n  (B) Hazard ratio per standard deviation, univariate Cox, adjusted for age and\n  stratified by stage. The lactate metabolic program (SLC16A3, LDHA, HK2) gives HR %.2f\n  (95%% CI %.2f-%.2f, p = %.3f). Its components are LDHA %.2f (p = %.3f), SLC16A3\n  %.2f (p = %.2f) and HK2 %.2f (p = %.2f). SPP1, MMP11, CA9 and VEGFA belong to the\n  fuller signature rather than the export program and are labelled separately; none\n  reaches significance alone. With %d events, a hazard ratio of %.2f would be\n  detected in about %.0f%% of studies this size, so the SLC16A3 result is\n  uninformative rather than negative.\n\n  (C) SLC16A3's hazard ratio before and after adjustment for myeloid content\n  (%.2f to %.2f), with myeloid content itself carrying no independent signal once\n  SLC16A3 is included (HR %.2f, p = %.2f).\n",
  LR_P, nrow(km_d), sum(km_d$OS), N, N_EVENTS, max(km_d$yrs), TMAX,
  ec$HR, ec$lo, ec$hi, ec$p,
  ld$HR, ld$p, s3$HR, s3$p, hk$HR, hk$p,
  N_EVENTS, s3$HR, 100*POWER_AT_SLC,
  adj$HR[adj$model=="SLC16A3 alone"], adj$HR[adj$model=="SLC16A3 | +myeloid"],
  adj$HR[adj$model=="myeloid | +SLC16A3"], adj$p[adj$model=="myeloid | +SLC16A3"]))

cat(sprintf(
"\n[SIGNATURES] Both are defined in 30_tcga_signature.R (lines 39-42):\n  export_core  = SLC16A3, LDHA, HK2                              HR %.2f  p = %.3f\n  sig_full     = SLC16A3, LDHA, HK2 + SPP1, MMP11 + CA9, VEGFA   HR %.2f  p = %.3f\n  Both are hand-specified gene sets rather than sets derived from the flux analysis.\n  The four genes unique to sig_full are individually null: SPP1 %.2f, MMP11 %.2f,\n  CA9 %.2f, VEGFA %.2f. sig_full therefore has the wider interval and weaker p\n  despite containing the same export genes.\n",
  ec$HR, ec$p, sf$HR, sf$p,
  uni$HR[uni$var=='SPP1'], uni$HR[uni$var=='MMP11'],
  uni$HR[uni$var=='CA9'],  uni$HR[uni$var=='VEGFA']))

cat("\n[NOTE] N_EVENTS = 184, confirmed against incremental_cox.txt (n = 417,\n       number of events = 184 in both fitted models). The 216 events reported by\n       30_tcga_signature.R are among the 514 analysable samples before the age\n       and stage filters.\n")

