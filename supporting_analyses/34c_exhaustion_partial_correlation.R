#!/usr/bin/env Rscript
# ==============================================================================
# 34c_exhaustion_partial_correlation.R
#   does myeloid export -> T-cell exhaustion survive adjustment for confounders?
#
#   Rscript 34c_exhaustion_partial_correlation.R --arm bulk_before
#
# THE THREAT. 34 reported a per-sample correlation between myeloid export
# allocation and T-cell exhaustion (rho = 0.24, n = 105). That correlation is
# satisfied TRIVIALLY if a third variable drives both. Three candidates, all
# plausible and all present in these data:
#
#   HYPOXIA  -- MCT4 is HIF-1a-driven, and hypoxia independently promotes T-cell
#               dysfunction. We separately found export tracks bulk hypoxia
#               signatures at rho = 0.60, so this is not hypothetical.
#   STAGE    -- export allocation rises with progression, and exhaustion rises
#               with advanced disease. A stage proxy would produce this
#               correlation with no direct coupling at all.
#   COHORT   -- five cohorts differ in platform, depth and composition.
#
# If the association is only a stage or hypoxia proxy, the cross-compartment
# framing collapses and the claim must be withdrawn. This must be settled before
# the introduction promises a specific coupling.
#
# TESTS
#   A. per-sample hypoxia score in MYELOID cells (HIF-1a targets), so the
#      confounder is measured in the same compartment as the exposure
#   B. partial Spearman: export ~ exhaustion controlling for hypoxia, stage,
#      cohort, and myeloid/T cell counts (proxy for purity/composition)
#   C. mixed model with cohort as random effect, stage and hypoxia as fixed
#   D. WITHIN-STAGE correlation -- if it holds inside carcinoma alone, it is not
#      a stage proxy
#
# DECISION: the claim survives only if the partial correlation stays clearly
# positive AND the within-carcinoma correlation is positive.
#
# Output -> results/tcell_dysfunction/partial/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(ppcor)
})

OUT <- file.path(ROOT, "results/tcell_dysfunction/partial")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

EXHAUSTION <- c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","ENTPD1","TOX")
# HIF-1a targets / hypoxia response, measured IN MYELOID cells
HYPOXIA    <- c("HIF1A","CA9","VEGFA","SLC2A1","PGK1","ALDOA","LDHA","ADM","NDRG1","BNIP3")
NORMAL_T <- c("NL","NL-contact","Normal"); PRE_T <- c("Pre","LP")
CANCER_T <- c("CA","CA-E","CA-A","CA-c","CA-p")

# ---- exposure: myeloid export allocation, per sample -------------------------
cw <- as.data.frame(readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))) %>%
  filter(n_cells >= GATES$min_cells)
mc <- mcols_of(cw); S <- as.matrix(cw[, mc]); S <- S / rowSums(S)
cw$export <- S[, "M_78"]
myl <- cw %>% filter(annot_4th == "Myeloid") %>%
  group_by(patient, tissue, study) %>%
  summarise(export = weighted.mean(export, n_cells), n_myeloid = sum(n_cells),
            .groups = "drop")

# ---- outcome + confounder: from the imputed expression used for grounding ----
read_ct <- function(ct, genes) {
  ep <- file.path(ROOT, "data/input/all_sites/bulk_before", ct, "train_expression.csv")
  mp <- file.path(ROOT, "metadata/all_sites/bulk_before", ct, "metadata_units.csv")
  if (!file.exists(ep)) return(NULL)
  e <- fread(ep); g <- as.character(e[[1]])
  keep <- intersect(genes, g); if (!length(keep)) return(NULL)
  E <- as.matrix(e[, -1]); rownames(E) <- g
  mu <- fread(mp); setnames(mu, 1, "id"); mu <- as.data.frame(mu)
  rownames(mu) <- as.character(mu$id); mu <- mu[colnames(E), ]
  Z <- t(scale(t(E[keep, , drop = FALSE])))
  tibble(patient = mu$patient, tissue = mu$tissue, study = mu$study,
         n = mu$n_cells, score = colMeans(Z, na.rm = TRUE))
}
tex <- read_ct("T", EXHAUSTION) %>% group_by(patient, tissue, study) %>%
  summarise(exhaustion = weighted.mean(score, n), n_tcell = sum(n), .groups="drop")
hyp <- read_ct("Myeloid", HYPOXIA) %>% group_by(patient, tissue, study) %>%
  summarise(hypoxia = weighted.mean(score, n), .groups="drop")

D <- myl %>% inner_join(tex, by = c("patient","tissue","study")) %>%
  left_join(hyp, by = c("patient","tissue","study")) %>%
  mutate(stage = case_when(tissue %in% NORMAL_T ~ 1, tissue %in% PRE_T ~ 2,
                           tissue %in% CANCER_T ~ 3, TRUE ~ NA_real_)) %>%
  filter(!is.na(stage), !is.na(hypoxia))
say("samples with export + exhaustion + hypoxia + stage: %d", nrow(D))
say("  stage distribution: %s",
    paste(names(table(D$stage)), table(D$stage), sep="=", collapse="  "))

# ---- A: how confounded is this? ---------------------------------------------
say("\n=== A: pairwise correlations among the variables ===")
vars <- c("export","exhaustion","hypoxia","stage")
M <- cor(D[, vars], method = "spearman", use = "complete.obs")
print(round(M, 3))
say("  READ: if hypoxia and stage both correlate with export AND exhaustion,")
say("  the raw association is a candidate proxy and must be adjusted.")

# ---- B: partial correlations -------------------------------------------------
say("\n=== B: partial Spearman, export ~ exhaustion ===")
raw <- cor.test(D$export, D$exhaustion, method = "spearman")
say("  unadjusted                      rho = %+.3f  p = %.3g", raw$estimate, raw$p.value)
pc <- function(ctrl, label) {
  d <- D[complete.cases(D[, c("export","exhaustion", ctrl)]), ]
  r <- ppcor::pcor.test(d$export, d$exhaustion, d[, ctrl, drop = FALSE],
                        method = "spearman")
  say("  + %-28s rho = %+.3f  p = %.3g  (n=%d)", label, r$estimate, r$p.value, nrow(d))
  r$estimate
}
r_h  <- pc("hypoxia", "hypoxia")
r_s  <- pc("stage", "stage")
r_hs <- pc(c("hypoxia","stage"), "hypoxia + stage")
r_all<- pc(c("hypoxia","stage","n_myeloid","n_tcell"), "hypoxia + stage + cell counts")

# ---- C: mixed model with cohort as random effect ----------------------------
say("\n=== C: cohort as random effect ===")
if (requireNamespace("lme4", quietly = TRUE)) {
  suppressPackageStartupMessages(library(lme4))
  D2 <- D %>% mutate(across(c(export, exhaustion, hypoxia), ~ as.numeric(scale(.x))))
  m <- lmer(exhaustion ~ export + hypoxia + stage + (1 | study), data = D2)
  cf <- summary(m)$coefficients
  say("  export coef = %+.3f (SE %.3f, t = %.2f)",
      cf["export","Estimate"], cf["export","Std. Error"], cf["export","t value"])
  say("  (|t| > 2 => export adds signal beyond hypoxia, stage and cohort)")
} else say("  lme4 not installed -- skipping")

# ---- D: within-stage, the cleanest test -------------------------------------
say("\n=== D: WITHIN-STAGE correlation (is it a stage proxy?) ===")
for (s in sort(unique(D$stage))) {
  d <- D %>% filter(stage == s)
  lbl <- c("normal","premalignant","carcinoma")[s]
  if (nrow(d) < 10) { say("  %-13s n=%d, too few", lbl, nrow(d)); next }
  ct <- suppressWarnings(cor.test(d$export, d$exhaustion, method="spearman"))
  say("  %-13s n=%3d  rho = %+.3f  p = %.3g", lbl, nrow(d), ct$estimate, ct$p.value)
}
say("  READ: a positive rho WITHIN carcinoma alone means the association is not")
say("  produced by stage. This is the single most important line in the script.")

# ---- verdict -----------------------------------------------------------------
ca <- D %>% filter(stage == 3)
r_ca <- if (nrow(ca) >= 10) cor(ca$export, ca$exhaustion, method="spearman",
                                use="complete.obs") else NA
say("\n=== VERDICT ===")
say("  unadjusted %.2f | adjusted (hyp+stage) %.2f | within-carcinoma %.2f",
    raw$estimate, r_hs, r_ca)
ok <- !is.na(r_hs) && r_hs > 0.10 && !is.na(r_ca) && r_ca > 0.10
say("  %s", if (ok)
  "SURVIVES -- report the adjusted and within-stage values alongside the raw one" else
  "DOES NOT SURVIVE -- the association is largely explained by stage and/or hypoxia.\n     Withdraw the cross-compartment coupling claim, or report it explicitly as\n     confounded with disease stage.")

write_csv(D, file.path(OUT, "per_sample_variables.csv"))
p <- ggplot(D, aes(export, exhaustion, colour = factor(stage))) +
  geom_point(size = 2, alpha = .8) + geom_smooth(method = "lm", se = FALSE) +
  scale_colour_manual(values = c(`1`="#92c5de",`2`="#f4a582",`3`="#b2182b"),
                      labels = c("normal","premalignant","carcinoma"), name = NULL) +
  labs(title = "Myeloid export vs T-cell exhaustion, by stage",
       subtitle = sprintf("unadjusted rho = %.2f; adjusted for hypoxia + stage = %.2f; within carcinoma = %.2f",
                          raw$estimate, r_hs, r_ca),
       x = "myeloid export allocation", y = "T-cell exhaustion score") +
  theme_bw(base_size = 11)
ggsave(file.path(OUT, "fig_partial.png"), p, width = 7.5, height = 5.5, dpi = 150, bg = "white")
say("\n-> %s", OUT)
