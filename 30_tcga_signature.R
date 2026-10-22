#!/usr/bin/env Rscript
# ==============================================================================
# 30_tcga_signature.R  --  does the lactate metabolic program matter in TCGA?
#
#   Rscript 30_tcga_signature.R
#
# WHAT BULK TCGA CAN AND CANNOT DO FOR THIS FINDING.
#   CANNOT: the isoform specificity (MCT4-in-myeloid vs MCT1-in-T) or per-cell
#           allocation -- bulk has no cell types. The progression axis -- TCGA-HNSC
#           is almost all carcinoma.
#   CAN:    cohort-scale survival, and HPV stratification with real sample sizes
#           (unlike the sparse HPV+ in the oral scRNA cohorts).
#
# WHETHER THE SIGNATURE IS PROGNOSTIC INDEPENDENTLY OF SPP1.
#   Meng et al. reported SPP1+ TAM infiltration -> worse OS in TCGA. So SPP1 alone
#   being prognostic is NOT new. The test that matters: does the EXPORT-SPECIFIC
#   gene SLC16A3 (MCT4) add prognostic value OVER SPP1? If a Cox model with
#   SLC16A3 beats SPP1-only (likelihood-ratio test), the export program is a
#   distinct signal. If SPP1 carries everything, the signature adds nothing beyond SPP1.
#
# STEPS
#   1. pull TCGA-HNSC expression + clinical (TCGAbiolinks) OR read a local matrix
#   2. build the export signature: SIG = EXPORT_CORE + NICHE_ID + HYPOXIA,
#      i.e. SLC16A3, LDHA, HK2, SPP1, MMP11, CA9, VEGFA (seven genes)
#   3. Cox OS: signature score, adjusted for age/stage
#   4. INCREMENTAL TEST: SLC16A3 (and the export core) over SPP1 alone (LRT)
#   5. HPV stratification: is the signature higher in HPV-? (TCGA has power here)
#   6. co-expression: does SLC16A3 track SPP1/MMP11 across tumors?
#
# Output -> results/tcga/
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(survival); library(ggplot2)
})
ROOT <- Sys.getenv("SCFEA_ROOT", unset = "/Users/Tommy/dechenli_lab_local/scfea26")
OUT  <- file.path(ROOT, "results", "tcga"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
say  <- function(...) message(sprintf(...))

# ---- the signature -----------------------------------------------------------
EXPORT_CORE <- c("SLC16A3", "LDHA", "HK2")           # the export-specific program
NICHE_ID    <- c("SPP1", "MMP11")                    # niche identity (SPP1 is the established marker, Meng et al.)
HYPOXIA     <- c("CA9", "VEGFA")                     # context (weak in our data; test anyway)
SIG <- unique(c(EXPORT_CORE, NICHE_ID, HYPOXIA))

# ==============================================================================
# 1. DATA -- cBioPortal "HNSC (TCGA, PanCancer Atlas)" local files.
#    expression: data_mrna_seq_v2_rsem.txt (Hugo_Symbol, Entrez, then samples; RSEM)
#    clinical  : hnsc_tcga_pan_can_atlas_2018_clinical_data.tsv
# ==============================================================================
TCGA_DIR  <- file.path(ROOT, "tcga")
EXPR_FILE <- file.path(TCGA_DIR, "hnsc_tcga_pan_can_atlas_2018", "data_mrna_seq_v2_rsem.txt")
CLIN_FILE <- file.path(TCGA_DIR, "hnsc_tcga_pan_can_atlas_2018_clinical_data.tsv")
stopifnot(file.exists(EXPR_FILE), file.exists(CLIN_FILE))

say("reading cBioPortal expression (RSEM)...")
ex <- data.table::fread(EXPR_FILE)
ex <- ex[Hugo_Symbol != "" & !is.na(Hugo_Symbol)]          # drop blank gene rows
ex <- ex[!duplicated(Hugo_Symbol)]                          # first entry per symbol
genes <- ex$Hugo_Symbol
M <- as.matrix(ex[, -(1:2)]); rownames(M) <- genes          # genes x samples (RSEM)
M[!is.finite(M)] <- 0
expr <- log2(M + 1)                                         # RSEM -> log2(x+1)
say("  expr: %d genes x %d samples (log2 RSEM)", nrow(expr), ncol(expr))

say("reading clinical...")
clin_raw <- readr::read_tsv(CLIN_FILE, show_col_types = FALSE)
clin <- clin_raw %>% transmute(
  sample  = `Sample ID`,
  OS.time = as.numeric(`Overall Survival (Months)`) * 30.4375,   # months -> days
  OS      = as.integer(grepl("DECEASED", `Overall Survival Status`)),
  age     = as.numeric(`Diagnosis Age`),
  stage   = `Neoplasm Disease Stage American Joint Committee on Cancer Code`,
  hpv     = `Subtype`,                                            # HNSC_HPV+ / HNSC_HPV-
  buffa_hypoxia  = as.numeric(`Buffa Hypoxia Score`),
  ragnum_hypoxia = as.numeric(`Ragnum Hypoxia Score`),
  winter_hypoxia = as.numeric(`Winter Hypoxia Score`),
  icd10   = `ICD-10 Classification`,                              # sub-site proxy
  icdo3_site = `International Classification of Diseases for Oncology, Third Edition ICD-O-3 Site Code`)
say("  clinical: %d rows | HPV: %s",
    nrow(clin), paste(names(table(clin$hpv)), table(clin$hpv), collapse="  "))
D <- list(expr = expr, clin = clin)
expr <- D$expr; clin <- D$clin
say("expr: %d genes x %d samples | clinical: %d rows",
    nrow(expr), ncol(expr), nrow(clin))

# Harmonize BOTH sides to the cBioPortal Sample ID (e.g. TCGA-4P-AA8J-01).
# The expression columns and the clinical `Sample ID` are already in this form;
# the previous truncation to a 12-char patient id dropped the join key.
norm_id <- function(x) sub("^(TCGA-[^-]+-[^-]+-[0-9]+).*", "\\1", x)  # keep through -01
colnames(expr) <- norm_id(colnames(expr))
expr <- expr[, !duplicated(colnames(expr)), drop=FALSE]
clin <- clin %>% mutate(sample = norm_id(sample)) %>% distinct(sample, .keep_all=TRUE)
say("  expr samples: %d | clinical samples: %d | overlap: %d",
    ncol(expr), nrow(clin), length(intersect(colnames(expr), clin$sample)))

have <- intersect(SIG, rownames(expr))
say("signature genes present: %s", paste(have, collapse=", "))
miss <- setdiff(SIG, rownames(expr)); if (length(miss)) say("  MISSING: %s", paste(miss, collapse=", "))

# ---- z-score each gene, build scores -----------------------------------------
Z <- t(scale(t(expr[have,,drop=FALSE])))                 # gene z across samples
score <- function(genes) colMeans(Z[intersect(genes, rownames(Z)),,drop=FALSE], na.rm=TRUE)
df <- tibble(sample = colnames(expr),
             sig_full   = score(SIG),
             export_core= score(EXPORT_CORE),
             SPP1       = Z["SPP1", ],
             SLC16A3    = if ("SLC16A3" %in% rownames(Z)) Z["SLC16A3", ] else NA_real_,
             MMP11      = if ("MMP11" %in% rownames(Z)) Z["MMP11", ] else NA_real_,
             LDHA       = if ("LDHA" %in% rownames(Z)) Z["LDHA", ] else NA_real_,
             HK2        = if ("HK2" %in% rownames(Z)) Z["HK2", ] else NA_real_,
             CA9        = if ("CA9" %in% rownames(Z)) Z["CA9", ] else NA_real_,
             VEGFA      = if ("VEGFA" %in% rownames(Z)) Z["VEGFA", ] else NA_real_) %>%
  inner_join(clin, by="sample") %>%
  filter(!is.na(OS.time), OS.time > 0, !is.na(OS)) %>%
  mutate(stage_simple = case_when(
    grepl("IV", stage) ~ "III-IV",
    grepl("III", stage) ~ "III-IV",
    grepl("II", stage) ~ "I-II",
    grepl("Stage I$|Stage IA|Stage IB|^I$", stage) ~ "I-II",
    TRUE ~ NA_character_))

# per-patient scores and survival, so downstream figures (Kaplan-Meier) can be
# rebuilt without re-running this script
write_csv(df %>% select(sample, export_core, sig_full, SLC16A3, LDHA, HK2,
                        SPP1, MMP11, CA9, VEGFA, OS, OS.time, age, stage_simple,
                        buffa_hypoxia, ragnum_hypoxia, winter_hypoxia),
          file.path(OUT, "per_patient_scores.csv"))

stage_ok <- mean(!is.na(df$stage_simple)) > 0.5
say("analysable samples (OS available): %d | events: %d | stage usable: %s (%.0f%% non-missing)",
    nrow(df), sum(df$OS), stage_ok, 100*mean(!is.na(df$stage_simple)))

# ==============================================================================
# 3. SURVIVAL: does the export signature predict OS?
# ==============================================================================
say("\n=== Cox OS ~ signature (continuous, z-scored) ===")
cox_one <- function(var, adj=TRUE) {
  base <- if (stage_ok) "age + strata(stage_simple)" else "age"
  f <- if (adj) sprintf("Surv(OS.time, OS) ~ %s + %s", var, base)
       else      sprintf("Surv(OS.time, OS) ~ %s", var)
  d <- df %>% filter(!is.na(.data[[var]]), !is.na(age))
  if (stage_ok) d <- d %>% filter(!is.na(stage_simple))
  m <- tryCatch(coxph(as.formula(f), data=d), error=function(e) NULL)
  if (is.null(m)) return(tibble(var=var, HR=NA, lo=NA, hi=NA, p=NA, n=nrow(d)))
  s <- summary(m)
  tibble(var=var, HR=s$coef[var,"exp(coef)"],
         lo=s$conf.int[var,"lower .95"], hi=s$conf.int[var,"upper .95"],
         p=s$coef[var,"Pr(>|z|)"], n=nrow(d))
}
cox_tab <- bind_rows(lapply(c("sig_full","export_core","SPP1","SLC16A3","MMP11","LDHA","HK2","CA9","VEGFA"),
                            cox_one, adj=TRUE))
print(as.data.frame(cox_tab %>% mutate(across(c(HR,lo,hi),~round(.x,2)), p=signif(p,3))),
      row.names=FALSE)
write_csv(cox_tab, file.path(OUT, "cox_univariate.csv"))

# ==============================================================================
# 4. THE KEY TEST: does SLC16A3 (export) add value OVER SPP1 alone?
# ==============================================================================
say("\n=== INCREMENTAL VALUE: SLC16A3 over SPP1 (the novelty test) ===")
d2 <- df %>% filter(!is.na(SPP1), !is.na(SLC16A3), !is.na(age))
adj <- if (stage_ok) "+ strata(stage_simple)" else ""
if (stage_ok) d2 <- d2 %>% filter(!is.na(stage_simple))
m_spp1  <- coxph(as.formula(sprintf("Surv(OS.time,OS) ~ SPP1 + age %s", adj)), data=d2)
m_both  <- coxph(as.formula(sprintf("Surv(OS.time,OS) ~ SPP1 + SLC16A3 + age %s", adj)), data=d2)
lrt <- anova(m_spp1, m_both)
say("  SPP1-only model:   AIC = %.1f", AIC(m_spp1))
say("  + SLC16A3 model:   AIC = %.1f", AIC(m_both))
say("  Likelihood-ratio test (does SLC16A3 add signal over SPP1?): p = %.4g",
    lrt$`P(>|Chi|)`[2])
print(as.data.frame(round(summary(m_both)$coef, 4)))
say("  READ: if LRT p < 0.05 and SLC16A3's own p < 0.05, the export program is")
say("  prognostic INDEPENDENT of SPP1 -- i.e. not explained by SPP1 alone.")
say("  If not, SPP1 carries the signal.")
capture.output(list(spp1=summary(m_spp1), both=summary(m_both), lrt=lrt),
               file=file.path(OUT, "incremental_cox.txt"))

# ==============================================================================
# 5. HPV STRATIFICATION -- done carefully.
#   Three traps this must avoid:
#   (a) the composite blends MCT4 with LDHA/HK2 -> test SLC16A3 ALONE too.
#   (b) HPV+ HNSC is enriched for OROPHARYNX, HPV- for oral cavity/larynx. Site
#       has its own stromal/immune composition, so a raw HPV difference may be an
#       ANATOMICAL SITE effect. Adjust for site (or restrict to oropharynx, which
#       has both arms).
#   (c) bulk cannot separate "more MCT4 per cell" from "more macrophages". A bulk
#       difference is CONSISTENT with the scRNA finding but does not prove the
#       cell-intrinsic part -- state that limitation.
# ==============================================================================
say("\n=== HPV: is lactate export (MCT4) higher in HPV-negative? ===")
hpv_col <- "hpv"
if (hpv_col %in% names(df) && any(!is.na(df$hpv))) {
  dh <- df %>% mutate(hpv2 = case_when(
      grepl("HPV\\+|HPV_pos|_pos", hpv, ignore.case=TRUE) ~ "Positive",
      grepl("HPV\\-|HPV_neg|_neg", hpv, ignore.case=TRUE) ~ "Negative",
      TRUE ~ NA_character_)) %>% filter(!is.na(hpv2))
  say("  HPV+ = %d | HPV- = %d", sum(dh$hpv2=="Positive"), sum(dh$hpv2=="Negative"))

  # ICD-10 as a sub-site proxy (C00-C06 oral cavity, C09-C10 oropharynx/tonsil,
  # C32 larynx). Not perfect, but lets us check the oropharynx confound.
  has_site <- "icd10" %in% names(dh) && any(!is.na(dh$icd10))
  if (has_site) {
    dh <- dh %>% mutate(site = case_when(
      grepl("C0[0-6]", icd10) ~ "oral_cavity",
      grepl("C09|C10|C01", icd10) ~ "oropharynx",
      grepl("C32", icd10) ~ "larynx",
      grepl("C12|C13", icd10) ~ "hypopharynx", TRUE ~ "other"))
    say("  ICD-10 sub-sites: %s",
        paste(names(table(dh$site)), table(dh$site), collapse="  "))
  } else say("  no ICD site proxy available")

  if (nrow(dh) > 20) {
    # (a) test each quantity separately
    say("\n  -- HPV- vs HPV+ (Wilcoxon), each quantity on its own --")
    for (v in c("SLC16A3","export_core","SPP1","MMP11","LDHA")) {
      if (!v %in% names(dh) || all(is.na(dh[[v]]))) next
      w <- wilcox.test(dh[[v]] ~ dh$hpv2)
      say("    %-11s HPV- %.2f vs HPV+ %.2f | p = %.3g", v,
          median(dh[[v]][dh$hpv2=="Negative"], na.rm=TRUE),
          median(dh[[v]][dh$hpv2=="Positive"], na.rm=TRUE), w$p.value)
    }

    # (b) does the SLC16A3 difference SURVIVE adjusting for anatomical site?
    if (has_site) {
      say("\n  -- SLC16A3 ~ HPV, ADJUSTED for ICD sub-site (the key test) --")
      dh2 <- dh %>% filter(!is.na(SLC16A3), site != "other")
      if (n_distinct(dh2$hpv2)==2 && nrow(dh2) > 30) {
        m_raw <- lm(SLC16A3 ~ hpv2, data=dh2)
        m_adj <- lm(SLC16A3 ~ hpv2 + site, data=dh2)
        say("    raw    HPV coef = %+.3f (p=%.3g)",
            coef(summary(m_raw))["hpv2Positive","Estimate"],
            coef(summary(m_raw))["hpv2Positive","Pr(>|t|)"])
        say("    +site  HPV coef = %+.3f (p=%.3g)",
            coef(summary(m_adj))["hpv2Positive","Estimate"],
            coef(summary(m_adj))["hpv2Positive","Pr(>|t|)"])
        say("    READ: if the HPV effect shrinks toward 0 after adding site, the")
        say("    'HPV difference' was largely the oropharynx/oral-cavity confound.")
      }

      # (c) WITHIN OROPHARYNX only -- same site, differ only by HPV
      oro <- dh %>% filter(site=="oropharynx", !is.na(SLC16A3))
      if (nrow(oro) > 15 && n_distinct(oro$hpv2)==2) {
        w <- wilcox.test(SLC16A3 ~ hpv2, data=oro)
        say("\n  -- WITHIN OROPHARYNX (n=%d): SLC16A3 HPV- %.2f vs HPV+ %.2f | p=%.3g",
            nrow(oro), median(oro$SLC16A3[oro$hpv2=="Negative"]),
            median(oro$SLC16A3[oro$hpv2=="Positive"]), w$p.value)
        say("     The confound-free HPV test: same site, differ only by HPV.")
      } else say("\n  (too few oropharynx samples with both arms for within-site test)")
    }

    # (d) BONUS: does MCT4/export track a VALIDATED bulk hypoxia score? The scRNA
    #     hypoxia test failed; the clinical file has Buffa/Ragnum/Winter scores, so
    #     we can ask the hypoxia question a second, independent way.
    say("\n  -- MCT4 & export vs validated hypoxia scores (scRNA hypoxia failed) --")
    for (hs in c("buffa_hypoxia","ragnum_hypoxia","winter_hypoxia")) {
      if (!hs %in% names(df) || all(is.na(df[[hs]]))) next
      r1 <- cor(df$SLC16A3, df[[hs]], use="complete.obs", method="spearman")
      r2 <- cor(df$export_core, df[[hs]], use="complete.obs", method="spearman")
      say("    %-16s SLC16A3 rho=%+.2f | export_core rho=%+.2f", hs, r1, r2)
    }
    say("    (positive rho = export tracks hypoxia at bulk level -> supports the")
    say("     hypoxia framing the single-cell data could not confirm)")

    dh$mct4 <- dh$SLC16A3
    ph <- ggplot(dh, aes(hpv2, mct4, fill=hpv2)) +
      geom_boxplot(width=.5, outlier.size=.6) +
      scale_fill_manual(values=c(Negative="#b2182b", Positive="#4393c3"), guide="none") +
      labs(title="MCT4 (SLC16A3) by HPV status, TCGA-HNSC",
           subtitle="bulk expression -- cannot separate per-cell export from macrophage abundance",
           x=NULL, y="SLC16A3 expression (z)") + theme_bw(base_size=12)
    ggsave(file.path(OUT,"fig_hpv_mct4.png"), ph, width=6, height=5, dpi=150, bg="white")
  }
} else say("  no HPV field -- set the column name in the clinical table (TCGA HPV is")
  say("  in a molecular-subtype or p16 field; name varies by release).")

# ==============================================================================
# 6. CO-EXPRESSION: does SLC16A3 track the niche markers across tumors?
# ==============================================================================
say("\n=== co-expression across tumors (bulk, suggestive) ===")
ce <- tibble(pair=c("SLC16A3~SPP1","SLC16A3~MMP11","SLC16A3~LDHA","SLC16A3~HK2"),
             rho=c(cor(df$SLC16A3, df$SPP1, use="complete.obs", method="spearman"),
                   cor(df$SLC16A3, df$MMP11, use="complete.obs", method="spearman"),
                   cor(df$SLC16A3, df$LDHA, use="complete.obs", method="spearman"),
                   suppressWarnings(cor(Z["SLC16A3",], Z["HK2",], use="complete.obs", method="spearman"))))
print(as.data.frame(ce %>% mutate(rho=round(rho,3))), row.names=FALSE)

# ---- KM figure for the full signature (median split) -------------------------
df <- df %>% mutate(grp = ifelse(sig_full > median(sig_full), "high", "low"))
fit <- survfit(Surv(OS.time, OS) ~ grp, data=df)
sd  <- survdiff(Surv(OS.time, OS) ~ grp, data=df)
pval <- 1 - pchisq(sd$chisq, length(sd$n)-1)
km <- broom::tidy(fit) %>% mutate(grp = sub("grp=","",strata))
kp <- ggplot(km, aes(time/365.25, estimate, colour=grp)) +
  geom_step(linewidth=1) +
  scale_colour_manual(values=c(high="#b2182b", low="#4393c3")) +
  labs(title="lactate metabolic signature and overall survival (TCGA-HNSC)",
       subtitle=sprintf("median split, log-rank p = %.3g", pval),
       x="years", y="overall survival", colour="signature") +
  coord_cartesian(ylim=c(0,1)) + theme_bw(base_size=12)
ggsave(file.path(OUT,"fig_km_signature.png"), kp, width=7, height=5, dpi=150, bg="white")
say("\nKM log-rank p = %.3g", pval)

say("\n-> %s", OUT)
say("HEADLINE to check: (1) export_core HR>1 and significant; (2) SLC16A3 LRT over")
say("SPP1 significant (= independent of SPP1); (3) signature higher in HPV-.")
