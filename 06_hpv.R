#!/usr/bin/env Rscript
# ==============================================================================
# 06_hpv.R  --  does HPV status shift metabolic ALLOCATION? (bulk_before)
#
#   Rscript 06_hpv.R --arm bulk_before
#
# DESIGN, informed by every mistake earlier in this pipeline:
#
#  * FLUX SHARE, never absolute. HPV+ and HPV- tumors differ in overall activity;
#    absolute flux would just re-report that (the 05 global-shift trap).
#
#  * WITHIN-COHORT, then meta. HPV+ is NOT evenly distributed: He and Koh have
#    ZERO HPV+ patients, so a pooled HPV+/- contrast is confounded with cohort
#    (HPV+ = Kim/Pittet/Tirosh, HPV- piles up in He/Koh). Fix: test HPV inside
#    each cohort that has BOTH arms (Kim, Pittet, Tirosh), then Stouffer-meta with
#    concordance required. He and Koh cannot inform the contrast and are dropped.
#
#  * EMPIRICAL REL_FLOOR. HPV effects are small in absolute terms (~0.002-0.01
#    share) but that is 3-15x a typical residual, so the scale-free REL_FLOOR
#    (|effect| >= 10% of the module's own share) is the right primary gate.
#
#  * TARGETED FIRST, scan second. The headline question is ONE line: is MCT4
#    lactate export (M_78) higher in HPV- than HPV+? the SPP1+ TAM / MMP11+ CAF niche
#    is the HPV- niche; if export is HPV--enriched, the metabolic finding and the
#    CoVarNet niche lock together.
#
#  * GLYCAN/MEVALONATE AUDIT. The old HPV claim was "concentrated in N-linked
#    glycan synthesis (8 of 15)". M_113 (mevalonate) and M_116 (CoQ) are FILED
#    under that supermodule but are NOT glycan. Every hit is re-labelled by its
#    actual top gene so that error cannot recur.
#
# Output -> results/all_sites/<ARM>/06_hpv/<level>/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(purrr); library(stringr) })

LEVELS      <- c("finaltype", "subtype")
MR          <- load_modref()
usable      <- readRDS(file.path(PATHS$results, "03_qc", "usable_modules.rds"))
HPV_COHORTS <- c("Kim", "Pittet", "Tirosh")     # cohorts with BOTH arms
MIN_PT_ARM  <- 5

MODULE_GENE <- file.path(ROOT, "reference", "module_gene_m168.csv")
if (!file.exists(MODULE_GENE)) MODULE_GENE <- file.path(ROOT, "module_gene_m168.csv")
mgraw <- read.csv(MODULE_GENE, row.names = 1, stringsAsFactors = FALSE)
top_gene <- setNames(vapply(rownames(mgraw), function(m){
  g <- as.character(unlist(mgraw[m,])); g <- g[!is.na(g)&g!=""&g!="nan"]
  if (length(g)) g[1] else NA_character_ }, character(1)), rownames(mgraw))

hpv_shift <- function(share, hpv) {
  neg <- share[hpv=="Negative"]; pos <- share[hpv=="Positive"]
  if (length(neg) < MIN_PT_ARM || length(pos) < MIN_PT_ARM) return(c(delta=NA, z=NA, p=NA))
  w <- suppressWarnings(wilcox.test(neg, pos))
  d <- median(neg) - median(pos)                 # + => higher in HPV-negative
  z <- qnorm(pmax(pmin(1 - w$p.value/2, 1-1e-16), 1e-16)) * sign(d)
  c(delta = d, z = z, p = w$p.value)
}
stouffer <- function(z, w) { ok <- is.finite(z)&is.finite(w); z<-z[ok]; w<-w[ok]
  if (!length(z)) return(NA_real_); sum(w*z)/sqrt(sum(w^2)) }

run_level <- function(level) {
  OUT <- file.path(PATHS$results, "06_hpv", level)
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  say("\n================ LEVEL: %s ================", level)

  cw <- readRDS(file.path(PATHS$cw_dir,
         sprintf("cw_%s.rds", if (level=="finaltype") "finaltype" else "subtype"))) %>%
    filter(n_cells >= GATES$min_cells,
           tissue %in% names(STAGE_MAP)[STAGE_MAP == "CA"],
           hpv %in% c("Negative","Positive"),
           study %in% HPV_COHORTS)
  unit_col <- if (level=="finaltype") "finaltype" else "annot_2nd"
  mcols <- intersect(mcols_of(cw), usable)

  pt <- cw %>% rename(celltype = annot_4th, unit = !!unit_col) %>% filter(!is.na(unit)) %>%
    group_by(study, celltype, unit, patient, hpv) %>%
    summarise(across(all_of(mcols), ~ weighted.mean(.x, n_cells)),
              n_cells = sum(n_cells), .groups="drop")
  X <- as.matrix(pt[, mcols]); Sh <- X / rowSums(X)
  ps <- bind_cols(pt %>% select(study, celltype, unit, patient, hpv), as_tibble(Sh))
  long <- ps %>% pivot_longer(all_of(mcols), names_to="module", values_to="share")

  pc <- long %>% group_by(study, celltype, unit, module) %>%
    summarise(n_neg = sum(hpv=="Negative"), n_pos = sum(hpv=="Positive"),
              tt = list(hpv_shift(share, hpv)),
              med_neg = median(share[hpv=="Negative"]),
              med_pos = median(share[hpv=="Positive"]), .groups="drop") %>%
    mutate(delta = map_dbl(tt,"delta"), z = map_dbl(tt,"z"), p = map_dbl(tt,"p"), tt=NULL) %>%
    filter(!is.na(z))

  meta <- pc %>% group_by(celltype, unit, module) %>%
    summarise(n_cohorts = n(),
              z_meta = stouffer(z, sqrt(pmin(n_neg, n_pos))),
              delta_meta = weighted.mean(delta, sqrt(pmin(n_neg,n_pos))),
              concordant = (all(z>0)|all(z<0)) & n()>=2,
              cohorts = paste(study, collapse="+"),
              med_neg = mean(med_neg), med_pos = mean(med_pos), .groups="drop") %>%
    filter(!is.na(z_meta)) %>%
    mutate(p_meta = 2*pnorm(-abs(z_meta)), fdr = p.adjust(p_meta,"BH"),
           rel_effect = abs(delta_meta) / pmax(med_neg, med_pos, 1e-9),
           higher_in = ifelse(delta_meta > 0, "HPV-negative", "HPV-positive"),
           gene = top_gene[module]) %>%
    left_join(MR$rxn, by="module") %>% left_join(MR$sm, by="module") %>%
    mutate(pass_fdr = fdr < FDR_THRESH,
           pass_rel = rel_effect >= REL_FLOOR,
           survivor = pass_fdr & pass_rel & concordant) %>%
    arrange(desc(survivor), fdr, desc(rel_effect))

  write_csv(meta, file.path(OUT, "hpv_all.csv"))
  write_csv(meta %>% filter(survivor), file.path(OUT, "hpv_survivors.csv"))

  funnel <- tibble(step=factor(c("tested","FDR<0.10","+ rel_floor 10%","+ concordant (survivors)"),
      levels=c("tested","FDR<0.10","+ rel_floor 10%","+ concordant (survivors)")),
    n=c(nrow(meta), sum(meta$pass_fdr), sum(meta$pass_fdr&meta$pass_rel), sum(meta$survivor)))
  say("cohorts used: %s | tested %s contrasts",
      paste(HPV_COHORTS, collapse="+"), format(nrow(meta), big.mark=","))
  say("funnel: %s", paste(sprintf("%s=%d", funnel$step, funnel$n), collapse=" -> "))
  say("direction of survivors: HPV- = %d | HPV+ = %d",
      sum(meta$survivor & meta$higher_in=="HPV-negative"),
      sum(meta$survivor & meta$higher_in=="HPV-positive"))

  say("\n-- supermodule label vs ACTUAL top gene (the old mislabel) --")
  aud <- meta %>% filter(survivor) %>% count(supermodule, gene) %>% arrange(desc(n))
  if (nrow(aud)) print(as.data.frame(head(aud, 12)), row.names=FALSE)
  say("  (M_113 top gene should read FDPS/mevalonate, M_116 COQ* -- NOT 'glycan')")

  if (sum(meta$survivor)) {
    say("\n-- HPV SURVIVORS (%s) --", level)
    print(as.data.frame(meta %>% filter(survivor) %>% slice_head(n=25) %>%
      transmute(celltype, unit, module, gene, reaction=str_trunc(reaction,26),
                higher_in, rel=round(rel_effect,2), fdr=signif(fdr,2),
                z=round(z_meta,1), coh=cohorts)), row.names=FALSE)
  } else say("\n-- no HPV survivors at %s --", level)

  say("\n================ TARGETED: is lactate export HPV--enriched? ================")
  KEY <- tibble(module=c("M_78","M_6","M_4"),
                what=c("MCT lactate export","LDH-assoc flux","glycolysis 3PD->pyr"))
  targ <- meta %>% inner_join(KEY, by="module") %>%
    filter(celltype %in% c("Myeloid","Fibroblast") |
           unit %in% c("MQ_C4_SPP1","FB_C1_MMP11","DC_C2_cDC2","MQ","cDC")) %>%
    transmute(what, celltype, unit, gene, higher_in,
              rel=round(rel_effect,3), delta=signif(delta_meta,3),
              z=round(z_meta,2), fdr=signif(fdr,3), cohorts, concordant)
  if (nrow(targ)) {
    print(as.data.frame(targ %>% arrange(match(what, KEY$what), desc(rel))), row.names=FALSE)
    say("\n  READ: MCT export (M_78) 'higher_in = HPV-negative' with concordant cohorts")
    say("  is the result that ties the lactate finding to the SPP1+ TAM / MMP11+ CAF niche (HPV-).")
  } else say("  no MCT/LDH/glycolysis contrasts passed coverage in myeloid/stroma")
  write_csv(targ, file.path(OUT, "targeted_lactate_hpv.csv"))

  pv <- ggplot(meta, aes(delta_meta, -log10(fdr), colour=survivor)) +
    geom_vline(xintercept=0, colour="grey50") +
    geom_hline(yintercept=-log10(FDR_THRESH), linetype=2, colour="grey60") +
    geom_point(alpha=.6, size=1.3) +
    scale_colour_manual(values=c(`FALSE`="grey75", `TRUE`="#b2182b")) +
    labs(title=sprintf("HPV allocation shifts (%s, %s)", ARM, level),
         subtitle="flux share, within-cohort meta (Kim+Pittet+Tirosh). + delta = higher in HPV-negative.",
         x="delta share (HPV-neg minus HPV-pos)", y="-log10 FDR") + theme_bw(base_size=11)
  ggsave(file.path(OUT,"fig_volcano.png"), pv, width=8, height=5.5, dpi=150, bg="white")

  meta
}

out <- lapply(LEVELS, run_level); names(out) <- LEVELS
say("\n=== SUMMARY (%s) ===", ARM)
for (lv in LEVELS) if (!is.null(out[[lv]]))
  say("  %-10s: %d HPV survivors", lv, sum(out[[lv]]$survivor))
say("\nNEXT: ground any survivor against its enzyme (07c-style) before believing it.")
