#!/usr/bin/env Rscript
# ==============================================================================
# 07b_grounding_scalefree.R  --  the HONEST version of the enzyme grounding.
#
#   Rscript 07b_grounding_scalefree.R --arm bulk_before
#
# 07 found real enzymes at median rho = 0.48 vs RANDOM genes at 0.20 (41% of
# random genes above |0.3|). Random genes should sit at ZERO. That positive
# baseline is a SHARED SCALE: units with more total expression have more total
# flux, so anything correlates with anything. It is the same global-offset
# confound that invented 132 fake progression hits in 05.
#
# 07 also reported a top-25 at rho ~= 0.999 -- every one of them a ONE-GENE
# module (M_12=FH, M_150=UMPS, M_170=GMPS...). For a single-gene module scFEA's
# flux is essentially a learned monotone transform OF THAT GENE, so correlating
# them recovers rho~1 BY CONSTRUCTION. That is circular, not validation, and it
# inflated the "real enzyme" median that the negative control was judged against.
#
# THIS SCRIPT FIXES BOTH:
#   1. SCALE-FREE. Correlate flux SHARE (module flux / unit total flux) against
#      expression SHARE (gene expr / unit total expr over the 661 scFEA genes).
#      Both sides lose the per-unit total, so a shared-scale artifact cannot
#      produce a correlation. The negative control is re-run on this basis: if
#      random genes now sit at ~0 and real enzymes stay high, the grounding is
#      REAL and every flux claim upstream is licensed.
#   2. NO CIRCULAR MODULES. Modules with < MIN_GENES genes are excluded from the
#      grounded set (reported separately, flagged as circular).
#
# ALSO: per-gene output, because module medians lie. M_6 (Pyruvate->Lactate) has
# median rho ~0.16 but LDHA alone correlates at ~0.93 -- the median is dragged
# down by LDHB, which catalyses the REVERSE reaction (lactate->pyruvate).
# Averaging a forward and a reverse isozyme is a meaningless summary. Claims must
# name the enzyme, not the gene set.
#
# Output -> results/all_sites/<ARM>/07b_grounding_scalefree/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(purrr); library(stringr) })

OUT <- file.path(PATHS$results, "07b_grounding_scalefree")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
MR <- load_modref()

MIN_N_COHORT <- 8     # units per cohort to compute a correlation
MIN_GENES    <- 2     # modules with 1 gene are CIRCULAR -> excluded from grounded
N_NEG        <- 300    # random non-module genes per module (negative control)
set.seed(1)

MODULE_GENE <- file.path(ROOT, "reference", "module_gene_m168.csv")
if (!file.exists(MODULE_GENE)) MODULE_GENE <- file.path(ROOT, "module_gene_m168.csv")
mg <- read.csv(MODULE_GENE, row.names = 1, stringsAsFactors = FALSE)
module_genes <- lapply(rownames(mg), function(m) {
  g <- as.character(unlist(mg[m, ])); unique(g[!is.na(g) & g != "" & g != "nan"])
})
names(module_genes) <- rownames(mg)

# ---- expression (same units as flux), then SHARE ------------------------------
say("building expression table...")
expr_list <- list()
for (ct in TYPES) {
# NOTE (arm-safety): both arms consume the SAME MAGIC-imputed expression -- the
# imputation ran once in shared/ and the arms differ only in whether pseudobulking
# happens before or after scFEA. Unit-level expression is therefore identical, so
# grounding reads the single bulk_before copy regardless of --arm. This also
# isolates the comparison: any arm difference in enzyme correlation reflects the
# FLUX method, not different input expression.
  ep <- file.path(ROOT, "data/input/all_sites", "bulk_before", ct, "train_expression.csv")
  mp <- file.path(file.path(ROOT, "metadata", "all_sites", "bulk_before"), ct, "metadata_units.csv")
  if (!file.exists(ep) || !file.exists(mp)) next
  e <- fread(ep); genes <- as.character(e[[1]])
  E <- as.matrix(e[, -1]); rownames(E) <- genes; ids <- colnames(E)
  mu <- fread(mp); setnames(mu, 1, "id"); mu[, id := as.character(id)]
  mu <- as.data.frame(mu); rownames(mu) <- mu$id; mu <- mu[ids, , drop = FALSE]
  expr_list[[ct]] <- as_tibble(t(E)) %>% mutate(id = ids) %>%
    left_join(mu %>% select(id, patient, tissue, finaltype, annot_4th, study, n_cells),
              by = "id") %>%
    group_by(patient, tissue, finaltype, annot_4th, study) %>%
    summarise(across(all_of(genes), ~ weighted.mean(.x, n_cells)),
              n_cells = sum(n_cells), .groups = "drop")
}
EX <- bind_rows(expr_list) %>% filter(n_cells >= GATES$min_cells)
gene_cols <- setdiff(names(EX), c("patient","tissue","finaltype","annot_4th","study","n_cells"))
key <- c("patient","tissue","finaltype","annot_4th","study")

FL <- readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds")) %>% filter(n_cells >= GATES$min_cells)
mcols <- mcols_of(FL)

dedup <- function(d) d %>% group_by(across(all_of(key))) %>% slice_head(n=1) %>% ungroup()
FLk <- dedup(FL %>% select(all_of(key), all_of(mcols)))
EXk <- dedup(EX %>% select(all_of(key), all_of(gene_cols)))
D <- inner_join(FLk, EXk, by = key)

# ---- THE TRANSFORM: both sides -> SHARE --------------------------------------
Fm <- as.matrix(D[, mcols]);     Fs <- Fm / rowSums(Fm)          # flux share
Em <- as.matrix(D[, gene_cols]); Es <- Em / rowSums(Em)          # expression share
DS <- bind_cols(D %>% select(all_of(key)), as_tibble(Fs), as_tibble(Es))
say("units: %s | modules: %d | genes: %d (all SHARE-transformed)",
    format(nrow(DS), big.mark=","), length(mcols), length(gene_cols))

fisher_meta <- function(rho, n) {
  ok <- is.finite(rho) & is.finite(n) & n > 3; rho <- rho[ok]; n <- n[ok]
  if (!length(rho) || sum(n-3) <= 0)
    return(tibble(rho_meta=NA_real_, lo=NA_real_, hi=NA_real_, k=0L))
  z <- atanh(pmin(pmax(rho,-0.999),0.999)); w <- n-3
  zm <- sum(w*z)/sum(w); se <- sqrt(1/sum(w))
  tibble(rho_meta = tanh(zm), lo = tanh(zm-1.96*se), hi = tanh(zm+1.96*se), k = length(z))
}
per_cohort_rho <- function(sub, a, b) {
  sub %>% group_by(study) %>% filter(n() >= MIN_N_COHORT) %>%
    summarise(n = n(), rho = suppressWarnings(cor(.data[[a]], .data[[b]],
              method="spearman", use="complete.obs")), .groups="drop") %>%
    filter(is.finite(rho))
}

# ---- pairs to ground ----------------------------------------------------------
grab <- function(p, src) {
  if (!file.exists(p)) return(tibble())
  x <- read_csv(p, show_col_types = FALSE)
  if (!all(c("celltype","module") %in% names(x))) return(tibble())
  x %>% distinct(celltype, module) %>% mutate(src = src)
}
hits <- bind_rows(
  grab(file.path(PATHS$results,"04_specialization/finaltype/specialization_survivors.csv"),"specialization"),
  grab(file.path(PATHS$results,"05c_progression_meta/finaltype/meta_survivors.csv"),"progression"),
  grab(file.path(PATHS$results,"05c_progression_meta/subtype/meta_survivors.csv"),"progression")) %>%
  group_by(celltype, module) %>%
  summarise(src = paste(sort(unique(src)), collapse="+"), .groups="drop") %>%
  filter(!is.na(module))
say("(celltype, module) pairs: %d", nrow(hits))

# ---- correlate ---------------------------------------------------------------
res <- list(); neg <- list()
for (i in seq_len(nrow(hits))) {
  ct <- hits$celltype[i]; mod <- hits$module[i]
  sub <- DS %>% filter(annot_4th == ct)
  if (nrow(sub) < MIN_N_COHORT) next
  gs <- intersect(module_genes[[mod]], gene_cols)
  if (!length(gs)) next

  for (g in gs) {
    pc <- per_cohort_rho(sub, mod, g); if (nrow(pc) < 1) next
    m <- fisher_meta(pc$rho, pc$n)
    lv <- numeric(0)
    if (nrow(pc) >= 2) for (s in pc$study) {
      k <- pc[pc$study != s, , drop=FALSE]
      if (nrow(k)) { mm <- fisher_meta(k$rho, k$n)
                     if (is.finite(mm$rho_meta)) lv <- c(lv, mm$rho_meta) }
    }
    res[[length(res)+1]] <- tibble(
      celltype=ct, module=mod, gene=g, src=hits$src[i],
      n_module_genes = length(gs),
      n_cohorts=m$k, rho=m$rho_meta, lo=m$lo, hi=m$hi,
      all_same_sign = all(pc$rho>0)|all(pc$rho<0),
      loso_min = if(length(lv)) min(lv) else NA_real_,
      loso_max = if(length(lv)) max(lv) else NA_real_)
  }
  ng <- sample(setdiff(gene_cols, module_genes[[mod]]), min(N_NEG, length(gene_cols)))
  for (g in ng) {
    pc <- per_cohort_rho(sub, mod, g); if (nrow(pc) < 1) next
    m <- fisher_meta(pc$rho, pc$n)
    neg[[length(neg)+1]] <- tibble(celltype=ct, module=mod, gene=g, rho=m$rho_meta,
                                   n_module_genes = length(gs))
  }
}
R <- bind_rows(res); NG <- bind_rows(neg)

R <- R %>% mutate(
  circular    = n_module_genes < MIN_GENES,     # 1-gene module: flux ~= that gene
  ci_excl_0   = (lo>0 & hi>0) | (lo<0 & hi<0),
  multi       = n_cohorts >= 2,
  loso_stable = !is.na(loso_min) & sign(loso_min)==sign(loso_max) &
                sign(loso_min)==sign(rho),
  grounded    = !circular & ci_excl_0 & all_same_sign & multi & loso_stable & rho > 0)
write_csv(R, file.path(OUT, "gene_level_scalefree.csv"))
write_csv(NG, file.path(OUT, "negative_control_scalefree.csv"))

# ================= THE DECISIVE COMPARISON ====================================
Rm  <- R  %>% filter(!circular)                 # multi-gene only
NGm <- NG %>% filter(n_module_genes >= MIN_GENES)
say("\n============================================================")
say("NEGATIVE CONTROL, SCALE-FREE  (multi-gene modules only)")
say("============================================================")
say("  real enzymes : median rho = %+.3f | %.0f%% |rho|>0.3 | n = %d",
    median(Rm$rho, na.rm=TRUE), 100*mean(abs(Rm$rho)>0.3, na.rm=TRUE), nrow(Rm))
say("  random genes : median rho = %+.3f | %.0f%% |rho|>0.3 | n = %d",
    median(NGm$rho, na.rm=TRUE), 100*mean(abs(NGm$rho)>0.3, na.rm=TRUE), nrow(NGm))
say("  [07, scale-CONFOUNDED, for reference: real 0.482 vs random 0.202]")
say("")
say("  If random genes are now ~0 and enzymes stay high -> the shared-scale")
say("  artifact is GONE and the flux genuinely tracks its own enzymes.")
say("  If BOTH collapse -> scFEA flux here is mostly a restatement of total")
say("  expression, and the upstream flux claims are not supported.")
if (nrow(Rm) && nrow(NGm)) {
  w <- suppressWarnings(wilcox.test(Rm$rho, NGm$rho))
  say("  Wilcoxon enzymes vs random: p = %.3g", w$p.value)
}
say("\n  circular (1-gene) modules excluded: %d correlations", sum(R$circular))

# ---- module-level (multi-gene only) ------------------------------------------
modsum <- Rm %>% group_by(celltype, module, src) %>%
  summarise(n_genes=n(), med_rho=median(rho,na.rm=TRUE), max_rho=max(rho,na.rm=TRUE),
            top_gene=gene[which.max(rho)], n_grounded=sum(grounded), .groups="drop") %>%
  left_join(MR$rxn, by="module") %>% arrange(desc(med_rho))
write_csv(modsum, file.path(OUT, "module_level_scalefree.csv"))

say("\n=== KEY MODULES (the ones carrying every claim) ===")
KEY <- c("M_2","M_3","M_4","M_6","M_78","M_113","M_116")
print(as.data.frame(modsum %>% filter(module %in% KEY) %>%
  transmute(celltype, module, reaction=str_trunc(reaction,30), src, n_genes,
            med_rho=round(med_rho,3), max_rho=round(max_rho,3), top_gene,
            grounded=sprintf("%d/%d", n_grounded, n_genes))), row.names=FALSE)

# ---- per-gene detail: module medians LIE (LDHA vs LDHB) ----------------------
say("\n=== PER-GENE: M_6 (Pyruvate->Lactate). LDHA is FORWARD; LDHB is REVERSE. ===")
say("    A module median that averages them is meaningless -- name the enzyme.")
print(as.data.frame(R %>% filter(module=="M_6") %>%
  transmute(celltype, gene, rho=round(rho,3), lo=round(lo,3), hi=round(hi,3),
            same_sign=all_same_sign, grounded) %>% arrange(celltype, desc(rho))),
  row.names=FALSE)

say("\n=== PER-GENE: M_78 (Lactate export) = MCT1/MCT4 ===")
print(as.data.frame(R %>% filter(module=="M_78") %>%
  transmute(celltype, gene, rho=round(rho,3), grounded) %>% arrange(desc(rho))),
  row.names=FALSE)

# ---- figures ------------------------------------------------------------------
cmp <- bind_rows(Rm %>% transmute(rho, set="module enzymes"),
                 NGm %>% transmute(rho, set="random genes"))
p1 <- ggplot(cmp, aes(rho, fill=set)) +
  geom_density(alpha=.55, colour=NA) + geom_vline(xintercept=0, colour="grey30") +
  scale_fill_manual(values=c(`module enzymes`="#b2182b", `random genes`="grey60")) +
  labs(title=sprintf("Scale-free grounding vs negative control (%s)", ARM),
       subtitle="flux SHARE vs expression SHARE; 1-gene (circular) modules excluded.\nRandom genes should now sit at zero.",
       x="meta Spearman rho", y="density", fill=NULL) +
  theme_bw(base_size=11)
ggsave(file.path(OUT,"fig_negative_control_scalefree.png"), p1, width=8, height=5, dpi=150, bg="white")

kg <- R %>% filter(module %in% KEY, !circular) %>%
  mutate(lab = paste0(celltype, " | ", module, " | ", gene))
if (nrow(kg)) {
  p2 <- ggplot(kg, aes(rho, reorder(lab, rho), colour=grounded)) +
    geom_vline(xintercept=0, colour="grey40") +
    geom_errorbarh(aes(xmin=lo, xmax=hi), height=.2, colour="grey65") +
    geom_point(size=2) +
    scale_colour_manual(values=c(`TRUE`="#b2182b", `FALSE`="grey65")) +
    labs(title="Enzyme grounding of the key modules (scale-free)",
         subtitle="glycolysis (M_2/3/4), lactate production (M_6), export (M_78), mevalonate (M_113), CoQ (M_116)",
         x="meta Spearman rho (flux share ~ expression share)", y=NULL) +
    theme_bw(base_size=8)
  ggsave(file.path(OUT,"fig_key_modules.png"), p2, width=10,
         height=max(6, 0.16*nrow(kg)), dpi=150, bg="white", limitsize=FALSE)
}
say("\n-> %s", OUT)
