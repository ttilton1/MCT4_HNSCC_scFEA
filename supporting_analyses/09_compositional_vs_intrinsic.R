#!/usr/bin/env Rscript
# ==============================================================================
# 09_compositional_vs_intrinsic.R
#
#   Rscript 09_compositional_vs_intrinsic.R --arm bulk_before
#
# The niche's lactate-export story has THREE separable parts, and a reviewer will
# insist we separate them. Total export by a cell state through progression =
#
#       (how MANY of those cells there are)      <- ABUNDANCE (compositional)
#     x (how much MARKER they express)           <- IDENTITY
#     x (what FRACTION of metabolism is export)  <- ALLOCATION (intrinsic)
#
# Our headline is the THIRD term (05c: MCT flux share rises). But if the first
# term drives everything -- if there are simply more SPP1+ TAM / MMP11+ CAF in
# carcinoma -- then the "metabolic reprogramming" framing collapses into "more
# TAMs", which is the infiltration result of Meng et al. The finding is only distinct if ALLOCATION
# rises INDEPENDENT of abundance.
#
# This is the compositional-vs-intrinsic split flagged in the old Axis 2 write-up
# (resident-macrophage reprogramming vs SPP1+ TAM infiltration). It tests all
# three terms on the same axis (oral NL->Pre->CA, He+Kim) for the niche cells, so
# we can state plainly which mechanism is operating.
#
# Output -> results/all_sites/<ARM>/09_compositional_vs_intrinsic/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(purrr); library(stringr) })

OUT <- file.path(PATHS$results, "09_compositional_vs_intrinsic")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# the niche cell states + their markers + the export module
NICHE <- tribble(
  ~celltype,    ~finaltype,     ~marker,  ~note,
  "Myeloid",    "MQ_C4_SPP1",   "SPP1",   "SPP1+ TAM",
  "Fibroblast", "FB_C1_MMP11",  "MMP11",  "MMP11+ CAF",
  "Myeloid",    "DC_C2_cDC2",   "CLEC10A","cDC2")
EXPORT_MODULE <- "M_78"

jt_test <- function(x, g) {
  ok <- is.finite(x) & !is.na(g); x <- x[ok]
  g <- droplevels(factor(g[ok], levels = STAGE_LEVELS))
  ns <- table(g); k <- length(ns); N <- sum(ns)
  if (k < 3 || any(ns < 2)) return(c(z=NA, p=NA))
  U <- 0; L <- levels(g)
  for (i in 1:(k-1)) for (j in (i+1):k) { a<-x[g==L[i]]; b<-x[g==L[j]]
    U <- U + sum(outer(a,b,"<")) + .5*sum(outer(a,b,"==")) }
  EU <- (N^2-sum(ns^2))/4; VU <- (N^2*(2*N+3)-sum(ns^2*(2*ns+3)))/72
  if (VU<=0) return(c(z=NA,p=NA)); z<-(U-EU)/sqrt(VU); c(z=z, p=2*pnorm(-abs(z)))
}
trend <- function(df, val) {   # per-patient value, JT per cohort, over stage
  d <- df %>% mutate(stage = factor(unname(STAGE_MAP[tissue]), levels=STAGE_LEVELS)) %>%
    filter(!is.na(stage))
  map_dfr(PROG_COHORTS, function(s) {
    x <- d %>% filter(study==s)
    if (n_distinct(x$stage) < 3) return(tibble(study=s, z=NA, p=NA,
        Normal=NA, Pre=NA, CA=NA, n_pt=n_distinct(x$patient)))
    tt <- jt_test(x[[val]], x$stage)
    tibble(study=s, z=tt["z"], p=tt["p"],
           Normal=median(x[[val]][x$stage=="Normal"], na.rm=TRUE),
           Pre   =median(x[[val]][x$stage=="Pre"],    na.rm=TRUE),
           CA    =median(x[[val]][x$stage=="CA"],     na.rm=TRUE),
           n_pt=n_distinct(x$patient))
  })
}

# ---- load flux (for allocation + abundance) and expression (for marker) ------
cw <- as.data.frame(readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))) %>%
  filter(!tissue %in% STAGE_EXCLUDED, study %in% PROG_COHORTS)
mcols <- mcols_of(cw)

# ABUNDANCE: fraction of the parent major type that is this finaltype, per sample.
# (n_cells of this unit / total n_cells of that major type in the same patient/tissue)
abund <- cw %>% group_by(patient, tissue, study, annot_4th) %>%
  mutate(parent_cells = sum(n_cells)) %>% ungroup() %>%
  mutate(frac_of_parent = n_cells / parent_cells)

# ALLOCATION: export flux SHARE, gated (needs a real pseudobulk)
alloc <- cw %>% filter(n_cells >= GATES$min_cells)
X <- as.matrix(alloc[, mcols]); alloc$export_share <- (X / rowSums(X))[, EXPORT_MODULE]

# MARKER: mean imputed marker expression from train_expression, per unit
say("loading marker expression...")
mk_list <- list()
for (ct in unique(NICHE$celltype)) {
  ep <- file.path(ROOT, "data/input/all_sites", ARM, ct, "train_expression.csv")
  mp <- file.path(PATHS$meta_dir, ct, "metadata_units.csv")
  if (!file.exists(ep)) next
  e <- fread(ep); g <- as.character(e[[1]]); E <- as.matrix(e[,-1]); rownames(E) <- g
  want <- intersect(NICHE$marker, g); if (!length(want)) next
  mu <- fread(mp); setnames(mu,1,"id"); mu[, id:=as.character(id)]
  mu <- as.data.frame(mu); rownames(mu) <- mu$id; mu <- mu[colnames(E),,drop=FALSE]
  mk_list[[ct]] <- as_tibble(t(E[want,,drop=FALSE])) %>%
    mutate(patient=mu$patient, tissue=mu$tissue, finaltype=mu$finaltype,
           study=mu$study, n_cells=mu$n_cells)
}
MK <- bind_rows(mk_list)

# ---- run all three terms for each niche cell state ---------------------------
say("\n============================================================")
say("THREE TERMS per niche cell state (oral axis, He + Kim)")
say("============================================================")
all_res <- list()
for (i in seq_len(nrow(NICHE))) {
  ct <- NICHE$celltype[i]; ft <- NICHE$finaltype[i]; mk <- NICHE$marker[i]
  say("\n--- %s (%s) | marker %s ---", ft, NICHE$note[i], mk)

  a <- abund %>% filter(annot_4th==ct, finaltype==ft) %>%
    group_by(patient, tissue, study) %>%
    summarise(frac_of_parent = sum(n_cells)/first(parent_cells), .groups="drop")
  ta <- trend(a, "frac_of_parent") %>% mutate(term="1. ABUNDANCE (% of parent)")

  al <- alloc %>% filter(annot_4th==ct, finaltype==ft) %>%
    group_by(patient, tissue, study) %>%
    summarise(export_share = weighted.mean(export_share, n_cells), .groups="drop")
  tl <- trend(al, "export_share") %>% mutate(term="3. ALLOCATION (export share)")

  tm <- tibble()
  if (mk %in% names(MK)) {
    mm <- MK %>% filter(finaltype==ft) %>%
      group_by(patient, tissue, study) %>%
      summarise(marker = weighted.mean(.data[[mk]], n_cells), .groups="drop")
    tm <- trend(mm, "marker") %>% mutate(term="2. MARKER (expression)")
  }

  r <- bind_rows(ta, tm, tl) %>% mutate(finaltype=ft, note=NICHE$note[i], marker=mk)
  all_res[[ft]] <- r
  print(as.data.frame(r %>% transmute(term, study, n_pt,
        Normal=signif(Normal,3), CA=signif(CA,3),
        CA_over_N = round(CA/Normal,2), z=round(z,2), p=signif(p,2))), row.names=FALSE)
}
RES <- bind_rows(all_res)
write_csv(RES, file.path(OUT, "three_terms_trends.csv"))

# ---- verdict per cell state --------------------------------------------------
say("\n============================================================")
say("VERDICT: is the story INTRINSIC (allocation up on its own) or")
say("just COMPOSITIONAL (more cells)?")
say("============================================================")
verdict <- RES %>% filter(!is.na(z)) %>%
  group_by(finaltype, note, term) %>%
  summarise(both_up = all(z > 0, na.rm=TRUE),
            both_sig = all(p < 0.10, na.rm=TRUE),
            median_CA_over_N = round(median(CA/Normal, na.rm=TRUE),2),
            .groups="drop")
print(as.data.frame(verdict), row.names=FALSE)
say("\nREAD:")
say("  If ALLOCATION rises (term 3) INDEPENDENT of ABUNDANCE (term 1), the finding")
say("  is INTRINSIC metabolic reprogramming -- not explained by increased TAM abundance alone.")
say("  If only ABUNDANCE rises, the story is compositional and the metabolic claim")
say("  is weak. If BOTH rise, state both: the niche expands AND commits more per cell.")

# ---- figure: three terms, one row per cell state -----------------------------
pdat <- RES %>% filter(!is.na(Normal)) %>%
  select(finaltype, note, term, study, Normal, Pre, CA) %>%
  pivot_longer(c(Normal, Pre, CA), names_to="stage", values_to="val") %>%
  mutate(stage=factor(stage, levels=STAGE_LEVELS),
         term=factor(term, levels=c("1. ABUNDANCE (% of parent)",
                                    "2. MARKER (expression)",
                                    "3. ALLOCATION (export share)")))
p <- ggplot(pdat, aes(stage, val, colour=study, group=study)) +
  geom_line(linewidth=.9) + geom_point(size=2.4) +
  facet_grid(note ~ term, scales="free_y", switch="y") +
  scale_colour_manual(values=c(He="#d6604d", Kim="#4393c3")) +
  labs(title="Compositional vs intrinsic: three terms of the lactate-export story",
       subtitle="Does the niche just get BIGGER (abundance), or do its cells commit MORE per cell (allocation)?\nThe headline requires ALLOCATION (rightmost column) to rise on its own.",
       x=NULL, y=NULL, colour="cohort") +
  theme_bw(base_size=10) + theme(strip.text.y.left=element_text(angle=0))
ggsave(file.path(OUT,"fig_three_terms.png"), p, width=11, height=7, dpi=150, bg="white")
say("\n-> %s", OUT)
