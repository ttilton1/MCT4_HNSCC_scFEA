#!/usr/bin/env Rscript
# ==============================================================================
# 01c_negative_flux_diagnostic.R
#
#   Rscript 01c_negative_flux_diagnostic.R --arm bulk_after
#
# 01b found ~20% of unit-level module means are NEGATIVE in bulk_after (bulk_before
# pseudobulk flux was strictly non-negative, range 0-4.66). Every scale-free result
# in this project uses share = flux / rowSums(flux), which is only interpretable
# for non-negative values. So before ANY bulk_after analysis runs we need to know:
#
#   Q1. WHERE are the negatives? If they sit in modules that are near-zero and
#       already unusable, flooring at 0 changes nothing that matters. If they
#       contaminate M_78 / M_6 / M_2 / M_3 / M_4 -- the modules carrying every
#       claim -- then bulk_after cannot be compared to bulk_before on share.
#
#   Q2. HOW BIG are they relative to the positive signal in the same module?
#       A -1e-9 next to a 4e-3 is numerical noise. A -0.03 next to a 3.8e-3
#       median is real mass and cannot be discarded silently.
#
#   Q3. Does FLOORING AT ZERO change the biology? scFEA modules are directional
#       (each is a specified reaction), so negative predicted flux is not reverse
#       flux -- it is the network's regression noise around zero for modules with
#       little enzyme expression in that cell. Flooring is therefore defensible,
#       but only if it does not reshuffle the module ranking within a unit.
#
# Output -> results/all_sites/bulk_after/01c_negative_diagnostic/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(ggplot2); library(stringr) })
stopifnot(ARM == "bulk_after")

OUT <- file.path(PATHS$results, "01c_negative_diagnostic")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
MR <- load_modref()

KEY <- c("M_2","M_3","M_4","M_6","M_78","M_113","M_116")   # carry every claim

CW <- readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))
mc <- mcols_of(CW)
X  <- as.matrix(CW[, mc, with = FALSE])
say("units %s | modules %d", format(nrow(X), big.mark=","), ncol(X))

# ---- Q1: per-module negativity ----------------------------------------------
permod <- tibble(
  module      = mc,
  frac_neg    = colMeans(X < 0),
  med_val     = apply(X, 2, median),
  med_abs     = apply(X, 2, function(v) median(abs(v))),
  min_val     = apply(X, 2, min),
  max_val     = apply(X, 2, max),
  # how big is the typical negative relative to the typical positive IN THIS MODULE
  neg_burden  = apply(X, 2, function(v) {
                  p <- v[v > 0]; n <- v[v < 0]
                  if (!length(p) || !length(n)) return(0)
                  median(abs(n)) / median(p) })) %>%
  left_join(MR$rxn, by = "module") %>%
  arrange(desc(frac_neg))

write_csv(permod, file.path(OUT, "per_module_negativity.csv"))

say("\n=== MOST NEGATIVE MODULES (are these the ones we care about?) ===")
print(as.data.frame(permod %>% slice_head(n = 15) %>%
  transmute(module, reaction = str_trunc(reaction, 30),
            pct_neg = round(100*frac_neg), med = signif(med_val,2),
            neg_burden = round(neg_burden,2))), row.names = FALSE)

say("\n=== THE KEY MODULES (these carry every claim) ===")
print(as.data.frame(permod %>% filter(module %in% KEY) %>%
  transmute(module, reaction = str_trunc(reaction, 30),
            pct_neg = round(100*frac_neg), med = signif(med_val,3),
            min = signif(min_val,2), neg_burden = round(neg_burden,2))),
  row.names = FALSE)
say("  READ: pct_neg near 0 for the key modules => the negatives live elsewhere")
say("  and flooring at zero cannot touch the headline results.")

# ---- Q2: is negativity a property of LOW-FLUX modules? -----------------------
say("\n=== is negativity concentrated in low-signal modules? ===")
cc <- cor(permod$frac_neg, permod$med_abs, method = "spearman")
say("  Spearman(frac_negative, median |flux|) = %+.2f", cc)
say("  (strongly negative => negatives are a near-zero-module phenomenon, i.e. noise)")

p1 <- ggplot(permod, aes(med_abs, frac_neg)) +
  geom_point(alpha=.6, size=1.6, colour="grey40") +
  geom_point(data = permod %>% filter(module %in% KEY),
             colour = "#b2182b", size = 3) +
  ggrepel::geom_text_repel(data = permod %>% filter(module %in% KEY),
             aes(label = module), colour = "#b2182b", size = 3.2) +
  scale_x_log10() +
  labs(title = "Negative flux is a low-signal phenomenon (bulk_after)",
       subtitle = "each point = one module. Red = the modules carrying our claims.",
       x = "median |flux| across units (log10)", y = "fraction of units negative") +
  theme_bw(base_size = 11)
ggsave(file.path(OUT,"fig_negativity_vs_signal.png"), p1, width=8, height=5.5,
       dpi=150, bg="white")

# ---- Q3: does flooring at zero change anything? ------------------------------
say("\n=== EFFECT OF FLOORING AT ZERO ===")
Xf <- pmax(X, 0)
say("  total flux mass removed: %.2f%%", 100*(sum(X) - sum(Xf))/sum(Xf))
say("  unit totals: raw median %.3f -> floored median %.3f",
    median(rowSums(X)), median(rowSums(Xf)))

S_raw <- X  / rowSums(X)
S_flr <- Xf / rowSums(Xf)
# does the WITHIN-UNIT ranking of modules survive flooring?
rank_cor <- vapply(seq_len(nrow(X)), function(i)
  suppressWarnings(cor(S_raw[i,], S_flr[i,], method = "spearman")), numeric(1))
say("  within-unit module rank correlation raw vs floored: median rho = %.4f",
    median(rank_cor, na.rm = TRUE))
say("  (near 1.0 => flooring does NOT reshuffle which modules dominate a unit)")

for (k in intersect(KEY, colnames(X))) {
  r <- suppressWarnings(cor(S_raw[,k], S_flr[,k], method = "spearman"))
  say("    %-6s share raw vs floored: rho = %.4f | %% units negative = %.1f%%",
      k, r, 100*mean(X[,k] < 0))
}

# ---- gate comparison across arms (the fair comparison) -----------------------
say("\n=== POST-GATE UNIT COUNTS (the sets we actually analyse) ===")
bb <- file.path(ROOT, "data/combined_wide/all_sites/bulk_before/cw_finaltype.rds")
if (file.exists(bb)) {
  B <- readRDS(bb)
  key3 <- c("patient","tissue","finaltype")
  Bg <- B  %>% filter(n_cells >= GATES$min_cells) %>% distinct(across(all_of(key3)))
  Ag <- CW %>% filter(n_cells >= GATES$min_cells) %>% distinct(across(all_of(key3)))
  say("  bulk_before gated: %s units", format(nrow(Bg), big.mark=","))
  say("  bulk_after  gated: %s units", format(nrow(Ag), big.mark=","))
  say("  overlap          : %s units", format(nrow(inner_join(Bg, Ag, by=key3)), big.mark=","))
  say("  (01b's 'shared units' line was wrong -- it joined on duplicated keys)")
}

say("\n=== RECOMMENDATION ===")
kn <- permod %>% filter(module %in% KEY) %>% pull(frac_neg)
if (max(kn) < 0.05 && median(rank_cor, na.rm=TRUE) > 0.99) {
  say("  SAFE TO FLOOR. The key modules are essentially never negative and")
  say("  flooring does not reshuffle within-unit module rankings. Add")
  say("  pmax(flux, 0) before the share transform in 05b/05c/06 and note it in")
  say("  methods as: negative per-cell predictions (regression noise around zero")
  say("  for low-expression modules) are floored at zero before normalisation.")
} else {
  say("  DO NOT FLOOR BLINDLY. Negativity touches the key modules (max %.1f%% of",
      100*max(kn))
  say("  units) and/or flooring changes module rankings. bulk_after may not be")
  say("  comparable to bulk_before on the share scale -- consider comparing the")
  say("  arms on rank-based measures, or restricting to modules with <5%% negatives.")
}
say("\n-> %s", OUT)
