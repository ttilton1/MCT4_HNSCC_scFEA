#!/usr/bin/env Rscript
# ==============================================================================
# 01d_floor_and_mark_comparable.R
#
#   Rscript 01d_floor_and_mark_comparable.R --arm bulk_after
#
# 01c established two things:
#
#  1. FLOORING AT ZERO IS SAFE. Removing negative predictions changes total flux
#     mass by 0.13%, preserves within-unit module ranking (median rho 0.9951),
#     and is EXACTLY neutral (rho = 1.0000) for every module carrying a claim:
#     M_78 (lactate export), M_6 (LDH), M_2, M_113, M_116 -- all 0% negative.
#     Negative per-cell scFEA output is regression noise around zero for modules
#     with little enzyme expression in that cell; pseudobulking averaged it away
#     in bulk_before, per-cell prediction exposes it. It is NOT reverse flux --
#     the modules are directional.
#
#  2. M_3 AND M_4 ARE NOT COMPARABLE ACROSS ARMS. In bulk_after they carry
#     ~4e-05 median flux, ~270x smaller than M_78, with 37-39% of units negative.
#     In bulk_before, M_4 was the single largest effect in the study. That is a
#     real arm disagreement, not a rounding issue, and any M_3/M_4 comparison
#     between arms would be meaningless. Note also that M_2 (G6P->G3P) carries
#     100x the flux of M_3 (G3P->3PD) in bulk_after despite being consecutive
#     steps -- scFEA's balance constraints evidently behave differently at
#     single-cell resolution than on pseudobulk. Worth stating in methods.
#
# This script writes the floored cw (keeping the raw as *_raw.rds) and saves a
# COMPARABLE MODULE list so downstream cross-arm claims can be restricted to
# modules that are well-behaved in BOTH arms.
#
# After this, 03/04/05c/06/07/09 run unchanged with --arm bulk_after.
# ==============================================================================
source("config_all_sites.R")
stopifnot(ARM == "bulk_after")
NEG_TOL <- 0.05      # a module must be <5% negative to be cross-arm comparable

for (lvl in c("finaltype", "subtype")) {
  f <- file.path(PATHS$cw_dir, sprintf("cw_%s.rds", lvl))
  if (!file.exists(f)) { say("missing %s -- skip", f); next }
  CW <- readRDS(f)
  mc <- mcols_of(CW)
  X  <- as.matrix(CW[, mc, with = FALSE])

  raw_path <- file.path(PATHS$cw_dir, sprintf("cw_%s_raw.rds", lvl))
  if (!file.exists(raw_path)) { saveRDS(CW, raw_path); say("kept raw -> %s", raw_path) }

  neg_frac <- colMeans(X < 0)
  Xf <- pmax(X, 0)
  say("[%s] units %s | mass change %+.3f%% | modules >=%.0f%% negative: %d",
      lvl, format(nrow(X), big.mark=","),
      100*(sum(Xf)-sum(X))/sum(X), 100*NEG_TOL, sum(neg_frac >= NEG_TOL))

  CW[, (mc) := as.data.table(Xf)]
  saveRDS(CW, f)
  say("  floored -> %s", f)

  if (lvl == "finaltype") {
    comparable <- names(neg_frac)[neg_frac < NEG_TOL]
    # also require the module to be usable in bulk_before (03 QC)
    ub <- file.path(ROOT, "results/all_sites/bulk_before/03_qc/usable_modules.rds")
    if (file.exists(ub)) {
      before_ok <- readRDS(ub)
      comparable <- intersect(comparable, before_ok)
      say("  usable in bulk_before: %d | comparable in BOTH arms: %d",
          length(before_ok), length(comparable))
    }
    saveRDS(comparable, file.path(PATHS$cw_dir, "comparable_modules.rds"))
    say("  -> comparable_modules.rds (%d modules)", length(comparable))

    KEY <- c("M_2","M_3","M_4","M_6","M_78","M_113","M_116")
    say("\n  KEY MODULE STATUS for cross-arm comparison:")
    for (k in KEY) {
      ok <- k %in% comparable
      say("    %-6s %-18s neg=%4.1f%%  %s", k,
          ifelse(ok, "COMPARABLE", "NOT comparable"),
          100*neg_frac[[k]],
          ifelse(ok, "", "<- exclude from arm-vs-arm claims"))
    }
  }
}

say("\nNEXT: 03_qc_flux.R --arm bulk_after, then 04 / 05c / 06 / 07c / 09.")
say("When comparing arms, intersect hits with comparable_modules.rds.")
