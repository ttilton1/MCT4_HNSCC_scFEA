#!/usr/bin/env Rscript
# 01d_floor_and_mark_comparable.R
#   Floor negative flux at zero and record which modules are comparable across
#   the two arms. Run with --arm bulk_after.
#   in   data/combined_wide/all_sites/<ARM>/cw_finaltype.rds
#   out  cw_finaltype.rds (floored; original kept as cw_finaltype_raw.rds)
#        comparable_modules.rds
#   note Flooring shifts 0.13% of total flux mass, preserves within-unit module
#        ranking (median rho 0.9951) and is exactly neutral for every module
#        carrying a claim: M_78, M_6, M_2, M_113 and M_116 have no negative
#        values. Negative per-cell output is regression noise around zero for
#        modules with little enzyme expression, not reverse flux; the modules
#        are directional.
#   note M_3 and M_4 are marked not comparable across arms. In bulk_after they
#        carry ~4e-05 median flux with 37-39% of units negative, whereas M_4 was
#        the largest single effect in bulk_before.
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
