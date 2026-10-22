#!/usr/bin/env Rscript
# ==============================================================================
# 01e_fix_unit_duplicates.R
#
#   Rscript 01e_fix_unit_duplicates.R --arm bulk_after
#
# PROBLEM. 03_qc reports 6,126 gated units for bulk_after, but only 5,305 distinct
# (patient, tissue, finaltype) triples exist. bulk_before has 5,459 rows / 5,431
# distinct -- essentially one row per unit. So bulk_after is SPLITTING units.
#
# CAUSE. 01b aggregated with
#     by = (patient, tissue, finaltype, annot_2nd, annot_4th, study, hpv, site,
#           site_group)
# If ANY of those annotation columns varies within a patient x tissue x finaltype
# -- e.g. two samples from the same patient/tissue carrying different `site` or
# `site_group` labels -- the unit fragments into several rows.
#
# WHY IT MATTERS. (1) the arms then have different unit counts and the comparison
# is not like-for-like; (2) fragments hold fewer cells each, so the n_cells >= 10
# gate removes different things in each arm; (3) near-duplicate rows double-count
# the same cells in every downstream test.
#
# FIX. Re-aggregate to the canonical unit -- patient x tissue x finaltype -- taking
# the n_cells-weighted mean of flux and, for annotation columns, the value from the
# largest constituent fragment. This reproduces bulk_before's unit definition
# exactly, so the arms become comparable.
#
# Run AFTER 01b/01d, BEFORE 03/04/05c/06/07/09.
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(data.table) })
stopifnot(ARM == "bulk_after")

KEY3 <- c("patient", "tissue", "finaltype")
ANNOT <- c("annot_2nd", "annot_4th", "study", "hpv", "site", "site_group")

for (lvl in c("finaltype", "subtype")) {
  f <- file.path(PATHS$cw_dir, sprintf("cw_%s.rds", lvl))
  if (!file.exists(f)) { say("missing %s -- skip", f); next }
  CW <- as.data.table(readRDS(f))
  mc <- mcols_of(CW)
  key <- if (lvl == "finaltype") KEY3 else c("patient","tissue","annot_2nd")
  ann <- setdiff(intersect(ANNOT, names(CW)), key)

  n0 <- nrow(CW); nd <- uniqueN(CW[, ..key])
  say("\n[%s] rows %s | distinct units %s | fragments %s",
      lvl, format(n0, big.mark=","), format(nd, big.mark=","),
      format(n0 - nd, big.mark=","))

  if (n0 == nd) { say("  already one row per unit -- nothing to do"); next }

  # ---- WHICH column is causing the split? ----
  say("  columns that vary within a unit:")
  for (a in ann) {
    v <- CW[, .(k = uniqueN(get(a))), by = key][, sum(k > 1)]
    if (v > 0) say("    %-12s varies in %s units", a, format(v, big.mark=","))
  }

  # ---- collapse: n_cells-weighted mean flux, modal annotation ----
  setorderv(CW, c(key, "n_cells"), c(rep(1, length(key)), -1))
  ANN <- CW[, lapply(.SD, first), by = key, .SDcols = ann]   # largest fragment wins
  FLX <- CW[, c(lapply(.SD, function(x) weighted.mean(x, n_cells)),
                .(n_cells = sum(n_cells))),
            by = key, .SDcols = mc]
  OUT <- merge(ANN, FLX, by = key)
  setcolorder(OUT, c(key, ann, "n_cells", mc))

  say("  collapsed: %s -> %s rows", format(n0, big.mark=","), format(nrow(OUT), big.mark=","))
  say("  gated (n_cells >= %d): %s units", GATES$min_cells,
      format(sum(OUT$n_cells >= GATES$min_cells), big.mark=","))

  bak <- file.path(PATHS$cw_dir, sprintf("cw_%s_fragmented.rds", lvl))
  if (!file.exists(bak)) saveRDS(CW, bak)
  saveRDS(OUT, f)
  say("  -> %s (pre-fix kept at %s)", f, basename(bak))
}

# ---- confirm the arms now line up -------------------------------------------
bb <- file.path(ROOT, "data/combined_wide/all_sites/bulk_before/cw_finaltype.rds")
if (file.exists(bb)) {
  B <- as.data.table(readRDS(bb)); A <- as.data.table(readRDS(
        file.path(PATHS$cw_dir, "cw_finaltype.rds")))
  Bg <- unique(B[n_cells >= GATES$min_cells, ..KEY3])
  Ag <- unique(A[n_cells >= GATES$min_cells, ..KEY3])
  say("\n=== ARM ALIGNMENT (gated, finaltype) ===")
  say("  bulk_before : %s units", format(nrow(Bg), big.mark=","))
  say("  bulk_after  : %s units", format(nrow(Ag), big.mark=","))
  say("  shared      : %s units", format(nrow(merge(Bg, Ag, by = KEY3)), big.mark=","))
  say("  after-only  : %s | before-only : %s",
      format(nrow(fsetdiff(Ag, Bg)), big.mark=","),
      format(nrow(fsetdiff(Bg, Ag)), big.mark=","))
  say("  (a large shared set is what makes the replication test meaningful --")
  say("   differences between arms are then METHOD, not different data)")
}
say("\nNEXT: re-run 03_qc_flux.R --arm bulk_after, then 04 / 05c / 07c / 09.")
