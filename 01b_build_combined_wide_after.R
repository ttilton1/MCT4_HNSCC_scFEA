#!/usr/bin/env Rscript
# ==============================================================================
# 01b_build_combined_wide_after.R  --  aggregate PER-CELL flux to units
#
#   Rscript 01b_build_combined_wide_after.R --arm bulk_after
#
# bulk_before = flux(mean expression)   -- pseudobulk FIRST, then scFEA
# bulk_after  = mean(flux per cell)     -- scFEA per cell, then average  <-- THIS
#
# Produces cw_finaltype.rds / cw_subtype.rds with the SAME columns and the SAME
# unit definition as the bulk_before arm, so every downstream script (03 QC,
# 04 specialization, 05c progression meta, 06 HPV, 07* grounding, 09) runs
# unchanged with --arm bulk_after.
#
# TWO THINGS THIS CHECKS BEFORE TRUSTING ANYTHING:
#
#  1. BARCODE ALIGNMENT. Flux rownames and cell_metadata rownames are both
#     "HN1_AAACGGGAGCCCAGCT-1". If they ever silently drift apart, every result
#     downstream would be wrong but not obviously broken -- so the overlap is
#     verified explicitly per cell type and the run ABORTS on a bad join rather
#     than proceeding with misaligned cells.
#
#  2. NEGATIVE FLUX. Per-cell scFEA output contains small negative values
#     (bulk_before pseudobulk flux was strictly non-negative). The scale-free
#     SHARE transform used everywhere downstream is flux / rowSums(flux), which
#     is only meaningful for non-negative values on a common scale. This script
#     quantifies how much negative mass there is per cell type and per unit AFTER
#     averaging, so we know whether share is safe here or needs a floor.
#
# Output -> data/combined_wide/all_sites/bulk_after/cw_{finaltype,subtype}.rds
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(data.table) })
stopifnot(ARM == "bulk_after")

FLUX_DIR <- file.path(ROOT, "data/flux/all_sites/bulk_after")
META_DIR <- file.path(ROOT, "data/input/all_sites/shared")
dir.create(PATHS$cw_dir, recursive = TRUE, showWarnings = FALSE)

META_COLS <- c("annot_2nd","annot_4th","finaltype","site","site_group",
               "hpv","study","tissue","patient","sample")

neg_report <- list(); unit_list <- list()

for (ct in TYPES) {
  fp <- file.path(FLUX_DIR, ct, "full_run", "flux.csv")
  mp <- file.path(META_DIR, ct, "cell_metadata.csv")
  if (!file.exists(fp)) { say("[%s] no flux.csv -- SKIP", ct); next }
  if (!file.exists(mp)) { say("[%s] no cell_metadata.csv -- SKIP", ct); next }

  say("\n[%s] reading per-cell flux...", ct)
  F <- fread(fp)
  setnames(F, 1, "barcode"); F[, barcode := as.character(barcode)]
  mcols <- setdiff(names(F), "barcode")

  M <- fread(mp)
  setnames(M, 1, "barcode"); M[, barcode := as.character(barcode)]

  # ---- 1. BARCODE VERIFICATION (abort rather than silently misalign) --------
  ov <- length(intersect(F$barcode, M$barcode))
  say("  cells: flux %s | meta %s | overlap %s (%.1f%%)",
      format(nrow(F), big.mark=","), format(nrow(M), big.mark=","),
      format(ov, big.mark=","), 100*ov/nrow(F))
  if (ov < 0.95 * nrow(F)) {
    say("  !! barcode overlap < 95%% -- flux and metadata are NOT aligned.")
    say("     flux head: %s", paste(head(F$barcode,2), collapse=" | "))
    say("     meta head: %s", paste(head(M$barcode,2), collapse=" | "))
    stop("aborting: refusing to aggregate misaligned cells for ", ct)
  }

  # ---- 2. NEGATIVE-FLUX DIAGNOSTIC (per cell, before averaging) -------------
  X <- as.matrix(F[, ..mcols])
  n_neg  <- sum(X < 0); frac_neg <- n_neg / length(X)
  worst  <- min(X)
  # how big are the negatives relative to the typical positive value?
  med_pos <- median(X[X > 0])
  neg_report[[ct]] <- data.table(
    celltype = ct, n_cells = nrow(F),
    frac_negative = frac_neg, most_negative = worst,
    median_positive = med_pos,
    neg_magnitude_ratio = abs(worst)/med_pos)
  say("  negative flux: %.3f%% of values | most negative %.2e | median positive %.2e",
      100*frac_neg, worst, med_pos)

  D <- merge(F, M[, c("barcode", META_COLS), with = FALSE], by = "barcode")
  rm(F, M, X); gc()

  # ---- aggregate: MEAN of per-cell flux within each unit --------------------
  # unit = patient x tissue x finaltype (identical to bulk_before)
  agg <- D[, c(lapply(.SD, mean, na.rm = TRUE), .(n_cells = .N)),
           by = .(patient, tissue, finaltype, annot_2nd, annot_4th,
                  study, hpv, site, site_group),
           .SDcols = mcols]
  say("  units: %s", format(nrow(agg), big.mark=","))
  unit_list[[ct]] <- agg
  rm(D); gc()
}

CW <- rbindlist(unit_list, use.names = TRUE, fill = TRUE)
mcols <- grep("^M_", names(CW), value = TRUE)
say("\n=== combined: %s units | %d modules | %d patients ===",
    format(nrow(CW), big.mark=","), length(mcols), uniqueN(CW$patient))

# ---- negative-flux summary AFTER averaging (this is what matters) -----------
NR <- rbindlist(neg_report)
say("\n=== NEGATIVE FLUX per cell type (per-cell, before averaging) ===")
print(as.data.frame(NR[, .(celltype, n_cells,
      pct_negative = round(100*frac_negative,3),
      most_negative = signif(most_negative,3),
      median_positive = signif(median_positive,3))]), row.names = FALSE)
write_csv(NR, file.path(PATHS$cw_dir, "negative_flux_report.csv"))

U <- as.matrix(CW[, ..mcols])
tot <- rowSums(U)
say("\n=== AFTER averaging to units: is the SHARE transform safe? ===")
say("  units with any negative module mean : %d of %d (%.1f%%)",
    sum(apply(U, 1, function(r) any(r < 0))), nrow(U),
    100*mean(apply(U, 1, function(r) any(r < 0))))
say("  units with negative TOTAL flux      : %d", sum(tot <= 0))
say("  min unit total = %.4g | median unit total = %.4g", min(tot), median(tot))
say("  fraction of all module means that are negative: %.3f%%",
    100*mean(U < 0))
if (sum(tot <= 0) > 0 || mean(U < 0) > 0.05) {
  say("  !! WARNING: share = flux/total is NOT safe as-is. Negative mass is")
  say("     material. Downstream 05b/05c/06 must floor at 0 or use a different")
  say("     scale-free transform. DO NOT run them until this is decided.")
} else {
  say("  OK: negatives are negligible; the share transform behaves as in bulk_before.")
}

# ---- write both levels, same schema as bulk_before --------------------------
saveRDS(CW, file.path(PATHS$cw_dir, "cw_finaltype.rds"))
say("\n-> %s", file.path(PATHS$cw_dir, "cw_finaltype.rds"))

# subtype level: roll UP from finaltype, n_cells-weighted (matches bulk_before)
CWs <- CW[, c(lapply(.SD, function(x) weighted.mean(x, n_cells)),
              .(n_cells = sum(n_cells))),
          by = .(patient, tissue, annot_2nd, annot_4th, study, hpv, site, site_group),
          .SDcols = mcols]
saveRDS(CWs, file.path(PATHS$cw_dir, "cw_subtype.rds"))
say("-> %s  (%s units)", file.path(PATHS$cw_dir, "cw_subtype.rds"),
    format(nrow(CWs), big.mark=","))

# ---- side-by-side with bulk_before, so the arms are comparable --------------
bb <- file.path(ROOT, "data/combined_wide/all_sites/bulk_before/cw_finaltype.rds")
if (file.exists(bb)) {
  B <- readRDS(bb)
  say("\n=== ARM COMPARISON (finaltype units) ===")
  say("  bulk_before: %s units | %d patients", format(nrow(B), big.mark=","), uniqueN(B$patient))
  say("  bulk_after : %s units | %d patients", format(nrow(CW), big.mark=","), uniqueN(CW$patient))
  key <- c("patient","tissue","finaltype")
  shared <- merge(as.data.table(B)[, ..key], CW[, ..key], by = key)
  say("  shared units: %s", format(nrow(shared), big.mark=","))
  say("  (differences are expected: per-cell prediction has no gate at this stage)")
}
say("\nNEXT: run 03_qc_flux.R --arm bulk_after, then 04/05c/06/07 unchanged.")
