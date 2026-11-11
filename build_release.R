#!/usr/bin/env Rscript
# build_release.R -- assemble the public data deposit.
#
#   Rscript build_release.R                # pseudobulk inputs only (approved)
#   Rscript build_release.R --with-flux    # also the scFEA flux output
#
# Takes the pseudobulk scFEA inputs, drops every unit from the He cohort
# (HRA004648, controlled access), and writes the remainder to release/.
#
# Inferred flux (--with-flux) is released for all five cohorts, He included.
# Flux is a derived model output rather than the dataset: 168 module values
# produced by MAGIC imputation and a graph neural network, aggregated over at
# least 10 cells per unit, and not invertible back to expression. The He
# expression inputs stay withheld. See FLUX_ALL below.
#
# in   data/input/all_sites/bulk_before/<Type>/train_expression.csv
#      metadata/all_sites/bulk_before/<Type>/metadata_units.csv
#      data/combined_wide/all_sites/bulk_before/cw_finaltype.rds   (--with-flux)
# out  release/pseudobulk/<Type>/train_expression.csv
#      release/pseudobulk/<Type>/metadata_units.csv
#      release/flux/cw_finaltype.csv, release/flux/cw_subtype.csv  (--with-flux)
#      release/flux/usable_modules.csv                             (--with-flux)
#      release/README.txt
#      release/unit_counts.csv
#
# Nothing is modified in place. Every He unit is dropped from pseudobulk/ before
# anything is written, and the per-type counts are printed so the exclusion is
# auditable.

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr)
})

ROOT <- Sys.getenv("SCFEA_ROOT", unset = "")
if (ROOT == "") stop("set SCFEA_ROOT to the data directory")

ARM       <- "bulk_before"
EXCLUDE   <- "He"                 # HRA004648: controlled access. Withheld from
                                  # pseudobulk/ (expression) only.
FLUX_ALL  <- TRUE                 # inferred flux is a derived model output, not
                                  # the dataset, and is released for every
                                  # cohort. Set FALSE to withhold He from flux
                                  # as well.
MIN_CELLS <- 10                   # see note below -- this gate is a disclosure
                                  # control, not just a quality filter

# WHY THE GATE IS APPLIED HERE. Units are ungated upstream: the 10-cell
# threshold is applied downstream of 01_build_combined_wide. An ungated unit
# with n_cells = 1 is a single cell's imputed profile carrying a subtype label.
# These files contain no barcodes, but a 661-gene profile can be matched back to
# its cell of origin in the public GEO data, which would recover the
# barcode-to-annotation mapping we are deliberately not releasing. Gating at 10
# removes every unit that could be matched to an individual cell, and matches
# what the paper analyses. The gate applies to the He flux rows too, so no
# released He unit is a single cell's profile.
TYPES     <- c("T","Myeloid","Fibroblast","Epithelial","Endothelial",
               "B","NK","Mast","Plasma")
WITH_FLUX <- "--with-flux" %in% commandArgs(TRUE)

OUT <- file.path(ROOT, "release")
dir.create(file.path(OUT, "pseudobulk"), recursive = TRUE, showWarnings = FALSE)
say <- function(...) message(sprintf(...))

say("building release from %s", ROOT)
say("excluding cohort from pseudobulk/: %s", EXCLUDE)
say("flux released for all cohorts: %s\n", FLUX_ALL)

counts <- list()

for (ct in TYPES) {
  ep <- file.path(ROOT, "data/input/all_sites", ARM, ct, "train_expression.csv")
  mp <- file.path(ROOT, "metadata/all_sites",   ARM, ct, "metadata_units.csv")
  if (!file.exists(ep) || !file.exists(mp)) { say("  skip %-12s (missing input)", ct); next }

  e  <- fread(ep)
  gene_col <- names(e)[1]
  mu <- as.data.frame(fread(mp))

  # the unit id column is whichever one matches the expression column names
  ids <- setdiff(names(e), gene_col)
  key <- names(mu)[which.max(sapply(mu, function(x) length(intersect(as.character(x), ids))))]
  rownames(mu) <- as.character(mu[[key]])

  if (!"study" %in% names(mu)) stop(sprintf("%s: metadata has no 'study' column", ct))

  if (!"n_cells" %in% names(mu)) stop(sprintf("%s: metadata has no 'n_cells' column", ct))

  keep_ids <- ids[ids %in% rownames(mu)]
  n_he     <- sum(mu[keep_ids, "study"] == EXCLUDE)
  keep_ids <- keep_ids[mu[keep_ids, "study"] != EXCLUDE]
  n_small  <- sum(mu[keep_ids, "n_cells"] < MIN_CELLS)
  keep_ids <- keep_ids[mu[keep_ids, "n_cells"] >= MIN_CELLS]
  n_drop   <- n_he

  if (!length(keep_ids)) { say("  skip %-12s (no units left after exclusion)", ct); next }

  dir.create(file.path(OUT, "pseudobulk", ct), recursive = TRUE, showWarnings = FALSE)
  fwrite(e[, c(gene_col, keep_ids), with = FALSE],
         file.path(OUT, "pseudobulk", ct, "train_expression.csv"))
  fwrite(mu[keep_ids, , drop = FALSE],
         file.path(OUT, "pseudobulk", ct, "metadata_units.csv"))

  say("  %-12s %4d kept | %3d He dropped | %4d under %d cells dropped",
      ct, length(keep_ids), n_he, n_small, MIN_CELLS)
  counts[[ct]] <- data.frame(celltype = ct, units_released = length(keep_ids),
                             withheld_He = n_he, withheld_small = n_small,
                             genes = nrow(e))
}

cnt <- bind_rows(counts)
write_csv(cnt, file.path(OUT, "unit_counts.csv"))
say("\npseudobulk total: %d units released, %d withheld (He + under-gate)",
    sum(cnt$units_released), sum(cnt$withheld_He) + sum(cnt$withheld_small))

# ---- optional: inferred flux for the released units --------------------------
if (WITH_FLUX) {
  say("\nflux output:")
  dir.create(file.path(OUT, "flux"), recursive = TRUE, showWarnings = FALSE)

  for (lv in c("finaltype", "subtype")) {
    p <- file.path(ROOT, "data/combined_wide/all_sites", ARM, sprintf("cw_%s.rds", lv))
    if (!file.exists(p)) { say("  skip cw_%s (not found)", lv); next }
    cw <- readRDS(p)
    if (!FLUX_ALL) cw <- cw %>% filter(study != EXCLUDE)
    n_small <- sum(cw$n_cells < MIN_CELLS)
    cw      <- cw %>% filter(n_cells >= MIN_CELLS)
    write_csv(cw, file.path(OUT, "flux", sprintf("cw_%s.csv", lv)))
    say("  cw_%-10s %5d kept (%d He) | %4d under %d cells dropped",
        lv, nrow(cw), sum(cw$study == EXCLUDE), n_small, MIN_CELLS)
  }

  # the module filter, without which the released flux cannot be filtered the
  # same way the paper filters it
  um <- file.path(ROOT, "results/all_sites", ARM, "03_qc", "usable_modules.rds")
  if (file.exists(um)) {
    u <- readRDS(um)
    write_csv(data.frame(module = u), file.path(OUT, "flux", "usable_modules.csv"))
    say("  usable_modules  %d modules", length(u))
  } else say("  usable_modules.rds not found")
}

# ---- deposit README ----------------------------------------------------------
writeLines(c(
"Pseudobulk scFEA inputs and inferred flux: MCT4-associated lactate export in",
"HNSCC myeloid cells",
"Tilton T, Park YM, Zhu H, Paolillo J, Zhang C, Lin D-C. bioRxiv 2026.",
"https://doi.org/10.64898/2026.10.13.758108",
"",
"CONTENTS",
"",
"pseudobulk/<CellType>/train_expression.csv",
"    661 scFEA metabolic genes (rows) by pseudobulk units (columns).",
"    Values are MAGIC-imputed, log-normalized means across the cells in each",
"    unit. Not counts, and not the full transcriptome.",
"",
"pseudobulk/<CellType>/metadata_units.csv",
"    One row per unit: patient, sample, study, tissue, site, HPV status,",
"    cell type annotation at three levels (annot_2nd, annot_4th, finaltype),",
"    cell count and median UMI.",
"",
"unit_counts.csv",
"    Units released and withheld per cell type, for pseudobulk/ only.",
"",
if (WITH_FLUX) paste(c(
"flux/cw_finaltype.csv, flux/cw_subtype.csv",
"    Inferred scFEA metabolic flux per unit across the 168 metabolic modules,",
"    with the same annotation columns as the metadata above. These are the flux",
"    values used in the paper, not a re-run. Running scFEA yourself on the",
"    released inputs alone will not reproduce them, because the published model",
"    was trained across all cohorts including He.",
"",
"    Released for all five cohorts, He included. Flux is a derived model output",
"    rather than the source data: module-level values produced by MAGIC",
"    imputation and a graph neural network, aggregated over at least 10 cells",
"    per unit, and not invertible back to expression.",
"",
"    Tissue labels are as the source cohorts provide them (NL, Pre, LP, CA-A,",
"    CA-E, CA, LN-in, LN-out, LN-normal, M, R). The mapping onto the four",
"    classes used in the paper, and the exclusion of lymph node, metastatic and",
"    recurrent tissue, is in config_all_sites.R in the code repository.",
"",
"flux/usable_modules.csv",
"    The modules passing quality control. Analyses in the paper are restricted",
"    to these.",
""), collapse = "\n") else "",
"UNITS",
"",
"Every unit here contains at least 10 cells. Units below that are withheld,",
"because a unit of one or two cells is close to a single-cell profile and could",
"be matched back to its cell of origin in the public source data.",
"",
"In pseudobulk/, a unit is one sample x cell subtype.",
"",
"In flux/, a unit is one patient x tissue x cell subtype. Samples from the same",
"patient and tissue are combined by n_cells-weighted mean, so there are fewer",
"rows than pseudobulk units.",
"",
"LEFT OUT OF THESE FILES",
"",
"Pseudobulk expression for the He cohort (HRA004648) is withheld. That dataset",
"is controlled access and readers should apply to the original repository",
"directly. Inferred flux for He IS included in flux/, so the paired progression",
"analyses reported in the paper can be reproduced from this deposit. Re-running",
"flux inference for He from expression requires applying to HRA004648.",
"",
"No single-cell data. These are pseudobulk aggregates of at least 10 cells each,",
"and no barcode-level annotation is included.",
"",
"SOURCE COHORTS REPRESENTED HERE",
"",
"pseudobulk/: GSE181919 (Kim), GSE198315 (Koh), GSE234933 (Pittet),",
"GSE182227 (Tirosh).",
"",
"flux/: the four above plus HRA004648 (He).",
"",
"CODE",
"",
"https://github.com/ttilton1/MCT4_HNSCC_scFEA"
), file.path(OUT, "README.txt"))

say("\n-> %s", OUT)
say("Check unit_counts.csv and README.txt before depositing.")
say("Confirm He appears in flux/ and nowhere under pseudobulk/.")
