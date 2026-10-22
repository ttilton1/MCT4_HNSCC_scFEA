#!/usr/bin/env Rscript
# ==============================================================================
# 01_build_combined_wide.R  --  Tier 0. Emits ONE schema for BOTH arms.
#
#   Rscript 01_build_combined_wide.R --arm bulk_before
#   Rscript 01_build_combined_wide.R --arm bulk_after
#
# bulk_before: flux is ALREADY per (sample x finaltype) unit -> read it, join
#              metadata_units.csv, roll sample -> patient (n_cells-weighted).
# bulk_after : flux is per CELL -> join cell_metadata.csv on barcode, pseudobulk
#              (mean flux) to patient x tissue x finaltype.
#
# BOTH then roll UP identically (n_cells-weighted mean of finaltype units) to
# subtype (annot_2nd) and major (annot_4th), so the arms differ ONLY in when
# pseudobulk happened relative to flux -- not in how levels were aggregated.
#
# NO GATE IS APPLIED HERE. n_cells travels with every row; gates are applied
# downstream (see 02_gate_sweep.R) so a threshold change never means a rebuild.
#
# Writes -> data/combined_wide/all_sites/<ARM>/
#             cw_finaltype.rds   cw_subtype.rds   cw_major.rds
# ==============================================================================
source("config_all_sites.R")

# ---- locate the flux file for a type (handles both naming conventions) --------
find_flux <- function(ct) {
  d <- file.path(PATHS$flux_dir, ct)
  cands <- c(file.path(d, "full_run", "flux.csv"),   # bulk_after (predict)
             file.path(d, "flux.csv"))               # pinned name
  hit <- cands[file.exists(cands)]
  if (length(hit)) return(hit[1])
  g <- list.files(d, pattern = "_SCimpute_[FT]_.*\\.csv$", full.names = TRUE,
                  recursive = TRUE)                  # bulk_before (timestamped)
  if (!length(g)) return(NA_character_)
  g[order(file.mtime(g), decreasing = TRUE)][1]      # newest
}

# n_cells-weighted mean of the module columns
wmean_modules <- function(d, keys, mcols) {
  d %>% group_by(across(all_of(keys))) %>%
    summarise(across(all_of(mcols), ~ weighted.mean(.x, n_cells, na.rm = TRUE)),
              n_cells = sum(n_cells), n_units = n(), .groups = "drop")
}

say("\n=== building combined_wide | ARM = %s ===", ARM)
ft_list <- list()

for (ct in TYPES) {
  fpath <- find_flux(ct)
  if (is.na(fpath)) { say("  [%s] no flux file -- skipping", ct); next }

  flux <- fread(fpath)
  setnames(flux, 1, "id")
  flux[, id := as.character(id)]
  mcols <- grep("^M_[0-9]+$", names(flux), value = TRUE)
  stopifnot(length(mcols) == 168)

  if (ARM == "bulk_before") {
    # flux rows = sample__finaltype units; metadata_units.csv keys on the same id
    mu <- fread(file.path(PATHS$meta_dir, ct, "metadata_units.csv"))
    setnames(mu, 1, "id"); mu[, id := as.character(id)]
    d <- merge(flux, mu, by = "id", all = FALSE)
    if (nrow(d) != nrow(flux))
      say("  [%s] WARNING: %d/%d flux rows matched metadata", ct, nrow(d), nrow(flux))
    # sample -> patient x tissue x finaltype (weighted by the cells behind each sample unit)
    ft <- wmean_modules(as_tibble(d), c(UNIT_KEYS, setdiff(COVARS, UNIT_KEYS)), mcols) %>%
      rename(n_samples = n_units)

  } else {
    # flux rows = barcodes; cell_metadata.csv (from 00a, shared/) keys on barcode
    cm_path <- file.path(PATHS$meta_dir, ct, "cell_metadata.csv")
    if (!file.exists(cm_path)) { say("  [%s] missing cell_metadata.csv -- skipping", ct); next }
    cm <- fread(cm_path)
    setnames(cm, 1, "id"); cm[, id := as.character(id)]
    if ("barcode" %in% names(cm)) cm[, barcode := NULL]        # dup of the index
    d <- merge(flux, cm, by = "id", all = FALSE)
    say("  [%s] joined %d / %d cells to metadata", ct, nrow(d), nrow(flux))
    d[, n_cells := 1L]
    ft <- wmean_modules(as_tibble(d), c(UNIT_KEYS, setdiff(COVARS, UNIT_KEYS)), mcols) %>%
      select(-n_units) %>% mutate(n_samples = NA_integer_)
  }

  say("  [%s] %s units (patient x tissue x finaltype) | n_cells med=%d min=%d max=%d",
      ct, format(nrow(ft), big.mark = ","),
      as.integer(median(ft$n_cells)), min(ft$n_cells), max(ft$n_cells))
  ft_list[[ct]] <- ft
}

cw_ft <- bind_rows(ft_list)
mcols <- mcols_of(cw_ft)
say("\nfinaltype level: %s units | %d major types | %d patients",
    format(nrow(cw_ft), big.mark=","), n_distinct(cw_ft$annot_4th), n_distinct(cw_ft$patient))

# ---- roll UP to subtype and major (identical logic in both arms) --------------
cw_st <- wmean_modules(cw_ft, c("patient","tissue","annot_2nd","annot_4th",
                                "site","site_group","hpv","study"), mcols) %>%
  rename(n_finaltypes = n_units)
cw_mj <- wmean_modules(cw_ft, c("patient","tissue","annot_4th",
                                "site","site_group","hpv","study"), mcols) %>%
  rename(n_finaltypes = n_units)
say("subtype level:  %s units", format(nrow(cw_st), big.mark=","))
say("major level:    %s units", format(nrow(cw_mj), big.mark=","))

saveRDS(cw_ft, file.path(PATHS$cw_dir, "cw_finaltype.rds"))
saveRDS(cw_st, file.path(PATHS$cw_dir, "cw_subtype.rds"))
saveRDS(cw_mj, file.path(PATHS$cw_dir, "cw_major.rds"))
say("\nwrote -> %s", PATHS$cw_dir)
say("NOTE: no cell gate applied. n_cells is carried on every row; gate downstream.")
