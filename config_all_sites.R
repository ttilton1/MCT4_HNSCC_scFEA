#!/usr/bin/env Rscript
# ==============================================================================
# config_all_sites.R -- shared config for the all_sites analysis pipeline.
# Sourced by every downstream script. Set ARM once; everything else follows.
#
#   ARM = "bulk_before"  pseudobulk expression -> scFEA   ( flux(mean) )
#   ARM = "bulk_after"   scFEA per cell -> pseudobulk     ( mean(flux) )
#
# KEY STRUCTURAL FACTS (differ from the old pooled/site-separated pipeline):
#  1. One model per MAJOR TYPE across ALL sites. So flux IS comparable across
#     Oral/LN/MR within a major type -- the old "Oral-only / within-model"
#     guardrail is obsolete. The guardrail that REMAINS: never compare flux
#     across major types (different trained networks).
#  2. bulk_before's native unit is (sample x finaltype). Subtype and major-type
#     levels are ROLLED UP from finaltype units (n_cells-weighted), identically
#     in both arms, so the arm comparison is not confounded by aggregation.
#  3. Effect floors from the old run (abs 1e-3, spec 0.01) were tuned to a
#     different flux scale and do NOT transfer. Absolute floors here are
#     DATA-DERIVED per arm (see derive_floor). rel_floor (fraction of module
#     flux) is scale-free and does transfer -- prefer it, esp. for HPV.
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(data.table); library(tibble)
})

# ---- ARM: override before sourcing, or pass --arm on the command line --------
if (!exists("ARM")) {
  .a <- commandArgs(trailingOnly = TRUE)
  ARM <- if (length(.a) && any(.a == "--arm")) .a[which(.a == "--arm") + 1] else "bulk_before"
}
stopifnot(ARM %in% c("bulk_before", "bulk_after"))

ROOT <- "/Users/Tommy/dechenli_lab_local/scfea26"

PATHS <- list(
  flux_dir   = file.path(ROOT, "output",   "all_sites", ARM),   # <Type>/ flux
  meta_dir   = file.path(ROOT, "metadata", "all_sites", ARM),   # <Type>/ metadata
  cw_dir     = file.path(ROOT, "data/combined_wide/all_sites", ARM),
  results    = file.path(ROOT, "results",  "all_sites", ARM),
  module_ref = file.path(ROOT, "module_ref.rds"),
  seurat     = file.path(ROOT, "data", "obj.final.2.0.rda")
)
invisible(lapply(PATHS[c("cw_dir","results")], dir.create, recursive = TRUE, showWarnings = FALSE))

TYPES <- c("T","Myeloid","Fibroblast","Epithelial","Endothelial","B","NK","Mast","Plasma")

# ---- UNIT ---------------------------------------------------------------------
# patient x tissue x finaltype. tissue (not site) because tissue is the finer,
# nested label (each tissue sits under exactly one site) and keeps M vs R and the
# normal subtypes separate. patient (not sample) avoids pseudoreplication when a
# patient contributes >1 sample of the same tissue.
UNIT_KEYS <- c("patient", "tissue", "finaltype")
COVARS    <- c("patient","tissue","site","site_group","hpv","study","annot_2nd","annot_4th")

# ---- GATES (SET FROM 02_gate_sweep.R on bulk_before) ---------------------------
# gate=10 chosen empirically: identical power to gate=8 (13 progression combos,
# 9 HPV-testable major types) while removing 350 more degenerate units, so it
# strictly dominates. gate=15 is where power breaks (loses 2 progression combos;
# Mast falls to 3 HPV+ patients). ONE gate for all analyses and BOTH arms, so the
# arm comparison is not confounded by different filtering.
#   Sensitivity to run on any surviving hit: 15 and 20.
#   KNOWN: NK is thin in bulk_before (median 2 cells/unit, 71% of units < 8 cells)
#          -> ~3/4 of NK units drop at gate 10. Treat NK from this arm as weak;
#          bulk_after does not have this problem (per-cell flux, no unit to fatten).
GATES <- list(
  min_cells             = 10,   # primary gate, all analyses
  min_cells_sens        = c(15, 20),
  min_patients          = 3,    # patients per group for a test to run
  min_per_group_hpv     = 5,    # patients per HPV arm
  min_units_dc          = 3     # units in a (patient x major type) to double-center
)

# ---- STATISTICS ---------------------------------------------------------------
FDR_THRESH <- 0.10
REL_FLOOR  <- 0.10     # |effect| >= 10% of the module's flux  (SCALE-FREE: transfers)
FLOOR_PCTL <- 0.90     # data-derived absolute floor = this pctl of |effect| distribution

# derive an absolute effect floor from the data itself (arm-specific scale)
derive_floor <- function(effects, pctl = FLOOR_PCTL) {
  e <- abs(effects[is.finite(effects)])
  if (!length(e)) return(0)
  unname(quantile(e, pctl, na.rm = TRUE))
}

# ---- PROGRESSION AXIS (within-cohort; never pool studies for this) -------------
# ORAL MUCOSA ONLY. Critical: LN-normal has site = "Distal_normal" and therefore
# site_group = "Oral", but it is a normal LYMPH NODE, not an early stage of oral
# carcinogenesis. Folding it into "Normal" would inject LN units (which carry a
# large global flux offset -- see QC Q2) into the Normal arm and manufacture a
# spurious trend. It is EXCLUDED. NL-contact (tumor-adjacent) is also excluded
# from the primary axis -- it is field-cancerized tissue, not distal normal --
# and is available as a separate check.
STAGE_MAP <- c(
  NL  = "Normal",                                   # distal normal oral mucosa
  Pre = "Pre", LP = "Pre",                          # premalignant
  CA  = "CA", `CA-E` = "CA", `CA-A` = "CA", `CA-c` = "CA", `CA-p` = "CA"
)
STAGE_LEVELS   <- c("Normal", "Pre", "CA")
STAGE_EXCLUDED <- c("LN-normal", "NL-contact", "LN", "LN-in", "LN-out", "M", "R")
PROG_COHORTS   <- c("He", "Kim")           # primary, replication
HE_EXTENSION   <- c("CA-E", "CA-A")        # He-only staged-cancer extension

# ---- module reference ---------------------------------------------------------
load_modref <- function() {
  mr <- as_tibble(readRDS(PATHS$module_ref))
  list(
    sm  = mr %>% transmute(module = module_id, supermodule = supermodule_name),
    rxn = mr %>% transmute(module = module_id,
                           reaction = paste0(compound_in, "\u2192", compound_out)),
    sm_order = mr %>% distinct(supermodule_name, supermodule_id) %>%
      mutate(supermodule_id = as.integer(supermodule_id)) %>%
      arrange(supermodule_id) %>% pull(supermodule_name))
}

say <- function(...) message(sprintf(...))
mcols_of <- function(d) grep("^M_[0-9]+$", names(d), value = TRUE)

say("config: ARM = %s | unit = %s", ARM, paste(UNIT_KEYS, collapse = " x "))
