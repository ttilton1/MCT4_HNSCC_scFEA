#!/usr/bin/env bash
# check_figure_scripts.sh - check the repo still works, without re-running the pipeline.
#
#   export SCFEA_ROOT=/Users/Tommy/dechenli_lab_local/scfea26
#   bash check_figure_scripts.sh
#
# Three tiers, cheapest first. Nothing is recomputed from raw data; the figure
# scripts read the existing intermediates and redraw.
#
#   1. PARSE    every R and Python file still parses
#   2. INPUTS   the ten files every figure depends on exist and are readable
#   3. REDRAW   re-run the figure and table scripts, compare against the
#               versions currently on disk
#
# Existing figures are copied to figures_backup/ first, so a failed
# redraw cannot destroy what the paper used.
set -o pipefail

cd "$(git rev-parse --show-toplevel 2>/dev/null)" || { echo "not in the repo"; exit 1; }
: "${SCFEA_ROOT:?set SCFEA_ROOT first}"

R="$SCFEA_ROOT/results/all_sites/bulk_before"
CW="$SCFEA_ROOT/data/combined_wide/all_sites/bulk_before"
fail=0

echo "=============================================================="
echo "1. PARSE"
echo "=============================================================="
for f in $(find . -name '*.R' -not -path './.git/*' | sort); do
  Rscript -e 'invisible(parse(commandArgs(TRUE)[1]))' "$f" >/dev/null 2>&1 \
    || { echo "  FAIL  $f"; fail=1; }
done
for f in $(find . -name '*.py' -not -path './.git/*' | sort); do
  python -m py_compile "$f" >/dev/null 2>&1 || { echo "  FAIL  $f"; fail=1; }
done
find . -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null
[ $fail -eq 0 ] && echo "  all files parse"
echo

echo "=============================================================="
echo "2. INPUTS the figures depend on"
echo "=============================================================="
check () {   # check <label> <path>
  if [ -f "$2" ]; then
    printf '  ok    %-34s %8s  %s\n' "$1" \
      "$(du -h "$2" | cut -f1)" "$(date -r "$2" '+%Y-%m-%d')"
  else
    printf '  MISS  %-34s %s\n' "$1" "$2"; fail=1
  fi
}
check "cw_finaltype.rds"            "$CW/cw_finaltype.rds"
check "usable_modules.rds"          "$R/03_qc/usable_modules.rds"
check "progression_all.csv"         "$R/05_progression/finaltype/progression_all.csv"
check "meta_all.csv"                "$R/05c_progression_meta/finaltype/meta_all.csv"
check "module_grounding_final.csv"  "$R/07c_grounding_final/module_grounding_final.csv"
check "gene_level_scalefree.csv"    "$R/07b_grounding_scalefree/gene_level_scalefree.csv"
check "hpv_all.csv"                 "$R/06_hpv/finaltype/hpv_all.csv"
check "tumor_vs_normal_protein.csv" "$SCFEA_ROOT/results/cptac/tumor_vs_normal_protein.csv"
check "cox_univariate.csv"          "$SCFEA_ROOT/results/tcga/cox_univariate.csv"
check "per_patient_scores.csv"      "$SCFEA_ROOT/results/tcga/per_patient_scores.csv"
echo

echo "  -- the matched-null patch must have run --"
if [ -f "$R/07c_grounding_final/module_grounding_final_singledraw_backup.csv" ]; then
  n=$(grep -c TRUE "$R/07c_grounding_final/module_grounding_final.csv")
  echo "  ok    backup present, $n TRUE rows (expect 44)"
  [ "$n" = "44" ] || { echo "  WARNING: expected 44 grounded rows, found $n"; fail=1; }
else
  echo "  MISS  no backup file: 07g_patch_matched_threshold.R may not have run"; fail=1
fi
echo

if [ $fail -ne 0 ]; then
  echo "Stopping before redraw -- fix the above first."
  exit 1
fi

echo "=============================================================="
echo "3. REDRAW figures and tables"
echo "=============================================================="
BK="$SCFEA_ROOT/figures_backup"
if [ -d "$SCFEA_ROOT/figures" ] && [ ! -d "$BK" ]; then
  cp -R "$SCFEA_ROOT/figures" "$BK"
  echo "  backed up figures/ -> $(basename "$BK")/"
fi
echo

for s in fig1_atlas_overview.R fig2_allocation_panels.R fig3_reorganisation_panels.R \
         fig4_progression_all_panels.R fig5_transporter_all_panels.R \
         fig6_clinical_panels.R figS_hypoxia.R \
         make_supplementary_tables.R table_S3_survivors.R make_table_S4_hpv.R \
         rank_modules_by_effect.R; do
  printf '  %-36s ' "$s"
  if Rscript "$s" > "check_logs/${s%.R}.log" 2>&1; then
    echo "ok"
  else
    echo "FAILED  (see check_logs/${s%.R}.log)"
    tail -3 "check_logs/${s%.R}.log" | sed 's/^/        /'
    fail=1
  fi
done

echo
echo "=============================================================="
if [ $fail -eq 0 ]; then
  echo "PASS. Every figure and table regenerated from the committed code."
  echo
  echo "Now compare by eye against the preprint. The logs carry the printed"
  echo "numbers, so these are worth reading:"
  echo "  grep -h 'M_78\\|0.715\\|survivor' check_logs/fig*.log"
  echo
  echo "If the figures match, delete the backup:"
  echo "  rm -rf $BK"
else
  echo "FAIL. Restore from the backup if a figure was overwritten badly:"
  echo "  rm -rf $SCFEA_ROOT/figures && mv $BK $SCFEA_ROOT/figures"
fi
exit $fail
