#!/usr/bin/env bash
# run_all.sh -- the local pipeline, in order.
#
#   export SCFEA_ROOT=/path/to/scfea_data
#   bash run_all.sh              # primary arm (bulk_before)
#   bash run_all.sh bulk_after   # robustness arm
#   bash run_all.sh --from 07b   # resume at a step
#   bash run_all.sh --list       # print the order and stop
#
# Stage 0 and Stage 1 run on a cluster and are NOT invoked here; see cloud_cluster/ and
# the README. This script starts from assembled scFEA flux output and carries
# through to the figures and tables.
#
# Each step is run with Rscript in its own process, so a failure stops the run
# rather than leaving a half-updated session.
set -o pipefail

ARM="bulk_before"
FROM=""
LIST=0
while [ $# -gt 0 ]; do
  case "$1" in
    --from)  FROM="$2"; shift 2 ;;
    --list)  LIST=1; shift ;;
    bulk_before|bulk_after) ARM="$1"; shift ;;
    *) echo "unknown argument: $1"; exit 1 ;;
  esac
done

cd "$(dirname "$0")" || exit 1

# step | script | pass --arm? | note
STEPS=(
  "01|01_build_combined_wide.R|yes|assemble flux into cw_finaltype.rds"
  "01d|01d_floor_and_mark_comparable.R|yes|floor negative flux, mark comparable modules"
  "01e|01e_fix_unit_duplicates.R|bulk_after|canonical units; bulk_after only"
  "03|03_qc_flux.R|yes|usable_modules.rds"
  "04|04_specialization.R|yes|specialization residuals (feeds 07b)"
  "05|05_progression.R|yes|progression_all.csv"
  "05c|05c_progression_meta.R|yes|meta_all.csv"
  "06|06_hpv.R|yes|hpv_all.csv"
  "07b|07b_grounding_scalefree.R|yes|flux-enzyme correlations + null pool"
  "07c|07c_grounding_final.R|yes|module_grounding_final.csv (single-draw null)"
  "07e|07e_null_max_matched.R|yes|matched null: 99.9th pctl of max over k genes"
  "07g|07g_patch_matched_threshold.R|yes|REQUIRED: writes matched thresholds into 07c output"
  "deg|make_progression_deg_edger.R|no|edgeR differential expression"
  "30|30_tcga_signature.R|no|TCGA survival"
  "33|33_tcga_deconvolution.R|no|TCGA adjusted for myeloid fraction"
  "32a|32a_build_cillo.R|no|GSE139324"
  "32b|32b_build_kurten.R|no|GSE164690"
  "32|32_independent_cohort.R|no|isoform specificity, three cohorts"
  "fig1|fig1_atlas_overview.R|no|Figure 1"
  "fig2|fig2_allocation_panels.R|no|Figure 2"
  "fig3|fig3_reorganisation_panels.R|no|Figure 3"
  "fig4|fig4_progression_all_panels.R|no|Figure 4"
  "fig5|fig5_transporter_all_panels.R|no|Figure 5"
  "fig6|fig6_clinical_panels.R|no|Figure 6"
  "figS|figS_hypoxia.R|no|Supplementary figure"
  "S12|make_supplementary_tables.R|no|Tables S1, S2"
  "S3|table_S3_survivors.R|no|Table S3"
  "S4|make_table_S4_hpv.R|no|Table S4"
  "rank|rank_modules_by_effect.R|no|module ranking by effect size"
)

if [ "$LIST" = 1 ]; then
  printf '%-6s %-36s %s\n' STEP SCRIPT NOTE
  for s in "${STEPS[@]}"; do
    IFS='|' read -r id scr arm note <<< "$s"
    printf '%-6s %-36s %s\n' "$id" "$scr" "$note"
  done
  exit 0
fi

# only a real run needs the data directory
: "${SCFEA_ROOT:?set SCFEA_ROOT to the data directory (see README)}"

echo "arm: $ARM | SCFEA_ROOT: $SCFEA_ROOT"
echo

started=0
[ -z "$FROM" ] && started=1

for s in "${STEPS[@]}"; do
  IFS='|' read -r id scr arm note <<< "$s"

  [ "$started" = 0 ] && [ "$id" = "$FROM" ] && started=1
  [ "$started" = 0 ] && continue

  # 31_cptac_protein.py is Python and needs its own environment -- not run here
  if [ ! -f "$scr" ]; then echo "MISSING  $scr"; exit 1; fi

  # steps that only apply to one arm
  if [ "$arm" = "bulk_after" ] && [ "$ARM" != "bulk_after" ]; then
    printf '  skip   %-36s (bulk_after only)\n' "$scr"; continue
  fi

  printf '>> %-6s %s\n' "$id" "$scr"
  if [ "$arm" = "no" ]; then
    Rscript "$scr" || { echo "FAILED at $id ($scr)"; exit 1; }
  else
    Rscript "$scr" --arm "$ARM" || { echo "FAILED at $id ($scr)"; exit 1; }
  fi
done

cat <<'EOF'

Done.

Not run by this script:
  cloud_cluster/              Stages 0-1 (MAGIC + scFEA) -- submit on the cluster first
  31_cptac_protein.py   needs the cptac 1.5.14 environment:
                          conda activate cptac_env && python 31_cptac_protein.py
                        run it before fig5, which reads its three CSVs
  supporting_analyses/  diagnostics and checks; see README
EOF
