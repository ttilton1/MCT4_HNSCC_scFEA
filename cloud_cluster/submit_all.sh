#!/usr/bin/env bash
# submit_all.sh -- submit Stages 0-1 on CARC in dependency order.
#
#   bash submit_all.sh            # print what WOULD be submitted, change nothing
#   bash submit_all.sh --go       # actually submit
#   bash submit_all.sh --go --arm before    # primary arm only
#   bash submit_all.sh --go --from 01       # skip imputation, start at export
#
# Each stage waits on the one before it via --dependency=afterok, so the whole
# chain can be queued at once and a failure stops everything downstream rather
# than feeding half-written input to the next step.
#
# Stage 2 onward runs locally -- see run_all.sh.
#
# NOTE ON MEMORY. The two heavy steps are scFEA predict for T in the secondary
# arm (700G, largemem: geneExprDf 168 x 661 x 244K plus a FloatTensor copy is
# ~108G each) and MAGIC for T. If largemem is queued solid,
# 03_predict_bulk_after_lowmem.slurm is the same work at 185G on main -- use one
# or the other, not both.
set -o pipefail

GO=0; ARM="both"; FROM=""
while [ $# -gt 0 ]; do
  case "$1" in
    --go)   GO=1; shift ;;
    --arm)  ARM="$2"; shift 2 ;;
    --from) FROM="$2"; shift 2 ;;
    *) echo "unknown argument: $1"; exit 1 ;;
  esac
done

cd "$(dirname "$0")/all_sites" || { echo "run this from cloud_cluster/"; exit 1; }

submit () {   # submit <label> <script> [dependency-job-id]
  local label="$1" script="$2" dep="${3:-}"
  local flags=""
  [ -n "$dep" ] && flags="--dependency=afterok:$dep"

  if [ "$GO" = 0 ]; then
    printf '  sbatch %-28s %s\n' "$flags" "$script"
    echo "DRYRUN"          # stand-in job id so the chain still prints
    return
  fi

  local out
  out=$(sbatch $flags "$script") || { echo "submit failed: $script" >&2; exit 1; }
  local id="${out##* }"
  printf '  %-22s job %s  %s\n' "$label" "$id" "$script" >&2
  echo "$id"
}

[ "$GO" = 0 ] && echo "DRY RUN -- nothing is submitted. Add --go to queue these." && echo

# ---- Stage 0: imputation ------------------------------------------------------
J00A=""; J00B=""; J00C=""
if [ -z "$FROM" ] || [ "$FROM" = "00" ]; then
  J00A=$(submit "00a export"      shared/00a_export_normalized_sparse.slurm)
  # T (array index 1) needs 6h rather than 1.5h, so it runs as its own job and
  # 00b_magic_impute.slurm covers indices 2-9. Both can run at once.
  J00B=$(submit "00b magic 2-9"   shared/00b_magic_impute.slurm "$J00A")
  J00T=$(submit "00b magic T"     shared/00b_magic_t.slurm "$J00A")
  J00C=$(submit "00c transpose"   shared/00c_transpose_imputed.slurm "$J00B:$J00T")
fi

DEP="$J00C"

# ---- Stage 1: pseudobulk + scFEA ----------------------------------------------
if [ "$ARM" = "both" ] || [ "$ARM" = "before" ]; then
  echo "  -- primary arm (bulk_before) --"
  J1B=$(submit "01 export before"  bulk_before/01_export_bulk_before.slurm "$DEP")
  J2B=$(submit "02 scfea before"   bulk_before/02_scfea_bulk_before.slurm "$J1B")
fi

if [ "$ARM" = "both" ] || [ "$ARM" = "after" ]; then
  echo "  -- secondary arm (bulk_after) --"
  J1A=$(submit "01 export after"   bulk_after/01_export_bulk_after.slurm "$DEP")
  J2A=$(submit "02 train after"    bulk_after/02_train_bulk_after.slurm "$J1A")
  J3A=$(submit "03 predict after"  bulk_after/03_predict_bulk_after.slurm "$J2A")
  J3T=$(submit "03 predict after T" bulk_after/03_predict_bulk_after_T.slurm "$J2A")
fi

echo
if [ "$GO" = 0 ]; then
  cat <<'EOF'
Nothing submitted. Re-run with --go when the order looks right.

Check first:
  - largemem queue depth, for 03_predict_bulk_after_T (700G, 12h)
  - that 00b_magic_impute.slurm says --array=2-9, not 1-9 (T is index 1 and
    runs separately; at 1.5h it times out)
EOF
else
  cat <<'EOF'
Queued. Watch with:
  squeue -u $USER
  sacct -u $USER --format=JobID,JobName%24,State,Elapsed,MaxRSS

Failures are appended to, per step:
  output/all_sites/<arm>/issues.txt, failed_train.txt, failed_predict.txt

When the chain finishes, sync the flux down and continue with run_all.sh locally.
EOF
fi
