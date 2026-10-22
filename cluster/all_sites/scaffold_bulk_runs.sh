#!/usr/bin/env bash
# scaffold_bulk_runs.sh
# Create the directory tree for the two all-sites pooled runs.
#
#   all_sites_bulk_after   Run A: per-cell scFEA, pseudobulk AFTER flux  -> mean(flux)
#   all_sites_bulk_first   Run B: pseudobulk expression FIRST, then scFEA -> flux(mean)
#
# Both runs use the log-normalized `data` layer (normalized by construction;
# no _datalayer suffix needed). Keyed on annot_4th, all sites (Oral+LN+MR) pooled.
#
# Usage:
#   bash scaffold_bulk_runs.sh carc
#   bash scaffold_bulk_runs.sh local
#
# Idempotent (mkdir -p): safe to re-run.
set -euo pipefail

TARGET="${1:-}"
RUNS=( all_sites_bulk_after all_sites_bulk_first )
TYPES=( T Myeloid Fibroblast Epithelial Endothelial B NK Mast Plasma )

case "$TARGET" in
  carc)
    BASE="/project2/dechenli_1477/ttilton/spring26/scfea"
    for RUN in "${RUNS[@]}"; do
      mkdir -p "$BASE/logs/$RUN"
      mkdir -p "$BASE/scripts/$RUN"
      mkdir -p "$BASE/models/$RUN"
      for T in "${TYPES[@]}"; do
        mkdir -p "$BASE/data/input/$RUN/$T"
        mkdir -p "$BASE/data/metadata/$RUN/$T"
        mkdir -p "$BASE/output/$RUN/$T"
      done
    done
    ;;

  local)
    BASE="/Users/Tommy/dechenli_lab_local/scfea26"
    for RUN in "${RUNS[@]}"; do
      mkdir -p "$BASE/data/combined_wide/$RUN"
      mkdir -p "$BASE/results/$RUN"
      mkdir -p "$BASE/scripts/$RUN"
      for T in "${TYPES[@]}"; do
        mkdir -p "$BASE/output/$RUN/$T"
        mkdir -p "$BASE/metadata/$RUN/$T"
      done
    done
    ;;

  *)
    echo "Usage: bash $0 [carc|local]" >&2
    exit 1
    ;;
esac

echo "Scaffold complete for target=$TARGET (base: $BASE)"
echo "----- created tree -----"
find "$BASE" -maxdepth 3 -type d -path "*all_sites_bulk*" | sort

