#!/bin/bash
#SBATCH --account=dechenli_1477
#SBATCH --job-name=magic_T
#SBATCH --output=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/00b_magic_T_%j.out
#SBATCH --error=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/00b_magic_T_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=main
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ttilton@usc.edu

# ==============================================================================
# 00b MAGIC for T ONLY -- rerun of the array task that hit the 1.5h wall.
# The failure was a TIMEOUT, not OOM: MaxRSS never approached 96G, and the .out
# showed PCA=283s then still in KNN at the wall. So exact is feasible on 244K,
# it's just slow. Only change vs the array: --time 1.5h -> 6h. Solver stays
# exact (matches scFEA + tasks 2-9). Memory stays 96G (proven fine). CPU/main.
#
# Writes to the SAME place the array would: overwrites nothing else, and if a
# partial imputed_expression.csv exists from a prior kill it is replaced.
# ==============================================================================
echo "Job started: $(date) on $(hostname) | T-only MAGIC rerun"

mkdir -p /project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites

module purge
module load conda
module load gcc/13.3.0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /home1/ttilton/.conda/envs/scfea_37
CONDA_ENV_PATH="/home1/ttilton/.conda/envs/scfea_37"
export CUDA_VISIBLE_DEVICES=""   # MAGIC is CPU

SCRIPTS="/project2/dechenli_1477/ttilton/spring26/scfea/scripts/all_sites/shared"
SHARED_BASE="/project2/dechenli_1477/ttilton/spring26/scfea/data/input/all_sites/shared"
MODULE_FILE="/project2/dechenli_1477/ttilton/spring26/scfea/scFEA_src/data/module_gene_m168.csv"

SOLVER="exact"
COMBO="T"
COMBO_DIR="${SHARED_BASE}/${COMBO}"
echo "MAGIC combo: ${COMBO} | solver=${SOLVER} | time=6h"

if [ ! -f "${COMBO_DIR}/normalized.mtx" ]; then
  echo "ERROR: missing ${COMBO_DIR}/normalized.mtx (run 00a first)"; exit 1
fi

start=$(date +%s)
"${CONDA_ENV_PATH}/bin/python" "${SCRIPTS}/00b_magic_impute.py" \
    --combo_dir   "${COMBO_DIR}" \
    --module_file "${MODULE_FILE}" \
    --solver      "${SOLVER}"
ec=$?; dur=$(( $(date +%s) - start ))
echo "Combo ${COMBO} | exit ${ec} | $((dur/60))m$((dur%60))s | -> ${COMBO_DIR}/imputed_expression.csv"
echo "Job finished: $(date)"
if [ ${ec} -ne 0 ]; then
  echo "FAILED(magic): ${COMBO} solver=${SOLVER}" >> "${SHARED_BASE}/failed_magic.txt"; exit 1
fi
echo "SUCCESS: ${COMBO}"

