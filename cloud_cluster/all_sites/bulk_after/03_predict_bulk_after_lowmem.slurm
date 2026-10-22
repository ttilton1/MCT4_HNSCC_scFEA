#!/bin/bash
#SBATCH --account=dechenli_1477
#SBATCH --job-name=bulkafter_pred_lm
#SBATCH --output=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/03_bulk_after_pred_lowmem_%A_%a.out
#SBATCH --error=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/03_bulk_after_pred_lowmem_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=185G
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --partition=main
#SBATCH --array=1-8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ttilton@usc.edu

# ==============================================================================
# bulk_after PREDICT, LOW-MEM / MAIN variant -- races the largemem version.
# Identical scFEA call to 03_predict_bulk_after.slurm EXCEPT:
#   - partition=main (no largemem queue wait -- the whole point)
#   - mem=185G (fits main's nodes: the 248G epyc-7513/7542 AND the 185G
#     xeon-4116; the 89G xeon can't host 185G so SLURM skips it)
#   - writes to <COMBO>/full_run_lowmem/  (SEPARATE dir, so it never clobbers
#     the largemem run's <COMBO>/full_run/ -- whichever finishes first wins)
#
# All but T (T stays on its own 700G largemem job). 185G is the empirical
# ceiling: our predict skips the 49107-gene dense load AND in-scFEA MAGIC that
# the pooled predict paid for, so peak is just geneExprDf (168x661xN) + the
# FloatTensor copy -- ~150-180G for Myeloid (99K), comfortably under 185G.
# If a combo OOMs on main, it'll still be covered by the largemem twin.
# ==============================================================================
echo "Job started: $(date) on $(hostname) | LOWMEM array task ${SLURM_ARRAY_TASK_ID}"

module purge
module load conda
module load gcc/13.3.0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /home1/ttilton/.conda/envs/scfea_37
CONDA_ENV_PATH="/home1/ttilton/.conda/envs/scfea_37"
export CUDA_VISIBLE_DEVICES=""   # force CPU

SCFEA_SRC="/project2/dechenli_1477/ttilton/spring26/scfea/scFEA_src/src"
SCFEA_DATA="/project2/dechenli_1477/ttilton/spring26/scfea/scFEA_src/data"
SHARED_BASE="/project2/dechenli_1477/ttilton/spring26/scfea/data/input/all_sites/shared"
OUTPUT_BASE="/project2/dechenli_1477/ttilton/spring26/scfea/output/all_sites/bulk_after"
MODEL_BASE="/project2/dechenli_1477/ttilton/spring26/scfea/models/all_sites/bulk_after"

# 8 combos (T handled separately), largest first
COMBOS=( Myeloid Fibroblast Epithelial Endothelial B NK Mast Plasma )
IDX=$((SLURM_ARRAY_TASK_ID - 1))
if [ ${IDX} -ge ${#COMBOS[@]} ]; then echo "Index out of range; nothing to do."; exit 0; fi
COMBO="${COMBOS[$IDX]}"
echo "Predicting combo: ${COMBO} (LOWMEM, full cohort, CPU main)"

INPUT_DIR="${SHARED_BASE}/${COMBO}"                      # full imputed matrix (wide; Python reads it fine)
OUTPUT_DIR="${OUTPUT_BASE}/${COMBO}/full_run_lowmem"     # SEPARATE from largemem's full_run/
MODEL_FILE="${MODEL_BASE}/${COMBO}_model.pt"
INPUT_FILE="${INPUT_DIR}/imputed_expression.csv"
mkdir -p "${OUTPUT_DIR}"

if [ ! -f "${INPUT_FILE}" ]; then echo "ERROR: missing ${INPUT_FILE}"; exit 1; fi
if [ ! -f "${MODEL_FILE}" ]; then echo "ERROR: missing ${MODEL_FILE}"; exit 1; fi
echo "Input: ${INPUT_FILE} ($(du -h ${INPUT_FILE} | cut -f1)) | Model: ${MODEL_FILE}"

start=$(date +%s)
"${CONDA_ENV_PATH}/bin/python" "${SCFEA_SRC}/scFEA.py" \
    --data_dir "${SCFEA_DATA}" \
    --input_dir "${INPUT_DIR}" \
    --test_file "imputed_expression.csv" \
    --moduleGene_file "module_gene_m168.csv" \
    --stoichiometry_matrix "cmMat_c70_m168.csv" \
    --cName_file "cName_c70_m168.csv" \
    --sc_imputation False \
    --output_flux_file "${OUTPUT_DIR}/flux.csv" \
    --output_balance_file "${OUTPUT_DIR}/balance.csv" \
    --mode predict \
    --model_file "${MODEL_FILE}" \
    --res_dir "${OUTPUT_DIR}"
ec=$?; dur=$(( $(date +%s) - start ))
echo "Combo ${COMBO} | exit ${ec} | $((dur/60))m$((dur%60))s | -> ${OUTPUT_DIR}"
echo "Job finished: $(date)"
if [ ${ec} -ne 0 ]; then echo "FAILED(predict-lowmem): ${COMBO}" >> "${OUTPUT_BASE}/failed_predict_lowmem.txt"; exit 1; fi
echo "SUCCESS: ${COMBO} (lowmem)"

