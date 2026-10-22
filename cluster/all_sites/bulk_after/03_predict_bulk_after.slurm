#!/bin/bash
#SBATCH --account=dechenli_1477
#SBATCH --job-name=bulkafter_pred
#SBATCH --output=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/03_bulk_after_pred_%A_%a.out
#SBATCH --error=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/03_bulk_after_pred_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=249G
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --partition=largemem
#SBATCH --constraint=epyc-9354
#SBATCH --array=1-8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ttilton@usc.edu

# ==============================================================================
# bulk_after PREDICT (per-cell flux), all major types EXCEPT T. Loads the saved
# model and predicts on the FULL shared imputed matrix (read directly from
# shared/<COMBO>/imputed_expression.csv -- NOT re-exported). --sc_imputation
# False (already imputed). CPU on epyc-9354 (PyTorch too old for A100/L40S;
# matches your pooled predict). Writes per-cell flux to <COMBO>/full_run/.
#
# 249G: memory driver is scFEA's geneExprDf expansion (168 modules x 661 genes
# x N cells) + the FloatTensor copy, NOT the 661-gene input. For the largest
# here (Myeloid 99K) that's ~44G + ~44G ~= 90-120G peak, so 249G is comfortable.
# This is far below your pooled predict's 700G because we skip the 49107-gene
# dense load and the in-scFEA MAGIC. T (244K) is heavier and runs separately.
# ==============================================================================
echo "Job started: $(date) on $(hostname) | array task ${SLURM_ARRAY_TASK_ID}"

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
echo "Predicting combo: ${COMBO} (full cohort, CPU epyc-9354)"

INPUT_DIR="${SHARED_BASE}/${COMBO}"                     # full imputed matrix lives here
OUTPUT_DIR="${OUTPUT_BASE}/${COMBO}/full_run"
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
if [ ${ec} -ne 0 ]; then echo "FAILED(predict): ${COMBO}" >> "${OUTPUT_BASE}/failed_predict.txt"; exit 1; fi
echo "SUCCESS: ${COMBO}"

