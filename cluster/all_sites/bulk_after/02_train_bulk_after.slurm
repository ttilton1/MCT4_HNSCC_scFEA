#!/bin/bash
#SBATCH --account=dechenli_1477
#SBATCH --job-name=bulkafter_train
#SBATCH --output=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/02_bulk_after_train_%A_%a.out
#SBATCH --error=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/02_bulk_after_train_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=gpu
#SBATCH --gpus-per-task=v100:1
#SBATCH --array=1-9
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ttilton@usc.edu

# ==============================================================================
# bulk_after TRAIN: one model per major type on its 10K top-UMI subsample.
# Same config as your pooled 02_train (v100, 96G, 100 epochs) EXCEPT
# --sc_imputation False: the subsample is already MAGIC-imputed (00b), so scFEA
# does no further imputation. Train writes the model + a byproduct flux of the
# 10K training cells (ignored for the arm; the real per-cell flux comes from the
# predict step on the FULL matrix). res_dir = output/.../<COMBO>/ for the train
# artifacts; predict writes to <COMBO>/full_run/.
# ==============================================================================
echo "Job started: $(date) | Node: $(hostname) | array ${SLURM_ARRAY_TASK_ID}"

module purge
module load conda
module load gcc/13.3.0
module load cuda/12.6.3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /home1/ttilton/.conda/envs/scfea_37
CONDA_ENV_PATH="/home1/ttilton/.conda/envs/scfea_37"
export PATH="${CONDA_ENV_PATH}/bin:$PATH"
export PYTHONPATH="${CONDA_ENV_PATH}/lib/python3.7/site-packages:$PYTHONPATH"

SCFEA_SRC="/project2/dechenli_1477/ttilton/spring26/scfea/scFEA_src/src"
SCFEA_DATA="/project2/dechenli_1477/ttilton/spring26/scfea/scFEA_src/data"
INPUT_BASE="/project2/dechenli_1477/ttilton/spring26/scfea/data/input/all_sites/bulk_after"
OUTPUT_BASE="/project2/dechenli_1477/ttilton/spring26/scfea/output/all_sites/bulk_after"
MODEL_BASE="/project2/dechenli_1477/ttilton/spring26/scfea/models/all_sites/bulk_after"
mkdir -p "${MODEL_BASE}"

COMBOS=( T Myeloid Fibroblast Epithelial Endothelial B NK Mast Plasma )
IDX=$((SLURM_ARRAY_TASK_ID - 1))
COMBO="${COMBOS[$IDX]}"
echo "Training combo: ${COMBO}"

INPUT_DIR="${INPUT_BASE}/${COMBO}"
OUTPUT_DIR="${OUTPUT_BASE}/${COMBO}"
MODEL_FILE="${MODEL_BASE}/${COMBO}_model.pt"
INPUT_FILE="${INPUT_DIR}/train_expression.csv"
mkdir -p "${OUTPUT_DIR}"

if [ ! -f "${INPUT_FILE}" ]; then echo "ERROR: missing ${INPUT_FILE}"; exit 1; fi
echo "Input: ${INPUT_FILE} ($(du -h ${INPUT_FILE} | cut -f1))"
"${CONDA_ENV_PATH}/bin/python" -c "import torch; print('CUDA:', torch.cuda.is_available())"

start=$(date +%s)
"${CONDA_ENV_PATH}/bin/python" "${SCFEA_SRC}/scFEA.py" \
    --data_dir "${SCFEA_DATA}" \
    --input_dir "${INPUT_DIR}" \
    --test_file "train_expression.csv" \
    --moduleGene_file "module_gene_m168.csv" \
    --stoichiometry_matrix "cmMat_c70_m168.csv" \
    --cName_file "cName_c70_m168.csv" \
    --sc_imputation False \
    --output_flux_file "${OUTPUT_DIR}/flux.csv" \
    --output_balance_file "${OUTPUT_DIR}/balance.csv" \
    --train_epoch 100 \
    --mode train \
    --model_file "${MODEL_FILE}" \
    --res_dir "${OUTPUT_DIR}"
ec=$?; dur=$(( $(date +%s) - start ))
echo "Combo ${COMBO} | exit ${ec} | $((dur/60))m$((dur%60))s | model -> ${MODEL_FILE}"
echo "Job finished: $(date)"
if [ ${ec} -ne 0 ]; then echo "FAILED(train): ${COMBO}" >> "${OUTPUT_BASE}/failed_train.txt"; exit 1; fi

