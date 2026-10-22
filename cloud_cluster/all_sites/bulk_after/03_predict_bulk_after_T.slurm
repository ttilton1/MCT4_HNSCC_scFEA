#!/bin/bash
#SBATCH --account=dechenli_1477
#SBATCH --job-name=bulkafter_pred_T
#SBATCH --output=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/03_bulk_after_pred_T_%j.out
#SBATCH --error=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/03_bulk_after_pred_T_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=700G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --partition=largemem
#SBATCH --constraint=epyc-9354
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ttilton@usc.edu

# ==============================================================================
# bulk_after PREDICT for T ONLY (244K cells). Same as the array but its own job
# because scFEA's geneExprDf (168 x 661 x 244K) + FloatTensor copy is ~108G each
# => ~220-250G peak, plus pandas-append overhead during construction. 700G is
# generous headroom (largemem node is 1498G); if MaxRSS comes back well under,
# you can trim on a rerun. CPU on epyc-9354, reads the full shared imputed
# matrix directly.
#
# If it ever OOMs even at 700G, the fallback is to chunk T's imputed matrix (as
# your pooled T predict did) and rbind the flux -- but with the 661-gene input
# (no 49107-gene load, no MAGIC) single-shot should hold.
# ==============================================================================
echo "Job started: $(date) on $(hostname) | T predict"

module purge
module load conda
module load gcc/13.3.0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /home1/ttilton/.conda/envs/scfea_37
CONDA_ENV_PATH="/home1/ttilton/.conda/envs/scfea_37"
export CUDA_VISIBLE_DEVICES=""

SCFEA_SRC="/project2/dechenli_1477/ttilton/spring26/scfea/scFEA_src/src"
SCFEA_DATA="/project2/dechenli_1477/ttilton/spring26/scfea/scFEA_src/data"
SHARED_BASE="/project2/dechenli_1477/ttilton/spring26/scfea/data/input/all_sites/shared"
OUTPUT_BASE="/project2/dechenli_1477/ttilton/spring26/scfea/output/all_sites/bulk_after"
MODEL_BASE="/project2/dechenli_1477/ttilton/spring26/scfea/models/all_sites/bulk_after"

COMBO="T"
INPUT_DIR="${SHARED_BASE}/${COMBO}"
OUTPUT_DIR="${OUTPUT_BASE}/${COMBO}/full_run"
MODEL_FILE="${MODEL_BASE}/${COMBO}_model.pt"
INPUT_FILE="${INPUT_DIR}/imputed_expression.csv"
mkdir -p "${OUTPUT_DIR}"

if [ ! -f "${INPUT_FILE}" ]; then echo "ERROR: missing ${INPUT_FILE}"; exit 1; fi
if [ ! -f "${MODEL_FILE}" ]; then echo "ERROR: missing ${MODEL_FILE}"; exit 1; fi
echo "Input: ${INPUT_FILE} ($(du -h ${INPUT_FILE} | cut -f1)) | mem 700G"

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
echo "T | exit ${ec} | $((dur/60))m$((dur%60))s | -> ${OUTPUT_DIR}"
echo "Job finished: $(date)"
if [ ${ec} -ne 0 ]; then echo "FAILED(predict): T" >> "${OUTPUT_BASE}/failed_predict.txt"; exit 1; fi
echo "SUCCESS: T"

