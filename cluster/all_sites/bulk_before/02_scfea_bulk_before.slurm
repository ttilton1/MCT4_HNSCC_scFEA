#!/bin/bash
#SBATCH --account=dechenli_1477
#SBATCH --job-name=bulkbefore_scfea
#SBATCH --output=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/02_bulk_before_%j.out
#SBATCH --error=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/02_bulk_before_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=main
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ttilton@usc.edu

# ==============================================================================
# bulk_before scFEA: train == predict on the SAME pseudobulk matrix per major
# type. scFEA runs its prediction block at the end of train mode, so the
# training-byproduct flux IS the result -- one row per (sample x finaltype) unit.
# No separate predict step, no subsample.
#
# --sc_imputation False: input is already MAGIC-imputed (00b) then pseudobulked,
# so no further imputation. Matrices are tiny (661 genes x a few-hundred units),
# so this is CPU and near-instant per combo; all 9 loop in one small job on main.
# ==============================================================================
echo "Job started: $(date) on $(hostname)"

mkdir -p /project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites

module purge
module load conda
module load gcc/13.3.0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /home1/ttilton/.conda/envs/scfea_37
PYENV="/home1/ttilton/.conda/envs/scfea_37"
export PATH="${PYENV}/bin:$PATH"
export PYTHONPATH="${PYENV}/lib/python3.7/site-packages:$PYTHONPATH"
export CUDA_VISIBLE_DEVICES=""    # force CPU
"${PYENV}/bin/python" -c "import torch; print('CUDA available:', torch.cuda.is_available(), '(expected False)')"

SCFEA_SRC="/project2/dechenli_1477/ttilton/spring26/scfea/scFEA_src/src"
SCFEA_DATA="/project2/dechenli_1477/ttilton/spring26/scfea/scFEA_src/data"
INPUT_BASE="/project2/dechenli_1477/ttilton/spring26/scfea/data/input/all_sites/bulk_before"
OUTPUT_BASE="/project2/dechenli_1477/ttilton/spring26/scfea/output/all_sites/bulk_before"
MODEL_BASE="/project2/dechenli_1477/ttilton/spring26/scfea/models/all_sites/bulk_before"
mkdir -p "${MODEL_BASE}" "${OUTPUT_BASE}"

COMBOS=( T Myeloid Fibroblast Epithelial Endothelial B NK Mast Plasma )

for COMBO in "${COMBOS[@]}"; do
  echo ""
  echo "---- ${COMBO} ($(date)) ----"
  INPUT_DIR="${INPUT_BASE}/${COMBO}"
  OUTPUT_DIR="${OUTPUT_BASE}/${COMBO}"
  MODEL_FILE="${MODEL_BASE}/${COMBO}_model.pt"
  INPUT_FILE="${INPUT_DIR}/train_expression.csv"
  mkdir -p "${OUTPUT_DIR}"

  if [ ! -f "${INPUT_FILE}" ]; then
    echo "  MISSING ${INPUT_FILE}; skipping"
    echo "MISSING_INPUT: ${COMBO}" >> "${OUTPUT_BASE}/issues.txt"; continue
  fi
  echo "  input: ${INPUT_FILE} ($(du -h ${INPUT_FILE} | cut -f1))"

  start=$(date +%s)
  "${PYENV}/bin/python" "${SCFEA_SRC}/scFEA.py" \
      --data_dir "${SCFEA_DATA}" \
      --input_dir "${INPUT_DIR}" \
      --test_file "train_expression.csv" \
      --moduleGene_file "module_gene_m168.csv" \
      --stoichiometry_matrix "cmMat_c70_m168.csv" \
      --cName_file "cName_c70_m168.csv" \
      --sc_imputation False \
      --train_epoch 100 \
      --mode train \
      --model_file "${MODEL_FILE}" \
      --res_dir "${OUTPUT_DIR}" || {
        ec=$?; echo "FAILED(bulk_before): ${COMBO} exit ${ec}" >> "${OUTPUT_BASE}/issues.txt"
        echo "  ${COMBO} FAILED (exit ${ec}), continuing"; continue; }
  dur=$(( $(date +%s) - start ))
  echo "  ${COMBO} done | $((dur/60))m$((dur%60))s | flux + model -> ${OUTPUT_DIR}"
done

echo ""
echo "Job finished: $(date)"
echo "Flux: output/all_sites/bulk_before/<combo>/*_module168_cell*_SCimpute_F_*.csv"
echo "Each row = one (sample x finaltype) unit. Next: local combined_wide for bulk_before."

