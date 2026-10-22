#!/bin/bash
#SBATCH --account=dechenli_1477
#SBATCH --job-name=magic_impute
#SBATCH --output=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/00b_magic_%A_%a.out
#SBATCH --error=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/00b_magic_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=96G
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=main
#SBATCH --array=2-9   # T is index 1, run separately via 00b_magic_t.slurm (6h; 1.5h times out)
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ttilton@usc.edu

# ==============================================================================
# 00b: MAGIC once per major type on the full normalized cell matrix. Shared by
# Run A (all_sites_bulk_after) and Run B (all_sites_bulk_first) so imputation is
# held identical across both.
#
# SOLVER: 'exact' (scFEA's default). Because 00b passes genes=<overlap> (~600
# module genes), the graph is still built on all 49107 genes but only the module
# genes are denoised, which cuts the exact-solver imputation ~80x -- exact should
# hold even on T (244K). If T OOMs or overruns, resubmit ONLY that index with
# SOLVER=approximate below (and, for uniform imputation, rerun all 9 approximate).
#
# Partition: main (CPU-only; MAGIC is CPU). At 96G largemem is no longer needed
# -- main's nodes are 185-248G (the 89G xeon-4116 can't host a 96G request, so
# SLURM just skips it). Left unconstrained to grab the first free main node;
# combos may land on mixed zen2/zen3/skylake, which is harmless for independent
# per-combo imputation. If you want all combos on ONE microarch (your original
# epyc-9354 consistency rationale), add:  #SBATCH --constraint=epyc-7513
#
# Sizing (array shares one --mem/--time, so both cover the worst case = T, 244K):
#   Benchmark: Plasma (6292 cells) did FULL-gene exact MAGIC in 15s. Ours
#   denoises only ~600 genes and takes sparse .mtx input, which graphtools
#   reduces via TruncatedSVD WITHOUT densifying the 244K x 49107 matrix. So T
#   peak is ~mmread (15-20G transient) + PCA embedding + ~600-gene output
#   (~1.2G) = ~30-40G, and runtime is dominated by I/O, not compute. 96G / 1.5h
#   is margin over that; the 8 smaller combos underuse it. CHECK T's MaxRSS on
#   the first run and trim if you want. If T somehow OOMs (older MAGIC densifying
#   during PCA), resubmit ONLY --array=1 at higher mem -- it fails fast.
# ==============================================================================
echo "Job started: $(date) on $(hostname) | array task ${SLURM_ARRAY_TASK_ID}"

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

SOLVER="exact"   # <- set to "approximate" only if a combo OOMs/overruns exact

COMBOS=(
  "T"            # 1  (244K -- the one to watch)
  "Myeloid"      # 2
  "Fibroblast"   # 3
  "Epithelial"   # 4
  "Endothelial"  # 5
  "B"            # 6
  "NK"           # 7
  "Mast"         # 8
  "Plasma"       # 9
)
IDX=$((SLURM_ARRAY_TASK_ID - 1))
COMBO="${COMBOS[$IDX]}"
COMBO_DIR="${SHARED_BASE}/${COMBO}"
echo "MAGIC combo: ${COMBO} | solver=${SOLVER}"

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

