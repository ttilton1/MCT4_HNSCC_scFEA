#!/bin/bash
#SBATCH --account=dechenli_1477
#SBATCH --job-name=shared_export
#SBATCH --output=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/00a_export_%j.out
#SBATCH --error=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/00a_export_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --partition=main
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ttilton@usc.edu

# ==============================================================================
# 00a: export the full normalized (data-layer) matrix per major type as sparse
# .mtx, plus genes/barcodes/cell_metadata (incl raw UMI), into
# data/input/all_sites/shared/<MajorType>/. This is 00b's (MAGIC) input and the
# shared upstream for both bulk_before and bulk_after.
#
# 128G on main: loads the 17GB Seurat object, extracts the data + counts layers
# once, then writes sparse .mtx per type (no dense materialization -- writeMM on
# the sparse matrix), so peak is ~obj + data_all + counts_all ~= 30-50G. main
# (no largemem, no constraint) since we never hold a dense matrix. R is CPU-only.
# ==============================================================================
echo "Job started: $(date) on $(hostname)"

mkdir -p /project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites

module load conda
source activate /project2/dechenli_1477/covarnet_analysis/spring26/envs/covarnet_r

Rscript /project2/dechenli_1477/ttilton/spring26/scfea/scripts/all_sites/shared/00a_export_normalized_sparse.R

ec=$?
echo "00a exit ${ec} | finished: $(date)"
exit ${ec}

