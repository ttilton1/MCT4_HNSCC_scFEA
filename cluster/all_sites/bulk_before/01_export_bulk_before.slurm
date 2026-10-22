#!/bin/bash
#SBATCH --account=dechenli_1477
#SBATCH --job-name=bulkbefore_export
#SBATCH --output=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/01_bulk_before_export_%j.out
#SBATCH --error=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/01_bulk_before_export_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --partition=main
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ttilton@usc.edu

# ==============================================================================
# bulk_before export: pseudobulk the shared imputed matrix to (sample x finaltype)
# means, one small matrix per major type. No Seurat object here -- reads the
# imputed CSVs from 00b. 32G covers T's 2.6 GB imputed matrix held dense
# (661 x 244153 ~ 1.3 GB) plus fread overhead; the other 8 are much smaller.
# CPU, main, ~10-15 min total.
# ==============================================================================
echo "Job started: $(date) on $(hostname)"

mkdir -p /project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites

module load conda
source activate /project2/dechenli_1477/covarnet_analysis/spring26/envs/covarnet_r

Rscript /project2/dechenli_1477/ttilton/spring26/scfea/scripts/all_sites/bulk_before/01_export_bulk_before.R

ec=$?
echo "01_export_bulk_before exit ${ec} | finished: $(date)"
exit ${ec}

