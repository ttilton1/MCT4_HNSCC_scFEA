#!/bin/bash
#SBATCH --account=dechenli_1477
#SBATCH --job-name=bulkafter_export
#SBATCH --output=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/01_bulk_after_export_%j.out
#SBATCH --error=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/01_bulk_after_export_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --partition=main
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ttilton@usc.edu

# ==============================================================================
# bulk_after export: build the top-UMI-within-(tissue,hpv) 10K training subsample
# per major type from the shared imputed matrix. Reads imputed CSVs (no Seurat).
# 32G covers T's 2.6 GB matrix held dense (~1.3 GB) + fread overhead. ~5-10 min.
# The predict step reads the FULL shared matrix directly, so nothing else here.
# ==============================================================================
echo "Job started: $(date) on $(hostname)"

mkdir -p /project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites

module load conda
source activate /project2/dechenli_1477/covarnet_analysis/spring26/envs/covarnet_r

Rscript /project2/dechenli_1477/ttilton/spring26/scfea/scripts/all_sites/bulk_after/01_export_bulk_after.R

ec=$?
echo "01_export_bulk_after exit ${ec} | finished: $(date)"
exit ${ec}

