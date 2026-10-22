#!/bin/bash
#SBATCH --account=dechenli_1477
#SBATCH --job-name=transpose_imputed
#SBATCH --output=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/00c_transpose_%j.out
#SBATCH --error=/project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites/00c_transpose_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --partition=main
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ttilton@usc.edu

# ==============================================================================
# Transpose the shared imputed matrices (genes x cells) to tall (cells x genes)
# so R's fread reads them reliably. Pure pandas reshape; values unchanged. 32G
# covers T's 1.3 GB dense frame + transpose. ~5-10 min for all 9.
# ==============================================================================
echo "Job started: $(date) on $(hostname)"

mkdir -p /project2/dechenli_1477/ttilton/spring26/scfea/logs/all_sites

module purge
module load conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /home1/ttilton/.conda/envs/scfea_37

python /project2/dechenli_1477/ttilton/spring26/scfea/scripts/all_sites/shared/00c_transpose_imputed.py

ec=$?
echo "00c exit ${ec} | finished: $(date)"
exit ${ec}

