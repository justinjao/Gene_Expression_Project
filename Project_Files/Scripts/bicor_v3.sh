#!/bin/bash
#
#SBATCH -c 8
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=network_analysis_bicor
#SBATCH --output=bicor_v3.out
#SBATCH --time=12:00:00

module load R

Rscript --vanilla ~/Analysis/V3/Bicor/Network_Construction_bicor.R

echo "FINISHED!"
