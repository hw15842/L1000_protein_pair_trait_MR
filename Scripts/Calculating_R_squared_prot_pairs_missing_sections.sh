#!/bin/bash

#SBATCH --job-name=rsq
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --array=2,6,31,44,45,46,83,90,91,92,136,139,145,149,151,192,231,245,285,286,289,300,340,342,347,348,359,362,369,370,373,379,381,385,396,397,398,401,412,426,427,444,446,458%1
#SBATCH --mem=10000M
#SBATCH --output=slurm-%A_%a.out
#SBATCH --partition=serial


echo "Running on ${HOSTNAME}"
module add languages/r/3.6.2
Rscript Calculating_R_squared_prot_pairs.R /mnt/storage/scratch/hw15842/repo/L1000_protein_pair_trait_MR/Data/ /mnt/storage/scratch/hw15842/protein_interactions/network_comparisons/Data/ --df_num=${SLURM_ARRAY_TASK_ID}



