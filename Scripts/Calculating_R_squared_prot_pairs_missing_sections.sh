#!/bin/bash

#SBATCH --job-name=rsq
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --array=2,31,44,90,145,151,231,340,347,359,362,370,373,379,385,401,412,426,444,446,458%1
#SBATCH --mem=10000M
#SBATCH --output=slurm-%A_%a.out
#SBATCH --partition=serial


echo "Running on ${HOSTNAME}"
module add languages/r/3.6.2
Rscript Calculating_R_squared_prot_pairs.R /mnt/storage/scratch/hw15842/repo/L1000_protein_pair_trait_MR/Data/ /mnt/storage/scratch/hw15842/protein_interactions/network_comparisons/Data/ --df_num=${SLURM_ARRAY_TASK_ID}



