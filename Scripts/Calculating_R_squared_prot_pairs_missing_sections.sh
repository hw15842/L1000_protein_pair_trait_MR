#!/bin/bash

#SBATCH --job-name=rsq
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --array=2,6,9,15,19,31,32,38,41,42,44,45,46,81,83,90,91,92,105,112,135,136,139,145,149,151,192,231,239,245,268,272,285,286,289,300,318,324,327,331,333,340,341,342,344,347,348,353,359,362,369,370,373,379,381,385,396,397,398,401,404,412,414,421,424,426,427,435,444,446,448,453,458,465%1
#SBATCH --mem=10000M
#SBATCH --output=slurm-%A_%a.out
#SBATCH --partition=serial


echo "Running on ${HOSTNAME}"
module add languages/r/3.6.2
Rscript Calculating_R_squared_prot_pairs.R /mnt/storage/scratch/hw15842/repo/L1000_protein_pair_trait_MR/Data/ /mnt/storage/scratch/hw15842/protein_interactions/network_comparisons/Data/ --df_num=${SLURM_ARRAY_TASK_ID}



