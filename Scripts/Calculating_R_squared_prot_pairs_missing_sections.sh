#!/bin/bash

#SBATCH --job-name=rsq
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --array=2,5,6,8,9,10,11,15,18,19,24,27,28,31,32,33,38,39,41,42,43,44,45,46,50,53,54,65,68,69,81,83,84,85,86,90,91,92,97,102,105,109,112,123,125,131,133,134,135,136,137,138,139,143,145,146,148,149,151,152,155,165,166,172,176,178,181,182,184,187,188,190,191,192,194,196,198,202,203,211,214,219,221,222,228,231,237,239,240,242,243,245,250,260,261,264,266,268,272,279,281,285,286,289,292,297,298,300,301,312,318,319,321,323,324,327,331,332,333,340,341,342,343,344,347,348,350,353,359,362,366,367,369,370,373,377,379,381,385,391,395,396,397,398,399,400,401,402,404,405,406,409,412,414,417,421,422,423,424,425,426,427,435,436,437,442,444,446,448,450,453,454,458,461,463,465%1
#SBATCH --mem=10000M
#SBATCH --output=slurm-%A_%a.out
#SBATCH --partition=serial


echo "Running on ${HOSTNAME}"
module add languages/r/3.6.2
Rscript Calculating_R_squared_prot_pairs.R /mnt/storage/scratch/hw15842/repo/L1000_protein_pair_trait_MR/Data/ /mnt/storage/scratch/hw15842/protein_interactions/network_comparisons/Data/ --df_num=${SLURM_ARRAY_TASK_ID}



