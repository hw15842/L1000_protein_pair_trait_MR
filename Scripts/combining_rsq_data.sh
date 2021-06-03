#!/bin/bash

#SBATCH --job-name=rsq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=12:00:00
#SBATCH --mem=100000M
#SBATCH --partition=cpu


echo "Running on ${HOSTNAME}"
module add languages/r/3.6.2
Rscript combining_rsq_data.R /mnt/storage/scratch/hw15842/repo/L1000_protein_pair_trait_MR/Data