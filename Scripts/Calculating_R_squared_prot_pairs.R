############################################
#### Calculating_R_squared_prot_pairs.R ####
############################################


args  <-  commandArgs(trailingOnly=TRUE)
results_location <- toString(args[1])
data_location <- toString(args[2])


library(dplyr)
library(TwoSampleMR)

setwd(paste0(results_location) ### saving these to the "data" folder in the git hub repository to use with the rmarkdown script 

getwd()

load("/mnt/storage/scratch/hw15842/protein_interactions/network_comparisons/Data/Zheng_pQTLs_with_file_names.rdata")
load("/mnt/storage/scratch/hw15842/protein_interactions/network_comparisons/Data/significant_MRs.rdata")

load(paste0(data_location, "Zheng_pQTLs_with_file_names.rdata")) ## not adding this file or the significant MRs to the Data folder as sig MRs too big so just extracting them from their original locations 
load(paste0(data_location, "significant_MRs.rdata"))


exposure_ids <- as.character(unique(Zheng_pQTLs_with_file_names$prot_file_name))
outcome_ids <- as.character(unique(significant_MRs$lower.outcome))

expdat <- lapply(exposure_ids, extract_instruments) %>% bind_rows()

save(expdat, file="expdat.rdata")

l <- list()
m <- list()
for(i in 1:length(outcome_ids))
{
  # get list of exposures to test on outcome
  exps <- subset(significant_MRs, id.outcome == outcome_ids[i])$id.exposure %>% unique()
  # Get instruments for the exposures that you need for this outcome
  snps <- subset(expdat, id.exposure  %in% exps)$SNP
  
  # extract
  l[[i]] <- extract_outcome_data(snps, outcome_ids[i])
  
  # harmonise
  e <- subset(expdat, exps %in% id.exposure)
  m[[i]] <- harmonise_data(e, l[[i]], action=1)
}

outdat <- bind_rows(l)

save(outdat, file="outdat.rdata")

exp_out_dat<- bind_rows(m)

save(exp_out_dat, file="exp_out_dat.rdata"

# get rsq

rsq_dat <- steiger_filtering(exp_out_dat)

save(rsq_dat, file="rsq_dat.rdata")



# save dat expdat and outdat so you don't need to extract again

