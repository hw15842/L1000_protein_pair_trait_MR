##############################
#### combining_rsq_data.R ####
##############################


setwd("/mnt/storage/scratch/hw15842/L1000/Data/")

library(dplyr)
library(plyr)



file_list_rsq_dat <- list.files(path="/mnt/storage/scratch/hw15842/L1000/Data/", pattern=glob2rx("rsq_dat_table_*.rdata"))

head(file_list_rsq_dat)

rsq_dat <- lapply(file_list_rsq_dat,function(x){
     df <- get(load(x))})


rsq_dat <- bind_rows(rsq_dat)

head(rsq_dat)
nrow(rsq_dat)


save(rsq_dat, file=paste0("/mnt/storage/scratch/hw15842/L1000/Data/rsq_dat.rdata"))