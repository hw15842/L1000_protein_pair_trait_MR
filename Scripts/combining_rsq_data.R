##############################
#### combining_rsq_data.R ####
##############################


args  <-  commandArgs(trailingOnly=TRUE)
data_location <- toString(args[1])

setwd(paste0(data_location))

library(dplyr)
library(plyr)



file_list_rsq_dat <- list.files(path=paste0(data_location), pattern=glob2rx("rsq_dat_table_*.rdata"))

head(file_list_rsq_dat)

rsq_dat <- lapply(file_list_rsq_dat,function(x){
     df <- get(load(x))})


rsq_dat <- bind_rows(rsq_dat)

head(rsq_dat)
nrow(rsq_dat)


save(rsq_dat, file=paste0(data_location, "rsq_dat.rdata"))