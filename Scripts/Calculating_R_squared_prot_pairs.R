############################################
#### Calculating_R_squared_prot_pairs.R ####
############################################


args  <-  commandArgs(trailingOnly=TRUE)
results_location <- toString(args[1])
data_location <- toString(args[2])


library(dplyr)
library(TwoSampleMR)
library(plyr)
library(data.table)

setwd(paste0(results_location)) ### saving these to the "data" folder in the git hub repository to use with the rmarkdown script 


load(paste0(data_location, "Zheng_pQTLs_with_file_names.rdata")) ## not adding this file or the significant MRs to the Data folder as sig MRs too big so just extracting them from their original locations 
load(paste0(data_location, "significant_MRs.rdata"))

#### Changing from this method from Gib as too hard to trouble shoot someone elses code


# exposure_ids <- as.character(unique(Zheng_pQTLs_with_file_names$prot_file_name))
# exposure_ids <- na.omit(exposure_ids)
# 
# outcome_ids <- as.character(unique(significant_MRs$lower.outcome))
# outcome_ids <- na.omit(outcome_ids)
# 
# 
# # expdat <- lapply(exposure_ids, extract_instruments) %>% bind_rows()
# load("expdat.rdata")
# 
# save(expdat, file="expdat.rdata")
# 
# l <- list()
# m <- list()
# for(i in 1:length(outcome_ids))
# {
#   # get list of exposures to test on outcome
#   exps <- subset(significant_MRs, id.outcome == outcome_ids[i])$id.exposure %>% unique()
#   # Get instruments for the exposures that you need for this outcome
#   snps <- subset(expdat, id.exposure  %in% exps)$SNP
  # 
#   # extract
#   l[[i]] <- extract_outcome_data(snps, outcome_ids[i])
  # 
#   # harmonise
#   e <- subset(expdat, exps %in% id.exposure)
#   m[[i]] <- harmonise_data(e, l[[i]], action=1)
# }
# 
# outdat <- bind_rows(l)
# 
# save(outdat, file="outdat.rdata")
# 
# exp_out_dat<- bind_rows(m)
# 
# save(exp_out_dat, file="exp_out_dat.rdata"
# 
# # get rsq
# 
# rsq_dat <- steiger_filtering(exp_out_dat)
# 
# save(rsq_dat, file="rsq_dat.rdata")





### My method for extracting the data 



df1 <- significant_MRs[c("exposure_platform_id", "lower.outcome","SNP")]
df2 <- Zheng_pQTLs_with_file_names[c("Platform_id", "prot_file_name")]
colnames(df2) <- c("exposure_platform_id", "exposure_file_name")

df3 <- merge(df1, df2, by="exposure_platform_id")



extract_dat_func <- function(x){

	df <- df3[x,]
	exp <- extract_instruments(df$exposure_file_name)
	out <- extract_outcome_data(df$SNP, df$lower.outcome)
	exp_out_harm <- harmonise_data(exp, out, action=1)
	rsq_dat <- steiger_filtering(exp_out_harm)
	return(rsq_dat)
}


rsq_dat <- lapply(1:nrow(df3), extract_dat_func)
rsq_dat <- ldply(rsq_dat, data.table)

save(rsq_dat, file="rsq_dat.rdata")






