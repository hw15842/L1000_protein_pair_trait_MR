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
# this file is wrong 
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

df3 <- na.omit(df3)


## The Suhre data is not the same as the mrbase ids so need to change it to same format

suhre_exp <- df3[grep("one.out.gz", df3$exposure_file_name),]
suhre_exp$exposure_platform_id <- gsub("-", "_", suhre_exp$exposure_platform_id)
suhre_exp$mrbase_id <- paste0("prot-c-", suhre_exp$exposure_platform_id)

df3 <- df3[-grep("one.out.gz", df3$exposure_file_name),]

df3$mrbase_id <- df3$exposure_file_name

df3 <- rbind(df3, suhre_exp)

## remove protb12 and 24 as they dont extract the same instruments as Chris used 

df3 <- df3[-grep("prot-b-12", df3$mrbase_id),]
df3 <- df3[-grep("prot-b-24", df3$mrbase_id),]


#### Split into sections of 100,000 to save sections as taking a long time and timing out / if errors looses everything 

df4 <- df3[1:100000,]
df5 <- df3[100001:200000,]
df6 <- df3[100001:200000,]
df7 <- df3[200001:300000,]
df8 <- df3[300001:400000,]
df9 <- df3[400001:466079,]




## extract the exposures first 

exposure_ids <- as.character(unique(df3$mrbase_id))

#expdat <- lapply(exposure_ids, extract_instruments) %>% bind_rows()
#save(expdat, file="expdat.rdata")
load("expdat.rdata")  ## loading in now as script ran out of time previously 


extract_dat_func <- function(x, df_num){

	df <- df_num[x,]
	exp <- subset(expdat, expdat$id.exposure == df$mrbase_id & expdat$SNP == df$SNP)
	out <- extract_outcome_data(df$SNP, df$lower.outcome, proxies=F)
	exp_out_harm <- harmonise_data(exp, out, action=1)
	rsq_dat <- steiger_filtering(exp_out_harm)
	return(rsq_dat)
}


rsq_dat_1 <- lapply(1:nrow(df4), extract_dat_func, df_num=df4)
save(rsq_dat_1, file="rsq_dat_1.rdata")
rsq_dat_table_1 <- ldply(rsq_dat_1, data.table)
save(rsq_dat_table_1, file="rsq_dat_table_1.rdata")


rsq_dat_2 <- lapply(1:nrow(df5), extract_dat_func, df_num=df5)
save(rsq_dat_2, file="rsq_dat_2.rdata")
rsq_dat_table_2 <- ldply(rsq_dat_2, data.table)
save(rsq_dat_table_2, file="rsq_dat_table_2.rdata")


rsq_dat_3 <- lapply(1:nrow(df6), extract_dat_func, df_num=df6)
save(rsq_dat_3, file="rsq_dat_3.rdata")
rsq_dat_table_3 <- ldply(rsq_dat_3, data.table)
save(rsq_dat_table_3, file="rsq_dat_table_3.rdata")


rsq_dat_4 <- lapply(1:nrow(df7), extract_dat_func, df_num=df7)
save(rsq_dat_4, file="rsq_dat_4.rdata")
rsq_dat_table_4 <- ldply(rsq_dat_4, data.table)
save(rsq_dat_table_4, file="rsq_dat_table_4.rdata")



rsq_dat_5 <- lapply(1:nrow(df8), extract_dat_func, df_num=df8)
save(rsq_dat_5, file="rsq_dat_5.rdata")
rsq_dat_table_5 <- ldply(rsq_dat_5, data.table)
save(rsq_dat_table_5, file="rsq_dat_table_5.rdata")



rsq_dat_6 <- lapply(1:nrow(df9), extract_dat_func, df_num=df9)
save(rsq_dat_6, file="rsq_dat_6.rdata")
rsq_dat_table_6 <- ldply(rsq_dat_6, data.table)
save(rsq_dat_table_6, file="rsq_dat_table_6.rdata")










