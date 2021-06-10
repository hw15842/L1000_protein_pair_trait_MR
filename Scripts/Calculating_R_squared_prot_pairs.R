############################################
#### Calculating_R_squared_prot_pairs.R ####
############################################


args  <-  commandArgs(trailingOnly=TRUE)
results_location <- toString(args[1])
data_location <- toString(args[2])
df_num <- toString(args[3])



library(dplyr)
library(TwoSampleMR)
library(plyr)
library(data.table)
library(ieugwasr)

setwd(paste0(results_location)) ### saving these to the "data" folder in the git hub repository to use with the rmarkdown script 


load(paste0(data_location, "Zheng_pQTLs_with_file_names.rdata")) ## not adding this file or the significant MRs to the Data folder as sig MRs too big so just extracting them from their original locations 
load(paste0(data_location, "significant_MRs.rdata"))


df_num <- unlist(regmatches(df_num, gregexpr('\\(?[0-9,.]+', df_num))) ### This stops it pasting the full "--section_of_subset=1" and just keeps the "1" 

df_num <- as.numeric(noquote(df_num))

df_num


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


df_split <- split(df3, (seq(nrow(df3))-1) %/% 1000) 


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




rsq_dat <- lapply(1:nrow(df_split[[df_num]]), extract_dat_func, df_num=df_split[[df_num]])
save(rsq_dat, file=paste0("rsq_dat_", df_num, ".rdata"))
rsq_dat_table <- ldply(rsq_dat, data.table)
save(rsq_dat_table, file=paste0("rsq_dat_table_", df_num, ".rdata"))



