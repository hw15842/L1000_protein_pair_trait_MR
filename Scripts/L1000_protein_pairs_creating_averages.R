

# Run interactively on BC4

library("data.table")
library("dplyr")
library("plyr")



## We know that in the sig_id codes, CGS001 = shRNA and all other codes are for cDNA

## we need to times all the z-scores for the shRNAs by -1 to make it comparable to the cDNA

prot_pairs_list_dfs <- get(load("prot_pairs_l1000_list_of_dataframes.rdata"))


# load(paste0(data_location,"prot_pairs_l1000_results.rdata"))
# 
# ### Make the Z- scores asjusted for shRNA
# 
# prot_pairs_l1000_results$perturb_type <- ifelse(grepl("CGS001",prot_pairs_l1000_results$sig_info), c("shRNA"), c("cDNA") )
# 
# prot_pairs_l1000_results$Z_score_adjusted <- ifelse(grepl("CGS001",prot_pairs_l1000_results$sig_info), (prot_pairs_l1000_results$Z * (-1)), prot_pairs_l1000_results$Z)



## for each dataframe (which is a protein pair) make a long dataframe for the average z-scores 

func1 <- function(dataframe){

	df <- prot_pairs_list_dfs[[dataframe]]
	df$Z_score_adjusted <- ifelse(grepl("CGS001",df$sig_info), (df$Z * (-1)), df$Z)

	df_new <- data.frame(	prot_1 = df[1,1], 
							prot_2 = df[1,2], 
							num_tests = nrow(df), 
							min = min(df$Z_score_adjusted),
							Q1 = quantile(df$Z_score_adjusted, 0.25),
							median = median(df$Z_score_adjusted),
							mean = mean(df$Z_score_adjusted),
							Q3 = quantile(df$Z_score_adjusted, 0.75),
							max = max(df$Z_score_adjusted)
							)

	row.names(df_new) <- NULL

	return(df_new)

}


L1000_prot_pair_averages <- lapply(1:length(prot_pairs_list_dfs), func1)
L1000_prot_pair_averages <- ldply(L1000_prot_pair_averages, data.table)


save(L1000_prot_pair_averages, file="L1000_prot_pair_averages.rdata")

