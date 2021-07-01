

## splitting_L1000_data_by_peturbagen.R ##


prot_pairs_list_dfs <- get(load("/mnt/storage/scratch/hw15842/L1000/Results/prot_pairs_l1000_list_of_dataframes.rdata"))


TOTAL_prot_pairs_list_dfs <- get(load("/mnt/storage/scratch/hw15842/repo/L1000_protein_pair_trait_MR/Data/TOTAL_prot_pairs_l1000_list_of_dataframes.rdata"))


# Split it based on the peturbagen type to check not any weird stuff going on 

func3 <- function(dataframe, total_or_sig, shRNA_T_OR_F){

	df <- total_or_sig[[dataframe]]
	df_peturb_type <- subset(df, grepl("CGS001",df$sig_info) == shRNA_T_OR_F)

	# getting the prot names from the original df so know what the prots are if they dont have shRNA or cDNA data 
	df_peturb_type <- data.frame(	prot_1 = df[1,1], 
							prot_2 = df[1,2], 
							num_tests = nrow(df_peturb_type), 
							min = min(df_peturb_type$Z),
							Q1 = quantile(df_peturb_type$Z, 0.25),
							median = median(df_peturb_type$Z),
							mean = mean(df_peturb_type$Z),
							Q3 = quantile(df_peturb_type$Z, 0.75),
							max = max(df_peturb_type$Z)
							)

	row.names(df_peturb_type) <- NULL

	return(df_peturb_type)

}



TOTAL_shRNA_L1000_prot_pair_averages <- lapply(1:length(TOTAL_prot_pairs_list_dfs), func3, total_or_sig=TOTAL_prot_pairs_list_dfs, shRNA_T_OR_F ="TRUE")
TOTAL_shRNA_L1000_prot_pair_averages <- ldply(TOTAL_shRNA_L1000_prot_pair_averages, data.table)
save(TOTAL_shRNA_L1000_prot_pair_averages, file="/mnt/storage/scratch/hw15842/repo/L1000_protein_pair_trait_MR/Data/TOTAL_shRNA_L1000_prot_pair_averages.rdata")



TOTAL_cDNA_L1000_prot_pair_averages <- lapply(1:length(TOTAL_prot_pairs_list_dfs), func3, total_or_sig=TOTAL_prot_pairs_list_dfs, shRNA_T_OR_F ="FALSE")
TOTAL_cDNA_L1000_prot_pair_averages <- ldply(TOTAL_cDNA_L1000_prot_pair_averages, data.table)
save(TOTAL_cDNA_L1000_prot_pair_averages, file="/mnt/storage/scratch/hw15842/repo/L1000_protein_pair_trait_MR/Data/TOTAL_cDNA_L1000_prot_pair_averages.rdata")


SIG_shRNA_L1000_prot_pair_averages <- lapply(1:length(prot_pairs_list_dfs), func3, total_or_sig=prot_pairs_list_dfs, shRNA_T_OR_F ="TRUE")
SIG_shRNA_L1000_prot_pair_averages <- ldply(SIG_shRNA_L1000_prot_pair_averages, data.table)
save(SIG_shRNA_L1000_prot_pair_averages, file="/mnt/storage/scratch/hw15842/repo/L1000_protein_pair_trait_MR/Data/SIG_shRNA_L1000_prot_pair_averages.rdata")



SIG_cDNA_L1000_prot_pair_averages <- lapply(1:length(prot_pairs_list_dfs), func3, total_or_sig=prot_pairs_list_dfs, shRNA_T_OR_F ="FALSE")
SIG_cDNA_L1000_prot_pair_averages <- ldply(SIG_cDNA_L1000_prot_pair_averages, data.table)
save(SIG_cDNA_L1000_prot_pair_averages, file="/mnt/storage/scratch/hw15842/repo/L1000_protein_pair_trait_MR/Data/SIG_cDNA_L1000_prot_pair_averages.rdata")



