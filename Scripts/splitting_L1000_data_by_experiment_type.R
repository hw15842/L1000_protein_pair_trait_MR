

## splitting_L1000_data_by_experiment_type.R ##


library("data.table")
library("dplyr")
library("plyr")



prot_pairs_list_dfs <- get(load("/mnt/storage/scratch/hw15842/L1000/Results/prot_pairs_l1000_list_of_dataframes.rdata"))

# for each datafrae get the experiment type 


func1 <- function(dataframe){

	df <- prot_pairs_list_dfs[[dataframe]]
	split_df <- data.frame(do.call('rbind', strsplit(as.character(df$sig_info), ':', fixed=TRUE )))
	df$experiment_name <- split_df$X1
	return(df)
}

prot_pairs_list_dfs <- lapply(1:length(prot_pairs_list_dfs), func1)



## for each prot pair (each dataframe) which experiment has the largest absolute Z score 

func2 <- function(dataframe){

	df <- prot_pairs_list_dfs[[dataframe]]
	largest_abs_Z_row <- df[which.max(abs(df$Z)),]
	largest_abs_Z_row$num_experiments <- nrow(df)
	return(largest_abs_Z_row)
}

max_abs_Z_experiment <- lapply(1:length(prot_pairs_list_dfs), func2)
max_abs_Z_experiment <- ldply(max_abs_Z_experiment, data.table)


save(max_abs_Z_experiment, file="/mnt/storage/scratch/hw15842/repo/L1000_protein_pair_trait_MR/Data/max_abs_Z_experiment.rdata")