

#######################################################
### extracting_all_possible_prot_pairs_L1000_data.R ###
#######################################################

## module add languages/r/4.0.3 
# has to be in R version 4

args  <-  commandArgs(trailingOnly=TRUE)
data_location <- toString(args[1])



setwd(paste0(data_location))

library("cmapR")
library("data.table")
library("dplyr")
library("plyr")
library("biomaRt")
library("tidyr")


## First part is exactly the same as Tom has done - only looking at certain perturbations as well

## read in the moderated Z scores file ##
# rid = either a vector of character or integer row indices or a path to a grp file containing character row indices. Only these indices will be parsed from the file.
# not sure why rid = 1 here .....
ds <- parse.gctx(paste0(data_location, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"), rid=1)

# mat - the data matrix
# rdesc - a data.frame of row annotations, with one row per matrix row
# cdesc - a data.frame of column annotations, with one row per matrix column
# rid - a character vector of unique row identifiers
# cid - a character vector of unique column identifiers
# src - a character string indicating the source (usually a file path) of the data


## In the ds file, extract the cids section
cids <- as.data.frame(ds@cid)

# call the column names of the cid dataframe "sig_name"
colnames(cids) <- c('sig_name')

# add a column to the cids dataframe called "index" which is just the row names (1 to 473647)
cids$index <- row.names(cids)

# read in the sig_info.txt file 
# "Metadata for each signature in the Level 5 matrix (metadata for the columns in the Level 5 data matrix)""
a <- fread(paste0(data_location, "GSE92742_Broad_LINCS_sig_info.txt"))

## filter the sig_infor.txt file based on the perturbagen type 
# here tom has kept:
# 	- trt_sh.cgs' = Consensus signature from shRNAs targeting the same gene
#	- 'trt_oe'  = cDNA for overexpression of wild-type gene
b <- filter(a, pert_type=='trt_sh.cgs'|pert_type=='trt_oe')

# subset the cid dataframe based on the sig_name column IDs that are in the filtered perturbagens file 
# takes it down to 58,925 rows (from 473,647 rows)
c <- cids[which(cids$sig_name %in% b$sig_id),]

# create a list of the index (row names) from datafrme c - i.e. the row names from the subsetsed pertrubagens list 
index <- as.numeric(c$index)

# select just the "sig_id" and "pert_iname" columns from the b dataset (the filtered pertrubagens dataset)
d <- b[,c('sig_id','pert_iname')]

# merge together d (sig_ id and pert_iname columns of the subsetsed pertrubagens dataset) and c (the )
d1 <- merge(d,c,by.x='sig_id',by.y='sig_name',all.x=T)
d2 <- d1[order(d1$index),]





##load data matrix of Z scores - re load the data but this time with just the "index" columns - the ones we have decided we want after filtering on the perturbations 
# again takes it from 473,647 data points to 58,925 columns
ds <- parse.gctx(paste0(data_location, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"), cid=index)



## map entrez genes to ENSG - think this is the more useful one according to Tom (ensemble # out in his script)
# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
load(paste0(data_location, "ensembl_database_downloaded_from_biomart.R")) # loading it rather than trying to extract it from site as doesnt always work from site - however if running at a later data may need to see if anything has changed. 
f <- getBM(attributes=c('hgnc_symbol','entrezgene_id'), mart=ensembl) # The entrezgene attribute was changed to entrezgene_id when Ensembl 97 was released

# in the new ds that just has the purturbagens Tom wanted, extract the row names from the matrix 
entrez <- as.data.frame(row.names(ds@mat))
colnames(entrez) <- c("entrez_id")

# add an index column which is the row names of the entrez data frame
entrez$index <- row.names(entrez)

# merge the entrez ids from the ds and the ones we got from biomart, keeping all the ones from the ds - gives us the entrez ID, the index and the hgnc_symbol
merge_entrez <- merge(entrez,f,by.x="entrez_id","entrezgene_id",all.x=T)
# keep just the unique entrez_ids
uniq_entrez <- merge_entrez[!duplicated(merge_entrez$entrez_id),]
# order them numerically based on the index value (which is the row names of the entrz dataframe made from the row names of the ds matrix)
sort_entrez <- uniq_entrez[order(as.numeric(uniq_entrez$index)),]


## load in the protein pairs - this data set is all the possible names for the proteins and their pairs, so more than there actually is
prot_pairs <- get(load(paste0(data_location,"MR_results_SNP_only.rdata")))
load(paste0(data_location, "protein_linker_file.rdata"))






##### Create the prot on prot data the same way we created it for the significant MRs

all_prot_data <- get(load(paste0(data_location, "MR_results_SNP_only.rdata")))

head(all_prot_data)

load(paste0(data_location,"protein_linker_file.rdata"))


prot_IDs <- data.frame(do.call('rbind', strsplit(as.character(all_prot_data$exposure), ':', fixed=TRUE )))

all_prot_data <- cbind(all_prot_data, prot_IDs$X1)
names(all_prot_data)[18] <- c("exposure_prot_name")

# library(stringr)
# all_prot_data$num_of_semicolons <- str_count(all_prot_data$exposure, ";")
# table(all_prot_data$num_of_semicolons)
# so max num of names is 7 as there are max 6 semi colons so the below is fine

all_prot_data <- separate(all_prot_data, "exposure_prot_name", paste("exposure_prot_name", 1:7, sep="_"), sep=";", extra="drop") 

protein_linker_file <- separate(protein_linker_file, "protein_name", paste("protein_name", 1:7, sep="_"), sep="\\.", extra="drop")

all_prot_data <- merge(all_prot_data, protein_linker_file, by.x="lower.outcome", by.y="protein_number")

colnames(all_prot_data)[26:32]  <- paste("outcome", colnames(all_prot_data)[26:32], sep = "_")

head(all_prot_data)

all_prot_data_2 <- all_prot_data[c(18:24,26:32)]

head(all_prot_data_2)

## want to take all the ones that have an extra name, then create extra rows for those exp out associations 


extra_names <- all_prot_data_2[grep("FALSE", is.na(all_prot_data_2$exposure_prot_name_2)),]


extra_names_function <- function(expo_col, outcome_col){

  new_df <- data.frame(cbind(extra_names[,expo_col], extra_names[,outcome_col]))
  return(new_df)
}

### need to run columns 1:7 against 8:14 all inclusive i.e. 1 against 8-14, 2 against 8-14 etc 


tester_df <- NULL

for(i in 1:7){

  df_new <- mapply(extra_names_function, i, 8:14)
  tester_df <- rbind(tester_df, df_new)
}

tester_new_2 <- lapply(1:7, function(x){data.frame(tester_df[,x])})
tester_new_3 <- ldply(tester_new_2, data.table)

tester_new_4 <- tapply(as.list(tester_new_3), gl(ncol(tester_new_3)/2, 2), as.data.frame)
colnames <- c("X1", "X2")
tester_new_5 <- lapply(tester_new_4, setNames, colnames)
tester_new_6 <- ldply(tester_new_5, data.frame)
tester_new_7 <- na.omit(tester_new_6)

head(tester_new_7)

#### now need to join up the data ###

not_extra_names <- all_prot_data_2[grep("TRUE", is.na(all_prot_data_2$exposure_prot_name_2)),]
not_extra_names <- data.frame(not_extra_names$exposure_prot_name_1, not_extra_names$outcome_protein_name_1)

head(not_extra_names)

names(not_extra_names)[1:2] <- c("X1", "X2")

TOTAL_all_prot_data_prot_on_prot_names <- data.frame(rbind(not_extra_names , tester_new_7[2:3]))

head(TOTAL_all_prot_data_prot_on_prot_names)

save(TOTAL_all_prot_data_prot_on_prot_names, file="TOTAL_all_prot_on_prot_data.rdata")  ### save all protein pairs 





























## load in the protein pairs - this data set is all the possible names for the proteins and their pairs, so more than there actually is
prot_pairs <- TOTAL_all_prot_data_prot_on_prot_names
prot_pairs_reverse <- prot_pairs[c("X2", "X1")]
names(prot_pairs_reverse) <- c("X1", "X2")

prot_pairs_both_directions <- rbind(prot_pairs, prot_pairs_reverse)
prot_pairs_both_directions <- unique(prot_pairs_both_directions)


# keep the proteins in collumn 1 (this is all of them as did the reverse of the pairs as well) that are in d2 dataset in the perturbations internal name column - d2 has been subset for the purturbations Tom wanted
prot_pairs_l1000 <- subset(prot_pairs_both_directions, X1 %in% d2$pert_iname) 
# from prot_pairs_l1000, keep just the proteins in column 2 that match the entrez dataset hgnc symbol - still all proteins as did the reverse
prot_pairs_l1000 <- subset(prot_pairs_l1000, X2 %in% sort_entrez$hgnc_symbol)  


# so prot_pairs_l1000 are the protein pairs that have a hgnc symbol in the d2 dataset and the entrez dataset
 

datalist=list()
# pull out Z scores
for (i in 1:nrow(prot_pairs_l1000)) {
  cols <- which(d2$pert_iname==as.character(prot_pairs_l1000[i,1])) # this is the row of the matrix and the first protein
  row <-which(sort_entrez$hgnc_symbol==as.character(prot_pairs_l1000[i,2])) # this is the column of the matrix and the second protein
  indices <- cbind(row,cols)
  #   if (nrow(d2[cols,])==1) {
  #     paste(prot_pairs_l1000[i,],d2[cols,],ds@mat[indices])  
  #   } else {
  datalist[[i]] <- cbind(prot_pairs_l1000[i,],d2[cols,][,1],ds@mat[indices])  # cbind the row of prot_pairs_l1000, the sig_id and the corresponding z score of the col and row of the matrix
}
prot_pairs_l1000_results = do.call(rbind, datalist)


colnames <- c("prot_1", "prot_2", "sig_info","Z") 
datalist <- lapply(datalist, setNames, colnames)

colnames(prot_pairs_l1000_results) <- c("prot_1", "prot_2", "sig_info","Z")

TOTAL_prot_pairs_l1000_results <- prot_pairs_l1000_results

# save(datalist, file="prot_pairs_l1000_list_of_dataframes.rdata")

save(TOTAL_prot_pairs_l1000_results, file="TOTAL_prot_pairs_l1000_results.rdata")

# datalist gives us a list of dataframes - one dataframe for each cpg snp pair - 3,390 dataframes
# prot_pairs_l1000_results gives us dataframe of all snp-cpg pairs and their associated sig_id and the z-score for that sig_id
# sig_id = (A CMap unique identification number assigned to each signature generated from L1000 data)


# Z-score   For a data point, the number of standard deviations that point is above or below the population mean is called its Z-score. 
# In the L1000 data processing pipeline, we compute a robust z-score for each gene in each sample. 
# The reference population used to compute the median and MAD is the expression of the given gene in every other well on the plate. 
# These z-score values correspond to level 4 data



### startby taking an average of all the different perturbations for each protein pair

ave_zscore_func <- function(prot_pair_dataset){

	df <- datalist[[prot_pair_dataset]]
	z_ave <- mean(df$Z)
	df <- df[1, 1:3]
	df <- as.data.frame(cbind(df, z_ave))
	names(df) <- c("prot_1", "prot_2", "sig_info", "Z_ave")
	return(df)
}

prot_pairs_l1000_ave_zscores <- lapply(1:length(datalist), ave_zscore_func)
TOTAL_prot_pairs_l1000_ave_zscores <- ldply(prot_pairs_l1000_ave_zscores, data.table)

save(TOTAL_prot_pairs_l1000_ave_zscores, file="TOTAL_prot_pairs_l1000_ave_zscores.rdata")



