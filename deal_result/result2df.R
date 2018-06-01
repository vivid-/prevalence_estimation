#!/usr/bin/R

###########################################################
# this script is used to deal with batches of results and #
# output a dataframe ######################################
###########################################################

args <- commandArgs()
result_dir <- args[6] # directory of the input estimated prevalence files
input_pattern <- args[7]
output_file <- args[8]
confidence <- as.numeric(args[9])

result_files <- list.files(result_dir,pattern=input_pattern)

result_df <- data.frame()
#result_df <- data.frame(gene="",population="",confidence_level =0, ) 
for(file in result_files){
	gene = strsplit(file,"_")[[1]][1]
	population = strsplit(file,"_")[[1]][2]
	
	
	test <- readLines(paste0(result_dir,file))
	tmp <- gsub("estimated prevalence:  ","",test[1])
	tmp.1 <- gsub(" ","",tmp)
	estimator <- 1e6 * as.numeric(tmp.1)
	tmp <- gsub("confidence interval with  95 % cofidence:  ","",test[2])
	lower <- as.numeric(strsplit(tmp," - ")[[1]][1]) * 1e6
	upper <- as.numeric(strsplit(tmp," - ")[[1]][2]) * 1e6

	tmp.df <- data.frame(gene=gene,population=population,confidence_level =confidence,
			     prevalence_estimated_perMillon = estimator, lowerBound_CI = lower, upperBound_CI = upper)
	if(length(result_df)==0){
		result_df <- tmp.df
	}else{
		result_df <- rbind(result_df,tmp.df)
	}

}

write.table(result_df,file= output_file,col.names=T,row.names=F,quote=F,sep="\t")


