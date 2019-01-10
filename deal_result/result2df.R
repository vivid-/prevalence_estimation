#!/usr/bin/R

###########################################################
# this script is used to deal with batches of results and #
# output a dataframe ######################################
###########################################################

args <- commandArgs()
result_dir <- args[6] # directory of the input estimated prevalence files
gene <- args[7]
output_file <- args[8]
confidence <- as.numeric(args[9])
#method <- args[10]

#input_pattern <- paste0(input_pattern,"_.*_",method)
#result_files <- list.files(result_dir,pattern=input_pattern)

result_df <- data.frame()
#result_df <- data.frame(gene="",population="",confidence_level =0, )

ethics <- c("AFR","ASJ","EAS","EUR","FIN","NFE","All")
for(ethic in ethics){
   tmp_file <- paste0(gene,"_",ethic,"_",confidence,".txt")
   test <- readLines(paste0(result_dir,tmp_file))
   tmp <- gsub("estimated prevalence:  ","",test[1])
   tmp.1 <- gsub(" ","",tmp)

   estimator <- 1e6 * as.numeric(tmp.1)
   tmp <- gsub("confidence interval with  95 % cofidence:  ","",test[2])
   lower <- as.numeric(strsplit(tmp," - ")[[1]][1]) * 1e6
   upper <- as.numeric(strsplit(tmp," - ")[[1]][2]) * 1e6

   # deal with results by direct method
   tmp_file <- paste0(gene,"_",ethic,"_direct_calculation.txt")
   test <- readLines(paste0(result_dir,tmp_file))
   tmp <- gsub("estimated prevalence using direct calculation:  ","",test[1])
   tmp.1 <- gsub(" ","",tmp)
   direct_estimator <- 1e6 * as.numeric(tmp.1)
   tmp.df <- data.frame(gene=gene,population=ethic,confidence_level =confidence,
                        prevalence_estimated_perMillon = estimator, lowerBound_CI = lower, upperBound_CI = upper,naive_estimator_perMillion=direct_estimator)
   result_df <- rbind(result_df,tmp.df)
}

write.table(result_df,file= output_file,col.names=T,row.names=F,quote=F,sep="\t")
