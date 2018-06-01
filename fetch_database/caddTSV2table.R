#!/usr/bin/R

##################################################
# this script is used to convert a CADD tsv file #
# to a table with column names ###################
##################################################

args <- commandArgs()
cadd_file <- args[6] # input file for format converting
output_file <- args[7]
cadd <- read.table(cadd_file,header=F,stringsAsFactors=F,sep="\t")
cadd_interested <- cadd[,c(1,2,3,5,111,112,113,114,115,116)]
colnames(cadd_interested) <- c("CHROM","POS","REF","ALT","PolyPhenCat","PolyPhenVal","SIFTcat","SIFTval","RawScore","PHRED")
write.table(cadd_interested,file= output_file,col.names=T,row.names=F,sep="\t",quote=F)


