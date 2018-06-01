#!/usr/bin/R

##############################################################
# this script is used to merge megered database_AF file with #
# CADD score file ############################################
##############################################################

args <- commandArgs()
mutation_af_file <- args[6]
cadd_file <- args[7]
output_file <- args[8]

mutation_af <- read.table(mutation_af_file,header=T,stringsAsFactors=F,sep="\t")
cadd <- read.table(cadd_file,header=T,stringsAsFactors=F,sep="\t")
#cadd_selected <- cadd[,c(1,2,5,111,112,113,114,115,116)]
all.join <- merge(x = mutation_af, y = cadd,
                  by = c("POS","REF","ALT"),
                  all = TRUE)

write.table(all.join,file=output_file,col.names=T,row.names=F,sep="\t",quote=F)

