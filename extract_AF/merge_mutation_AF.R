#!/usr/bin/R

##########################################
# this script is used to merge a AF file #
# and mutation database files ############
##########################################

args <- commandArgs()
mutation_file <- args[6]
af_file <- args[7]
output_file <- args[8]
#cat(mutation_file_file)

options(stringsAsFactors=F)
# read the mutation database annotation file
mutation <- read.table(mutation_file,header=T,sep="\t")
# read the AF file
af <- read.table(af_file,header=T,sep="\t")
# merge the mutation database annotation file and AF file
all.join <- merge(x = mutation, y = af,
                  by = c("POS","REF","ALT"),
                  all = TRUE)

write.table(all.join,file=output_file,col.names=T,row.names=F,sep="\t",quote=F)
