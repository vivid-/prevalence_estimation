#!/usr/bin/R

#############################################
# this script is used to merge the mutation #
# info from two databases ###################
#############################################

args <- commandArgs()
clinvar_file <- args[6]
egl_file <- args[7]
output_file <- args[8]

options(stringsAsFactors=F)
# load mutation databases
clinvar <- read.table(clinvar_file,header=T,sep="\t")
egl <- read.table(egl_file,header=T,sep="\t")
colnames(clinvar) <- c("Chrom","POS","ID","REF","ALT","QUAL","FILTER","CLNDN.clinvar","MC.clinvar","CLINDISDB.clinvar","sig.clinvar")
colnames(egl) <- c("Order","Gene","Exon","Nucleotide_Change.emory","Protein_Change.emory","Alias_Listing.emory","sig.emory","Last_Reviewed.emory","POS","REF","ALT")

# merge two dataframes
all.join <- merge(x=clinvar, y=egl,
		  by=c("POS","REF","ALT"),
                  all=TRUE)

# output the result dataframe
write.table(all.join,file=output_file,col.names=T,row.names=F,quote=F,sep="\t")

