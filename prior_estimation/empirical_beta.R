#!/usr/bin/R

########################################################
## this script is used to categorize variants in ExAC ##
########################################################
options(stingsAsFactos=F)
args <- commandArgs()

categorize_varaints <- function(tb,type){
    ind <- which(dat$Consequence==type)
    return(tb[ind,])
}


input <- args[6]
file <- gzfile(input,"r")
#file <- gzfile("/ysm-gpfs/scratch60/wl382/WES/Data/ExAC/ExAC.r0.3.1.sites.vep.canonical.table.gz","r")
dat <- read.table(file,header=T,sep="\t")
### first to filter varaints by PASS ###
dat.passed <- dat[dat$FILTER=="PASS",]

result_predix <- args[7]
#result_predix <- "/ysm-gpfs/scratch60/wl382/WES/Data/ExAC/ExAC.r0.3.1.sites.vep.canonical"
types <- c("frameshift_variant","splice_acceptor_variant","splice_donor_variant","stop_gained","missense_variant")
#types <- c("missense_variant")
for(type in types){
	dat.tmp <- categorize_varaints(dat.passed,type)
	# extract AC and AN #
        dat.1 <- dat.tmp[,c('CHROM','POS','REF','ALT','ID','AN','AN_Adj','AC','AC_Adj')]
	result_file <- paste0(result_predix,"_",type,".txt")
	write.table(dat.1,result_file,col.names=T,row.names=F,quote=F,sep="\t")	
}
