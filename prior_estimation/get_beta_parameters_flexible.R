#!/usr/bin/R

######################################################
# this script is used to estimate parameters in beta #
# distribution of allele frequency for each category #
# variants ###########################################
######################################################
args <- commandArgs()
options(stringsAsFactors=F)

empiri.parameters<- function(ACs,ANs){
	qs <- as.numeric(ACs)/as.numeric(ANs)
	estimated.mean = mean(na.omit(qs))
	eastimated.var = var(na.omit(qs))
	A = estimated.mean
	B = eastimated.var
	w = A-1 + sqrt((A*(1-A)^3)/B)
	v = A/(1-A)*w
	params = c(v,w)
	return(params)
}
result <- data.frame()
#types <- c("frameshift_variant","splice_acceptor_variant","splice_donor_variant","stop_gained","missense_variant")
#types <-  c("frameshift_variant","splice_acceptor_variant","splice_donor_variant","stop_gained","missense_variant","exon_variant","UTR_variant","other_variant")
input_file <- args[6]
output <- args[7]
input <- read.table(input_file,header=T,stringsAsFactors=F,sep="\t")
types <- unique(input$type)
for(type in types){
  dat <- input[which(input$type==type),]
  ACs = as.numeric(na.omit(dat$AC_Adj))
  ANs = as.numeric(na.omit(dat$AN_Adj))
  # get the minor allele counts
  minor.indx = which(ACs>ANs-ACs)
  ACs[minor.indx] = ANs[minor.indx] - ACs[minor.indx]
  
  estimated.params <- empiri.parameters(ACs,ANs)
  tmp <- data.frame(type = type, v=estimated.params[1],w = estimated.params[2])
  if(dim(result)[1]==0){
    result <- tmp
  }else{
    result <- rbind(result,tmp)
  }
}
#for(type in types){
#        dat <- read.table(paste0(input_predix,type,".txt"),header=T)
#	#dat <- read.table(paste0("/ysm-gpfs/scratch60/wl382/WES/Data/ExAC/ExAC.r0.3.1.sites.vep.canonical_",type,".txt"),header=T)
#	ACs = as.numeric(na.omit(dat$AC_Adj))
#	ANs = as.numeric(na.omit(dat$AN_Adj))
#	# get the minor allele counts 
#	minor.indx = which(ACs>ANs-ACs)
#	ACs[minor.indx] = ANs[minor.indx] - ACs[minor.indx]

#	estimated.params <- empiri.parameters(ACs,ANs)
#	tmp <- data.frame(type = type, v=estimated.params[1],w = estimated.params[2])
#	if(dim(result)[1]==0){
#		result <- tmp
#	}else{
#		result <- rbind(result,tmp)
#	}
#}

write.table(result,file = output,col.names=T,row.names=F,quote=F,sep="\t")
