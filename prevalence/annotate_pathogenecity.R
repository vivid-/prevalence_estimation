#!/usr/bin/R

###########################################################################
######## this script is used to identify possible pathogenic SNPs  ########
###########################################################################
args <- commandArgs()

######################################################################
## find variants who have been annotated as pathogenic in X (argument)
## mutation databases
######################################################################
variant.annotated <- function(df, idc.mutation.db,databases,X){
  ### find the indices of variants who have been annotated in different
  ### mutation databases
  annotated.times <- rep(0,dim(df)[1])
  for(i in 1:length(idc.mutation.db)){
    
      #idc.tmp <- grep("Pathogenic",df[,idc.mutation.db[i]],ignore.case = T)
      idc.tmp <- which(df[,idc.mutation.db[i]] == "Pathogenic")
      if(X==3){
        idc.tmp.1 <- grep("likely Pathogenic",df[,idc.mutation.db[i]],ignore.case = T)
        idc.tmp <- setdiff(idc.tmp,idc.tmp.1)
      }
      
      idc.no <- grep("conflicting",df[,idc.mutation.db[i]],ignore.case = T)
      
      idc.tmp <- setdiff(idc.tmp,idc.no)
      if(databases[i] =="HGMD.pathogenicity"){
        tmp.idc = which(df[,idc.mutation.db[i]]!="-")        
        idc.tmp = union(idc.tmp,tmp.idc)
      }
      if(databases[i] == "UMD-Predictor.prediction"){
        idx.info <- which(colnames(df) == "UMD.Score")
        tmp.idc = which(as.numeric(df[,idx.info])>80)
        idc.tmp = intersect(idc.tmp,tmp.idc)
      }
      
    annotated.times[idc.tmp] <- annotated.times[idc.tmp]+1
  }
  
  
  return(which(annotated.times>=X))
}

######################################################################
## get rid of variants whose frequency is larger than 0.0005 and no ##
## pathogenicity in any mutation database ##
######################################################################
filter.variant <- function(df,databases,alle.frq.coln){
  ### get the index for the allele frequency column ###
  idx.alle.frq = which(colnames(df)==alle.frq.coln)
  
  ### get the indicies for mutation databases ###
  idc.mutation.db = match(databases,colnames(df))
  
  ### first get the variants whose frq is less than 0.0005 ###
  idc.rare = which(df[,idx.alle.frq]<=0.0005)
  # df.rare = df[which(df[,idx.alle.frq]<=0.0005),]
  
  ### get the variants whose frq is larger than 0.0005 but found pathogenic in 
  ### mutation databases
  idc.common = which(df[,idx.alle.frq]>0.0005)
  
  idc.patheo = variant.annotated(df,idc.mutation.db,databases,1)
  idc.common.annotated = intersect(idc.common,idc.patheo)
  # df.annotated = df[idc.filtered,]
  # cat(83%in%idc.patheo)
  
  return(union(idc.rare,idc.common.annotated))
}

######################################################################
## annotate variants based on their allele frq and mutation database
## annotation information
######################################################################
annot.var <- function(df,databases,alle.frq.coln,annotation.coln){
  
  ### first filter variants based their frq ###
  idc.filtered <- filter.variant(df,databases,alle.frq.coln)
  df.filtered <- df[idc.filtered,]
  
  ## extract nonsense, frameshift, splicing acceptor and splicing donor ##
  annotat.idx <- which(colnames(df.filtered)==annotation.coln)
  nonss.idx <- grep("stop_gained",df.filtered[,annotat.idx])
  splc.acp.idx <- grep("splice_acceptor",df.filtered[,annotat.idx])
  splc.dn.idx <- grep("splice_donor",df.filtered[,annotat.idx])
  frms.idx <- grep("frameshift",df.filtered[,annotat.idx])
  patheo.idx.1 <- c(nonss.idx,splc.acp.idx,splc.dn.idx,frms.idx) # indices in the idc.filtered
  
  if(length(patheo.idx.1)==0){
	df.left = df.filtered
	idx.left <- idc.filtered
  }else{
  	df.left = df.filtered[-patheo.idx.1,]
 	idx.left <- idc.filtered[-patheo.idx.1]
  }
  idc.mutation.db <- match(databases,colnames(df))
  ## extract variants with more than 2 mutation databases annotated ##
  annot.idx <- variant.annotated(df.left, idc.mutation.db,databases,1)
    
  ## indices of identified patheogenic variants in the original df
  patheo.idc <- c(idx.left[annot.idx],idc.filtered[patheo.idx.1])
  
  patheogenic <- rep(0,dim(df)[1])
  patheogenic[patheo.idc] <- 1
  
  return(patheogenic)
  
}





## read the input file 
input_file <- args[6]
dat <- read.table(input_file,header=T,stringsAsFactors=F,sep="\t")
# calculated the intergrated Allele Frequency based on the data
allele.freq = dat$Allele.Count/dat$Allele.Number
# remove NA's
dat.noNa <- dat[-which(is.na(allele.freq)),]
AF = allele.freq[-which(is.na(allele.freq))]
dat.noNa <- cbind(dat.noNa,AF)
# get the annotation based on the WGS and WES data
annotation <- dat$Annotation.genome
annotation[which(is.na(annotation))] <- dat$Annotation.exome[which(is.na(annotation))]
# remove rows with NA AF
Annotation <- annotation[-which(is.na(allele.freq))]
dat.noNa <- cbind(dat.noNa,Annotation)

#### extract patheogenic variants ########
databases = c("sig.emory","sig.clinvar")
frq.coln = "AF"
annotation.coln = "Annotation"
patheo.info <- annot.var(dat.noNa,databases = databases,alle.frq.coln = frq.coln,annotation.coln = annotation.coln)
result <- cbind(dat.noNa,patheo.info)

# output
output_file <- args[7]
write.table(result,file=args[7],col.names=T,row.names=F,sep="\t",quote=F)
