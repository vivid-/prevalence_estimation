#!/usr/bin/R

#######################################################
## this script is used to extract clinVar annotation ##
## from the database ##################################
#######################################################

options(stringsAsFactors=F)
args <- commandArgs()

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
extrac.info <- function(info.col,values){
	info.split <- strsplit(info.col,";")
        df <- data.frame()
	for(i in 1:length(info.split)){
                info.tmp <- info.split[[i]]
                values.tmp <- data.frame()

 		for(j in 1:length(values)){
                        idx <- grep(paste0("^",values[j],"="),info.tmp,perl=T)
			
                        if(length(idx)== 0){
                                val.tmp <- NA
                        }else{
				#val.tmp <- rep("",length(index_pre))			
                                val.tmp <- (gsub(paste0(values[j],"="),"",info.tmp[idx]))
				val.tmp <- getmode(val.tmp)
                        }
                        if(j == 1){
                                values.tmp <- data.frame(val.tmp)
                        }else{
                                values.tmp <- cbind(values.tmp,val.tmp)
                        }
                #       values.tmp <- c(values.tmp,val.tmp)
                }
		               colnames(values.tmp) <- values
                if(i==1){
                        df <- values.tmp
                }else{

                        df <- rbind(df,values.tmp)
                }

	}
 return(df)
}

info.file <- args[6]
output <- args[7]
#info.file <- "/ysm-gpfs/scratch60/wl382/WES/Data/clinvar/clinvar_20180401_15.42696977-42704515_CANP3.vcf"
dat <- read.table(info.file,header=F,sep="\t",comment.char="@",stringsAsFactors=F)
types = c("CLNDN","MC","CLNDISDB","CLNSIG")
df.tmp <- extrac.info(dat$V8,types)
result <- data.frame(dat[,-8],df.tmp)
colnames(result) <- c("Chrom","Position","ID","Reference","Alternate","QUAL","FILTER",types)
write.table(result,file=output,col.names=T,row.names=F,quote=F,sep="\t")
