#!/usr/bin/R

#################################################
# this script is used to transform VCF files to #
# readable table files ##########################
#################################################

args <- commandArgs()

# deal with the information column and get the important information
extract.infor <- function(info.col,values,gene_name){
	info.split <- strsplit(info.col,";")
	df <- data.frame()
	#for(i in length(values)){
	for(i in 1:length(info.split)){
		info.tmp <- info.split[[i]]
		values.tmp <- data.frame()
		for(j in 1:length(values)){
			if(values[j]=="Annotation"){
				#cat(i,"\n")
				idx <- grep(paste0("^","CSQ","="),info.tmp,perl=T)
				
				if(length(idx)==0){
					val.tmp <- ""
					
				}else{
					tmp <- strsplit(info.tmp[idx],",")
					test_id <- grep(gene_name,tmp[[1]])
			#		variant_idx <- grep("variant",tmp[[1]])
					tmp_annotations <- c()
					for(id_i in test_id){
						test_tmp <- strsplit(tmp[[1]][id_i],"\\|")
						tmp_annotations <- c(tmp_annotations,test_tmp[[1]][2])
					}
                                        # categorize variant annotations
					if(length(grep("frame",tmp_annotations))!=0){
						val.tmp <- "frameshift_variant"
}else if (length(grep("frame",tmp_annotations))!=0){
						val.tmp <- "missense_variant"
}else if (length(grep("splice_acceptor_variant",tmp_annotations))!=0){
                                                val.tmp <- "splice_acceptor_variant"
}else if (length(grep("splice_donor_variant",tmp_annotations))!=0){
                                                val.tmp <- "splice_donor_variant"
}else if (length(grep("stop_gained",tmp_annotations))!=0){
                                                val.tmp <- "stop_gained"
}else{
						val.tmp <- tmp_annotations[1]
}
					
				}
			}else{
				idx <- grep(paste0("^",values[j],"="),info.tmp,perl=T)
				val.tmp <- as.numeric(gsub(paste0(values[j],"="),"",info.tmp[idx]))
			}
			if(j == 1){
				values.tmp <- data.frame(val.tmp)
			}else{
				values.tmp <- cbind(values.tmp,val.tmp)
			}
		#	values.tmp <- c(values.tmp,val.tmp)
		}

		colnames(values.tmp) <- values
		if(i==1){
			df <- values.tmp
		}else{

			df <- rbind(df,values.tmp)
		}

	}
	#}
	return(df)
}


data <- read.table(args[6],header=F,stringsAsFactors=F)
gene_symbol <- args[7]
# information that we want to extract
values <- c("AC","AN","AF","AC_FIN","AC_NFE","AN_FIN","AN_NFE","Annotation")
info <- extract.infor(data[,8],values,gene_symbol)
result <- data.frame(data[,-8],info)
colnames(result) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","Allele.Count","Allele.Number","Allele.Frequency","AC_FIN","AC_NFE","AN_FIN","AN_NFE","Annotation")

# output the result
write.table(result,file=args[8],col.names=T,row.names=F,quote=F,sep="\t")

