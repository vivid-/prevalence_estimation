#!/usr/bin/R

####################################
# this script is used to merge WGS #
# and WES data tables ##############
####################################

args <- commandArgs()

options(stringsAsFactors=F)
merge.alleleFreq <- function(dat1,dat2){
	# outer join two dataframes
        all.join <- merge(x = dat1, y = dat2,
                        by = c("CHROM","POS","REF","ALT"),
			all = TRUE)
        
	afs <- rep(0.0,length=dim(all.join)[1])
	# for variants existed in one database
	afs[which(!is.na(all.join$Allele.Number.exome))] = all.join$Allele.Frequency.exome[which(!is.na(all.join$Allele.Number.exome))]
	afs[which(!is.na(all.join$Allele.Number.genome))] = all.join$Allele.Frequency.genome[which(!is.na(all.join$Allele.Number.genome))]
	# for intersected variants
	ind.intersect = intersect(which(!is.na(all.join$Allele.Number.exome)),which(!is.na(all.join$Allele.Number.genome)))
	afs[ind.intersect] = (all.join$Allele.Count.genome[ind.intersect] + all.join$Allele.Count.exome[ind.intersect])/(all.join$Allele.Number.genome[ind.intersect] + all.join$Allele.Number.exome[ind.intersect])
	tmp <- colnames(all.join)
	
	
	result <- cbind(all.join,afs)
	colnames(result) <- c(tmp,"Allele Frequency integrated")
	
	acs <- rep(0.0,length=dim(all.join)[1])
	ans <- rep(0.0,length=dim(all.join)[1])

	# afs[which(!is.na(all.join$Allele.Number.y))] = all.join$Allele.Count.y[which(!is.na(all.join$Allele.Number.y))] + 
	all.join[which(is.na(all.join$Allele.Number.exome)),]$Allele.Number.exome = 0
	all.join[which(is.na(all.join$Allele.Number.genome)),]$Allele.Number.genome = 0	
	all.join[which(is.na(all.join$Allele.Count.genome)),]$Allele.Count.genome = 0
	all.join[which(is.na(all.join$Allele.Count.exome)),]$Allele.Count.exome = 0
	acs <- all.join$Allele.Count.genome + all.join$Allele.Count.exome
	ans <- all.join$Allele.Number.genome + all.join$Allele.Number.exome
	tmp <- colnames(all.join)
	result <- cbind(all.join,acs,ans)
	
	colnames(result) <- c(tmp,"Allele.Count","Allele.Number")
	return(result)
}


genome_file <- args[6]
exome_file <- args[7]
output_file <- args[8]
genome <- read.table(genome_file,header=T,sep="\t")
exome <- read.table(exome_file,header=T,sep="\t")
colnames(genome) <- c("CHROM","POS","ID.genome","REF","ALT","QUAL.genome","FILTER.genome","Allele.Count.genome","Allele.Number.genome","Allele.Frequency.genome","AC_FIN.genome","AC_NFE.genome","AN_FIN.genome","AN_NFE.genome","Annotation.genome")
colnames(exome) <- c("CHROM","POS","ID.exome","REF","ALT","QUAL.exome","FILTER.exome","Allele.Count.exome","Allele.Number.exome","Allele.Frequency.exome","AC_FIN.exome","AC_NFE.exome","AN_FIN.exome","AN_NFE.exome","Annotation.exome")
result <- merge.alleleFreq(genome,exome)

write.table(result,file=output_file,col.names=T,row.names=F,sep="\t",quote=F)
