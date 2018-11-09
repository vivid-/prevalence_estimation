####################################
# this script is used to calculate #
# disease prevalence in the direct #
# way ##############################
####################################

library(plyr)

args <- commandArgs()
input.file <- args[6] # input file for AF and also patheogenic annotation
population <- args[7] # specify the interested population, NFE (Non_Finish_European),FIN (Finish),European,All
output_file <- args[8] # specify the output path

# read the input file
df <- read.table(input.file,header=T,sep="\t",quote="")
# only focused on the variants annotated as pathogenic
all.df <- df[which(df$patheo.info==1),]


# deal with NAs for AC or AN columns
all.df$AC_FIN.genome[which(is.na(all.df$AC_FIN.genome))] = 0
all.df$AN_FIN.genome[which(is.na(all.df$AN_FIN.genome))] = 0
all.df$AC_FIN.exome[which(is.na(all.df$AC_FIN.exome))] = 0
all.df$AN_FIN.exome[which(is.na(all.df$AN_FIN.exome))] = 0

all.df$AC_NFE.genome[which(is.na(all.df$AC_NFE.genome))] = 0
all.df$AN_NFE.genome[which(is.na(all.df$AN_NFE.genome))] = 0
all.df$AC_NFE.exome[which(is.na(all.df$AC_NFE.exome))] = 0
all.df$AN_NFE.exome[which(is.na(all.df$AN_NFE.exome))] = 0

##########

##########
# for different population, find the corresponding AF and AN column
if(population == "FIN"){
        AC_interested <- all.df$AC_FIN.exome + all.df$AC_FIN.genome
        AN_interested <- all.df$AN_FIN.exome + all.df$AN_FIN.genome
}
if(population == "NFE"){
        AC_interested <- all.df$AC_NFE.exome + all.df$AC_NFE.genome
        AN_interested <- all.df$AN_NFE.exome + all.df$AN_NFE.genome
}
if(population == "EUR"){
        AC_interested <- all.df$AC_FIN.exome + all.df$AC_FIN.genome + all.df$AC_NFE.exome + all.df$AC_NFE.genome
        AN_interested <- all.df$AN_FIN.exome + all.df$AN_FIN.genome + all.df$AN_NFE.exome + all.df$AN_NFE.genome
}
if(population == "All"){
        AC_interested <- all.df$Allele.Count
        AN_interested <- all.df$Allele.Number
}
if(population == "ASJ"){
        AC_interested <- all.df$AC_ASJ.exome + all.df$AC_ASJ.genome
        AN_interested <- all.df$AN_ASJ.exome + all.df$AN_ASJ.genome
}
if(population == "AFR"){
        AC_interested <- all.df$AC_AFR.exome + all.df$AC_AFR.genome
        AN_interested <- all.df$AN_AFR.exome + all.df$AN_AFR.genome
}
if(population == "EAS"){
        AC_interested <- all.df$AC_EAS.exome + all.df$AC_EAS.genome
        AN_interested <- all.df$AN_EAS.exome + all.df$AN_EAS.genome
}

afs_original <- AC_interested/AN_interested
# remove NAs
afs_original  <- na.omit(afs_original)
direct_calculation <- (1-prod(1-afs_original))^2

sink(output_file)
cat("estimated prevalence using direct calculation: ", direct_calculation,"\n")
sink()



