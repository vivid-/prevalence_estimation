#!/usr/bin/R

###################################
# this script is used to estimate #
# disease prevalence ##############
###################################

library(plyr)
#library("openxlsx")

# estimate the mean of AF posterior distribution
pre_w_v_plugged <- function(v,w,n,x){
 # tmp <- ((x+v-1)^2)/(2*n+v+w-1)^2
  tmp.1 <- (x+v)/(2*n+v+w)
  
  return(tmp.1)
}


### get the number of variables belonging to each category ###
get_mj <- function(vs,ws){
  tmp.df  <- data.frame(v=vs,w=ws)
  
  tmp=aggregate(list(numdup=rep(1,nrow(tmp.df))), tmp.df, length)
  return(tmp)
}


#######################################################################################
############# use normal distribution to approximate the sum of beta r.v.s ############
#######################################################################################

#### get the mean of the approximated normal distribution ####
get_mean_normal <- function(vs,ws){
  tmp <- vs/(vs+ws)
  return(sum(tmp))
}

#### get the variance of the approximated normal distribution ####
get_var_normal <- function(vs,ws){
  vars <- vs*ws/((vs+ws+1)*(vs+ws)^2)
  return(sum(vars))
}

###################################################################################################
############# use chi square distribution to approximate the squared sum of beta r.v.s ############
###################################################################################################


# get the noncentral parameter for- the approximated chi-square distribution
get_noncentral_param <- function(mu,sigma){
  return((mu/sigma)^2)
}

# get the posterior params 
get_posterior_w_v <- function(x,n,w,v){
  a <- x+v
  b <- 2*n-x+w
  return(list(a,b))
}


options(stringsAsFactors = F)
args <- commandArgs()
input.file <- args[6] # input file for AF and also patheogenic annotation
param.file <- args[7]
population <- args[8] # specify the interested population, NFE (Non_Finish_European),FIN (Finish),European,All
params <- read.table(param.file,header=T)

# read the input file
df <- read.table(input.file,header=T,sep="\t")


types = c("frameshift_variant","splice_acceptor_variant","splice_donor_variant","missense_variant","stop_gained","exon_variant","UTR_variant","other_variant")
#type_names = c("frameshift","splice acceptor","splice donor","missense","stop gained")

# only focused on the variants annotated as pathogenic
all.df <- df[which(df$patheo.info==1),]

# deal with NAs for AC or AN columns
na2num <- function(dat,colnm){
    ind <- which(colnames(dat)==colnm)
    dat[which(is.na(dat[,ind])),ind] = 0
    return(dat)
}


# deal with NAs for AC or AN columns
all.df$AC_FIN.genome[which(is.na(all.df$AC_FIN.genome))] = 0
all.df$AN_FIN.genome[which(is.na(all.df$AN_FIN.genome))] = 0
all.df$AC_FIN.exome[which(is.na(all.df$AC_FIN.exome))] = 0
all.df$AN_FIN.exome[which(is.na(all.df$AN_FIN.exome))] = 0

all.df$AC_NFE.genome[which(is.na(all.df$AC_NFE.genome))] = 0
all.df$AN_NFE.genome[which(is.na(all.df$AN_NFE.genome))] = 0
all.df$AC_NFE.exome[which(is.na(all.df$AC_NFE.exome))] = 0
all.df$AN_NFE.exome[which(is.na(all.df$AN_NFE.exome))] = 0
all.df <- na2num(all.df,"AC_ASJ.genome")
all.df <- na2num(all.df,"AC_AFR.genome")
all.df <- na2num(all.df,"AC_EAS.genome")
all.df <- na2num(all.df,"AN_ASJ.genome")
all.df <- na2num(all.df,"AN_AFR.genome")
all.df <- na2num(all.df,"AN_EAS.genome")
all.df <- na2num(all.df,"AC_ASJ.exome")
all.df <- na2num(all.df,"AC_AFR.exome")
all.df <- na2num(all.df,"AC_EAS.exome")
all.df <- na2num(all.df,"AN_ASJ.exome")
all.df <- na2num(all.df,"AN_AFR.exome")
all.df <- na2num(all.df,"AN_EAS.exome")


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
	all.df <- all.df[which(AC_interested!=0),]
        AC_interested <- AC_interested[which(AC_interested!=0)]
        AN_interested <- all.df$AN_ASJ.exome + all.df$AN_ASJ.genome
}
if(population == "AFR"){
        AC_interested <- all.df$AC_AFR.exome + all.df$AC_AFR.genome
	all.df <- all.df[which(AC_interested!=0),]
        AC_interested <- AC_interested[which(AC_interested!=0)]
        AN_interested <- all.df$AN_AFR.exome + all.df$AN_AFR.genome
}
if(population == "EAS"){
        AC_interested <- all.df$AC_EAS.exome + all.df$AC_EAS.genome
	all.df <- all.df[which(AC_interested!=0),]
        AC_interested <- AC_interested[which(AC_interested!=0)]
        AN_interested <- all.df$AN_EAS.exome + all.df$AN_EAS.genome
}




af_changed <- rep(0,dim(all.df)[1])
posterior_param <- data.frame(v=rep(0,dim(all.df)[1]),w=rep(0,dim(all.df)[1]))
inds_updated <- c()
for(i in 1:length(types)){
  ind.tmp <- grep(types[i],all.df$Annotation)
  inds_updated <- c(inds_updated,ind.tmp)
  v = params$v[which(params$type==types[i])]
  w = params$w[which(params$type==types[i])]
  ac = AC_interested[ind.tmp]
  an = AN_interested[ind.tmp]

  post_params <- get_posterior_w_v(ac,an/2,w,v)
  posterior_param[ind.tmp,1]<- post_params[[1]]
  posterior_param[ind.tmp,2]<- post_params[[2]]
  af_changed_tmp <- pre_w_v_plugged(v,w,an/2,ac)
  af_changed[ind.tmp] <- af_changed_tmp
}

# consider the other variant for varint type not included
ind.tmp = setdiff(1:dim(all.df)[1],inds_updated)
v = params$v[which(params$type==types[length(types)])]
w = params$w[which(params$type==types[length(types)])]
ac = AC_interested[ind.tmp]
an = AN_interested[ind.tmp]
post_params <- get_posterior_w_v(ac,an/2,w,v)
posterior_param[ind.tmp,1]<- post_params[[1]]
posterior_param[ind.tmp,2]<- post_params[[2]]
af_changed_tmp <- pre_w_v_plugged(v,w,an/2,ac)
af_changed[ind.tmp] <- af_changed_tmp



mu <- get_mean_normal(posterior_param$v[which(posterior_param$w!=0)],posterior_param$w[which(posterior_param$w!=0)])
sigma.2 <- get_var_normal(posterior_param$v[which(posterior_param$w!=0)],posterior_param$w[which(posterior_param$w!=0)])
# get the estimates
E.s.2 = mu^2 + sigma.2
prev.estimated =  E.s.2

# get the confidence interval
confidence <- args[9]
alpha = 1 - as.numeric(confidence)
# lower bound
lb <- qchisq(alpha/2,1,ncp = mu^2/sigma.2)*sigma.2
# upper bound
up <- qchisq(1-alpha/2,1,ncp = mu^2/sigma.2)*sigma.2

# output the result
output_file <- args[10] # destination for output
sink(output_file)
cat("estimated prevalence: ",prev.estimated,"\n")
cat("confidence interval with ",as.numeric(confidence)*100,"% cofidence: ",lb,"-",up,"\n")


