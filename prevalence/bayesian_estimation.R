#!/usr/bin/R

###################################
# this script is used to estimate #
# disease prevalence ##############
###################################

library(plyr)
library("openxlsx")

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

# get the value of Marcum Q-function for CDF
get_marcum_q_func <- function(M,a,b,thresh){
  tmp.1 = exp(-0.5*a^2)
  max_num = 1000
  tmp.sum = 0
  for(k in 0:max_num){
    tmp.2 = (0.5*a^2)^k
    tmp.3 = gamma(M+k) -gamma(0.5*b^2)
    tmp.4 = gamma(M+k)
    tmp.5 = factorial(k)
    
    tmp.now = tmp.1*tmp.2/tmp.4 *(tmp.3/tmp.5)
    cat(tmp.now,"\n")
    if(abs(tmp.now)>thresh){
      tmp.sum = tmp.sum + tmp.now
    }else{
      return(1-tmp.sum)
    }
    
  }
  
  return(1-tmp.sum)
}

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
input.file <- args[6]
param.file <- args[7]
params <- read.table(param.file,header=T)
af_file <- args[8] # input file for the allele frequency
patho_file <- args[9] # input file for patheogenecity annotation
#library("openxlsx")
af.table <- read.table(af_file,header = T,stringsAsFactors=F,fill = T,sep="\t")
types = c("frameshift_variant","splice_acceptor_variant","splice_donor_variant","missense_variant","stop_gained")
type_names = c("frameshift","splice acceptor","splice donor","missense","stop gained")

dat.1 <- read.table(patho_file,header = T,stringsAsFactors=F,fill = T,sep="\t")
dat.patheo <- dat.1[which(patheo.info==1),]
all.df <- merge(x=dat.patheo,y=af.table,
                by.x = c("Chrom","Position","Reference","Alternate"),
                by.y = c("CHROM","POS","REF","ALT"),all.x = T)

# deal with NAs for AC or AN columns
all.df$AC_FIN.x[which(is.na(all.df$AC_FIN.x))] = 0
all.df$AN_FIN.x[which(is.na(all.df$AN_FIN.x))] = 0
all.df$AC_FIN.y[which(is.na(all.df$AC_FIN.y))] = 0
all.df$AN_FIN.y[which(is.na(all.df$AN_FIN.y))] = 0

all.df$AC_NFE.x[which(is.na(all.df$AC_NFE.x))] = 0
all.df$AN_NFE.x[which(is.na(all.df$AN_NFE.x))] = 0
all.df$AC_NFE.y[which(is.na(all.df$AC_NFE.y))] = 0
all.df$AN_NFE.y[which(is.na(all.df$AN_NFE.y))] = 0

##########
# remove the NAs in the matched af.table
all.df <- all.df[-which(is.na(all.df[,47])),]

##########
# ! should add an input argument here to see which column/population they want to estimate

af_changed <- rep(0,dim(all.df)[1])
posterior_param <- data.frame(v=rep(0,dim(all.df)[1]),w=rep(0,dim(all.df)[1]))
for(i in 1:length(types)){
  ind.tmp <- which(all.df$Annotation==type_names[i])
  v = params$v[which(params$type==types[i])]
  w = params$w[which(params$type==types[i])]
 # ac= all.df[ind.tmp,45]
#  an = all.df[ind.tmp,46]
  ac = all.df$AC_NFE.x[ind.tmp] + all.df$AC_NFE.y[ind.tmp]+all.df$AC_FIN.x[ind.tmp] + all.df$AC_FIN.y[ind.tmp]
  an = all.df$AN_NFE.x[ind.tmp] + all.df$AN_NFE.y[ind.tmp]+ all.df$AN_FIN.x[ind.tmp] + all.df$AN_FIN.y[ind.tmp]
#  ac = all.df$AC_FIN.x[ind.tmp] + all.df$AC_FIN.y[ind.tmp]
#  an = all.df$AN_FIN.x[ind.tmp] + all.df$AN_FIN.y[ind.tmp]
#  ac = all.df$Allele.Count.European..Non.Finnish.[ind.tmp]#+all.df$Allele.Count.European..Finnish.[ind.tmp]
#  an = all.df$Allele.Number.European..Non.Finnish.[ind.tmp]#+all.df$Allele.Number.European..Finnish.[ind.tmp]
  post_params <- get_posterior_w_v(ac,an/2,w,v)
  posterior_param[ind.tmp,1]<- post_params[[1]]
  posterior_param[ind.tmp,2]<- post_params[[2]]
  af_changed_tmp <- pre_w_v_plugged(v,w,an/2,ac)
  af_changed[ind.tmp] <- af_changed_tmp
}




