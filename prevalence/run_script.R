#!/usr/bin/R
#genes <- c("DYSF")
genes <- c("FKRP","TCAP","ANO5","SGCA","CAPN3","SGCB","SGCD","SGCG","SGCA","DYSF")
populations <- c("All","EUR","NFE","FIN")

sink("run_test_allMissense.sh")
for(g in genes){
	for (p in populations){
		test <- paste0("Rscript bayesian_estimation.R /ysm-gpfs/scratch60/wl382/WES/prevalence_estimation/prevalence_estimation/data/",g,"_mutation_AF_patho_allMissense.txt /ysm-gpfs/scratch60/wl382/WES/data/beta_parameter_prior_ExAC.txt ",p," 0.95 /ysm-gpfs/scratch60/wl382/WES/prevalence_estimation/prevalence_estimation/result/",g,"_",p,"_0.95_allMissense.txt")
		cat(test,"\n")
	}
}
