#!/usr/bin/R

genes <- c("FKRP","TCAP","ANO5","SGCA","CAPN3","SGCB","SGCD","SGCG","SGCA")
populations <- c("All","EUR","NFE","FIN")

sink("run_test.sh")
for(g in genes){
	for (p in populations){
		test <- paste0("Rscript bayesian_estimation.R ./data/",g,"_mutation_AF_patho.txt ../beta_parameter_prior_ExAC.txt ",p," 0.95 ./data/",g,"_",p,"_0.95.txt")
		cat(test,"\n")
	}
}
