genes <- c("FKRP","TCAP","ANO5","SGCA","CAPN3","SGCB","SGCD","SGCG","SGCA")
sink("run_annot.sh")
for(g in genes){
	test <- paste0("Rscript annotate_pathogenecity.R /ysm-gpfs/scratch60/wl382/WES/prevalence_estimation/prevalence_estimation/extract_AF/data/",g,"_mutation_AF.txt ./data/",g,"_mutation_AF_patho.txt")
	cat(test,"\n")
}
