genes <- c("FKRP","TCAP","ANO5","SGCA","CAPN3","SGCB","SGCD","SGCG","SGCA")
sink("run_annot_allMis.sh")
for(g in genes){
	test <- paste0("Rscript annotate_pathogenecity_allMissense.R /ysm-gpfs/scratch60/wl382/WES/data/data/",g,"_mutation_AF.txt /ysm-gpfs/scratch60/wl382/WES/data/data/",g,"_mutation_AF_patho_allMissense.txt")
	cat(test,"\n")
}
