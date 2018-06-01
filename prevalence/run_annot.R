genes <- c("FKRP","TCAP","ANO5","SGCA","CAPN3","SGCB","SGCD","SGCG","DYSF")
sink("run_annot_allMis.sh")
for(g in genes){
	test <- paste0("Rscript annotate_pathogenecity.R /ysm-gpfs/scratch60/wl382/WES/Data/CADD/",g,"_mutation_AF_CADD.txt /ysm-gpfs/scratch60/wl382/WES/Data/CADD/",g,"_mutation_AF_CADD_patho_20.txt")
	cat(test,"\n")
}
