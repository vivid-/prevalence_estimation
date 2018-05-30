genes <- c("FKRP","TCAP","ANO5","SGCA","CAPN3","SGCB","SGCD","SGCG","DYSF")
sink("run_annot_allMis.sh")
for(g in genes){
	test <- paste0("Rscript annotate_pathogenecity_allMissense.R ../data/",g,"_mutation_AF.txt ../data/",g,"_mutation_AF_patho_allMissense.txt")
	cat(test,"\n")
}
