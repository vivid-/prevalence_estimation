# extract variants based on their locations
chr="$1"
sta="$2"
end="$3"
gene="$4"

tabix /ysm-gpfs/scratch60/wl382/WES/Data/CADD/ExAC_r0.3_inclAnno.tsv.gz ${chr}:${sta}-${end} > /ysm-gpfs/scratch60/wl382/WES/Data/CADD/CADD_${chr}_${sta}-${end}_${gene}.tsv


Rscript caddTSV2table.R /ysm-gpfs/scratch60/wl382/WES/Data/CADD/CADD_${chr}_${sta}-${end}_${gene}.tsv /ysm-gpfs/scratch60/wl382/WES/Data/CADD/CADD_${chr}_${sta}-${end}_${gene}.txt

python ../extract_AF/get_minimal_representation.py --INPUT /ysm-gpfs/scratch60/wl382/WES/Data/CADD/CADD_${chr}_${sta}-${end}_${gene}.txt --OUTPUT /ysm-gpfs/scratch60/wl382/WES/Data/CADD/CADD_${chr}_${sta}-${end}_${gene}_miniRepresented.txt

Rscript merge_CADD_mutation_AF.R /ysm-gpfs/scratch60/wl382/WES/prevalence_estimation/prevalence_estimation/data/${gene}_mutation_AF.txt /ysm-gpfs/scratch60/wl382/WES/Data/CADD/CADD_${chr}_${sta}-${end}_${gene}_miniRepresented.txt /ysm-gpfs/scratch60/wl382/WES/Data/CADD/${gene}_mutation_AF_CADD.txt
