# steps for extracting the variant AFs from gnomAD
chr=$1 # no "chr"
sta=$2
end=$3
gene=$4

Rscript /ysm-gpfs/scratch60/wl382/WES/prevalence_estimation/prevalence_estimation/extract_AF/vcf2table.R gnomad.exomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted.vcf ${gene} gnomad.exomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted.table

Rscript /ysm-gpfs/scratch60/wl382/WES/prevalence_estimation/prevalence_estimation/extract_AF/vcf2table.R gnomad.genomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted.vcf ${gene} gnomad.genomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted.table

python /ysm-gpfs/scratch60/wl382/WES/prevalence_estimation/prevalence_estimation/extract_AF/get_minimal_representation.py --INPUT gnomad.exomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted.table --OUTPUT gnomad.exomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted_miniRepresented.table

python /ysm-gpfs/scratch60/wl382/WES/prevalence_estimation/prevalence_estimation/extract_AF/get_minimal_representation.py --INPUT gnomad.genomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted.table --OUTPUT gnomad.genomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted_miniRepresented.table

Rscript /ysm-gpfs/scratch60/wl382/WES/prevalence_estimation/prevalence_estimation/extract_AF/merge_genome_exome.R gnomad.genomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted_miniRepresented.table gnomad.exomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted_miniRepresented.table gnomad.genomes_exomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted_miniRepresented.table
