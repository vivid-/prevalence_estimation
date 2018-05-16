# steps for extracting the variant AFs from gnomAD
chr=$1 # no "chr"
stat=$2
end=$3
gene=$4

# 1. download the WGS and WES files from gnomAD
wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz.tbi
wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr${chr}.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr${chr}.vcf.bgz.tbi

# 2. retrieve WGS and WES based on their genomic coordinates
tabix -fh gnomad.exomes.r2.0.2.sites.vcf.bgz ${chr}:${stat}-${end} > gnomad.exomes.r2.0.2.sites_${chr}.${stat}-${end}_${gene}.vcf
tabix -fh gnomad.genomes.r2.0.2.sites.chr${chr}.vcf.bgz  ${chr}:${sta}-${end} > gnomad.genomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}.vcf

# 3. parse mutiple ATL alleles into different lines with just one ALT allele
vcf_parser gnomad.exomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}.vcf --split -o gnomad.exomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted.vcf
vcf_parser gnomad.genomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}.vcf --split -o gnomad.genomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted.vcf

# 4. vcf2table
Rscript vcf2table.R gnomad.exomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted.vcf gnomad.exomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted.table
Rscript vcf2table.R gnomad.genomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted.vcf gnomad.genomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted.table

# 5. mini represent vcf information
python get_minimal_representation.py --INPUT gnomad.exomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted.table --OUTPUT gnomad.exomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted_miniRepresented.table
python get_minimal_representation.py --INPUT gnomad.genomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted.table --OUTPUT gnomad.genomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted_miniRepresented.table

# 6. merge WES and WGS result
Rscript merge_exome_genome.R gnomad.genomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted_miniRepresented.table gnomad.exomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted_miniRepresented.table gnomad.genomes_exomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted_miniRepresented.table



