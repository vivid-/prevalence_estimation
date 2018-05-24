# fecth clinvar database file
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20180429.vcf.gz
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20180429.vcf.gz.tbi

# extract variants based on their locations
chr="$1"
sta="$2"
end="$3"
gene="$4"
tabix clinvar_20180429.vcf.gz ${chr}:${sta}-${end} > clinvar_${chr}_${sta}-${end}_${gene}.vcf
Rscript extract_patho_annotation.R clinvar_${chr}_${sta}-${end}_${gene}.vcf clinvar_${chr}_${sta}-${end}_${gene}_info.table

