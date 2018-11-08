chr=$1
sta=$2
end=$3
gene=$4
# confidence score
cfs=$5
clinvar=$6

# Downlad pathogenicity annotation from mutation database
mkdir ./result
mkdir ./data
## For EGL database
### download the database for this gene
python ./fetch_database/extract_EGL.py --gene ${gene} --out ./data/${gene}_EGL.txt
### to extract genomic coordinates for annotated variants in the EGL database.
python ./fetch_database/convert_HGVS.py --inpu ./data/${gene}_EGL.txt --out ./data/${gene}_EGL_loc.txt

## For ClinVar databse
### download the databse from clinvar
cd ./fetch_database/;sh extract_clinvar.sh ${chr} ${sta} ${end} ${gene} ${clinvar}

cd ..
## merge two databases
Rscript ./fetch_database/merge_database.R ./data/clinvar_${chr}_${sta}-${end}_${gene}_info.table ./data/${gene}_EGL_loc.txt ./data/${gene}_mutation.txt
### getting the minirepresentational reference alleles and alternative alleles
python ./extract_AF/get_minimal_representation.py --INPUT ./data/${gene}_mutation.txt --OUTPUT ./data/${gene}_mutation_miniR.txt


##########################################################################################################################################

# Download allele frequency files for genes of interests
## download AF files from gnomAD
cd ./extract_AF ;sh bash_pre.sh ${chr} ${sta} ${end} ${gene}
cd ..
## merge AF files with mutation files
Rscript ./extract_AF/merge_mutation_AF.R ./data/${gene}_mutation_miniR.txt ./data/gnomad.genomes_exomes.r2.0.2.sites_${chr}.${sta}-${end}_${gene}_splitted_miniRepresented.table  ./data/${gene}_mutation_AF.txt


##########################################################################################################################################

# Estimate beta paramters for AF priors
## categorize variants
#wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/manuscript_data/ExAC.r1.sites.vep.canonical.table.gz
#wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/manuscript_data/ExAC.r1.sites.vep.canonical.table.gz.tbi
#Rscript ./prior_estimation/empirical_beta.R ExAC.r0.3.1.sites.vep.canonical.table.gz ExAC.r0.3.1.sites.vep.canonical
## get the estimated beta priors for each category
#Rscript ./prior_estimation/get_beta_parameters.R ExAC.r0.3.1.sites.vep.canonical_ ./data/beta_parameter_prior_ExAC.txt


##########################################################################################################################################

# Estimate prevalences for diseases of interests
## annotate variants
Rscript ./prevalence/annotate_pathogenecity.R ./data/${gene}_mutation_AF.txt ./data/${gene}_mutation_AF_patho.txt

## calculate the prevalence
Rscript ./prevalence/bayesian_estimation.R ./data/${gene}_mutation_AF_patho.txt ./data/beta_parameter_prior_ExAC.txt All ${cfs} ./result/${gene}_All_${cfs}.txt
Rscript ./prevalence/bayesian_estimation.R ./data/${gene}_mutation_AF_patho.txt ./data/beta_parameter_prior_ExAC.txt NFE ${cfs} ./result/${gene}_NFE_${cfs}.txt
Rscript ./prevalence/bayesian_estimation.R ./data/${gene}_mutation_AF_patho.txt ./data/beta_parameter_prior_ExAC.txt FIN ${cfs} ./result/${gene}_FIN_${cfs}.txt
Rscript ./prevalence/bayesian_estimation.R ./data/${gene}_mutation_AF_patho.txt ./data/beta_parameter_prior_ExAC.txt EUR ${cfs} ./result/${gene}_EUR_${cfs}.txt
Rscript ./prevalence/bayesian_estimation.R ./data/${gene}_mutation_AF_patho.txt ./data/beta_parameter_prior_ExAC.txt EAS ${cfs} ./result/${gene}_EAS_${cfs}.txt
Rscript ./prevalence/bayesian_estimation.R ./data/${gene}_mutation_AF_patho.txt ./data/beta_parameter_prior_ExAC.txt ASJ ${cfs} ./result/${gene}_ASJ_${cfs}.txt
Rscript ./prevalence/bayesian_estimation.R ./data/${gene}_mutation_AF_patho.txt ./data/beta_parameter_prior_ExAC.txt AFR ${cfs} ./result/${gene}_AFR_${cfs}.txt


