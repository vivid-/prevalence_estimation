# Estimating prevalences for monogenic autosomic recessive diseases
Complete version of prevalence estimation

## Install required libraries
We are using `R version 3.4.1` and `python 2.7.11` here. 
For R, `plyr` library is required, which can be installed via 
``` R
install.packages("plyr")
```

For python, you can install all required libraries using the command line
``` bash
pip install -r requirements.txt
```

Besides, you need to have [tabix](http://wiki.wubrowse.org/How_to_install_tabix) installed beforehand.  


## Downlad pathogenicity annotation from mutation database
For EGL database, the python script `extract_EGL.py` would be used to fetch information from the website.
``` bash
python ./fetch_database/extract_EGL.sh --gene [GENE] --out [OUT]
## for example: 
python ./fetch_database/extract_EGL.sh --gene DYSF --out ../DYSF_EGL.txt
```
Another python script `convert_HGVS.py` would be furtherly used to extract genomic coordinates for annotated variants in the EGL database.
```bash
python ./fetch_database/convert_HGVS.py --inpu [INPU] --out [OUT]
## for example:
python ./fetch_database/convert_HGVS.py --inpu ../DYSF_EGL.txt --out DYSF_EGL_loc.txt
```

For ClinVar database, we would extract the annotation information from the ClinVar database
```bash
# download the databse
sh ./fetch_database/extract_clinvar.sh [CHR] [START] [STOP] [GENE]
## for example:
sh ./fetch_database/extract_clinvar.sh 2 71680852 71913898 DYSF
```

Merge two database files
```bash
Rscript ./fetch_database/merge_database.R [clinvar_file] [EGL_file] [database_output]
## for example:
Rscript ./fetch_database/merge_database.R clinvar2_71680852-71913898_DYSF_info.table DYSF_EGL_loc.txt DYSF_mutation.txt

# getting the minirepresentational reference alleles and alternative alleles
python ./extract_AF/get_minimal_representation.py --INPUT [database_output] --OUTPUT [database_miniRepresented_output]
## for example:
python ./extract_AF/get_minimal_representation.py --INPUT DYSF_mutation.txt DYSF_mutation_miniR.txt
```

## Download allele frequency files for genes of interests
```bash
# download AF files from gnomAD
sh ./extract_AF/bash_pre.sh [CHR] [START] [STOP] [GENE]
## for example:
sh ./extract_AF/bash_pre.sh 2 71680852 71913898 DYSF

# merge AF files with mutation files
Rscript ./extract_AF/merge_mutation_AF.R [database_output] [AF_output] [merged_database_AF]
## for example:
Rscript ./extract_AF/merge_mutation_AF.R DYSF_mutation_miniR.txt gnomad.genomes_exomes.r2.0.2.sites_2.71680852-71913898_DYSF_splitted_miniRepresented.table DYSF_mutation_AF.txt
```

## Estimate beta paramters for AF priors

Here, we are using whole exome sequencing data from [ExAC](ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/manuscript_data/) to estimate beta distribution parameters for categorized variants. We are only considering frameshift, splice acceptor/donor, stop gained and missense variants here.

```bash
# categorize variants
Rscript ./prior_estimation/empirical_beta.R [WES_INPUT] [OUTPUT_prefix]
## for example:
Rscript ./prior_estimation/empirical_beta.R ExAC.r0.3.1.sites.vep.canonical.table.gz xAC.r0.3.1.sites.vep.canonical

# get the estimated beta priors for each category
Rscript ./prior_estimation/get_beta_parameters.R [INPUT_prefix] [prior_output]
## for example:
Rscript ./prior_estimation/get_beta_parameters.R xAC.r0.3.1.sites.vep.canonical ExAC_beta_priors.txt 
``` 

## Estimate prevalences for diseases of interests

In this part, result files from previous steps would be used to estimate the disease prevalence. Besides, we can estimate population-specific disease prevalences seperately, like European (EUR), Finnish (FIN), non-Finnish European (NFE) and all populations (All).

```bash
# annotate the variants
Rscript ./prevalence/annotate_pathogenecity.R [INPUT AF_MUTATION_FILE] [OUTPUT ANNOTATED_AF_MUTATION_FILE]
## for example:
Rscript ./prevalence/annotate_pathogenecity.R DYSF_mutation_AF.txt DYSF_mutation_AF_patho.txt

# calculate the prevalence
Rscript ./prevalence/bayesian_estimation.R [ANNOTATED_AF_MUTATION_FILE] [PARAMETER_FILE] [POPULATION] [CONFIDENCE] [OUTPUT PREVELANCE_RESULT_FILE]
## for example:
Rscript ./prevalence/bayesian_estimation.R DYSF_mutation_AF_patho.txt ExAC_beta_priors.txt All 0.95 DYSF_All_0.95.txt
```

