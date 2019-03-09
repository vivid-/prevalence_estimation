# Estimating prevalences for monogenic autosomic recessive diseases
Estimating prevalence for limb-girdle muscular dystrophy based on public sequencing databases, see https://www.biorxiv.org/content/early/2018/12/20/502708 for more method details.

To replicate results presented in the manuscript above (which used the clinvar database updated in 20180429), please run the script after installation of all required libraries 
```bash
$ sh run_replicate.sh [chr] [start] [end] [gene_symbol] [confidence level]
# for example
$ sh run_replicate.sh 2 71680852 71913898 DYSF 0.95
```

## Install required libraries
We are using `R version 3.4.1` and `python 2.7.11` here. 
For R, `plyr` library is required, which can be installed via 
``` R
> install.packages("plyr")
```

For python, you can install all required libraries using the command line
``` bash
$ pip install -r requirements.txt
```

Besides, you need to have [tabix](http://wiki.wubrowse.org/How_to_install_tabix) installed beforehand.  


## Disease prevalence estimation
For a autosomal recessive disease whose causal gene is known, we can use the following script to estimate its disease prevalence in different populations. Only European (EUR), Finnish (FIN), Non-Finnish European (NFE) and all population in gnomAD data are supported now.
```bash
$ sh run.sh [chr] [start] [end] [gene_symbol] [confidence level] [ the most updated clinvar version]
# please go to ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/ for the most updated clinvar database version or the archived version in 2018 (ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_*/20XX/)
# for example
$ sh run.sh 2 71680852 71913898 DYSF 0.95 clinvar_20190108
or
$ sh run.sh 2 71680852 71913898 DYSF 0.95 archive_2.0/2018 clinvar_20180429
```
For more detailed computation steps, please go to `details.md` file.

## Utilize self-defined variant categories for estimation
To integrate more information for categorizing variants and get more specific allele frequency priors for each self-defined variant category, we provide a script  `prior_estimation/get_beta_parameters_flexible.R` for prior caculation. The script accept input file formatted as:
```
CHROM   POS     REF     ALT     ID      AN      AN_Adj  AC      AC_Adj  type
1       69428   T       G       rs140739101     99358   80618   2141    1985    exon_variant_20larger
1       69590   T       A       rs141776804     93836   83862   110     103     exon_variant_20larger
```
The detailed calculation is as follows:
```bash
# calculate the distribution priors for each category
$ Rscript prior_estimation/get_beta_parameters_flexible.R [category input file] [prior output file]
# for example
$ Rscript prior_estimation/get_beta_parameters_flexible.R ExAC_cates.txt ExAC_priors.txt

# calculate disease prevalence using the different priors
$ Rscript prevalence/bayesian_estimation_flexible.R [gene_alleleFrequeny_pathogeneicityAnnotation file] [prior file] [population] [category input file] [gene chromosome] [confidence level] [output result file]
$ Rscript prevalence/bayesian_estimation_flexible.R ./data/CAPN3_mutation_AF_patho.txt ExAC_priors.txt All ExAC_cates.txt 7 0.95 ./result/CAPN3_All_0.95.txt
```

## Annotate pathogenicity considering unseen variants in gnomAD
To consider variants annotated as pathogenic ones in mutation databases like ClinVar but not found in gnomAD database, we provide a script `prevalence/annotate_pathogenecity_includeUnseenVariant.R`. After re-annotating the variants, we can use the script `prevalence/bayesian_estimation_unseen.R` to recalculate the prevalence estimators。 

The detailed process is as follows:
```bash
# annotate the variants
$ Rscript prevalence/annotate_pathogenecity_includeUnseenVariant.R [mutation AF input file] [annotated output file]
# for example
$ Rscript prevalence/annotate_pathogenecity_includeUnseenVariant.R CAPN3_mutation_AF.txt CAPN3_mutation_AF_unseen_patho.txt
# re-calculate disease prevalence 
$ Rscript prevalence/bayesian_estimation_unseen.R [annotated file] [prior file] [population] [confidence score] [output file]
# for example
$ Rscript prevalence/bayesian_estimation_unseen.R CAPN3_mutation_AF_unseen_patho.txt ./data/beta_parameter_prior_ExAC.txt AFR 0.95 CAPN3_AFR_0.95_unseen.txt 
```



## Pre-computed results
For convinience, you can find the pre-computed files under the `data` directory and directly run the step `Estimate prevalences for diseases of interests` in the `details.md` file.


