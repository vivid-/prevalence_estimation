# Estimating prevalences for monogenic autosomic recessive diseases
Complete version of prevalence estimation

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
$ sh run.sh [chr] [start] [end] [gene_symbol] [confidence level]
# for example
$ sh run.sh 2 71680852 71913898 DYSF 0.95
```
For more detailed computation steps, please go to `details.md` file.
