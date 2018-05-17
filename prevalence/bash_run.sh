#!/usr/bin/sh

# annotate the variants
Rscript annotate_pathogenecity.R [INPUT AF_MUTATION_FILE] [OUTPUT ANNOTATED_AF_MUTATION_FILE]

# calculate the prevalence
Rscript bayesian_estimation.R [ANNOTATED_AF_MUTATION_FILE] [PARAMETER_FILE] [POPULATION] [CONFIDENCE] [OUTPUT PREVELANCE_RESULT_FILE]
