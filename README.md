## Introduction

This repository provides information regarding the construction of a
polygenic risk score (PRS) for Alzhaimer (AD) that we developed (AD-PRS)
in manuscript [A polygenic risk score for Alzheimer’s disease
constructed using APOE-region variants has stronger association than
APOE alleles with mild cognitive impairment in Hispanic/Latino adults in
the
U.S.](https://alzres.biomedcentral.com/articles/10.1186/s13195-023-01298-3 "A polygenic risk score for Alzheimer’s disease constructed using APOE-region variants has stronger association than APOE alleles with mild cognitive impairment in Hispanic/Latino adults in the")

First, it provides instructions for constructing the AD-PRS based on
summary statistics from GWAS. We provide the relevant summary statistics
(see folder “Variant\_weight”). Additionally, we provide instructions on
generating unweighted variant weights. This step simplifies the process
of generating the PRS by combining multiple variant weights into a
single overall variant weight.

## Required packages

We used [PRSice
2.3.1.e](https://choishingwan.github.io/PRSice/ "PRSice 2.3.1.e") to
generate PRS. We provide example code that also uses PRSice to construct
PRS based on the provided summary statistics in folder
“Variant\_weight”.

Other software and packages that we used, but may not be necessary for
others to construct the PRS, are as follows:

We performed the analysis using R version 4.0.2.

We used the following packages from CRAN: dplyr, tidyverse, data.table,
purrr, pROC.

install.packages(“dplyr”)
