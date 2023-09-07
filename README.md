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

## Un-weighted PRSsum construction

To create an unweighted PRSsum, please follow the instructions in our
repository: <https://github.com/nkurniansyah/Hypertension_PRS>. The
details on how to create the selected PRS can be found in the paper. We
offer an alternative approach for simplification by providing a function
to combine variant weights into a single overall variant weight. This
approach have slightly different when we compute PRSsum like in we
explain in the paper. Below are the instructions:

1.  Clone
    [PRSsum\_Simple](https://github.com/nkurniansyah/PRSsum_Simple "PRSsum_Simple")
2.  in terminal : git clone
    <https://github.com/nkurniansyah/PRSsum_Simple>
3.  then you need to source the function in ../PRSsum\_Simple/Code/\*

<!-- -->

    source("../PRSsum_Simple/Code/create_PRSsum.R")
    source("../PRSsum_Simple/Code/match_allele.R")

### Step 1: Prepare the input files.

To prevent errors, please ensure that you follow the instructions
provided in PRSsum\_Simple. We do not include rsIDs in the variant
weight files since some of them lack rsIDs.

    library(data.table)
    library(dplyr)
    library(tidyverse)

    studies<-c("FINNGEN","Jun","Bellenguez","Kunkle","Kunkle_AFR")
    study_list<-list()
    for(study in studies){
      variant_weight_file<- paste0("Variant_weight/",study,".csv")
      
      variant_weight<-fread(variant_weight_file, data.table = F)
      colnames(variant_weight)<-c("Chromosome","Position", "Allele1","Allele2", "BETA")

      study_list[[study]]<- variant_weight
    }

    #example of variants weight
    head(study_list[["Jun"]])

    ##   Chromosome  Position Allele1 Allele2    BETA
    ## 1          1 207577223       T       C  0.1707
    ## 2          2 127090354       A       G  0.1384
    ## 3          2 127115459       T       C -0.1198
    ## 4          2 127132356       A       G  0.1359
    ## 5          2 127135234       T       C  0.1563
    ## 6          8  27362470       T       C  0.1011

Next, we will align all the alleles with the reference SNP. In this
project, we have utilized Bellenguez variant weights as references. To
expedite computations for this project, we have pre-processed the
reference SNP data, including only the SNPs used to create the PRS.

    Bellenguez_ref_snp<- fread("Variant_weight/Bellenguez_ref_snp.csv", data.table = F)
    ad_variants_weight_clean<- match_allele(refrence_snp=Bellenguez_ref_snp, 
                                         list_variants_weight=study_list,
                                         match_by_position=TRUE)

    ## variant weight will match using chromosome and position

    ## 23 alleles from FINNGEN are being flipped to align them with the reference SNP

    ## variant weight will match using chromosome and position

    ## 1 alleles from Jun are being flipped to align them with the reference SNP

    ## variant weight will match using chromosome and position

    ##  All the allele from Bellenguez are match with the reference SNP

    ## variant weight will match using chromosome and position

    ## 8 alleles from Kunkle are being flipped to align them with the reference SNP

    ## variant weight will match using chromosome and position

    ##  All the allele from Kunkle_AFR are match with the reference SNP

    head(ad_variants_weight_clean)

    ##   Chromosome  Position Allele1 Allele2            rsID   FINNGEN    Jun
    ## 1          2 127121925       A       G 2:127121925:A:G  0.120350     NA
    ## 2          2 127135234       T       C 2:127135234:T:C  0.126170 0.1563
    ## 3         11  60177107       T       C 11:60177107:T:C -0.112623     NA
    ## 4         19  43481304       A       G 19:43481304:A:G -0.134282     NA
    ## 5         19  43620216       T       G 19:43620216:T:G  0.230152     NA
    ## 6         19  43638619       A       G 19:43638619:A:G  0.127756     NA
    ##   Bellenguez Kunkle Kunkle_AFR
    ## 1     0.0797 0.0837         NA
    ## 2     0.1686 0.1693         NA
    ## 3         NA     NA         NA
    ## 4         NA     NA         NA
    ## 5         NA     NA         NA
    ## 6         NA     NA         NA

### Step 2: Create summation variant weight

we use the TOPMed mean and SD values of each variants weight. Using the
same scaling throughout guarantees that effect size estimates are
similarly interpreted across all datasets and individuals who use this
PRS.

    PRSsum_Scaling<- fread("Variant_weight/PRSsum_Scaling.csv")
    PRSsum_Scaling

    ##         Study         Mean          SD N_variants
    ## 1: Bellenguez  0.000415514 0.000337105       1937
    ## 2:    FINNGEN  0.002860917 0.007873218         81
    ## 3:        Jun  0.014803928 0.007789044         85
    ## 4:     Kunkle -0.005668489 0.000503101      12002
    ## 5: Kunkle_AFR -0.037393299 0.006307699        157

After prepare for all the input, now we can generate create summation
variant weight. follow this command below:

    ad_unweigted_variants_weight<- create_prsum(variant_weights=ad_variants_weight_clean, 
                                           PRSsum_scaling=PRSsum_Scaling, 
                                           weight_file=NULL,
                                           chr_col_name="Chromosome", 
                                           pos_col_name="Position", 
                                           effect_allele_col_name="Allele1",
                                           other_allele_col_name="Allele2",
                                           rsID_col_name=NULL)

    ## Run unweighted PRSsum

    head(ad_unweigted_variants_weight)

    ##   chr_name chr_position effect_allele other_allele effect_weight
    ## 1        2    127121925             A            G    0.16231748
    ## 2        2    127135234             T            C    0.36008115
    ## 3       11     60177107             T            C   -0.08829982
    ## 4       19     43481304             A            G   -0.10528112
    ## 5       19     43620216             T            G    0.18044608
    ## 6       19     43638619             A            G    0.10016454

    write.table(ad_unweigted_variants_weight, file = "Variant_weight/Unweighted_PRSsum.txt", row.names = F, quote = F, sep="\t")

### Step 3: Construct PRS and perform the asscoiation test

In this paper, we utilized PRSice2 to create the PRS and employed
HCHS/SOL as the reference panel (for further details, please refer to
the paper). To construct the PRS and conduct the association test,
please follow the instructions provided in our repository at
<https://github.com/nkurniansyah/Hypertension_PRS>.
