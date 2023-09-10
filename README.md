## Introduction

This repository provides information regarding the construction of a
polygenic risk score (PRS) for Alzhaimer (AD) that we developed (AD-PRS)
in manuscript [A polygenic risk score for Alzheimer’s disease
constructed using APOE-region variants has stronger association than
APOE alleles with mild cognitive impairment in Hispanic/Latino adults in
the
U.S.](https://alzres.biomedcentral.com/articles/10.1186/s13195-023-01298-3 "A polygenic risk score for Alzheimer’s disease constructed using APOE-region variants has stronger association than APOE alleles with mild cognitive impairment in Hispanic/Latino adults in the")

First, it provides instructions for constructing the AD-PRS from GWAS
summary statistics, which are available in the “Variant\_weight” folder
and can also be downloaded from the PGS catalog (Link to be added).
Additionally, we provide an instructions for creating an unweighted
PRSsum (AD-PRS) and an alternative method that simplifies the process by
combining multiple variant weights into a single overall weight before
constructing the PRS.

In this repository, we also include code and examples on how to perform
association analysis, which can be found in the “./Code/\*” directory.

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

## PRS construction

Our AD-PRS is a sum of multiple AD PRS (FINNGEN, Jun et al., Kunkle et
al 2019.,Bellenguez et al., and kunkle et al. 2021). Summary statistics
to create the each of the AD-PRS are provided here in
subfolder(./Variant\_weight/\*). The variant positions are in Hg38 and
Alelle1 is the effect alelle.

Specific GWAS used:

    ##             GAS name    Race/ethnicity     References
    ## 1            FinnGen  Finnish European              -
    ## 2         Jun et al.      Multi-ethnic PMID: 28183528
    ## 3 Kunkle et al. 2019 European ancestry PMID: 30820047
    ## 4  Bellenguez et al. European ancestry PMID: 35379992
    ## 5 Kunkle et al. 2021  African ancestry PMID: 33074286
    ##                                                   Site to download
    ## 1 https://www.FINNGEN.fi/en/access_results R8 accessed in Dec 2022
    ## 2                         https://www.niagads.org/datasets/ng00056
    ## 3                         https://www.niagads.org/datasets/ng00075
    ## 4                  https://www.ebi.ac.uk/gwas/studies/GCST90027158
    ## 5                         https://www.niagads.org/datasets/ng00100

The variants weight provided were selected from those in the complete
GWAS based on clumping parameter below, where we used the HCHS/SOL
TOPMed imputed Freeze5 dataset used in the paper as an LD reference
panel. To select the specific tuning parameter (LD parameters, p-value
threshold) for each AD PRS, we applied an approach where we optimized
the coefficient of variation (CV) computed over the estimated effect
sizes of the candidate PRS across 4 independent subset of the training
dataset. See manuscript for more detail.

The table below provides, for each trait-specific GWAS used, the
following information:

1.  GWAS\_name: GWAS name (which cohort/study the GWAS summary
    statistics are from)  
2.  Threshold: p-value threshold for selecting SNPs into the PRS  
3.  Distance: distance in kilo base-pairs used for clumping (SNPs were
    removed from consideration based on LD with other SNPs within a
    window of this distance)  
4.  R2: maximum LD for inclusion of SNPs within the distance-based
    window of another SNP that was already selected into the PRS.  
5.  SOL\_mean: the mean of the PRS after it was constructed in the SOL
    population. That is, each of the SOL participants had a PRS value.
    This is the mean of these values.  
6.  SOL\_sd: the standard deviation (SD) of the PRS after it was
    constructed in the SOL population. That is, each of the SOL
    participants had a PRS value. This is the SD of these values.  
7.  No\_SNP: Number of SNP selected from Clumping and Thresholding  

<!-- -->

    ##             GWAS_name Threshold Distance  R2     SOL_mean      SOL_sd No_SNP
    ## 1   Bellenguez et al.     1e-03   1000kb 0.1  0.000415514 0.000337105   1937
    ## 2             FinnGen     1e-06    250kb 0.1  0.002860917 0.007873218     81
    ## 3          Jun et al.     1e-05    250kb 0.2  0.014803928 0.007789044     85
    ## 4  Kunkle et al. 2019     1e-02   1000kb 0.3 -0.005668489 0.000503101  12002
    ## 5 Kunkle  et al. 2021     1e-04    250kb 0.3 -0.037393299 0.006307699    157

## PRSice command for PRS construction

This command is to construct PRS using the summary statistics that we
provide. No clumping is needed and no selection of SNPs. The summary
statistics are already based on the specific set of SNPs selected after
clumping and setting a p-value threshold. Note that genetic data files
need to be specified in the –target argument.


    Rscript ./PRSice.R \
     --dir ./PRS_Output \
     --prsice ./PRSice_linux/PRSice_linux \
     --base ./Variant_weight/. \
     --target ./Genotype \
     --thread 2 \
     --chr Chromosome 
     --bp Position 
     --A1 Allele1 
     --A2 Allele2 
     --pvalue PValue \
     --bar-levels Threshold \
     --stat BETA 
     --all-score T \
     --out ./out_prs \
     --no-clump T
     --print-snp T \
     --ignore-fid T 
     --no-regress T 
     --fastscore T 
     --model add 
     --no-full T 
     --chr-id c:l:a:b

## Constructing PRSsum based on multiple AD-PRS

After constructing multiple AD PRS, the AD-PRS is derived using the
PRSsum approach. This involves taking an unweighted sum of the scaled
trait-specific PRS. We utilize the mean and standard deviation (SD)
values from the SOL for each AD-specific PRS for scaling. Additionally,
we include the SOL mean and SD of the AD-PRS for the final scaling step.
Using the same scaling throughout guarantees that effect size estimates
are similarly interpreted across all datasets and individuals who use
this PRS.

    ##        SOL_mean   SOL_sd
    ## 1 -7.879378e-16 2.993377

See code below to construct PRSsum

    studies <- c("FINNGEN", "Jun","Kunkle","Kunkle_AFR","Bellenguez")
    out<-list()
    for(study in studies){
      
      
      prs_output <-paste0("./", study,".all_score")
      prs_df <-fread(prs_output, data.table=F)
      prs_df <- prs_df %>% dplyr::select(-IID)
      colnames(prs_df)<- c("sample.id", study)
      
      #standardize trait-prs using the mean and sd from TOPMed 
      cur_SOL_mean <- scaling_info$SOL_mean[which(scaling_info$GWAS_name == study)]
      cur_SOL_sd <- scaling_info$SOL_sd[which(scaling_info$GWAS_name == study)]
      
      prs_df[, study]<- (prs_df[, study] - cur_SOL_mean)/cur_SOL_sd
      out[[study]]<- prs_df
      
      
    }

    combine_prs <- purrr::reduce(out, full_join , by="sample.id")



    prssum<- data.frame(sample.id=prssum$sample.id, 
                        PRSsum=apply(prssum[,], 1, sum))

    prssum[,"PRSsum"]<- (prssum[,"PRSsum"] - SOL_AD_PRS_mean_sd$SOL_mean))/SOL_AD_PRS_mean_sd$SOL_sd

## Alternative to construct PRSsum

We offer an alternative approach for simplification by providing a
function to combine variant weights into a single overall variant
weight. This approach have slightly different when we compute PRSsum
like in we explain in the paper. Below are the instructions:

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

    ## variant weight will be matched using chromosome and position

    ## 23 alleles from FINNGEN are being flipped to align them with the reference SNP

    ## variant weight will be matched using chromosome and position

    ## 1 alleles from Jun are being flipped to align them with the reference SNP

    ## variant weight will be matched using chromosome and position

    ##  All the allele from Bellenguez are matched with the reference SNP

    ## variant weight will be matched using chromosome and position

    ## 8 alleles from Kunkle are being flipped to align them with the reference SNP

    ## variant weight will be matched using chromosome and position

    ##  All the allele from Kunkle_AFR are matched with the reference SNP

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

    #write.table(ad_unweigted_variants_weight, file = "Variant_weight/Unweighted_PRSsum.txt", row.names = F, quote = F, #sep="\t")

### Step 3: Construct PRS

Below is PRSice2 command to construct PRSsum

    Rscript ./PRSice.R \
     --dir ./PRS_Output \
     --prsice ./PRSice_linux/PRSice_linux \
     --base ./Variant_weight/Unweighted_PRSsum.txt \
     --target ./Genotype \
     --thread 2 \
     --chr Chromosome 
     --bp Position 
     --A1 Allele1 
     --A2 Allele2 
     --pvalue PValue \
     --bar-levels Threshold \
     --stat BETA 
     --all-score T \
     --out ./out_prs \
     --no-clump T
     --print-snp T \
     --ignore-fid T \
     --no-regress T \ 
     --fastscore T \ 
     --model add \
     --no-full T \
     --score sum \
     --chr-id c:l:a:b

next, we can scale the PRSsum using SOL final scaling.

    prssum<- fread("../PRS/Unweighted_PRSsum.all_score", data.table = F)
    prssum<- prssum %>% dplyr::select(-IID)
    colnames(prssum)<- c("sample.id", "PRSsum")

    prssum[, "PRSsum"]<- (prssum[, "PRSsum"]-SOL_AD_PRS_mean_sd$SOL_mean)/SOL_AD_PRS_mean_sd$SOL_sd

## Example code for association analysis

“We conducted an association analysis using mixed models, which were
implemented using the GENESIS R package, for all individuals.
Additionally, we offer an alternative method for conducting association
tests using linear regression specifically designed for unrelated
individuals. Below, you will find an example code for this alternative
approach, which utilizes functions provided in the ‘Code’ folder

    library(GENESIS)
    library(GWASTools)
    library(pROC)


    source("./Code/*")

    # add name of file of IDs of unrelated individuals as needed

    pheno<- fread(phenotype_file, data.table=F)


    # merge PRSsum with phenotype

    pheno_df<-left_join(pheno,prssum, by="person_id")


    covarites_prs<- c("age","sex","site","ethnic_background","education","apoe4","apoe2",paste0("PC_",1:5),"PRSsum")

    outcome<-"MCI"

    ## Kinship matrix

    covMatlist<-getobj(covMatlist)


    assoc_df<- run_assoc_mixmodel(pheno=pheno,
                                  outcome=outcome,
                                  covars_prs=covarites_prs, 
                                  covmat=covMatlist,
                                  group.var=NULL)

    assoc_df
