# Estimate PRS Uncertainty for Biobank
Here we demonstrate how to estimate PRS uncertainty (full posterior distribution of genetic value) as described in [Large uncertainty in individual PRS estimation impacts PRS-based risk stratification](https://doi.org/10.1101/2020.11.30.403188). 

## Data
Prepare the following data:
- `train_bfile`: Genotype file for training individuals. 
- `train_sumstats`: GWAS association summary statistics output from plink. 
- `val_bfile`: Genotype file for validation individuals.
- `val_pheno`: Phenotype file for validation individuals.
- `test_bfile`: Genotype file for testing individuals

## Estimate PRS Uncertainty
Install the following packages:
```{r}
install.packages('optparse', 'bigsnpr', 'tibble', 'readr', 'dplyr')
```
Download the `prs_uncertainty.R` and `summary.R` scripts:
              
Compute the PRS uncertainty:
```{shell}
#!/bin/sh
#$ -cwd
#$ -j y
#$ -l h_data=20G,h_rt=10:00:00 -pe shared 6
#$ -o ./job_out
#$ -t 1-22

chr_i=$SGE_TASK_ID

Rscript prs_uncertainty.R
    --chr_i=$chr_i
    --train_bfile=train_bfile \
    --train_sumstats=all.assoc.linear \
    --val_bfile=val_bfile \
    --val_pheno=height.val.regressed_pheno \
    --test_bfile=test_bfile \
    --out_dir=out_dir  #store outputfiles \
    --cache_dir=cache_file #store LD cache\ 
    --n_cores=6 \
    --num_burn_in=100 \
    --num_iter=500 \
```
To reduce runtime and memory required for training PRS models, we train PRS models separtely for each of the twentytwo chromosomes and then merge them to obtain a final PRS model for the whole genome using the script `summary.R`

```{shell}
#!/bin/sh
#$ -cwd
#$ -j y
#$ -l h_data=10G,h_rt=1:00:00 
#$ -o ./job_out

Rscript summary.R --out_dir=out_dir
```

