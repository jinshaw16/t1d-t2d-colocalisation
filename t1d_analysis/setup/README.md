# T1D analysis set-up

# Scripts and running order

## Collating genotype data:
```
creating_genotype_data_anuhba.R
```
Reads in the genotype data from Illu and Affy collections & alligns to genome build 37

```
generate_outcome.R
```
Regenerates the outcome file for Illu and Affy so can trace back all the way.


## Imputation then performed on the Michigan Imputation server using the genotype data aligned to the same strand as the 1000 genomes strand.

## PCA

```
rename_variants_Affy.sh	
```
Renames the variants in the VCF file to chrom:pos:alleleA:alleleB for Affy collection


```
rename_variants_Illumina.sh
```
Renames the variants in the VCF file to chrom:pos:alleleA:alleleB for Illumina collection


```
convert_hrc_to_plink_Affy.sh
```
Generates PLINK files with variants filtered at r2<0.3 and maf>0.01 for Affy cohort


```
convert_hrc_to_plink_illu.sh
```
Generates PLINK files with variants filtered at r2<0.3 and maf>0.01 for Illumina cohort


```
pca_both.R (uses dropdups_pca.R)
```
Combines chromosomes, prunes variants, performs PCA in both cohorts.




