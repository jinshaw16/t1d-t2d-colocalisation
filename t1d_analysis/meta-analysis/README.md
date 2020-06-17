# T1D meta-analysis

# Scripts and running order

## UK association analyses:

```
assoctest_newml_all.R
```
Runs SNPTEST on every variant for Affy and Illu, adjusting for top 3 PCs.


```
meta-analysis_all_dups_drop.R
```
Reads in the SNPTEST results and meta-analyses between collections. Light QC to remove variants and reduce file size


```
qc_post_imputation_all_dups_drop.R
```
More in-depth SNP QC


## UK + Sardinia meta-analysis
```
meta_sardinia_all_dups_drop.R
```
Carrys out the meta analysis with the Sardinia data.

