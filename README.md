# t1d-t2d-colocalisation
repository with code to perform analysis of t1d and t2d overlap analysis

# Scripts and running order

Refer to the following READMEs for guidance on the running order for each part of the analysis.

## t1d analysis
```
./t1d_analysis/setup/README.md
```
Analysis set up.

```
./t1d_analysis/meta-analysis/README.md
```
Peforms the association analysis and meta-analysis across collections, including UK and Sardinian cohort.

```
./t1d_analysis/conditional/README.md
```
Performs conditional analyses for all T1D FDR<0.01 associated regions


## t2d analysis
Summary stats used (unadjusted for BMI) from Mahajan et. al, Nature Genetics 2018  doi: 10.1038/s41588-018-0241-6


## t1d t2d colocalisation analyses
```
./t1d_t2d_colocalisation/README.md
```
Identifies regions of overlap (FDR <0.01 to both diseases) and examines colocalisations between diseases in those regions.
