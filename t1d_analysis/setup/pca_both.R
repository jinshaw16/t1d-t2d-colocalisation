#pca_both.R

#producing the .sh file here and running it for both cohorts...
library(ggplot2)

pcarun<-function(cohortname){
sink(file=paste0("/users/todd/jinshaw/programs/t1d_risk/pca_calc_",cohortname,".sh"))
cat(paste0("#!/bin/bash
#pca_calc_",cohortname,".sh
           
#runs SNP QC and then pca for the independent individuals in this cohort combined for an initial gwas scan.
           
echo \"Analysis started at \"
date

for i in {1..22}
do
/users/todd/jinshaw/programs/t1d_risk/dropdups_pca.R ",cohortname," $i
plink --bfile /well/todd/users/jinshaw/t1d_risk/processed/",cohortname,
"_chr$i --exclude /well/todd/users/jinshaw/t1d_risk/processed/dropdup_",cohortname,
"_$i.txt --make-bed --out /well/todd/users/jinshaw/t1d_risk/processed/",cohortname,"_bi_chr$i
done

rm /well/todd/users/jinshaw/t1d_risk/processed/combinefile_",
cohortname,".txt

echo \"Combining all chromosomes:\"
for i in {2..22}
do
echo /well/todd/users/jinshaw/t1d_risk/processed/",cohortname,
"_bi_chr$i >> /well/todd/users/jinshaw/t1d_risk/processed/combinefile_",
cohortname,".txt
done
           
plink --bfile /well/todd/users/jinshaw/t1d_risk/processed/",cohortname,
"_bi_chr1 --merge-list /well/todd/users/jinshaw/t1d_risk/processed/combinefile_",cohortname,
".txt --out /well/todd/users/jinshaw/t1d_risk/processed/",cohortname,"_all


#remove MHC region 
plink --bfile /well/todd/users/jinshaw/t1d_risk/processed/",
cohortname,"_all --chr 6 --from-bp 25000000 --to-bp 35000000 --make-just-bim --out /well/todd/users/jinshaw/t1d_risk/processed/",
cohortname,"_independents --threads 32
plink --bfile /well/todd/users/jinshaw/t1d_risk/processed/",
cohortname,"_all --exclude /well/todd/users/jinshaw/t1d_risk/processed/",
cohortname,"_independents.bim --make-bed --out /well/todd/users/jinshaw/t1d_risk/processed/",
cohortname,"_nomhc
plink --bfile /well/todd/users/jinshaw/t1d_risk/processed/",
cohortname,"_nomhc --indep-pairwise 50 5 0.2 --autosome --threads 32 --out /well/todd/users/jinshaw/t1d_risk/processed/",
cohortname,"_independents_unfiltered
           
#and actually remove the pruned SNPs:
plink --bfile /well/todd/users/jinshaw/t1d_risk/processed/",cohortname,
"_nomhc --extract /well/todd/users/jinshaw/t1d_risk/processed/",
cohortname,"_independents_unfiltered.prune.in --out /well/todd/users/jinshaw/t1d_risk/processed/",
cohortname,"_independents_filtered --make-bed
           
#now perform PCA:
plink --bfile /well/todd/users/jinshaw/t1d_risk/processed/",cohortname,
"_independents_filtered --pca --out /well/todd/users/jinshaw/t1d_risk/processed/pcafiles/",cohortname,
"pcs_independents_filtered
           
echo \"Analysis finished at \"
date"))
sink()


system(paste0("chmod a=rwx /users/todd/jinshaw/programs/t1d_risk/pca_calc_",cohortname,".sh"))
system(paste0("bash /users/todd/jinshaw/programs/t1d_risk/pca_calc_",cohortname,".sh"))



co<-read.table(file=paste0("/well/todd/users/jinshaw/t1d_risk/processed/",cohortname,"_independents.fam"))
colnames(co)<-c("pedigree","member","mother","father","sex","affected")
vecs<-read.table(file=paste0("/well/todd/users/jinshaw/t1d_risk/processed/pcafiles/",cohortname,"pcs_independents_filtered.eigenvec"), header=F, as.is=T)
labs<-paste0("PC",1:20)
colnames(vecs)<-c("pedigree","member",labs)
co<-merge(co, vecs, by=c("pedigree","member"))

png(file=paste0("~/output/t1d_risk/hrc/pca/one_vs_two_",cohortname,".png"),
    res=500, units="cm", height=20, width=20)
ggplot(data=co, aes(PC1, PC2)) + geom_point()
dev.off()
}

pcarun("affy")
pcarun("illu")
