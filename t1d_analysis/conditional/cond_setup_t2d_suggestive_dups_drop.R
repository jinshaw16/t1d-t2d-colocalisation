#cond_setup_t2d_suggestive_dups_drop.R
#read all files in to R and convert to dosage format for regression analysis
library(vcfR)
library(jimisc)
library(parallel)

#since this is going to computationally intensive, just read in chromosomes that are significantly asosicated 
#index variant plus 2Mb region around the index variant

load(file="/well/todd/users/jinshaw/t2d_risk/regionsoverlapping_dups_drop.RData")
signals<-l
signals$position<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\2",signals$idmin)
signals$position<-as.numeric(signals$position)
#keeping only those not genomewide associated with T1D since have done the genomewide ones in cond_setup_all_dups_drop.R
signals<-signals[signals$t1dminp>5*10^-8,]
#and remove other known T1D genes:
signals<-signals[!signals$gene %in% c("CENPW"),]

load(file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/resinit_post_imp_qc_adjusted_newml_5pc_diff_dups_drop.RData")

#use BCF tools to keep only the SNPs of interest and also to include SNPs from the 1000 genomes project that aren't on the HRC panel (indels + some other SNPs):
f<-function(snp, cohort){
s<-signals[signals$IDmin==snp,]
min<-s$position-1000000
max<-s$position+1000000
chr<-s$chrom
pos<-s$position
sink(file=paste0("~/programs/t1d_risk/vcf/conditional/",snp,"_get_",cohort,"_all.sh"))
cat(paste0("#",snp,"_get_",cohort,"_all.sh\n\n"))
cat(paste0("/apps/well/bcftools/1.4.1/bin/bcftools view -r ",chr,
":",min,"-",max," /well/todd/users/jinshaw/t1d_risk/T1D.",cohort,".chr",chr,
".vcf.gz > /well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/",snp,
"_out",cohort,".vcf \nbgzip -f /well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/",
snp,"_out",cohort,".vcf \ntabix -f /well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/",snp,
"_out",cohort,".vcf.gz \n"))
sink()
system(paste0("chmod a=rwx ~/programs/t1d_risk/vcf/conditional/",snp,"_get_",cohort,"_all.sh"))
create_job(path=paste0("~/programs/t1d_risk/vcf/conditional/"),subname=paste0(snp,"_get_",cohort,"_all"),jobname=paste0(snp,"_get_",cohort,"_all"),R=F,
projectletter="c",qletter="c",qlength="short")
system(paste0("qsub ~/programs/t1d_risk/vcf/conditional/",snp,"_get_",cohort,"_all"))
}
lapply(signals$IDmin, f, cohort="Affy")
lapply(signals$IDmin, f, cohort="Illumina")


