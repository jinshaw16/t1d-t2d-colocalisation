#assoctest_newml_all.R
library(jimisc)

#generating snptest scripts to run on HPC (adjusted test):
f1<-function(cohort, chromosome){
sink(file=paste0("~/programs/t1d_risk/vcf/snptest_scripts/",cohort,"/chr_",chromosome,"_adjusted_all.sh"))
cat(paste0("/apps/well/snptest/2.5.4-beta3_CentOS6.6_x86_64_dynamic/snptest_v2.5.4-beta3 -data /well/todd/users/jinshaw/t1d_risk/T1D.",cohort,".chr",
chromosome,".vcf.gz /well/todd/users/jinshaw/t1d_risk/processed/vcf/",cohort,"/",cohort,
".sample -filetype vcf -genotype_field GP -pheno phenotype -frequentist 1 -method newml -cov_all_continuous "))q
if(cohort=="Affy"){
cat(paste0("-include_samples /well/todd/users/jinshaw/t1d_risk/processed/vcf/affy/include_affy "))
}
cat(paste0("-o /well/todd/users/jinshaw/t1d_risk/results/vcf_all/",
cohort,"/chr_",chromosome,"_adjusted"))
sink()

system(paste0("chmod a=rwx ~/programs/t1d_risk/vcf/snptest_scripts/",cohort,"/chr_",chromosome,"_adjusted_all.sh"))
create_job(path=paste0("/users/todd/jinshaw/programs/t1d_risk/vcf/snptest_scripts/",cohort,"/"),
subname=paste0("chr_",chromosome,"_adjusted_all"), jobname=paste0("c",chromosome,"_adj_all"),qletter="c", projectletter="c", qlength="short")
system(paste0("qsub ~/programs/t1d_risk/vcf/snptest_scripts/",cohort,"/chr_",chromosome,"_adjusted_all"))
}

junk<-lapply(c(1:22),f1,cohort="Affy")
junk1<-lapply(c(1:22),f1,cohort="Illumina")

