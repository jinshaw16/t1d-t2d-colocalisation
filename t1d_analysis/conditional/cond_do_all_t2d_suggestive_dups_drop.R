#cond_do_all_t2d_suggestive_dups_drop.R

#creating scripts to run in SNPTEST and get conditional regression results:
library(jimisc)
library(parallel)

load(file="/well/todd/users/jinshaw/t2d_risk/regionsoverlapping.RData")
signals<-l
signals$position<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\2",signals$idmin)
signals$position<-as.numeric(signals$position)
#keeping only those not genomewide associated with T1D since have done the genomewide ones in cond_setup_5pc_diff.R
signals<-signals[signals$t1dminp>5*10^-8,]
#and remove other known T1D genes:
signals<-signals[!signals$gene %in% c("CENPW"),]
signals$ID<-signals$IDmin
signals$chromosome<-signals$chrom

load(file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/resinit_post_imp_qc_adjusted_newml_5pc_diff.RData")
metaorig<-meta
metaorig$rsid<-metaorig$id

origaf<-metaorig[!is.na(metaorig$est_Affy),c("ID","rsid")]
origil<-metaorig[!is.na(metaorig$est_Illumina),c("ID","rsid")]
orig<-metaorig[,c("ID","rsid")]


signals$rsid<-paste0(signals$chromosome,":",signals$position)
signals$conditional=""
signals$index<-signals$rsid

#generate the snptest scripts:
snptestscripts<-function(snp,cohort, adj=NULL){
s<-signals[signals$ID==snp,]
id<-s$id
if (is.null(adj)){
sink(file=paste0("/users/todd/jinshaw/programs/t1d_risk/vcf/conditional/snptest_scripts/",cohort,"/",snp,"_do_all.sh"))
cat(paste0("/apps/well/snptest/2.5.4-beta3_CentOS6.6_x86_64_dynamic/snptest_v2.5.4-beta3 -data ", 
"/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/",snp,"_out",cohort,".vcf.gz ",
"/well/todd/users/jinshaw/t1d_risk/processed/vcf/",cohort,"/",cohort,".sample -filetype vcf -genotype_field GP -pheno phenotype ",
"-frequentist 1 -method newml -condition_on ",id," "))
if(cohort=="Affy"){
cat(paste0("-include_samples /well/todd/users/jinshaw/t1d_risk/processed/vcf/affy/include_affy "))
}

cat(paste0("-o /well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_",cohort,"_all \n"))
sink()
system(paste0("chmod a=rwx /users/todd/jinshaw/programs/t1d_risk/vcf/conditional/snptest_scripts/",cohort,"/",snp,"_do_all.sh"))
create_job(path=paste0("~/programs/t1d_risk/vcf/conditional/snptest_scripts/",cohort,"/"),subname=paste0(snp,"_do_all"),
jobname=paste0(snp,"_do_all"),R=F,projectletter="c",qletter="c",qlength="short")
system(paste0("qsub ~/programs/t1d_risk/vcf/conditional/snptest_scripts/",cohort,"/",snp,"_do_all"))
}
if (!is.null(adj)){
j<-strsplit(adj, split="+", fixed=T)
n<-length(j[[1]])
adjusted<-paste0(j[[1]]," add")
adjusted<-paste(adjusted, collapse=" ")
sink(file=paste0("/users/todd/jinshaw/programs/t1d_risk/vcf/conditional/snptest_scripts/",cohort,"/",snp,"_do_all_",n,".sh"))
cat(paste0("/apps/well/snptest/2.5.4-beta3_CentOS6.6_x86_64_dynamic/snptest_v2.5.4-beta3 -data ",
"/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/",snp,"_out",cohort,".vcf.gz ",
"/well/todd/users/jinshaw/t1d_risk/processed/vcf/",cohort,"/",cohort,".sample -filetype vcf -genotype_field GP -pheno phenotype ",
"-frequentist 1 -method newml -condition_on ",adjusted," "))
if(cohort=="Affy"){
cat(paste0("-include_samples /well/todd/users/jinshaw/t1d_risk/processed/vcf/affy/include_affy "))
}
cat(paste0("-o /well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_",cohort,"_do_all_",n,"\n"))
sink()
system(paste0("chmod a=rwx /users/todd/jinshaw/programs/t1d_risk/vcf/conditional/snptest_scripts/",cohort,"/",snp,"_do_all_",n,".sh"))
create_job(path=paste0("~/programs/t1d_risk/vcf/conditional/snptest_scripts/",cohort,"/"),subname=paste0(snp,"_do_all_",n),
jobname=paste0(snp,"_do_all_",n),R=F,
projectletter="c",qletter="c",qlength="short")
system(paste0("qsub ~/programs/t1d_risk/vcf/conditional/snptest_scripts/",cohort,"/",snp,"_do_all_",n))
}
}

#second signals?
lapply(signals$ID,snptestscripts, cohort="Affy")
lapply(signals$ID,snptestscripts, cohort="Illumina")


metathem<-function(snp, adj=NULL){
if(is.null(adj)){
if(file.exists(paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Affy_all"))){
af<-read.table(file=paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Affy_all"), skip=10, header=T, as.is=T)
afin=T
}
if(!file.exists(paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Affy_all"))){
afin=F
}

if(file.exists(paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Illumina_all"))){
il<-read.table(file=paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Illumina_all"), skip=10, header=T, as.is=T)
ilin=T
}
if(!file.exists(paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Illumina_all"))){
ilin=F
}
}
if(!is.null(adj)){
j<-strsplit(adj, split="+", fixed=T)
n<-length(j[[1]])
if(file.exists(paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Affy_do_all_",n))){
af<-read.table(file=paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Affy_do_all_",n), skip=10, header=T, as.is=T)
afin=T
}
if(!file.exists(paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Affy_do_all_",n))){
afin=F
}
if(file.exists(paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Illumina_do_all_",n))){
il<-read.table(file=paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Illumina_do_all_",n), skip=10, header=T, as.is=T)
ilin=T
}
if(!file.exists(paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Illumina_do_all_",n))){
ilin=F
}
}

if(afin==T & ilin==T){
#remove	SNPs that were dropped from the	original meta analysis due to QC:
af<-af[af$rsid %in% origaf$rsid,]
il<-il[il$rsid %in% origil$rsid,]

meta<-merge(af,il,by=c("alternate_ids","rsid","chromosome","position","alleleA","alleleB"))
meta$diff.x<-abs(meta$cases_info.x-meta$controls_info.x)
meta$diff.y<-abs(meta$cases_info.y-meta$controls_info.y)
meta<-meta[meta$diff.x<0.05 & meta$diff.y<0.05,]

meta$weight.x<-1/(meta$frequentist_add_se_1.x^2)
meta$betawt.x<-meta$weight.x*meta$frequentist_add_beta_1.add.phenotype.1.x

meta$weight.y<-1/(meta$frequentist_add_se_1.y^2)
meta$betawt.y<-meta$weight.y*meta$frequentist_add_beta_1.add.phenotype.1.y

meta$weightsum<-rowSums(meta[,c("weight.x", "weight.y")], na.rm=TRUE)
meta$seall<-sqrt(1/meta$weightsum)
meta$sumbeta<-rowSums(meta[,c("betawt.x", "betawt.y")], na.rm=TRUE)
meta$beta<-meta$sumbeta/meta$weightsum
meta$z<-meta$beta/meta$seall
meta$pmeta<-2*pnorm(-abs(meta$z))
m<-meta[meta$pmeta==min(meta$pmeta,na.rm=T),]
m<-m[!is.na(m$rsid),]
}
if(afin==F & ilin==T){
meta=il[il$rsid %in% orig$rsid,]
meta$diff<-abs(meta$cases_info-meta$controls_info)
meta<-meta[meta$diff<0.05,]
meta$beta<-meta$frequentist_add_beta_1.add.phenotype.1
meta$z<-meta$beta/meta$frequentist_add_se_1
meta$pmeta<-2*pnorm(-abs(meta$z))
m<-meta[meta$pmeta==min(meta$pmeta,na.rm=T),]
m<-m[!is.na(m$rsid),]
}
if(afin==T & ilin==F){
meta=af[af$rsid %in% orig$rsid,]
meta$diff<-abs(meta$cases_info-meta$controls_info)
meta<-meta[meta$diff<0.05,]
meta$beta<-meta$frequentist_add_beta_1.add.phenotype.1
meta$z<-meta$beta/meta$frequentist_add_se_1
meta$pmeta<-2*pnorm(-abs(meta$z))
m<-meta[meta$pmeta==min(meta$pmeta,na.rm=T),]
m<-m[!is.na(m$rsid),]
}

m<-merge(m,orig, by="rsid")
m$or<-exp(m$beta)
m$p_uk<-m$pmeta
id<-orig[orig$ID==snp,"rsid"]
if(is.null(adj)){
m$conditional<-paste0(id)
}
if(!is.null(adj)){
m$conditional<-paste0(adj)
}
m$index<-id
m<-m[,c("ID","rsid","chromosome","position","alleleA","alleleB","or","p_uk","conditional","index")]
message(paste0("Done ",snp))
return(m)
}

secondsigs<-mclapply(signals$ID, metathem, mc.cores=6)
secondsigs<-do.call("rbind", secondsigs)

signals$conditional=""
signals$alleleA<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\3",signals$idmin)
signals$alleleB<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\4",signals$idmin)
signals$or<-exp(signals$beta)
signals$p_uk<-signals$t1dminp
signals$index<-signals$idmin
signals<-signals[,colnames(signals) %in% colnames(secondsigs)]

signalsall<-rbind(signals,secondsigs)
secondsigs$newsig<-paste0(secondsigs$conditional,"+",secondsigs$rsid)
secondsigs<-secondsigs[secondsigs$p_uk<5*10^-5,]
secondsigs<-secondsigs[secondsigs$conditional!=secondsigs$rsid,]
#third signals?
s<-signals[signals$index %in% secondsigs$index,]

mapply(snptestscripts,snp=s$ID, adj=secondsigs$newsig, cohort="Affy")
mapply(snptestscripts,snp=s$ID, adj=secondsigs$newsig, cohort="Illumina")


thirdsigs<-mcmapply(metathem, snp=s$ID, adj=secondsigs$newsig, SIMPLIFY=FALSE, mc.cores=6)
thirdsigs<-do.call("rbind", thirdsigs)

signalsall<-rbind(signalsall,thirdsigs)
thirdsigs<-thirdsigs[thirdsigs$p_uk<5*10^-5,]
signalsall<-signalsall[signalsall$p_uk<5*10^-5,]

save(signalsall, file="/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/all_signals_t2d_suggestive_dups_drop.RData")
write.table(signalsall, file="/well/todd/users/jinshaw/output/t1d_risk/hrc/vcf/conditional/all_signals_t1d_suggesgtive_dups_drop.txt",
col.names=T, row.names=F, quote=F, sep="\t")
