#meta_sardinia_all_dups_dropped.R

library(qqman)
library(ggplot2)
library(jimisc)
library(humarray)
library(gridExtra)
library(ggbio)
library(gridExtra)
#combining the sardinian imputed results wit hthe HRC imputed results from UK case-control.

#load sardinian summary stats:
readitin<-function(chrom){
l<-read.table(file=paste0("/well/todd/users/jinshaw/t1d_risk/sardinia/TD.chr",chrom,".epacts.gz"),header=T,as.is=T, comment.char="")
l$dummy<-gsub("/.*","",l$MARKER_ID)
l$dummy1<-gsub(".*/","",l$MARKER_ID)
l$allele1<-gsub(".*_","",l$dummy)
l$allele2<-gsub("_.*","",l$dummy1)
l$id<-gsub("_.*","",l$MARKER_ID)
l$id<-paste0(l$id,":",l$allele1,":",l$allele2)
l$chromosome<-l$X.CHROM
l$position<-l$START
l<-l[,c("id","allele1","allele2","CALLRATE","MAF","PVALUE","BETA","SEBETA","CHISQ","NS.CASE",
"NS.CTRL","AF.CASE","AF.CTRL")]
l<-l[l$MAF>0.01,]
message(paste0("DONE CHROMOSOME ",chrom))
return(l)
}

#sard<-lapply(c(1:22),readitin)
#sard<-do.call("rbind",sard)
#save(sard, file="/well/todd/users/jinshaw/t1d_risk/sardinia/results_init.RData")
load(file="/well/todd/users/jinshaw/t1d_risk/sardinia/results_init.RData")

load(file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/resinit_post_imp_qc_adjusted_newml_5pc_diff_dups_dropped.RData")
meta$alleleA<-meta$alleleA_Illumina
meta$alleleB<-meta$alleleB_Illumina
s<-merge(meta,sard, by="id", all.x=T)
#now for the few SNPs where the different reference panels used different strands:

notin<-meta[!meta$id %in% sard$id,]
notin$did<-paste0(notin$chromosome,":",notin$position)
sard$did<-sub("(^.*:.*)[:](.*)[:](.*)","\\1",sard$id)
notin<-notin[notin$did %in% sard$did,]
notinsard<-sard[sard$did %in% notin$did,]
look<-merge(notin,notinsard,by="did")

look$switch<-ifelse(look$alleleA=="A" & look$alleleB=="C" &
look$allele1=="T" & look$alleleB=="G",1,
ifelse(look$alleleA=="A" & look$alleleB=="G" &
look$allele1=="T" & look$alleleB=="C",2,
ifelse(look$alleleA=="T" & look$alleleB=="C" &
look$allele1=="A" & look$alleleB=="G",3,
ifelse(look$alleleA=="T" & look$alleleB=="G" &
look$allele1=="A" & look$alleleB=="C",4,
ifelse(look$alleleA=="G" & look$alleleB=="A" &
look$allele1=="C" & look$alleleB=="T",5,
ifelse(look$alleleA=="G" &	look$alleleB=="T" &
look$allele1=="C" & look$alleleB=="A",6,
ifelse(look$alleleA=="C" & look$alleleB=="A" &
look$allele1=="G" & look$alleleB=="T",7,
ifelse(look$alleleA=="C" &	look$alleleB=="T" &
look$allele1=="G" & look$alleleB=="A",8,0))))))))
table(look$switch)
#can't get any more SNPs here - looks like all were done on the same strand.



#merge the datasets together:
comb<-s
message(paste0(nrow(comb)," SNPs included"))
p<-nrow(comb)
comb$dummy<-paste0(comb$alleleA,":",comb$alleleB)
comb$dummy1<-paste0(comb$allele1,":",comb$allele2)
comb$strand<-ifelse(comb$alleleA!=comb$allele1 & comb$alleleA!=comb$allele2,1,0)
look<-comb[comb$strand==1,]
table(look$dummy,look$dummy1)
#all appear to be allgined to the same strand - nice and easy to proceed with meta analysis:

comb$var_uk<-comb$seall^2
comb$var_sard<-comb$SEBETA^2
comb$est_uk<-comb$beta
comb$est_sard<-comb$BETA

weightcalc<-function(cohort,name){
  cohort[,paste0("weight_",name)]<-1/cohort[,paste0("var_",name)]
  cohort[,paste0("betawt_", name)]<-cohort[,paste0("est_",name)]*cohort[,paste0("weight_",name)]
  return(cohort)
}
comb<-weightcalc(comb, "uk")
comb<-weightcalc(comb, "sard")


comb$weightsum<-rowSums(comb[,c("weight_uk", "weight_sard")], na.rm=TRUE)
comb$sefinal<-sqrt(1/comb$weightsum)
comb$sumfinal<-rowSums(comb[,c("betawt_uk", "betawt_sard")], na.rm=TRUE)
comb$betafinal<-comb$sumfinal/comb$weightsum
comb$zfinal<-comb$betafinal/comb$sefinal
comb$pfinal<-2*pnorm(-abs(comb$zfinal))
comb<-comb[order(comb$chromosome, comb$position),]
comb$p_uk<-comb$pmeta
comb$p_sard<-comb$PVALUE
comb$all_maf_sard<-comb$MAF
comb$shouldbe<-((comb$NS.CASE/(comb$NS.CASE+comb$NS.CTRL))*comb$AF.CASE) + ((comb$NS.CTRL/(comb$NS.CASE+comb$NS.CTRL))*comb$AF.CTRL)
comb$cases_maf_sard<-ifelse(abs(comb$shouldbe-comb$all_maf_sard)<0.01,comb$AF.CASE, 1-comb$AF.CASE)
comb$controls_maf_sard<-ifelse(abs(comb$shouldbe-comb$all_maf_sard)<0.01,comb$AF.CTRL, 1-comb$AF.CTRL)
comb$minor_sard<-ifelse(abs(comb$shouldbe-comb$all_maf_sard)<0.01,comb$allele2, comb$allele1)
comb$se_uk<-comb$seall
comb$se_sard<-comb$SEBETA
#get the weighted allele frequency of the B allele over UK and Sardinian cohorts:
af<-read.table(file="/well/todd/users/jinshaw/t1d_risk/processed/vcf/affy/affy.sample", skip=2, as.is=T)
il<-read.table(file="/well/todd/users/jinshaw/t1d_risk/processed/vcf/illu/illu.sample",	skip=2,	as.is=T)
nuk<-sum(nrow(af[af$V4==0,]), nrow(il[il$V4==0,]))
nsard<-max(comb$NS.CTRL,na.rm=T)
weight_uk<-nuk/(nuk+nsard)
weight_sard<-nsard/(nuk+nsard)
comb$controls_AF_B_all<-ifelse(!is.na(comb$controls_AF_B) & !is.na(comb$AF.CTRL),
weight_uk*comb$controls_AF_B + weight_sard*comb$AF.CTRL,
ifelse(!is.na(comb$controls_AF_B) & is.na(comb$AF.CTRL), comb$controls_AF_B,
ifelse(is.na(comb$controls_AF_B) & !is.na(comb$AF.CTRL), comb$AF.CTRL,NA)))
comb$controls_AF_B<-comb$controls_AF_B_all
comb$all_info_Illumina<-ifelse(is.na(comb$est_Illumina),NA,comb$all_info_Illumina)
comb$all_info_Affy<-ifelse(is.na(comb$est_Affy),NA,comb$all_info_Affy)
comb$meanr2<-rowMeans(comb[,c("all_info_Illumina","all_info_Affy")],na.rm=T)

final<-comb[,c("id","ID","chromosome","position","alleleA","alleleB","est_Affy","est_Illumina","se_Affy","se_Illumina","cases_maf_Affy",
"controls_maf_Affy","cases_maf_Illumina","controls_maf_Illumina","all_maf_Affy","all_maf_Illumina",
"cases_info_Affy","controls_info_Affy","all_info_Affy","cases_info_Illumina","controls_info_Illumina","all_info_Illumina",
"est_uk","se_uk","p_uk","est_sard","se_sard","p_sard","all_maf_sard","cases_maf_sard", "controls_maf_sard","betafinal",
"sefinal", "pfinal","meanr2", "minor_sard","controls_AF_B")]
save(final, file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/meta_sardinia_adjusted_newml_5pc_diff_dups_dropped.RData")

load(file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/meta_sardinia_adjusted_newml_5pc_diff_dups_dropped.RData")
#export to a text file:
write.table(final,file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/summary_stats_af_sard_newml_5pc_diff_dups_dropped.txt", col.names=T, quote=F, row.names=F, sep="\t")
meta<-final
meta<-meta[!is.na(meta$pfinal) & meta$pfinal!=0,]
meta$Chisq <- qchisq(meta$pfinal,1, lower=F)
lam<-median(meta$Chisq)/qchisq(0.5,1)
max=-1*log10(min(meta$pfinal))

png(file="/well/todd/users/jinshaw/output/t1d_risk/hrc/vcf_all/all_pic_post_qc_adjusted_5pc_diff_dups_dropped.png",
    res=200, width=15, height=18, units="cm")
par(mfrow=c(2,1))
manhattan(meta, p="pfinal", bp="position", chr="chromosome",ylim=c(0,20))
qq(meta$pfinal)
text(x=1,y=max*0.8, labels=expression(paste(lambda,"=")))
text(x=1.5, y=max*0.8, labels=paste0(round(lam,2)))
dev.off()



#Export Table of all associated (index) SNPs
#New regions defined as anything >1/2 Megabase from first index SNP. Excluding MHC:
load(file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/meta_sardinia_adjusted_newml_5pc_diff_dups_dropped.RData")

final<-final[!(final$chromosome==6 & final$position>25000000 & final$position<35000000),]

getmin<-function(frame){
min<-frame[frame$pfinal==min(frame$pfinal),]
return(min)
}

dropnext<-function(min, results){
droplow<-min$position-500000
drophigh<-min$position+500000
chrom<-min$chromosome
m<-results[!(results$chromosome==chrom & results$position>droplow & results$position<drophigh),]
return(m)
}

n<-getmin(final)
dropped<-dropnext(n, final)
both<-n

while(n$pfinal<5*10^-8){
n<-getmin(dropped)
dropped<-dropnext(n,dropped)
both<-rbind(both, n)
}

both<-both[order(both$chromosome, both$position),]
both$or<-exp(both$beta)
both$or_sard<-exp(both$est_sard)
both$or_uk<-exp(both$est_uk)
b<-both[,c("ID","chromosome", "position","alleleA","alleleB","or_sard", "or_uk", "or", "p_sard", "p_uk", "pfinal")]
b$or<-round(b$or, digits=2)
b$or_sard<-round(b$or_sard, digits=2)
b$or_uk<-round(b$or_uk, digits=2)
b<-b[order(b$chromosome,b$position),]
b$pfinal<-format(b$pfinal,digits=3, scientific=T)
b$p_sard<-format(b$p_sard,digits=3,scientific=T)
b$p_uk<-format(b$p_uk,digits=3,scientific=T)

#there are a number of SNPs that are likely only here due to LD with the main signla: dropping these now (though formal testing might be required in the future - tricky with the sardinians as only summary stats avaialble)
#b<-b[!b$ID %in% c("rs72685677","rs12828640","rs10850001"),]


write.table(b,file="/well/todd/users/jinshaw/output/t1d_risk/hrc/vcf_all/sigregions_5pc_diff_dups_dropped.txt", col.names=T, row.names=F, quote=F, sep=";")


