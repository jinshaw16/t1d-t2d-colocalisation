#meta-analysis_all_dups_dropped.R
#reads in the results from SNPtest and performs the meta analysis to get the final results and saves
library(qqman)
library(ggplot2)
library(jimisc)
library(humarray)
library(gridExtra)
library(ggbio)


#read in data from both cohorts:


getdat<-function(cohort, chromosome){
message(paste0("reading in chromosome ",chromosome," - ", cohort))
r<-read.table(file=paste0("/well/todd/users/jinshaw/t1d_risk/results/vcf_all/",cohort,"/chr_",chromosome,"_adjusted"), skip=13, header=T, as.is=T)

r$a0<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\3",r$rsid)
r$a1<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\4",r$rsid)
w<-which(r$a0==r$alleleB & r$a1==r$alleleA)

message(paste0(length(w), " allele switches during SNPTEST procedure - do nothing"))
r$id<-r$alternate_ids
r[,paste0("alleleA_",cohort)]<-r$alleleA
r[,paste0("alleleB_",cohort)]<-r$alleleB
r[,paste0("all_info_",cohort)]<-r$all_info
r[,paste0("cases_info_",cohort)]<-r$cases_info
r[,paste0("controls_info_",cohort)]<-r$controls_info
r[,paste0("all_maf_",cohort)]<-r$all_maf
r[,paste0("cases_maf_",cohort)]<-r$cases_maf
r[,paste0("controls_maf_",cohort)]<-r$controls_maf

r[,paste0("est_",cohort)]<-r$frequentist_add_beta_1.add.phenotype.1
r[,paste0("se_",cohort)]<-r$frequentist_add_se_1
r[,paste0("p_",cohort)]<-r$frequentist_add_wald_pvalue_1
r[,paste0("var_",cohort)]<-r$frequentist_add_se_1^2
r[,paste0("controls_AA_",cohort)]<-r$controls_AA
r[,paste0("controls_AB_",cohort)]<-r$controls_AB
r[,paste0("controls_BB_",cohort)]<-r$controls_BB
r[,paste0("weight_",cohort)]<-1/r[,paste0("var_",cohort)]
r[,paste0("betawt_", cohort)]<-r[,paste0("est_",cohort)]*r[,paste0("weight_",cohort)]

#remove SNPs with an info difference of >5% between cases and controls:
r[,paste0(cohort,"diff")]<-abs(r[,paste0("cases_info_",cohort)]-r[,paste0("controls_info_",cohort)])

w<-which(r[,paste0(cohort,"diff")]>0.05 & r$alternate_ids!="11:2182224:A:T")
r[w,paste0("est_",cohort)]<-NA
r[w,paste0("se_",cohort)]<-NA
r[w,paste0("betawt_",cohort)]<-NA
r[w,paste0("weight_",cohort)]<-NA
r[w,paste0("var_",cohort)]<-NA
r[w,paste0("p_",cohort)]<-NA
#and any SNPs where either the cases or the controls have an info of <0.3:
w<-which(r[,paste0("cases_info_",cohort)]<0.3 | r[,paste0("controls_info_",cohort)]<0.3)
r[w,paste0("est_",cohort)]<-NA
r[w,paste0("se_",cohort)]<-NA
r[w,paste0("betawt_",cohort)]<-NA
r[w,paste0("weight_",cohort)]<-NA
r[w,paste0("var_",cohort)]<-NA
r[w,paste0("p_",cohort)]<-NA

r<-r[,c("id","chromosome","position",paste0("alleleA_",cohort),paste0("alleleB_",cohort),
paste0("all_info_",cohort),paste0("all_maf_",cohort),paste0("cases_maf_",cohort),
paste0("controls_maf_",cohort),paste0("est_",cohort),
paste0("se_",cohort),paste0("p_",cohort),paste0("var_",cohort),
paste0("cases_info_",cohort), paste0("controls_info_",cohort),
paste0("weight_",cohort), paste0("betawt_", cohort),
paste0("controls_AA_",cohort),paste0("controls_AB_",cohort),paste0("controls_BB_",cohort))]
message(paste0("Done chromosome ",chromosome))
return(r)
}

illu<-lapply(c(1:22),getdat, cohort="Illumina")
illu<-do.call("rbind",illu)

affy<-lapply(c(1:22),getdat, cohort="Affy")
affy<-do.call("rbind",affy)
#carry out the meta-analysis:
meta<-merge(illu, affy, by=c("id","chromosome","position"), all=TRUE)

meta$controls_AF_B<-ifelse(!is.na(meta$controls_AA_Illumina) & !is.na(meta$controls_AA_Affy),
(meta$controls_AB_Illumina + (2*meta$controls_BB_Illumina) + meta$controls_AB_Affy + (2*meta$controls_BB_Affy)) / 
((2*meta$controls_AA_Illumina) + (2*meta$controls_AB_Illumina) + (2*meta$controls_BB_Illumina) + 
(2*meta$controls_AA_Affy) + (2*meta$controls_AB_Affy) + (2*meta$controls_BB_Affy)),
ifelse(is.na(meta$controls_AA_Illumina) & !is.na(meta$controls_AA_Affy), (meta$controls_AB_Affy + 
(2*meta$controls_BB_Affy))/((2*meta$controls_AA_Affy) + (2*meta$controls_AB_Affy) + (2*meta$controls_BB_Affy)),
ifelse(is.na(meta$controls_AA_Affy) & !is.na(meta$controls_AA_Illumina), (meta$controls_AB_Illumina + 
(2*meta$controls_BB_Illumina))/((2*meta$controls_AA_Illumina) + (2*meta$controls_AB_Illumina)	+ (2*meta$controls_BB_Illumina)), NA))) 
meta$weightsum<-rowSums(meta[,c("weight_Illumina", "weight_Affy")], na.rm=TRUE)
meta$seall<-sqrt(1/meta$weightsum)
meta$sumbeta<-rowSums(meta[,c("betawt_Illumina", "betawt_Affy")], na.rm=TRUE)
meta$beta<-meta$sumbeta/meta$weightsum
meta$z<-meta$beta/meta$seall
meta$pmeta<-2*pnorm(-abs(meta$z))
meta<-meta[order(meta$chromosome, meta$position),]
meta<-meta[!is.na(meta$pmeta),]
save(meta, file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/resinit_adjusted_newml_5pc_diff_dups_dropped.RData")

