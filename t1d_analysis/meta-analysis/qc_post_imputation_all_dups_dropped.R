#qc_post_imputation_all.R

#artificial associations can occur under certain circumstances in imputation.
#to minimise the probability of observing false positive associations, I am copying Nick's exclusion criteria with a slight modification:
#1)Info score must be >0.3 for inclusion (Nick used 0.25)
#2)Difference in MAF of >0.05 between Illumina and Affymetrix controls
#3)Difference in MAF of >0.05 between Illumina and Affymetrix cases
#4)Difference in MAF of 0.05 between all controls and HRC european MAF
#5)Difference in effect size of >0.5 
#6) MAF<0.01

#load data :
#1 : already done this in meta-analysis_all.R

#2 : Difference in MAF of >0.05 between Illumina and Affymetrix controls
load(file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/resinit_adjusted_newml_5pc_diff_dups_dropped.RData")

meta$contdiff<-ifelse((!is.na(meta$est_Affy) & !is.na(meta$est_Illumina)),abs(meta$controls_maf_Affy-meta$controls_maf_Illumina),NA)
message(paste0("Removing ",nrow(meta[meta$contdiff>0.05 & !is.na(meta$contdiff),]), " SNPs due to difference in MAF between Affy and Illu controls of >0.05"))

drop<-meta[meta$contdiff>0.05 & !is.na(meta$contdiff),]
drop<-drop[,c("id","contdiff")]

#3 : Difference in MAF of >0.05 between Illumina and Affymetrix cases
meta$casediff<-ifelse((!is.na(meta$est_Affy) & !is.na(meta$est_Illumina)),abs(meta$cases_maf_Affy-meta$cases_maf_Illumina),NA)
message(paste0("Removing ",nrow(meta[meta$casediff>0.05 & !is.na(meta$casediff),]), " SNPs due to difference in MAF between Affy and Illu cases of >0.05"))


drop1<-meta[meta$casediff>0.05 & !is.na(meta$casediff) & meta$id!="11:2182224:A:T",]
drop1<-drop1[,c("id","casediff")]

#4)Difference in MAF of	0.05 between all controls and HRC european MAF
ref<-read.table(file="/well/todd/users/jinshaw/t1d_risk/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz", comment.char="", sep="\t",as.is=T, header=T)
ref$id<-paste0(ref$X.CHROM,":",ref$POS,":",ref$REF,":",ref$ALT)

#combined controls allele frequencies (weighted average by sample size):
csco<-merge(meta, ref, by="id",all.x=T)
csco$diff_controls<-abs(csco$AF_EXCLUDING_1000G-csco$controls_AF_B)

message(paste0("Removing ",nrow(csco[csco$diff_controls>0.05 & !is.na(csco$diff_controls),]), " SNPs due  to difference in MAF between population and HRC panel"))
drop2<-csco[csco$diff_controls>0.05 & !is.na(csco$diff_controls),]
drop2<-drop2[,c("id","diff_controls")]

d<-merge(drop, drop1, by="id",all=T)
d<-merge(d,drop2,by="id",all=T)

csco$effectdiff<-abs(csco$est_Affy-csco$est_Illumina)
drop3<-csco[csco$effectdiff>0.5 & !is.na(csco$effectdiff),]

d<-merge(d, drop3,by="id",all=T)

drop4<-csco[csco$AF<0.01 | csco$AF>0.99,]
drop4$MAF<-T
drop4<-drop4[,c("id","MAF")]
d<-merge(d,drop4,by="id",all=T)

message(paste0("In total, ",nrow(d)," SNPs removed prior to association testing due to MAF or effect size differences"))

save(d,file="/well/todd/users/jinshaw/t1d_risk/processed/vcf_all/exclusions/other_excl_newml_5pc_diff_dups_dropped.RData")
meta<-csco[!csco$id %in% d$id,]

save(meta, file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/resinit_post_imp_qc_adjusted_newml_5pc_diff_dups_dropped.RData")

load(file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/resinit_post_imp_qc_adjusted_newml_5pc_diff_dups_dropped.RData")
final<-meta[!(meta$chromosome==6 & meta$position>25000000 & meta$position<35000000),]

getmin<-function(frame){
min<-frame[frame$pmeta==min(frame$pmeta),]
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

while(n$pmeta<5*10^-8){
n<-getmin(dropped)
dropped<-dropnext(n,dropped)
both<-rbind(both, n)
}

both<-both[order(both$chromosome, both$position),]

write.table(both,file="/well/todd/users/jinshaw/output/t1d_risk/hrc/vcf_all/meta_combined/sigregions_5pc_diff_dups_dropped.txt",col.names=T, sep=";",row.names=F,quote=F)
