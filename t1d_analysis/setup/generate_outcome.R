#generate_outcome.R
#this genertaes the outcome (t1d risk) from the raw data files obtains on Ipswich:

system("/apps/well/bcftools/1.4.1/bin/bcftools query -l /well/todd/users/jinshaw/t1d_risk/T1D.Affy.chr1.vcf.gz >  /well/todd/users/jinshaw/t1d_risk/affysamp")

system("/apps/well/bcftools/1.4.1/bin/bcftools query -l /well/todd/users/jinshaw/t1d_risk/T1D.Illumina.chr1.vcf.gz >  /well/todd/users/jinshaw/t1d_risk/illusamp")


a<-read.table(file="/well/todd/users/jinshaw/t1d_risk/affysamp")
i<-read.table(file="/well/todd/users/jinshaw/t1d_risk/illusamp")
a$order<-c(1:nrow(a))
i$order<-c(1:nrow(i))

#got these files from /fs3/projects/todd/ipswich/wtccc/support_data/
one<-read.table(file="/well/todd/users/jinshaw/dil_samples/lookups/wtccc/Affx_20061121fa1_sample_58C.txt", sep="\t")
one$V1<-toupper(one$V1)
two<-read.table(file="/well/todd/users/jinshaw/dil_samples/lookups/wtccc/Affx_20061121fa1_sample_BD.txt", sep="\t")
two$V1<-toupper(two$V1)
three<-read.table(file="/well/todd/users/jinshaw/dil_samples/lookups/wtccc/Affx_20061121fa1_sample_NBS.txt",sep="\t")
three$V1<-toupper(three$V1)
four<-read.table(file="/well/todd/users/jinshaw/dil_samples/lookups/wtccc/Affx_20061121fa1_sample_T1D.txt",sep="\t")
four$V1<-toupper(four$V1)


a$phenotype<-ifelse(a$V1 %in% one$V1 | a$V1 %in% two$V1 | a$V1 %in% three$V1,0,ifelse(a$V1 %in% four$V1,1,NA))
adum<-a$V1
a<-cbind(adum,a)
a$missing<-0
a$phenotype<-as.character(a$phenotype)
vecs<-read.table(file=paste0("/well/todd/users/jinshaw/t1d_risk/processed/pcafiles/affypcs_independents_filtered.eigenvec.gz"), header=F, as.is=T)
labs<-paste0("PC",1:20)
colnames(vecs)[1:5]<-c("ID_1","ID_2","PC1","PC2","PC3")
vecs<-vecs[,c(1:5)]
a<-a[order(a$order),]
a<-a[,c("adum","V1","missing","phenotype","order")]
colnames(a)[1:2]<-c("ID_1","ID_2")
vecs<-merge(vecs,a,by=c("ID_1","ID_2"))
vecs<-vecs[order(vecs$order),]
a<-vecs[,c("ID_1","ID_2","missing","phenotype","PC1","PC2","PC3")]
write.table(a[,c(1,4)], file="/well/todd/users/jinshaw/t1d_risk/processed/vcf/affy/affy.sample_arc", row.names=F, col.names=T, quote=F)
a$PC1<-as.character(a$PC1)
a$PC2<-as.character(a$PC2)
a$PC3<-as.character(a$PC3)
dummy<-data.frame(ID_1=0,ID_2=0,missing=0,phenotype="B",PC1="C",PC2="C",PC3="C")
as<-rbind(dummy,a)
write.table(as, file="/well/todd/users/jinshaw/t1d_risk/processed/vcf/affy/affy.sample", row.names=F, col.names=T, quote=F)

#however, there are some individuals included in both affy and illu that need excluding...
dups<-read.table(file="/well/todd/users/jinshaw/nick_genotypes/related.kin0",header=T,as.is=T)
inc<-a[!a$ID_1 %in% dups$FID1,]
inc<-inc$ID_1
write.table(inc, file="/well/todd/users/jinshaw/t1d_risk/processed/vcf/affy/include_affy", quote=F, col.names=F, row.names=F)

#got these files from /fs3/projects/todd/ipswich/T1DGC/R-objects/support

load(file="/well/todd/users/jinshaw/dil_samples/lookups/t1dgc/T1DGC-sample-support.RData")
t1dgc.sample.support$V1<-t1dgc.sample.support$sample
rownames(t1dgc.sample.support)<-t1dgc.sample.support$V1
t1dgc.sample.support<-t1dgc.sample.support[rownames(t1dgc.sample.support) %in% i$V1,c("V1","t1d")]
t1dgc.sample.support<-t1dgc.sample.support[i$V1,]
i<-merge(i,t1dgc.sample.support,by="V1",all.x=T)
load(file="/well/todd/users/jinshaw/dil_samples/lookups/t1dgc/WTCCC-sample-support.RData")
wtccc.sample.support$V1<-wtccc.sample.support$sanger.id
wtccc.sample.support<-wtccc.sample.support[wtccc.sample.support$V1 %in% i$V1,]
wtccc.sample.support<-wtccc.sample.support[!duplicated(wtccc.sample.support$V1),]
rownames(wtccc.sample.support)<-wtccc.sample.support$V1
wtccc.sample.support<-wtccc.sample.support[rownames(wtccc.sample.support) %in% i$V1,c("V1","t1d")]
wtccc.sample.support<-wtccc.sample.support[i$V1,]
i<-merge(i,wtccc.sample.support,by="V1",all.x=T)
i$t1d<-ifelse(is.na(i$t1d.x),i$t1d.y,
ifelse(is.na(i$t1d.y),i$t1d.x,NA))
i$phenotype<-ifelse(i$t1d==2,1,ifelse(i$t1d==1,0,NA))
i<-i[order(i$order),]
idum<-i$V1
i<-i[,c("V1","phenotype","order")]
i$phenotype<-as.character(i$phenotype)
i<-cbind(idum,i)
i$missing<-0
vecs1<-read.table(file=paste0("/well/todd/users/jinshaw/t1d_risk/processed/pcafiles/illupcs_independents_filtered.eigenvec.gz"), header=F, as.is=T)
labs<-paste0("PC",1:20)
colnames(vecs1)[1:5]<-c("ID_1","ID_2","PC1","PC2","PC3")
colnames(i)[1:2]<-c("ID_1","ID_2")
vecs1$ID_1<-vecs1$ID_2
i$ID_1<-as.character(i$ID_1)
i$ID_2<-as.character(i$ID_2)
i<-merge(i,vecs1,by=c("ID_1","ID_2"))
i<-i[order(i$order),]
i<-i[,c("ID_1","ID_2","missing","phenotype","PC1","PC2","PC3")]
write.table(i[,c(1,4)], file="/well/todd/users/jinshaw/t1d_risk/processed/vcf/illu/illu.sample_arc", row.names=F, col.names=T, quote=F)

dummy<-data.frame(ID_1=0,ID_2=0,missing=0,phenotype="B",PC1="C",PC2="C",PC3="C")
i$PC1<-as.character(i$PC1)
i$PC2<-as.character(i$PC2)
i$PC3<-as.character(i$PC3)
i<-rbind(dummy,i)
write.table(i, file="/well/todd/users/jinshaw/t1d_risk/processed/vcf/illu/illu.sample", row.names=F, col.names=T, quote=F)
