#creating_genotype_data_anuhba.R
library(snpStats)
library(jimisc)
library(stringr)
#Illumina first:
illusamps<-read.table(file="/well/todd/users/jinshaw/nick_cooper_map/impute-9-99000001-100000001.out_samples.gz", header=T, as.is=T)
illusamps<-illusamps[-1,]

readinillu<-function(chr){
case<-get(load(file=paste0("/well/todd/users/jinshaw/nick_genotypes/illu/case.snp.chr",chr,".RData")))
cont<-get(load(file=paste0("/well/todd/users/jinshaw/nick_genotypes/illu/control.snp.chr",chr,".RData")))
sang<-get(load(file=paste0("/well/todd/users/jinshaw/nick_genotypes/illu/sanger.snp.chr",chr,".RData")))

all<-rbind(case, cont)
all<-rbind(all, sang)

all<-as(all,"SnpMatrix")

all<-all[rownames(all) %in% illusamps$ID_1,]
message(paste0("DONE CHROMOSOME ",chr))
return(all)
}
illu<-lapply(c(1:22),readinillu)
illu<-do.call(cbind, illu)

getsupillu<-function(chr){
sup<-get(load(file=paste0("/well/todd/users/jinshaw/nick_genotypes/illu/snp-support-",chr,".RData")))
sup<-sup[rownames(sup) %in% colnames(illu),]
sup$alleles<-str_count(sup$ilmn_alleles,"/")
sup<-sup[sup$alleles<2,]
sup$allele.1<-ifelse(sup$ilmn_top_or_bot=="T",gsub("/.*", "",sup$ilmn_alleles),
                     ifelse(sup$ilmn_top_or_bot=="B", gsub(".*/","",sup$ilmn_alleles),
                            ifelse(sup$ilmn_top_or_bot=="" & sup$b36_strand=="+",gsub(".*/","",sup$b36_alleles),
                                   ifelse(sup$ilmn_top_or_bot=="" & sup$b36_strand=="-", gsub("/.*","",sup$b36_alleles),NA))))
sup$allele.2<-ifelse(sup$ilmn_top_or_bot=="T",gsub(".*/", "",sup$ilmn_alleles),
                     ifelse(sup$ilmn_top_or_bot=="B", gsub("/.*","",sup$ilmn_alleles),
                            ifelse(sup$ilmn_top_or_bot=="" & sup$b36_strand=="+",gsub("/.*","",sup$b36_alleles),
                                   ifelse(sup$ilmn_top_or_bot=="" & sup$b36_strand=="-", gsub(".*/","",sup$b36_alleles),NA))))
sup$cM<-NA
sup$snp.name<-rownames(sup)
sup$chromosome<-sup$b36_chr
#sup<-sup[,c("snp.name", "cM", "chromosome", "b36_start", "allele.1", "allele.2")]
message(paste0("DONE CHROMOSOME ", chr))
return(sup)
}
support<-lapply(c(1:22), getsupillu)
support<-do.call("rbind", support)

illu<-illu[,colnames(illu) %in% support$snp.name]

#update coordinates to build 37:
sup<-support
support<-support[,c("snp.name", "cM", "chromosome", "b36_start", "allele.1", "allele.2")]
support<-liftthem(support, chain="hg18ToHg19.over.chain",
         position="b36_start", updateto="37")
rownames(support)<-support$snp.name
support<-support[colnames(illu),]

write.plink(file.base="/well/todd/users/jinshaw/nick_genotypes/illu/illudat_update_alleles",
snps=illu, chromosome=as.numeric(support$chromosome), genetic.distance=support$cM,
position=support$position37, allele.1=support$allele.1, allele.2=support$allele.2)

#export a snp list for anubha of snps whcih may be wrong (those with missing ilmn_strand data).
dodgy<-sup[sup$ilmn_top_or_bot=="",]
dodgy$thinkok<-ifelse(dodgy$b36_strand=="+",1,0)
dodgy<-dodgy[,c("snp.name","thinkok")]
write.table(dodgy, file="/well/todd/users/jinshaw/nick_genotypes/illu/illudat_unsure_alleles.tab",
            quote=F, sep="\t", col.names=T, row.names=F)

#Affymetrix next:
affysamps<-read.table(file='/well/todd/users/jinshaw/nick_cooper_map/affy.sample.txt.gz', header=T, as.is=T)
affysamps<-affysamps[-1,]

readinaffy<-function(chr){
af<-get(load(file=paste0("/well/todd/users/jinshaw/nick_genotypes/affy/converted_genotypes_ALL_",chr,"_90.RData"))) 
af<-af[rownames(af) %in% affysamps$ID_1,]
message(paste0("DONE CHROMOSOME ", chr))
return(af)
}
affy<-lapply(c("01","02","03","04","05","06","07","08","09",10:22), readinaffy)
affy<-do.call("cbind",affy)

getsupaff<-function(chr){
sup<-get(load(paste0(file="/well/todd/users/jinshaw/nick_genotypes/affy/col_support_",chr,".RData")))

sup$snp.name=sup$rs.number
sup$b35<-sup$position
sup$cM<-NA
sup<-sup[,c("snp.name", "chromosome", "b35", "cM", "allele.1", "allele.2", "assay.id")]
return(sup)
}

support<-lapply(c("01","02","03","04","05","06","07","08","09",10:22), getsupaff)
support<-do.call("rbind", support)

#liftover to gr37 
chr<-paste0("chr", support$chromosome)
start<-support$b35
snp.name<-support$snp.name

chain.file.path<-"/well/todd/users/jinshaw/liftover_allign/hg17ToHg19.over.chain"

# convert to Genomic Ranges
rownames(support)<-support$assay.id
input<-GRanges(
  seqname=Rle(chr),
  ranges=IRanges(start=start,end=start),
  snp.name=snp.name)


c<-import.chain(chain.file.path) ## hg17ToHg19.over.chain
alligned<-unlist(liftOver(input,c))
names(mcols(alligned))<-'snp.name'

alligned<-data.frame(snp.name=alligned@elementMetadata@listData$snp.name, position=alligned@ranges@start)
support<-merge(support, alligned, by="snp.name")
rownames(support)<-support$assay.id
support<-support[!is.na(support$snp.name),]
affy<-affy[,colnames(affy) %in% support$assay.id]
support<-support[colnames(affy),]
colnames(affy)<-support$snp.name
rownames(support)<-support$snp.name

write.plink(file.base="/well/todd/users/jinshaw/nick_genotypes/affy/affydat",
snps=affy, chromosome=as.numeric(support$chromosome), genetic.distance=support$cM,
position=support$position, allele.1=support$allele.1, allele.2=support$allele.2)




