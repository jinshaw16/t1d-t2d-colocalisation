#relatedness_test.R
#are there any relatives or duplicates in the analysis set? I think there might be approx 1000 in the 1958 BC genotypes using Affymetrix and Illumina.
#checking this now then will exclude if there are duplicates


library(snpStats)
library(jimisc)
library(stringr)
library(snpStatsWriter)
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

afs<-read.table(file="/well/todd/users/jinshaw/t1d_risk/affysamp",as.is=T)
affy<-affy[rownames(affy) %in% afs$V1,]

getsupaff<-function(chr){
sup<-get(load(paste0(file="/well/todd/users/jinshaw/nick_genotypes/affy/col_support_",chr,".RData")))

sup$snp.name=sup$rs.number
sup$b35<-sup$position
sup$cM<-NA
sup<-sup[,c("snp.name", "chromosome", "b35", "cM", "allele.1", "allele.2", "assay.id")]
return(sup)
}

supportaf<-lapply(c("01","02","03","04","05","06","07","08","09",10:22), getsupaff)
supportaf<-do.call("rbind", supportaf)

#liftover to gr37
chr<-paste0("chr", supportaf$chromosome)
start<-supportaf$b35
snp.name<-supportaf$snp.name

chain.file.path<-"/well/todd/users/jinshaw/liftover_allign/hg17ToHg19.over.chain"


# convert to Genomic Ranges
rownames(supportaf)<-supportaf$assay.id
input<-GRanges(
  seqname=Rle(chr),
  ranges=IRanges(start=start,end=start),
  snp.name=snp.name)


c<-import.chain(chain.file.path) ## hg17ToHg19.over.chain
alligned<-unlist(liftOver(input,c))
names(mcols(alligned))<-'snp.name'

alligned<-data.frame(snp.name=alligned@elementMetadata@listData$snp.name, position=alligned@ranges@start)
supportaf<-merge(supportaf, alligned, by="snp.name")
rownames(supportaf)<-supportaf$assay.id
supportaf<-supportaf[!is.na(supportaf$snp.name),]
affy<-affy[,colnames(affy) %in% supportaf$assay.id]
supportaf<-supportaf[colnames(affy),]
colnames(affy)<-supportaf$snp.name
rownames(supportaf)<-supportaf$snp.name


#take intersect & allign, or at least make sure same strand used:

af<-affy[,colnames(affy) %in% colnames(illu)]
il<-illu[,colnames(illu) %in% colnames(affy)]
sup<-support[colnames(il),]
supaf<-supportaf[colnames(af),]

af<-af[,colnames(il)]
supaf<-supaf[rownames(sup),]

csaf<-col.summary(af)
csil<-col.summary(il)
csaf$allele.1<-supaf$allele.1
csaf$allele.2<-supaf$allele.2
csil$allele.1<-sup$allele.1
csil$allele.2<-sup$allele.2

w<-which(csaf$allele.1==csil$allele.2 & csaf$allele.2==csil$allele.1)
af<-switch.alleles(af, snps=w)
csaf<-col.summary(af)
csaf$allele.1<-supaf$allele.1
csaf$allele.2<-supaf$allele.2
csaf$allele.1[w]<-supaf$allele.2[w]
csaf$allele.2[w]<-supaf$allele.1[w]
csaf$snp<-rownames(csaf)
csil$snp<-rownames(csil)

keep<-csaf[csaf$allele.1==csil$allele.1 & csaf$allele.2==csil$allele.2,]

both<-merge(keep, csil,by="snp")
both$diff<-abs(both$RAF.x-both$RAF.y)
both<-both[both$diff<0.01,]
library(ggplot2)
ggplot(data=both, aes(RAF.x,RAF.y)) + geom_point()

csaf<-csaf[rownames(csaf) %in% both$snp,]
csil<-csil[rownames(csil) %in% both$snp,]
af<-af[,rownames(csaf)]
il<-il[,rownames(csil)]
support<-support[colnames(il),]
supportaf<-supportaf[colnames(af),]
supportaf$allele.1<-csaf$allele.1
supportaf$allele.2<-csaf$allele.2

#######################################################
#keep intersect, combine, prune, then test relatedness#
#######################################################

write.plink(file.base="/well/todd/users/jinshaw/nick_genotypes/illu/illudat_intersect",
snps=il, chromosome=as.numeric(support$chromosome), genetic.distance=support$cM,
position=support$position37, allele.1=support$allele.1, allele.2=support$allele.2)


write.plink(file.base="/well/todd/users/jinshaw/nick_genotypes/affy/affydat_intersect",
snps=af, chromosome=as.numeric(supportaf$chromosome), genetic.distance=supportaf$cM,
position=supportaf$position, allele.1=supportaf$allele.1, allele.2=supportaf$allele.2)



system("plink --bfile /well/todd/users/jinshaw/nick_genotypes/illu/illudat_intersect --bmerge /well/todd/users/jinshaw/nick_genotypes/affy/affydat_intersect --make-bed --out /well/todd/users/jinshaw/nick_genotypes/both_intersect")
system("plink --bfile /well/todd/users/jinshaw/nick_genotypes/both_intersect --maf 0.01 --hwe 0.0005 --geno 0.05 --make-bed --out /well/todd/users/jinshaw/nick_genotypes/both_intersect_postqc")
system("plink --bfile /well/todd/users/jinshaw/nick_genotypes/both_intersect_postqc --indep-pairwise 1000 50 0.2 --out /well/todd/users/jinshaw/nick_genotypes/pruned_all")
system("plink --bfile /well/todd/users/jinshaw/nick_genotypes/both_intersect_postqc --exclude /well/todd/users/jinshaw/nick_genotypes/pruned_all.prune.out --make-bed --out /well/todd/users/jinshaw/nick_genotypes/pruned_intersect")
system(paste0("~/software/plink2 --bfile /well/todd/users/jinshaw/nick_genotypes/pruned_intersect --exclude range /well/todd/users/jinshaw/aad/under_7/mhc.txt --make-bed --out /well/todd/users/jinshaw/nick_genotypes/pruned_intersect_nomhc"))

system("king -b /well/todd/users/jinshaw/nick_genotypes/pruned_intersect_nomhc.bed --related --prefix /well/todd/users/jinshaw/nick_genotypes/related") 

#as expected, there are 1380 duplicates that need to be removed from the Affymetrix analysis.

