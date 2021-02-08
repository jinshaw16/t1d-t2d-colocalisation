#coloc_t2d_all_dups_drop_fdr_all.R
#looking at colocalisation of disease association and T2D
library(stringr)
library(snpStats)
library(coloc)
library(ggplot2)
library(ggbio)
library(annotSnpStats)
library(GenomicRanges)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(plyr)
library(dplyr)
library(cowplot)
library(scales)
#load the disease associations with the <7s:

mydir <-"/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/"
outdir<-"/well/todd/users/jinshaw/output/t1d_risk/hrc/vcf/meta_sardinia/newml/coloc/fdr/up/"

#load T1D summary stats:
t1d<-read.table(file=paste0(mydir,"summary_stats_af_sard_newml_5pc_diff_dups_dropped.txt.gz"),header=T, as.is=T)
#load T2D summary stats:
t2d<-read.table(file="/well/todd/users/jinshaw/t2d_risk/Mahajan.NatGenet2018b.T2D.European.txt.gz",header=T, as.is=T)

#allign:
t2d$id<-paste0(t2d$Chr,":",t2d$Pos,":",t2d$NEA,":",t2d$EA)
t2d$alt<-paste0(t2d$Chr,":",t2d$Pos,":",t2d$EA,":",t2d$NEA)
w<-which(t2d$alt %in% t1d$id)
t2d[w,"dum1"]<-t2d[w,"NEA"]
t2d[w,"dum2"]<-t2d[w,"EA"]
t2d[w,"NEA"]<-t2d[w,"dum2"]
t2d[w,"EA"]<-t2d[w,"dum1"]
t2d[w,"Beta"]<-t2d[w,"Beta"]*-1
t2d[w,"EAF"]<-1-t2d[w,"EAF"]
t2d$id<-paste0(t2d$Chr,":",t2d$Pos,":",t2d$NEA,":",t2d$EA)
table(t2d$id %in% t1d$id)

#find each region that overlaps (reproduce Anubhas graph):
comb<-merge(t1d,t2d,by="id")
comb<-comb[order(comb$chromosome,comb$position),]
library(plyr)
library(dplyr)
l<-ddply(comb,.(chromosome),summarise, value=max(position))
l$sum<-cumsum(l$value)
l[,"next"]<-c(0,l$sum[1:(nrow(l)-1)])
l<-l[,c("next","chromosome")]

comb1<-merge(comb,l,by="chromosome")
comb1<-comb1[order(comb1$chromosome, comb1$position),]
comb1$x<-comb1[,"next"]+comb1$position


comb1$t1d_logp<-log10(comb1$pfinal)*-1
comb1$t2d_logp<-log10(comb1$Pvalue)


#examine the regions that look like they have some overlap (FDR<0.01 in both):
comb1<-comb1[!(comb1$chromosome==6 & comb1$position>25000000 & comb1$position<35000000),]
comb1$pfdr_t1d<-p.adjust(comb1$pfinal,method="BH")
comb1$pfdr_t2d<-p.adjust(comb1$Pvalue,method="BH")
save(comb1, file="/well/todd/users/jinshaw/t2d_risk/both_sumstats.RData")

#get all T1D FDR 1% regions:
getmin<-function(frame,p){
min<-frame[frame[,p]==min(frame[,p]),]
return(min)
}

dropnext<-function(min, results){
droplow<-min$position-500000
drophigh<-min$position+500000
chrom<-min$chromosome
m<-results %>%
filter(!(chromosome==chrom & position>droplow & position<drophigh))
m<-as.data.frame(m) 
#[!(results$chromosome==chrom & results$position>droplow & results$position<drophigh),]
return(m)
}

n<-getmin(comb1, "pfinal")
dropped<-dropnext(n, comb1)
both<-n

while(n$pfdr_t1d<0.01){
n<-getmin(dropped,"pfinal")
dropped<-dropnext(n,dropped)
both<-rbind(both, n)
}

lesst1d<-both[order(both$chromosome, both$position),]

comb1$zt2d<-abs(comb1$Beta/comb1$SE)
getmax<-function(frame,z){
min<-frame[frame[,z]==max(frame[,z]),]
return(min)
}

n<-getmax(comb1,"zt2d")
n<-n[1,]
dropped<-dropnext(n, comb1)
both<-n
while(n$pfdr_t2d<0.01){
n<-getmax(dropped,"zt2d")
n<-n[1,]
dropped<-dropnext(n,dropped)
both<-rbind(both, n)
}

lesst2d<-both[order(both$chromosome, both$position),]
lesst2d$min<-lesst2d$position-250000
lesst2d$max<-lesst2d$position+250000
lesst1d$min<-lesst1d$position-250000
lesst1d$max<-lesst1d$position+250000

#save(lesst1d, lesst2d, file="/well/todd/users/jinshaw/t2d_risk/regionsoverlapping_dups_drop_fdr_all.RData")
#remove those in LD with each other, likely consistuting the same genetic signal (done using LD link):


load(file="/well/todd/users/jinshaw/t2d_risk/regionsoverlapping_dups_drop_fdr_all.RData")




#get r2 between the variants (break into chromosome + disease):

dolds<-function(disease, chr, disframe){

j<-disframe[disframe$chromosome==chr,]
if(nrow(j)>1){
j<-j[,c("chromosome","position")]
write.table(j,file=paste0("/well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_",chr,".txt"),col.names=F,row.names=F,quote=F, sep="\t")

sink(file=paste0("/users/todd/jinshaw/programs/t1d_risk/vcf/coloc/ld/",disease,"_",chr))
cat(paste0("#!/bin/bash

#$ -cwd -V
#$ -N ",disease,"_",chr," -j y
#$ -P todd.prjc -q short.qc
bcftools view /well/todd/users/jinshaw/t1d_risk/T1D.Illumina.chr",chr,".vcf.gz ",
"--regions-file /well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_",chr,".txt -o /well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_",chr,
".vcf -O vcf
bcftools view /well/todd/users/jinshaw/t1d_risk/T1D.Affy.chr",chr,".vcf.gz ",
"--regions-file /well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_",chr,".txt -o /well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_",chr,
"_affy.vcf -O vcf
plink -vcf /well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_",chr,
".vcf --make-bed --const-fid --out /well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_chr",chr," 
plink -vcf /well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_",chr,
".vcf --make-bed --const-fid --out /well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_chr",chr,"_affy

#now in an appropriate format, using ldstore to extract LD info and SNP info:
#LD info:
plink --bfile /well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_chr",chr,
" --r2 --ld-window-r2 0.1 --out /well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_",chr,"
plink --bfile /well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_chr",chr,
"_affy --r2 --ld-window-r2 0.1 --out /well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_",chr,"_affy\n"))
sink()
system(paste0("qsub /users/todd/jinshaw/programs/t1d_risk/vcf/coloc/ld/",disease,"_",chr))
}
}
lapply(c(1:22),dolds, disframe=lesst1d, disease="t1d")
lapply(c(1:22),dolds, disframe=lesst2d, disease="t2d")


#read the LD information in:
getld<-function(chrom,disease, disframe){
if(!file.exists(paste0("/well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_chr",chrom,".bim"))){
ldall<-NULL
}
if(file.exists(paste0("/well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_chr",chrom,".bim"))){
snps<-read.table(file=paste0("/well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_chr",chrom,".bim"),header=F,as.is=T)
snps1<-read.table(file=paste0("/well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_chr",chrom,"_affy.bim"),header=F,as.is=T)

d<-disframe[disframe$chromosome==chrom,]
lds<-read.table(file=paste0("/well/todd/users/jinshaw/t1d_risk/processed/lds/",disease,"_",chrom,".ld"),header=T,as.is=T)
if(nrow(lds)==0){
ldall<-NULL
}
if(nrow(lds)>0){
gethem<-function(i){
k<-lds[i,]
d1<-d[d$id %in% c(k$SNP_A,k$SNP_B),]
ldal<-d1[d1$pfinal==max(d1$pfinal),]
return(ldal)
}
ldall<-lapply(c(1:nrow(lds)),gethem)
ldall<-do.call("rbind",ldall)
}
message(paste0("Done ",chrom))
}
return(ldall)
}
t1dld<-lapply(c(1:22),getld,disease="t1d",disframe=lesst1d)
t1dld<-do.call("rbind",t1dld)

t2dld<-lapply(c(1:22),getld,disease="t2d",disframe=lesst2d)
t2dld<-do.call("rbind",t2dld)


#so can remove these from the list of associated regions:
t1dhits<-lesst1d[!lesst1d$id %in% t1dld$id,]
t2dhits<-lesst2d[!lesst2d$id %in% t2dld$id,]
#save(t1dhits, t2dhits, file="/well/todd/users/jinshaw/t2d_risk/regionsoverlapping_dups_drop_fdr_all_nolds.RData")

load(file="/well/todd/users/jinshaw/t2d_risk/regionsoverlapping_dups_drop_fdr_all_nolds.RData")
#use GRanges to identify any of these 0.5Mb regions that overlap across diseases:
library(GenomicRanges)
t2dh<-GRanges(seqnames=paste0("chr",t2dhits$chromosome),
IRanges(t2dhits$min, end=t2dhits$max))
t1dh<-GRanges(seqnames=paste0("chr",t1dhits$chromosome),
IRanges(t1dhits$min,end=t1dhits$max))

both<-mergeByOverlaps(t2dh,t1dh)
b<-as.data.frame(both)
b$dup1<-duplicated(b$t2dh.start)
b$dup2<-duplicated(b$t1dh.start)
b$anydups<-ifelse(b$dup1==TRUE | b$dup2==TRUE,0,1)
b$group<-cumsum(b$anydups)

mergethem<-function(group){
gs<-b[b$group==group,]
gs$mins<-min(gs$t2dh.start,gs$t1dh.start)
gs$maxs<-max(gs$t2dh.end,gs$t1dh.end)
gs$start<-gs$mins
gs$end<-gs$maxs
gs<-gs[1,]
return(gs)
}
bs<-lapply(c(unique(b$group)),mergethem)
bs<-do.call("rbind",bs)



colocranges<-GRanges(bs$t2dh.seqnames,
IRanges(start=bs$start, end=bs$end),
line=bs$group)

comb1$chr<-paste0("chr",comb1$chromosome)

#get univariable summary stats to use in the conditional analyses later:
getstats<-function(line){
colo<-as.data.frame(colocranges)
co<-colo[colo$line==line,]
cs<-comb1[comb1$chr==co$seqnames & comb1$position>co$start & comb1$position<co$end,]
return(cs)
}
resboth<-lapply(c(1:nrow(as.data.frame(colocranges))),getstats)
names(resboth)<-c(1:nrow(as.data.frame(colocranges)))
save(colocranges,resboth,file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/conditional/univariable_all.RData")

load(file="/well/todd/users/jinshaw/t2d_risk/both_sumstats.RData")

ranges<-as.data.frame(colocranges)

testcoloc<-function(line,p,se,beta, plot=TRUE,frame){
message(paste0("Testing colocalisation of diseases in line ",line,"\n"))
r<-ranges[line,]
chrom<-as.numeric(gsub("chr","",r$seqnames))
mini<-r$start
maxi<-r$end
bd<-frame[frame$chromosome==chrom & frame$position>mini & frame$position<maxi,]


bd$diff<-abs(bd$EAF-bd$controls_AF_B)
bd<-bd[bd$diff<0.05,]

t1<-bd
t1$pvalues<-t1[,p]
t1$type="cc"
t1$beta=t1[,beta]
t1$varbeta=t1[,se]^2
t1$MAF<-ifelse(t1$controls_AF_B>0.5,
1-t1$controls_AF_B, t1$controls_AF_B)
t1$s<-0.39
t1$N<-5268+7977

t2<-bd
t2$pvalues<-t2$Pvalue
t2$type="cc"
t2$beta=t2$Beta
t2$varbeta=t2$SE^2
t2$MAF<-ifelse(t2$EAF>0.5,
1-t2$EAF, t2$EAF)
t2$s<-0.09
t2$N<-898114
df1=t1[,c("id","beta","varbeta","type","MAF","s","N")]
df2=t2[,c("id","beta","varbeta","type","MAF","s","N")]
rownames(df1)<-df1$id
rownames(df2)<-df2$id
df1<-df1[rownames(df2),]

coloc<-coloc.abf(df1, df2, p1 = 1e-04, p2 = 1e-04,
      p12 = 5e-06)
probs<-coloc[[1]]
probs<-data.frame(line=line, P0=probs[2],P1=probs[3],P2=probs[4],P3=probs[5],P4=probs[6])

bd$logp_t2d<-log10(bd$Pvalue)*-1
bd$logp_t1d<-log10(bd[,p])*-1
top<-bd[bd$logp_t2d==max(bd$logp_t2d),]
top<-top[top$logp_t1d==max(top$logp_t1d),]
probs$id<-top$id
probs$ID<-top$ID
probs$position<-top$position
probs$t1d_beta<-top[,beta]
probs$t1d_se<-top[,se]
probs$t1d_p<-top[,p]
probs$t2d_beta<-top$Beta
probs$t2d_se<-top$SE
probs$t2d_p<-top$Pvalue
probs$alleleA<-top$alleleA
probs$alleleB<-top$alleleB
probs$minpos<-mini
probs$maxpos<-maxi
probs$chromosome<-chrom
if(plot==TRUE){
system(paste0("tabix /well/todd/users/jinshaw/t1d_risk/T1D.Affy.chr",chrom,".vcf.gz ",chrom,":",ceiling(mini),"-",ceiling(maxi), 
" -h > /well/todd/users/jinshaw/t1d_risk/regs/",line,"_region_affy.vcf"))
system(paste0("plink --vcf /well/todd/users/jinshaw/t1d_risk/regs/",line,"_region_affy.vcf --double-id --make-bed --out /well/todd/users/jinshaw/t1d_risk/regs/",line))


init<-read.plink(fam=paste0("/well/todd/users/jinshaw/t1d_risk/regs/",line,".fam"),
bim=paste0("/well/todd/users/jinshaw/t1d_risk/regs/",line,".bim"),
bed=paste0("/well/todd/users/jinshaw/t1d_risk/regs/",line,".bed"))
init<-annot.plink(init)
init<-init[,colnames(init) %in% bd$id]
init<-init[,bd$id]

tsnp<-top$id
lds<-ld(init[,!colnames(init) %in% tsnp],init[,tsnp],stats="R.squared")
next1<-data.frame(tsnp=1)
colnames(next1)=tsnp
lds<-rbind(lds,next1)
rownames(lds)[nrow(lds)]<-tsnp
lds$id=rownames(lds)
dfb1<-merge(bd,lds,by="id")
colnames(dfb1)[ncol(dfb1)]<-"ld"
o1<-ggplot(data=dfb1, aes(position, logp_t2d, colour=ld)) + geom_point() +
scale_y_continuous(name=bquote('T2D '*-log[10]*~italic(p)~'value'))  +
scale_color_gradient2(name=paste0("LD with T2D \nindex variant"),low = "black",
midpoint = 0.5,
mid = "orange",
high = "red",
space="Lab") +
scale_x_continuous(labels=comma) #+
#annotate("text", x=min+75000, y=max(b$logp_t2d), label=paste0("Colocalisation \nposterior probability=",round(probs[6],digits=3))) #+
#annotate("text", x=max-75000, y=max(b$logp_t2d), label=paste0(gene))
o2<-ggplot(data=dfb1, aes(position, logp_t1d, colour=ld)) + geom_point() +
scale_y_continuous(name=bquote('T1D '*-log[10]*~italic(p)~'value')) +
scale_color_gradient2(name=paste0("LD with T2D \nindex variant"),low = "black",
midpoint = 0.5,
mid = "orange",
high = "red",
space="Lab") +
scale_x_continuous(labels=comma)


data(genesymbol, package = "biovizBase")
g <- genesymbol[seqnames(genesymbol) == paste0('chr',chrom)]
gr <- GRanges(
      seqnames = Rle(c(paste0("chr",chrom)), c(1)),
      ranges = IRanges(mini, end = maxi))

o<-as.matrix(findOverlaps(g,gr))
g<-g[o[,1],]

posi<-top$position
gs<-data.frame(g)
gs$diff<-ceiling((gs$end-gs$start)/2)
gs$mid<-gs$start+gs$diff
gs$diff1<-abs(gs$mid-posi)
ks<-gs[gs$diff1==min(gs$diff1),]
ks<-ks[1,]
gene<-ks$symbol

t1<-autoplot(Homo.sapiens, which = g) +
scale_x_continuous(labels=comma)

t<-tracks(t1,o1,o2, xlab=paste0("Position along chromosome ",chrom, " (bp)"))
png(file=paste0(outdir,gene,"_",line,"_tracks_5pc_diff_dups_drop.png"), height=20, width=20, units="cm", res=400)
print({t})
dev.off()
#ggsave(t,file=paste0(outdir,gene,"_",line,"_tracks_5pc_diff_dups_drop.png"), height=20, width=20, units="cm", dpi=400)
}
if(plot==FALSE){
data(genesymbol, package = "biovizBase")
g <- genesymbol[seqnames(genesymbol) == paste0('chr',chrom)]
gr <- GRanges(
      seqnames = Rle(c(paste0("chr",chrom)), c(1)),
      ranges = IRanges(mini, end = maxi))

o<-as.matrix(findOverlaps(g,gr))
g<-g[o[,1],]

posi<-top$position
gs<-data.frame(g)
gs$diff<-ceiling((gs$end-gs$start)/2)
gs$mid<-gs$start+gs$diff
gs$diff1<-abs(gs$mid-posi)
ks<-gs[gs$diff1==min(gs$diff1),]
ks<-ks[1,]
gene<-ks$symbol
}
probs$gene=gene
return(probs)
}

res1<-mclapply(c(1:nrow(ranges)),testcoloc, p="pfinal",se="sefinal",beta="betafinal",frame=comb1,mc.cores=4)
res1<-do.call("rbind",res1)
res1$t1dsig=1
res1$t2dsig=1
res1<-res1[,c("line","id","ID","chromosome","position","gene","alleleA","alleleB","t1d_beta","t1d_se","t1d_p",
"t2d_beta","t2d_se","t2d_p","t1dsig","t2dsig","P0","P1","P2","P3","P4","minpos","maxpos")]

write.table(res1, file=paste0(outdir,"all_probs.txt"),col.names=T,row.names=F, sep="\t",quote=F)

load(file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/resinit_adjusted_newml_5pc_diff_dups_dropped.RData")
meta$diff1<-abs(meta$controls_info_Illumina-meta$cases_info_Illumina)
meta$diff2<-abs(meta$controls_info_Affy-meta$cases_info_Affy)
meta<-meta[meta$diff1<0.01 | meta$diff2<0.01 & !is.na(meta$diff1) & !is.na(meta$diff2),]
meta$rsid<-meta$id
t1d<-meta
t1d<-t1d[!is.na(t1d$chromosome),]

#t2d<-read.table(file="/well/todd/users/jinshaw/t2d_risk/Mahajan.NatGenet2018b.T2D.European.txt.gz",header=T, as.is=T)
#allign:
#t2d$id<-paste0(t2d$Chr,":",t2d$Pos,":",t2d$NEA,":",t2d$EA)
#t2d$alt<-paste0(t2d$Chr,":",t2d$Pos,":",t2d$EA,":",t2d$NEA)
#w<-which(t2d$alt %in% t1d$id)
#t2d[w,"dum1"]<-t2d[w,"NEA"]
#t2d[w,"dum2"]<-t2d[w,"EA"]
#t2d[w,"NEA"]<-t2d[w,"dum2"]
#t2d[w,"EA"]<-t2d[w,"dum1"]
#t2d[w,"Beta"]<-t2d[w,"Beta"]*-1
#t2d[w,"EAF"]<-1-t2d[w,"EAF"]
#t2d$id<-paste0(t2d$Chr,":",t2d$Pos,":",t2d$NEA,":",t2d$EA)
#table(t2d$id %in% t1d$id)

#find each region that overlaps (reproduce Anubhas graph):
#comb<-merge(meta,t2d,by="id")
#comb<-comb[order(comb$chromosome,comb$position),]


#save(comb,file="/well/todd/users/jinshaw/t2d_risk/both_sumstats_uk.RData")
load(file="/well/todd/users/jinshaw/t2d_risk/both_sumstats_uk.RData")
comb$controls_AF_B<-comb$controls_maf_Illumina
res2<-mclapply(c(1:nrow(ranges)),testcoloc, p="pmeta",se="seall",beta="beta",frame=comb, plot=FALSE,mc.cores=4)
res2<-do.call("rbind",res2)
