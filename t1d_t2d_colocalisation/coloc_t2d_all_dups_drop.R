#coloc_t2d_all_dups_drop.R
#looking at colocalisation of disease association and T2D
library(stringr)
library(snpStats)
library(coloc)
library(ggplot2)
library(ggbio)
library(annotSnpStats)
library(Sushi)
library(GenomicRanges)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#load the disease associations with the <7s:

mydir <-"/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/"
outdir<-"/well/todd/users/jinshaw/output/t1d_risk/hrc/vcf/meta_sardinia/newml/coloc/"

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
comb1$col<-ifelse(comb1$chromosome %in% c(1,3,5,7,9,11,13,15,17,19,21),1,0)

jpeg(file=paste0(outdir,"t1d_t2d_man_dups_dropped.jpeg"), height=16, width=16, units="cm", res=300)
ggplot(data=comb1, aes(t2d_logp, x,colour=as.factor(col))) + geom_point() +
geom_point(data=comb1, aes(t1d_logp,x,colour=as.factor(col))) + geom_vline(xintercept=8,colour="red") +
geom_vline(xintercept=-8,colour="red")	+ 
geom_vline(xintercept=0,colour="black") +
coord_cartesian(xlim=c(-50,50)) + 
scale_colour_manual(values=c("black","blue")) +
theme(axis.text.y=element_blank(),
axis.title.y=element_blank(),
axis.ticks.y=element_blank(),
axis.text.x=element_blank(),
axis.title.x=element_blank(),
axis.ticks.x=element_blank(),
legend.title = element_blank()) +
guides(colour=FALSE)
dev.off()

#examine the regions that look like they have some overlap (5*10-5 in one and 5*10^-8 in the other):

lesst1d<-comb1[comb1$pfinal<5*10^-5,]
lesst1d<-lesst1d[lesst1d$Pvalue<5*10^-8,]
lesst2d<-comb1[comb1$Pvalue<5*10^-5,]
lesst2d<-lesst2d[lesst2d$pfinal<5*10^-8,]
lesst2d<-lesst2d[!lesst2d$id %in% lesst1d$id,]
lessall<-rbind(lesst1d,lesst2d)
#remove MHC:
lessall<-lessall[!(lessall$chromosome==6 & lessall$position>25000000 & lessall$position<35000000),]
lessall<-lessall[order(lessall$chromosome,lessall$position),]
lessall$dum<-c(lessall$x[2:nrow(lessall)],lessall$x[nrow(lessall)])
lessall$diff<-lessall$dum-lessall$x
lessall$new<-ifelse(lessall$diff>150000,1,0)
lessall$new1<-c(0,lessall[1:(nrow(lessall)-1),"new"])
lessall$group<-cumsum(lessall$new1)

#what is the range in each region?
l<-ddply(lessall,.(group),summarise, ma=max(position), mi=min(position), chrom=min(chromosome))
l$range<-l$ma-l$mi
l1<- lessall %>%
group_by(group) %>%
filter(t1d_logp==max(t1d_logp))
l$idmin<-l1$id
l$IDmin<-l1$ID
l$t1dminp<-l1$pfinal
l$beta<-l1$betafinal
l$se<-l1$sefinal
l$t2dminp<-l1$Pvalue
#all less than 0.5Mb, so chosing this as the standard width (from the median position):

lessall$gene<-ifelse(lessall$group==0,"MACF1",
ifelse(lessall$group==1,"MAEA",
ifelse(lessall$group==2,"CENPW",
ifelse(lessall$group==3,"JAZF1",
ifelse(lessall$group==4,"GLIS3",
ifelse(lessall$group==5,"INS",
ifelse(lessall$group==6,"BCAR1",
ifelse(lessall$group==7,"BPTF",
ifelse(lessall$group==8,"MTMR3-ASCC2",NA)))))))))
l$gene<-ifelse(l$group==0,"MACF1",
ifelse(l$group==1,"MAEA",
ifelse(l$group==2,"CENPW",
ifelse(l$group==3,"JAZF1",
ifelse(l$group==4,"GLIS3",
ifelse(l$group==5,"INS",
ifelse(l$group==6,"BCAR1",
ifelse(l$group==7,"BPTF",
ifelse(l$group==8,"MTMR3-ASCC2",NA)))))))))
save(l,lessall, file="/well/todd/users/jinshaw/t2d_risk/regionsoverlapping_dups_drop.RData")

testcoloc<-function(t1set,gene, p, beta, se, plot=TRUE){
message(paste0("Testing colocalisation of disease and ",gene,"\n"))
chrom<-lessall[lessall$gene==gene,"chromosome"][1]
po<-l[l$gene==gene,]
mid<-po$mi+((po$ma-po$mi)/2)
min<-mid-250000
max<-mid+250000
t1<-t1set[t1set$chromosome==chrom & t1set$position>min & t1set$position<max,]
t2<-t2d[t2d$Chr==chrom & t2d$Pos>min & t2d$Pos<max,]

b<-merge(t1,t2,by="id")
b$diff<-abs(b$EAF-b$controls_AF_B)
b<-b[b$diff<0.05,]

t1<-t1[t1$id %in% b$id,]
t2<-t2[t2$id %in% b$id,]

t1$pvalues<-t1[,p]
t1$type="cc"
t1$beta=t1[,beta]
t1$varbeta=t1[,se]^2
t1$MAF<-ifelse(t1$controls_AF_B>0.5,
1-t1$controls_AF_B, t1$controls_AF_B)
t1$s<-0.39

t2$pvalues<-t2$Pvalue
t2$type="cc"
t2$beta=t2$Beta
t2$varbeta=t2$SE^2
t2$MAF<-ifelse(t2$EAF>0.5,
1-t2$EAF, t2$EAF)
t2$s<-0.09
df1=t1[,c("id","beta","varbeta","type","MAF","s")]
df2=t2[,c("id","beta","varbeta","type","MAF","s")]
rownames(df1)<-df1$id
rownames(df2)<-df2$id
df1<-df1[rownames(df2),]

coloc<-coloc.abf(df1, df2, p1 = 1e-04, p2 = 1e-04,
      p12 = 5e-06)
probs<-coloc[[1]]

b$logp_t2d<-log10(b$Pvalue)*-1
b$logp_t1d<-log10(b[,p])*-1
top<-b[b$logp_t2d==max(b$logp_t2d),]
chrom=top$chromosome
#if(gene!="BPTF"){
system(paste0("tabix /well/todd/users/jinshaw/t1d_risk/T1D.Affy.chr",chrom,".vcf.gz ",chrom,":",ceiling(min),"-",ceiling(max), " -h > /well/todd/users/jinshaw/t1d_risk/",gene,"_region_affy.vcf"))
#}
#if(gene=="BPTF"){
#system(paste0("tabix /well/todd/users/jinshaw/t1d_risk/T1D.Illumina.chr",chrom,".vcf.gz ",chrom,":",ceiling(min),"-",ceiling(max), " -h > /well/todd/users/jinshaw/t1d_risk/",gene,"_region_affy.vcf"))
#}
system(paste0("plink --vcf /well/todd/users/jinshaw/t1d_risk/",gene,"_region_affy.vcf --double-id --make-bed --out /well/todd/users/jinshaw/t1d_risk/",gene))

init<-read.plink(fam=paste0("/well/todd/users/jinshaw/t1d_risk/",gene,".fam"),
bim=paste0("/well/todd/users/jinshaw/t1d_risk/",gene,".bim"),
bed=paste0("/well/todd/users/jinshaw/t1d_risk/",gene,".bed"))
init<-annot.plink(init)
init<-init[,colnames(init) %in% b$id]
init<-init[,b$id]

tsnp<-top$id
lds<-ld(init[,!colnames(init) %in% tsnp],init[,tsnp],stats="R.squared")
next1<-data.frame(tsnp=1)
colnames(next1)=tsnp
lds<-rbind(lds,next1)
rownames(lds)[nrow(lds)]<-tsnp
lds$id=rownames(lds)
dfb1<-merge(b,lds,by="id")
colnames(dfb1)[ncol(dfb1)]<-"ld"
o1<-ggplot(data=dfb1, aes(position, logp_t2d, colour=ld)) + geom_point() +
scale_y_continuous(name=bquote('T2D '*-log[10] *' p-value'))  +
scale_color_gradient2(name=paste0("LD with T2D \nindex variant"),low = "black",
midpoint = 0.5,
mid = "orange",
high = "red",
space="Lab") + 
annotate("text", x=min+75000, y=max(b$logp_t2d), label=paste0("Colocalisation \nposterior probability=",round(probs[6],digits=3))) #+
#annotate("text", x=max-75000, y=max(b$logp_t2d), label=paste0(gene))
o2<-ggplot(data=dfb1, aes(position, logp_t1d, colour=ld)) + geom_point() +
scale_y_continuous(name=bquote('T1D '*-log[10] *' p-value')) +
scale_color_gradient2(name=paste0("LD with T2D \nindex variant"),low = "black",
midpoint = 0.5,
mid = "orange",
high = "red",
space="Lab")

data(genesymbol, package = "biovizBase")
g <- genesymbol[seqnames(genesymbol) == paste0('chr',chrom)]
gr <- GRanges(
      seqnames = Rle(c(paste0("chr",chrom)), c(1)),
      ranges = IRanges(min, end = max))

o<-as.matrix(findOverlaps(g,gr))
g<-g[o[,1],]


t1<-autoplot(Homo.sapiens, which = g)

t<-tracks(t1,o1,o2, xlab=paste0("Position along chromosome ",chrom))

if (plot==TRUE){
ggsave(t,file=paste0(outdir,gene,"_tracks_5pc_diff_dups_drop.png"), height=20, width=20, units="cm", dpi=400)

}
return(probs)
}
testcoloc("MACF1",p="pfinal", beta="betafinal",se="sefinal",t1set=t1d)
testcoloc("MAEA",p="pfinal", beta="betafinal",se="sefinal",t1set=t1d)
testcoloc("CENPW",p="pfinal", beta="betafinal",se="sefinal",t1set=t1d)
testcoloc("JAZF1",p="pfinal", beta="betafinal",se="sefinal",t1set=t1d)
testcoloc("GLIS3",p="pfinal", beta="betafinal",se="sefinal",t1set=t1d)
testcoloc("INS",p="pfinal", beta="betafinal",se="sefinal",t1set=t1d)
testcoloc("BCAR1",p="pfinal", beta="betafinal",se="sefinal",t1set=t1d)
testcoloc("BPTF",p="pfinal", beta="betafinal",se="sefinal",t1set=t1d)
testcoloc("MTMR3-ASCC2",p="pfinal", beta="betafinal",se="sefinal",t1set=t1d)

#do these again but removing variants where the difference in info score betwen cases and controls >0.01 just to check they all hold up to this level of scrutiny:
load(file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/resinit_adjusted_newml_5pc_diff_dups_dropped.RData")
meta$diff1<-abs(meta$controls_info_Illumina-meta$cases_info_Illumina)
meta$diff2<-abs(meta$controls_info_Affy-meta$cases_info_Affy)
meta<-meta[meta$diff1<0.01 | meta$diff2<0.01 & !is.na(meta$diff1) & !is.na(meta$diff2),]
meta$rsid<-meta$id
testcoloc("MACF1",p="pmeta",plot=F,beta="beta",se="seall",t1set=meta)
testcoloc("MAEA",p="pmeta",plot=F,beta="beta",se="seall",t1set=meta)
testcoloc("CENPW",p="pmeta",plot=F,beta="beta",se="seall",t1set=meta)
testcoloc("JAZF1",p="pmeta",plot=F,beta="beta",se="seall",t1set=meta)
testcoloc("GLIS3",p="pmeta",plot=F,beta="beta",se="seall",t1set=meta)
testcoloc("INS",p="pmeta",plot=F,beta="beta",se="seall",t1set=meta)
testcoloc("BCAR1",p="pmeta",plot=F,beta="beta",se="seall",t1set=meta)
testcoloc("BPTF",p="pmeta",plot=F,beta="beta",se="seall",t1set=meta)
testcoloc("MTMR3-ASCC2",p="pmeta",plot=F,beta="beta",se="seall",t1set=meta)



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
if(file.exists(paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Affy_all_",adj))){
af<-read.table(file=paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Affy_all_",adj), skip=10, header=T, as.is=T)
afin=T
}
if(!file.exists(paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Affy_all_",adj))){
afin=F
}
if(file.exists(paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Illumina_all_",adj))){
il<-read.table(file=paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Illumina_all_",adj), skip=10, header=T, as.is=T)
ilin=T
}
if(!file.exists(paste0("/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/",snp,"_Illumina_all_",adj))){
ilin=F
}
}

if(afin==T & ilin==T){
#remove SNPs that were dropped from the original meta analysis due to QC:
af<-af[af$rsid %in% origaf$rsid,]
il<-il[il$rsid %in% origil$rsid,]

meta<-merge(af,il,by=c("alternate_ids","rsid","chromosome","position","alleleA","alleleB"))
meta$diff.x<-abs(meta$cases_info.x-meta$controls_info.x)
meta$diff.y<-abs(meta$cases_info.y-meta$controls_info.y)
if(!is.null(adj)){
if(snp=="rs689" & adj=="1_cond_on_2"){
meta$diff.x<-ifelse(meta$rsid=="11:2182224:A:T",NA,meta$diff.x)
meta$diff.y<-ifelse(meta$rsid=="11:2182224:A:T",NA,meta$diff.y)
}
}
meta<-meta[meta$diff.x<0.05 & meta$diff.y<0.05 | is.na(meta$diff.x) | is.na(meta$diff.y),]

meta$effectdiff<-abs(meta$frequentist_add_beta_1.add.phenotype.1.x-meta$frequentist_add_beta_1.add.phenotype.1.y)
meta<-meta[meta$effectdiff<0.5 & !is.na(meta$effectdiff),]

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
}
if(afin==F & ilin==T){
meta=il[il$rsid %in% orig$rsid,]
meta$diff<-abs(meta$cases_info-meta$controls_info)
meta$diff<-abs(meta$cases_info-meta$controls_info)
meta<-meta[meta$diff<0.05 | is.na(meta$diff),]
meta$beta<-meta$frequentist_add_beta_1.add.phenotype.1
meta$seall<-meta$frequentist_add_se_1
meta$z<-meta$beta/meta$frequentist_add_se_1
meta$pmeta<-2*pnorm(-abs(meta$z))
}
if(afin==T & ilin==F){
meta=af[af$rsid %in% orig$rsid,]
meta$diff<-abs(meta$cases_info-meta$controls_info)
meta$diff<-abs(meta$cases_info-meta$controls_info)
meta<-meta[meta$diff<0.05 | is.na(meta$diff),]
meta$beta<-meta$frequentist_add_beta_1.add.phenotype.1
meta$seall<-meta$frequentist_add_se_1
meta$z<-meta$beta/meta$frequentist_add_se_1
meta$pmeta<-2*pnorm(-abs(meta$z))
}
return(meta)
message(paste0("Done ",snp))
}




#and now conditional signals on other associations in the region (T1D + T2D):
load(file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/resinit_adjusted_newml_5pc_diff_dups_dropped.RData")
metaorig<-meta
metaorig$rsid<-metaorig$id

origaf<-metaorig[!is.na(metaorig$est_Affy),c("rsid","position")]
origil<-metaorig[!is.na(metaorig$est_Illumina),c("rsid","position")]
orig<-metaorig[,c("rsid","position")]


#which of the regions above have secondary association signals?
load(file="/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/all_signals_dups_drop.RData")
#from the T1D perspective: INS has 3 signals, whilst there is evidence of 2 signals at  ZC3H4
#as for T2D, there are BCAR1, BPTF, GLIS3, HHEX-IDE, INS, (KCNQ1?), PAM

#check with conditional on these:
l1<-metathem("rs689", adj="1_cond_on_2")
l1<-l1[!is.na(l1$beta),]
l2<-metathem("rs689", adj="2_cond_on_1")
l2<-l2[!is.na(l2$beta),]
now<-metaorig[,c("rsid","controls_AF_B")]

l1<-merge(now, l1,by="rsid")
l2<-merge(now, l2,by="rsid")
t1s<-l1
t1s$id<-t1s$rsid

gett2dallign<-function(t1d, stub){
t21<-read.table(file=paste0("/well/todd/users/jinshaw/t2d_risk/conditional/coloc.input.locus.",stub,".txt"),header=T,as.is=T)
#allign:
t21$id<-paste0(t21$CHR,":",t21$POS,":",t21$NEA,":",t21$EA)
t21$alt<-paste0(t21$CHR,":",t21$POS,":",t21$EA,":",t21$NEA)
w<-which(t21$alt %in% t1d$id)
t21[w,"dum1"]<-t21[w,"NEA"]
t21[w,"dum2"]<-t21[w,"EA"]
t21[w,"NEA"]<-t21[w,"dum2"]
t21[w,"EA"]<-t21[w,"dum1"]
t21[w,"BETA"]<-t21[w,"BETA"]*-1
t21[w,"EAF"]<-1-t21[w,"EAF"]
t21$id<-paste0(t21$CHR,":",t21$POS,":",t21$NEA,":",t21$EA)
t21<-t21[t21$id %in% t1d$id,]
return(t21)
}
ts<-lapply(c("INS.IGF2_1","INS.IGF2_2",
"INS.IGF2_3","INS.IGF2_4","INS.IGF2_5"),gett2dallign, t1d=t1d)

t1s<-list(l1,l2)


testcoloc<-function(t1, t2, gene, stub, plot=TRUE){

t1$id<-t1$rsid
b<-merge(t1,t2,by="id")
b$diff<-abs(b$EAF-b$controls_AF_B)
b<-b[b$diff<0.05,]

t1<-t1[t1$id %in% b$id,]
t2<-t2[t2$id %in% b$id,]

t1$pvalues<-t1$pmeta
t1$type="cc"
t1$varbeta=t1$seall^2
t1$MAF<-ifelse(t1$controls_AF_B>0.5,
1-t1$controls_AF_B, t1$controls_AF_B)
t1$s<-0.39

t2$pvalues<-t2$P
t2$type="cc"
t2$beta=t2$BETA
t2$varbeta=t2$SE^2
t2$MAF<-ifelse(t2$EAF>0.5,
1-t2$EAF, t2$EAF)
t2$s<-0.09
df1=t1[,c("id","beta","varbeta","type","MAF","s")]
df2=t2[,c("id","beta","varbeta","type","MAF","s")]
rownames(df1)<-df1$id
rownames(df2)<-df2$id
df1<-df1[rownames(df2),]

coloc<-coloc.abf(df1, df2, p1 = 1e-04, p2 = 1e-04,
      p12 = 5e-06)
probs<-coloc[[1]]

b$logp_t2d<-log10(b$P)*-1
b$logp_t1d<-log10(b$pmeta)*-1
top<-b[b$logp_t2d==max(b$logp_t2d),]
chrom=top$chromosome
min<-min(b$position)
max<-max(b$position)
system(paste0("tabix /well/todd/users/jinshaw/t1d_risk/T1D.Affy.chr",chrom,".vcf.gz ",chrom,":",ceiling(min),"-",ceiling(max), " -h > /well/todd/users/jinshaw/t1d_risk/",gene,"_region_affy.vcf"))
system(paste0("plink --vcf /well/todd/users/jinshaw/t1d_risk/",gene,"_region_affy.vcf --double-id --make-bed --out /well/todd/users/jinshaw/t1d_risk/",gene))

init<-read.plink(fam=paste0("/well/todd/users/jinshaw/t1d_risk/",gene,".fam"),
bim=paste0("/well/todd/users/jinshaw/t1d_risk/",gene,".bim"),
bed=paste0("/well/todd/users/jinshaw/t1d_risk/",gene,".bed"))
init<-annot.plink(init)
init<-init[,colnames(init) %in% b$id]
init<-init[,b$id]

tsnp<-top$id
lds<-ld(init[,!colnames(init) %in% tsnp],init[,tsnp],stats="R.squared")
next1<-data.frame(tsnp=1)
colnames(next1)=tsnp
lds<-rbind(lds,next1)
rownames(lds)[nrow(lds)]<-tsnp
lds$id=rownames(lds)
dfb1<-merge(b,lds,by="id")
colnames(dfb1)[ncol(dfb1)]<-"ld"
o1<-ggplot(data=dfb1, aes(position, logp_t2d, colour=ld)) + geom_point() +
scale_y_continuous(name=bquote('T2D '*-log[10] *' p-value')) +
scale_color_gradient2(name=paste0("LD with T2D \nindex variant"),low = "black",
midpoint = 0.5,
mid = "orange",
high = "red",
space="Lab") +
annotate("text", x=min+75000, y=max(b$logp_t2d), label=paste0("Colocalisation \nposterior probability=",round(probs[6],digits=3))) #+
#annotate("text", x=max-75000, y=max(b$logp_t2d), label=paste0(gene))
o2<-ggplot(data=dfb1, aes(position, logp_t1d, colour=ld)) + geom_point() +
scale_y_continuous(name=bquote('T1D '*-log[10] *' p-value')) +
scale_color_gradient2(name=paste0("LD with T2D \nindex variant"),low = "black",
midpoint = 0.5,
mid = "orange",
high = "red",
space="Lab")
t<-tracks(o1,o2, xlab=paste0("Position along chromosome ",chrom))


if (plot==TRUE){
ggsave(t,file=paste0(outdir,gene,"_tracks_5pc_diff_",stub,"_dups_drop.png"), height=20, width=20, units="cm", dpi=400)

}
return(probs)
}
#examing the colocalisation in the INS region, all signals with all other signals for t1d and t2d:
testcoloc(t1=t1s[[1]],t2=ts[[1]],plot=T,stub="t1dtop_t2dtop", gene="INS")
testcoloc(t1=t1s[[1]],t2=ts[[2]],plot=T,stub="t1dtop_t2dsecond", gene="INS")
testcoloc(t1=t1s[[1]],t2=ts[[3]],plot=T,stub="t1dtop_t2dthird", gene="INS")
testcoloc(t1=t1s[[1]],t2=ts[[4]],plot=T,stub="t1dtop_t2dfourth", gene="INS")
testcoloc(t1=t1s[[1]],t2=ts[[5]],plot=T,stub="t1dtop_t2dfifth", gene="INS")
testcoloc(t1=t1s[[2]],t2=ts[[1]],plot=T,stub="t1dsecond_t2dtop", gene="INS")
testcoloc(t1=t1s[[2]],t2=ts[[2]],plot=T,stub="t1dsecond_t2dsecond", gene="INS")
testcoloc(t1=t1s[[2]],t2=ts[[3]],plot=T,stub="t1dsecond_t2dthird", gene="INS")
testcoloc(t1=t1s[[2]],t2=ts[[4]],plot=T,stub="t1dsecond_t2dfourth", gene="INS")
testcoloc(t1=t1s[[2]],t2=ts[[5]],plot=T,stub="t1dsecond_t2dfifth", gene="INS")


#try the colocalising one again but remove all variants with >0.01 difference in info score:
try<-t1s[[2]]
try$diff1<-abs(try$controls_info.x-try$cases_info.x)
try$diff2<-abs(try$controls_info.y-try$cases_info.y)
try<-try[try$diff1<0.01 | try$diff2<0.01 & !is.na(try$diff1) & !is.na(try$diff2),]
testcoloc(t1=try,t2=ts[[1]],plot=F,stub="t1dsecond_t2dtop", gene="INS")




#now examining the seconadry associations at BCAR1, BPTF, GLIS3, HHEX-IDE, INS, (KCNQ1?), PAM:
#using the univariable summary stats for T1D
t1d$rsid<-t1d$id
t1d$pmeta<-t1d$pfinal
t1d$beta<-t1d$betafinal
t1d$seall<-t1d$sefinal
load(file="/well/todd/users/jinshaw/t2d_risk/regionsoverlapping_dups_drop.RData")
#BCAR1:

gett1d<-function(gene, number,plot=T){
message(paste0("Testing colocalisation of disease and ",gene,"\n"))
chrom<-lessall[lessall$gene==gene,"chromosome"][1]
po<-l[l$gene==gene,]
mid<-po$mi+((po$ma-po$mi)/2)
min<-mid-250000
max<-mid+250000
t1<-t1d[t1d$chromosome==chrom & t1d$position>min & t1d$position<max,]

nums<-c(1:number)
examine<-paste0(gene,"_",nums)

ts<-lapply(examine,gett2dallign, t1d=t1)

docoloc<-function(num,stub,gene){
testcoloc(t1=t1,t2=ts[[num]],plot=plot,stub=stub, gene=gene)
}
orders<-c("top","second","third","fourth")
stubs<-paste0("t1dtop_t2d",orders)
stubs<-stubs[c(1:number)]
do<-mapply(docoloc,num=c(1:number),stub=stubs,gene=c(rep(gene,number)))
}

gett1d("BCAR1",2)
gett1d("GLIS3",3)

#do these again but removing variants where the difference in info score betwen cases and controls >0.01 just to check they all hold up to this level of scrutiny:
load(file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/meta/resinit_post_imp_qc_adjusted_newml_5pc_diff_dups_dropped.RData")
t1dorig<-t1d
t1d<-meta
t1d$diff1<-abs(meta$controls_info_Illumina-meta$cases_info_Illumina)
t1d$diff2<-abs(meta$controls_info_Affy-meta$cases_info_Affy)
t1d<-t1d[t1d$diff1<0.01 | t1d$diff2<0.01 & !is.na(t1d$diff1) & !is.na(t1d$diff2),]
t1d$rsid<-t1d$id

gett1d("BCAR1",2,plot=F)
gett1d("GLIS3",3, plot=F)
 
