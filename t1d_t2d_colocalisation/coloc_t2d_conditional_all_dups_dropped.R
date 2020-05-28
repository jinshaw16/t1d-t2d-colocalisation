#coloc_t2d_conditional_all_dups_dropped.R

library(snpStats)
library(annotSnpStats)
library(coloc)
library(ggplot2)
library(GenomicRanges)
library(ggbio)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(purrr)

#receieved conditional results from the 42 identified overlapping regions to FDR significance from Anubha.
#reading these in and assessing colocalisation between all the possible T1D and T2D signals.

outdir<-"/well/todd/users/jinshaw/output/t1d_risk/hrc/vcf/meta_sardinia/newml/coloc/fdr/conditional/"
#load T1D signals:
load(file="/well/todd/users/jinshaw/t1d_risk/vcf/processed/conditional/results/all_signals_t2d_suggestive.RData")


#read in the T2D signals:
lines<-system(paste0("ls /well/todd/users/jinshaw/t2d_risk/conditional/new/"),intern=T)
lines<-gsub("_.*","",lines)
lines<-lines[!duplicated(lines)]

gett2d<-function(line){
l<-length(system(paste0("ls /well/todd/users/jinshaw/t2d_risk/conditional/new/ | grep ",line,"_ "),intern=T))
getsums<-function(line, sig){
t2<-read.table(file=paste0("/well/todd/users/jinshaw/t2d_risk/conditional/new/",line,"_",sig,".cma.cojo.gz"),header=T, as.is=T)
return(t2)
}
out<-lapply(c(1:l),getsums, line=line)
return(out)
}

t2dres<-lapply(lines,gett2d)
names(t2dres)<-lines

#read in the T1D signals:
gett1d<-function(line){
l<-length(system(paste0("ls /well/todd/users/jinshaw/t1d_risk/results/vcf_all/conditional/ | grep ^",line,"_"),intern=T))
getts<-function(line, sig){
g<-get(load(file=paste0("/well/todd/users/jinshaw/t1d_risk/results/vcf_all/conditional/",line,"_",sig,".RData")))
return(g)
}
if(l==1)
sigs<-"0"
if(l==3)
sigs<-c("0_nonorig","0")
out<-lapply(sigs,getts,line=line)
return(out)
}
linest1<-signalsall[duplicated(signalsall$line),"line"]
t1dres<-lapply(linest1, gett1d)
names(t1dres)<-linest1

load(file="/well/todd/users/jinshaw/t1d_risk/results/vcf_all/conditional/univariable_all.RData")


testcoloc<-function(line,frame1, frame2,beta,se,p, t1dsig, t2dsig, plot=TRUE){
message(paste0("Testing colocalisation of diseases in line ",line,"\n"))

bd<-merge(frame1, frame2, by="SNP")
bd<-bd[!is.na(bd$Beta),]
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
t1$N<-t1$Nt1

t2<-bd
t2$pvalues<-t2$Pvalue
t2$type="cc"
t2$beta=t2$Beta
t2$varbeta=t2$SE^2
t2$MAF<-ifelse(t2$EAF>0.5,
1-t2$EAF, t2$EAF)
t2$s<-0.09
t2$N<-t2$n
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
probs$alleleA<-top$alleleA.x
probs$alleleB<-top$alleleB.x
probs$minpos<-min(bd$position)
probs$maxpos<-max(bd$position)
probs$chromosome<-max(bd$chromosome)
if(plot==TRUE){
system(paste0("tabix /well/todd/users/jinshaw/t1d_risk/T1D.Affy.chr",max(bd$chromosome),".vcf.gz ",max(bd$chromosome),":",ceiling(min(bd$position)),"-",ceiling(max(bd$position)), " -h > /well/todd/users/jinshaw/t1d_risk/regs/",line,"_region_affy.vcf"))
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
scale_y_continuous(name=bquote('T2D '*-log[10] *' p-value'))  +
scale_color_gradient2(name=paste0("LD with T2D \nindex variant"),low = "black",
midpoint = 0.5,
mid = "orange",
high = "red",
space="Lab") #+
#annotate("text", x=min+75000, y=max(b$logp_t2d), label=paste0("Colocalisation \nposterior probability=",round(probs[6],digits=3))) #+
#annotate("text", x=max-75000, y=max(b$logp_t2d), label=paste0(gene))
o2<-ggplot(data=dfb1, aes(position, logp_t1d, colour=ld)) + geom_point() +
scale_y_continuous(name=bquote('T1D '*-log[10] *' p-value')) +
scale_color_gradient2(name=paste0("LD with T2D \nindex variant"),low = "black",
midpoint = 0.5,
mid = "orange",
high = "red",
space="Lab")

data(genesymbol, package = "biovizBase")
g <- genesymbol[seqnames(genesymbol) == paste0('chr',max(bd$chromosome))]
gr <- GRanges(
      seqnames = Rle(c(paste0("chr",max(bd$chromosome))), c(1)),
      ranges = IRanges(min(bd$position), end = max(bd$position)))

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

t1<-autoplot(Homo.sapiens, which = g)

t<-tracks(t1,o1,o2, xlab=paste0("Position along chromosome ",max(bd$chromosome)))
ggsave(t,file=paste0(outdir,gene,"_",line,"_tracks_5pc_diff_dups_drop_t1dsig_",t1dsig,"_t2dsig_",t2dsig,".png"), height=20, width=20, units="cm", dpi=400)
}
if(plot==FALSE){
data(genesymbol, package = "biovizBase")
g <- genesymbol[seqnames(genesymbol) == paste0('chr',max(bd$chromosome))]
gr <- GRanges(
      seqnames = Rle(c(paste0("chr",max(bd$chromosome))), c(1)),
      ranges = IRanges(min(bd$position), end = max(bd$position)))

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







#now want to carry out all the colocalisations across all the conditional signals.
#1) Conditional T2D signals on univariable T1D associations:
t2dcondassoc<-function(line){
t2dassocs<-t2dres[[as.character(line)]]
t1s<-resboth[[as.character(line)]]

frame1<-t1s
frame1$SNP<-frame1$id
frame1$Nt1<-ifelse(is.na(frame1$est_Illumina), 5268, ifelse(is.na(frame1$est_Affy),7977, 
ifelse(!is.na(frame1$est_Illumina) & !is.na(frame1$est_Affy),5268+7977,NA)))
frame1<-frame1[,c("SNP","id","ID","alleleA","alleleB","est_uk","se_uk","controls_AF_B","p_uk","Nt1","chromosome","position")]
colocthem<-function(t2dassocs,sig, frame1){
frame2<-t2dassocs[[sig]]
if(sig>1){
getmins<-function(signals){
getrid<-t2dassocs[[signals]]
getrid<-getrid[!is.na(getrid$pC),]
getrid<-getrid[getrid$pC==min(getrid$pC) & !is.na(getrid$pC),]
return(getrid)
}
getrids<-lapply(c(2:sig),getmins)
getrids<-do.call("rbind",getrids)
frame2<-frame2[!frame2$SNP %in% getrids$SNP,]
}
frame2$chrom<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\1",frame2$SNP)
frame2$alleleA<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\3",frame2$SNP)
frame2$alleleB<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\4",frame2$SNP)
frame2$SNPalt<-paste0(frame2$chrom,":",frame2$bp,":",frame2$alleleB,":",frame2$alleleA)


frame2$SNP<-ifelse(frame2$SNPalt %in% frame1$SNP,frame2$SNPalt, frame2$SNP)
frame2$alleleA<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\3",frame2$SNP)
frame2$alleleB<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\4",frame2$SNP)
w<-which(frame2$refA != frame2$alleleB)
frame2[w,"b"]<-frame2[w,"b"]*-1
frame2[w,"freq"]<-1-frame2[w,"freq"]
frame2$Pvalue<-frame2$pC
frame2$Beta<-frame2$bC
frame2$SE<-frame2$bC_se
frame2$EAF<-frame2$freq

out<-testcoloc(line=line,frame1=frame1, frame2=frame2, beta="est_uk",se="se_uk", p="p_uk",t1dsig=1, t2dsig=sig)
return(out)
}
probs<-lapply(c(1:length(t2dassocs)), colocthem, t2dassocs=t2dassocs, frame1=frame1)
probs<-do.call("rbind",probs)
probs$t1dsig=1
probs$t2dsig=c(1:length(t2dassocs))
probs$line=line
message(paste0("Done ",line))
return(probs)
}
allconds1<-lapply(lines, t2dcondassoc)
allconds1<-do.call("rbind",allconds1)







#2) Conditional T1D signals on the univariable T2D associations:
t1dcondassoc<-function(line){
t1dassocs<-t1dres[[as.character(line)]]
t2s<-resboth[[as.character(line)]]
t2s$n<-t2s$Neff
t2s$SNP<-t2s$id
t2s$chrom<-t2s$chromosome
t2s$bp<-t2s$position
frame2<-t2s[,c("SNP","n","alleleA","alleleB","Beta","SE","EAF","Pvalue","chrom","bp")]
colocthem2<-function(t1dassocs,sig){
frame1<-t1dassocs[[sig]]
frame1$SNP<-frame1$id
frame1$Nt1<-5268+7977
p1<-resboth[[as.character(line)]]
p1<-p1[,c("id","alleleA","alleleB","controls_AF_B")]
frame1$ord<-c(1:nrow(frame1))
frame1<-merge(frame1,p1,by=c("id","alleleA","alleleB"))
frame1<-frame1[order(frame1$ord),]
frame1<-frame1[,c("SNP","id","ID","alleleA","alleleB","est_uk","se_uk","controls_AF_B","p_uk","Nt1","chromosome","position")]

out<-testcoloc(line=line,frame1=frame1, frame2=frame2, beta="est_uk",se="se_uk", p="p_uk",t1dsig=sig, t2dsig=1)
return(out)
}
probs<-lapply(c(1:length(t1dassocs)), colocthem2, t1dassocs=t1dassocs)
probs<-do.call("rbind",probs)
probs$t2dsig=1
probs$t1dsig=c(1:length(t1dassocs))
probs$line=line
return(probs)
}
allconds2<-lapply(linest1, t1dcondassoc)
allconds2<-do.call("rbind",allconds2)




#3) conditional T1D signals on conditional T2D signals:
condcond<-function(line){
t2dassocs<-t2dres[[as.character(line)]]
t1dassocs<-t1dres[[as.character(line)]]
l<-c(1:length(t2dassocs))
l1<-c(1:length(t1dassocs))
cb<-cross2(l1, l)
cb<-do.call("rbind",cb)

docomb<-function(sig1,sig2){
frame1<-t1dassocs[[sig1]]
frame1$SNP<-frame1$id
frame1$Nt1<-5268+7977
p1<-resboth[[as.character(line)]]
p1<-p1[,c("id","alleleA","alleleB","controls_AF_B")]
frame1$ord<-c(1:nrow(frame1))
frame1<-merge(frame1,p1,by=c("id","alleleA","alleleB"))
frame1<-frame1[order(frame1$ord),]
frame1<-frame1[,c("SNP","id","ID","alleleA","alleleB","est_uk","se_uk","controls_AF_B","p_uk","Nt1","chromosome","position")]


frame2<-t2dassocs[[sig2]]
frame2$chrom<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\1",frame2$SNP)
frame2$alleleA<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\3",frame2$SNP)
frame2$alleleB<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\4",frame2$SNP)
frame2$SNPalt<-paste0(frame2$chrom,":",frame2$bp,":",frame2$alleleB,":",frame2$alleleA)


frame2$SNP<-ifelse(frame2$SNPalt %in% frame1$SNP,frame2$SNPalt, frame2$SNP)
frame2$alleleA<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\3",frame2$SNP)
frame2$alleleB<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\4",frame2$SNP)
w<-which(frame2$refA != frame2$alleleB)
frame2[w,"b"]<-frame2[w,"b"]*-1
frame2[w,"freq"]<-1-frame2[w,"freq"]
frame2$Pvalue<-frame2$pC
frame2$Beta<-frame2$bC
frame2$SE<-frame2$bC_se
frame2$EAF<-frame2$freq

out<-testcoloc(line=line,frame1=frame1, frame2=frame2, beta="est_uk",se="se_uk", p="p_uk",t1dsig=sig1, t2dsig=sig2)
out$t1dsig=sig1
out$t2dsig=sig2
return(out)
}
probs<-mapply(docomb,sig1=do.call("rbind",cb[,1])[,1],sig2=do.call("rbind",cb[,2])[,1],SIMPLIFY=F)
probs<-do.call("rbind",probs)
probs$line=line
return(probs)
}

lboth<-lines[lines %in% linest1]
allconds3<-lapply(lboth,condcond)
allconds3<-do.call("rbind",allconds3)



#Finally, checking line 21 T2D with line 20 T1D (both very close to each other):
t2dassocs<-t2dres[[as.character("21")]]
t1dassocs<-t1dres[[as.character("20")]]
l<-c(1:length(t2dassocs))
l1<-c(1:length(t1dassocs))
cb<-cross2(l1, l)
cb<-do.call("rbind",cb)

docomb<-function(sig1,sig2){
frame1<-t1dassocs[[sig1]]
frame1$SNP<-frame1$id
frame1$Nt1<-5268+7977
p1<-resboth[[as.character("20")]]
p1<-p1[,c("id","alleleA","alleleB","controls_AF_B")]
frame1$ord<-c(1:nrow(frame1))
frame1<-merge(frame1,p1,by=c("id","alleleA","alleleB"))
frame1<-frame1[order(frame1$ord),]
frame1<-frame1[,c("SNP","id","ID","alleleA","alleleB","est_uk","se_uk","controls_AF_B","p_uk","Nt1","chromosome","position")]


frame2<-t2dassocs[[sig2]]
frame2$chrom<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\1",frame2$SNP)
frame2$alleleA<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\3",frame2$SNP)
frame2$alleleB<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\4",frame2$SNP)
frame2$SNPalt<-paste0(frame2$chrom,":",frame2$bp,":",frame2$alleleB,":",frame2$alleleA)


frame2$SNP<-ifelse(frame2$SNPalt %in% frame1$SNP,frame2$SNPalt, frame2$SNP)
frame2$alleleA<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\3",frame2$SNP)
frame2$alleleB<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\4",frame2$SNP)
w<-which(frame2$refA != frame2$alleleB)
frame2[w,"b"]<-frame2[w,"b"]*-1
frame2[w,"freq"]<-1-frame2[w,"freq"]
frame2$Pvalue<-frame2$pC
frame2$Beta<-frame2$bC
frame2$SE<-frame2$bC_se
frame2$EAF<-frame2$freq
out<-testcoloc(line=2021,frame1=frame1, frame2=frame2, beta="est_uk",se="se_uk", p="p_uk",t1dsig=sig1, t2dsig=sig2, plot=FALSE)
out$t1dsig=sig1
out$t2dsig=sig2
return(out)
}
probs<-mcmapply(docomb,sig1=do.call("rbind",cb[,1])[,1],sig2=do.call("rbind",cb[,2])[,1],SIMPLIFY=F,mc.cores=4)
probs<-do.call("rbind",probs)



allresults<-rbind(allconds1,allconds2)
allresults<-rbind(allresults, allconds3)


allresults<-allresults[,c("line","id","ID","chromosome","position","gene","alleleA","alleleB","t1d_beta","t1d_se","t1d_p","t2d_beta","t2d_se","t2d_p","t1dsig","t2dsig","P0","P1","P2","P3","P4","minpos","maxpos")]
write.table(allresults, file=paste0(outdir,"allres_fdr.txt"),col.names=T,row.names=F,quote=F, sep="\t")
