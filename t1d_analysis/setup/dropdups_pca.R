#!/usr/bin/env Rscript
#dropping any variants with duplicate positions.
args = commandArgs(trailingOnly=TRUE)


r<-read.table(file=paste0("/well/todd/users/jinshaw/t1d_risk/processed/",args[1],"_chr",args[2],".bim"),as.is=T)

drop<-r[duplicated(r$V4),'V2']
write.table(drop,file=paste0("/well/todd/users/jinshaw/t1d_risk/processed/dropdup_",args[1],"_",args[2],".txt"), quote=F, col.names=F, row.names=F)

