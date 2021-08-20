library(stringr)
ff=dir('.','*.annot.gz$')
df=data.frame(ff,str_extract(ff,'\\.\\d+\\.'))
df=split(df,df[,2])
dd=lapply(df,function(dd) {
x=apply(dd,1,function(d) read.table(gzfile(d[1]),header=T))
x=do.call(cbind,x)
colnames(x)=data.frame(do.call(rbind, strsplit(dd[,1], '.', fixed=T)))[,1]
cbind(read.table(gzfile(paste0('../baseline/baseline',dd[1,2],'annot.gz')),header=T)[,1:5],x)})
lapply(seq_along(df), function(i) write.table(dd[[i]],file=paste0('cortex',names(df)[i],'annot'),sep='\t',quote=F,row.names=F))
