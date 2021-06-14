Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
library(ggplot2)
library(umap)
set.seed(123)

args = commandArgs(trailingOnly=TRUE)
canvas.size=function(d) (max(d[,1])-min(d[,1])) * (max(d[,2])-min(d[,2]))
lw=function(d) length(which(d))
rr=3:15
df=do.call(rbind,lapply(rr,function(r){
x=read.table(gzfile(paste0(args[1],r,'.gz')),row.names=1)
i=apply(x,1,function(d) lw(d>0))
j=apply(x,2,function(d) lw(d>0))
x=x[i>=3,j>=3]
if (class(x)!="data.frame"||any(dim(x)<20))
return(data.frame())
print(r)
y=umap(x)
d=unique(round(y$layout*10))
print(nrow(d)/canvas.size(d))
data.frame(y$layout,r)
}))
p=ggplot(df,aes(X1,X2))+geom_point(size = 0.1)+facet_wrap(~r)
ggsave(p,file=paste0(args[1],'.pdf'))

