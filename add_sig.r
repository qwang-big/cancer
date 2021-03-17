Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
library(ggplot2)
library(reshape2)

add_sig=function(x,y,d,s=1){
x=x[x$geneID %in% y[,1],]
x[,2:3]=round(x[,2:3]/10)*10
x1=dcast(x,x~y, value.var = "MIDCounts", fun.aggregate = sum)
x <- x1[,-1]
rownames(x) <- x1[,1]
d1=melt(as.matrix(x)*s)
colnames(d1)=colnames(d)
d1
}

x=read.table('Sample30-T-FB13-DP8400015282BL_B5.bin1.Lasso.gene.txt',header=T)
y=read.table('tm1.txt')
d=unique(round(x[,2:3]/10)*10)
d$MIDCounts=1
d1=add_sig(x,y,d)
y=read.table('tm2.txt')
d2=add_sig(x,y,d,-1)
d=rbind(d,d1)
d=rbind(d,d2)
#d[d$MIDCounts>20,3]=20
ggplot(d[d$MIDCounts!=0,], aes(x,y, fill=MIDCounts)) + geom_point(size=1, shape=23, stroke=0) + scale_fill_gradient2(low='blue',mid='white',high='red',midpoint = 0, limits = c(-10,10))
