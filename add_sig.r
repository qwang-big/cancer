Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
library(ggplot2)
library(reshape2)

theme_set(theme_bw())
theme_update(panel.grid.minor = element_line(colour = NA),
panel.grid.major = element_line(colour = NA))

colMax <- function(X) apply(X, 2, max)
colMin <- function(X) apply(X, 2, min)

add_sig_rect=function(x,y,X,Y,s=1){
x=x[x$geneID %in% y,]
x=x[x$x %in% X,]
x=x[x$y %in% Y,]
x[,-1]
}

add_sig=function(x,y,d,s=1){
x=x[x$geneID %in% y,]
x[,J]=round(x[,J]/K)*K
x1=dcast(x,x~y, value.var = "MIDCounts", fun.aggregate = sum)
x <- x1[,-1]
rownames(x) <- x1[,1]
d1=melt(as.matrix(x)*s)
colnames(d1)=colnames(d)
d1
}

J=2:3
K=20
x=read.table('Sample31-T-FJ13_DP8400016195TL_F4.bin1.Lasso.gene.txt',header=T)
y=unique(x$geneID[grep('^MT-',x$geneID)])
d=unique(round(x[,J]/K)*K)
d$MIDCounts=0
m1=add_sig(x,y,d)
m1$g=1
m1$MIDCounts = m1$MIDCounts/10
v=matrix(c(colMin(x[,J]),max=colMax(x[,J])),nrow=2,byrow=T)
v1=diff(v)
vf=v1[1]/v1[2]
N=20
M=round(v1[2]/v1[1]*N)
x1=seq(v[1,1],v[2,1])
y1=seq(v[1,2],v[2,2])
X=split(x1,cut_number(x1,N))
Y=split(y1,cut_number(y1,M))
for(i in 1:N){
dx=X[[i]]
for(j in 1:M){
dy=Y[[j]]
m= add_sig_rect(x,y,dx,dy)
if(nrow(m)<10) next
m$g=2
den=round(nrow(m)/(max(dx)-min(dx))/(max(dy)-min(dy)),3)
p=ggplot() + geom_point(data=m, aes(x,y, fill=MIDCounts), size=1, shape=22, stroke=0) + geom_point(data=m1, aes(x,y, fill=MIDCounts), size=1, shape=22, stroke=0) +geom_rect(data=m1, aes(xmin=min(dx), xmax=max(dx), ymin=min(dy), ymax=max(dy)), color="black", fill=NA)+annotate(geom="text", x=min(dx)+10, y=min(dy), label=den)+ facet_wrap(~g, scales = "free") + scale_fill_gradient2(low='blue',mid='white',high='red',midpoint = 0, limits = c(-10,10), oob=scales::squish)
ggsave(p,filename=paste0('figs/',i,'-',j,'.jpg'),width=13.6,height=7)
}}
