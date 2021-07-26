options(stringsAsFactors = FALSE)
colorpanel <- function (n, low, high) {
    low <- col2rgb(low)
    high <- col2rgb(high)
    red <- seq(low[1, 1], high[1, 1], length = n)/255
    green <- seq(low[3, 1], high[3, 1], length = n)/255
    blue <- seq(low[2, 1], high[2, 1], length = n)/255
    rgb(red, blue, green)
}
ff=dir()
lapply(ff[grep('txt$',ff)],function(f){
print(f)
x=read.table(f,as.is=T,sep='\t')
x=unique(x)
df=do.call(rbind,strsplit(substr(x$V1,1,nchar(x$V1)-4),'_'))
df=data.frame(df)
x=cbind(df,x)
x[,1]=as.integer(x[,1])
x[,2]=as.integer(x[,2])
df=read.table(text=do.call(rbind,strsplit(substr(x$V3,2,nchar(x$V3)-1),'_')))
x=cbind(x[,1:2],df)
cl=colorpanel(10,'white','red')
x[,3:6]=round(10*x[,3:6])+1
p=function(d,cl,tt) plot(d,pch=22,bg=cl,cex=1.1,main=tt, xaxt='n', yaxt='n', xlab='', ylab='')
d=x
#d[,2]=max(d[,2])-d[,2]
png(paste0(f,'.png'),width=7*max(d[,1]),height=8*max(d[,2]))
p(d[,1:2], cl[d$V3], 'LUAD')
dev.off()

})
