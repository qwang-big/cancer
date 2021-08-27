x=read.table('down.xls')
x$V3=seq_len(nrow(x))
ggplot(x,aes(V3,V1,fill=V4))+geom_bar(stat='identity')+facet_wrap(.~V2,scale='free_x')

ff=dir()
ff=ff[grep('res',ff)]
df=do.call(cbind,lapply(ff,function(f) read.table(f,header=T)$Enrichment_p))
x=read.table(ff[1],header=T)
rownames(df)=x[,1]
colnames(df)=unlist(lapply(strsplit(gsub('_','.',ff), '.', fixed=T),function(d)d[2]))
df=df[2:5,]
rownames(df)=c('Geriatric','Middle','Old','Young')
colnames(df)[3]='Smoked'
df=-log(df)
df=melt(df)
ggplot(data=df, aes(x=Var2, y=value, fill=Var1)) + geom_bar(stat="identity", position=position_dodge()) + theme_minimal() + labs(x='',y='Enrichment -logP')
