x=read.table(text='0	658	Epithelial cells	FB11final_50
1	605	Epithelial cells	FB11final_50
2	536	Macrophages	FB11final_50
3	520	Epithelial cells	FB11final_50
4	446	Epithelial cells	FB11final_50
5	430	Epithelial cells	FB11final_50
6	412	Macrophages	FB11final_50
7	341	Epithelial cells	FB11final_50
0	658	Fibroblasts	FJ13final_50_2
1	605	Epithelial cells	FJ13final_50_2
2	536	Fibroblasts	FJ13final_50_2
3	520	Epithelial cells	FJ13final_50_2
4	446	Fibroblasts	FJ13final_50_2
5	430	Epithelial cells	FJ13final_50_2
6	412	Epithelial cells	FJ13final_50_2
7	341	Macrophages	FJ13final_50_2
8	240	Myofibroblasts	FJ13final_50_2
9	22	Epithelial cells	FJ13final_50_2
0	658	Epithelial cells	FD16final_50
1	605	Epithelial cells	FD16final_50
2	536	Epithelial cells	FD16final_50
4	446	Macrophages	FD16final_50
5	430	Myofibroblasts	FD16final_50
6	412	Macrophages	FD16final_50
7	341	Epithelial cells	FD16final_50
0	658	Epithelial cells	FB12final_50
1	605	Epithelial cells	FB12final_50
2	536	Macrophages	FB12final_50
3	520	Epithelial cells	FB12final_50
4	446	Epithelial cells	FB12final_50
5	430	Macrophages	FB12final_50
6	412	Epithelial cells	FB12final_50
7	341	Epithelial cells	FB12final_50
8	240	Epithelial cells	FB12final_50
0	658	Epithelial cells	FE11final_50
1	605	Epithelial cells	FE11final_50
2	536	Epithelial cells	FE11final_50
3	520	Epithelial cells	FE11final_50
4	446	Macrophages	FE11final_50
5	430	Myofibroblasts	FE11final_50
0	658	Epithelial cells	FD13final_50
1	605	Epithelial cells	FD13final_50
2	536	Epithelial cells	FD13final_50
3	520	Epithelial cells	FD13final_50
4	446	Epithelial cells	FD13final_50
4	430	B cells	FD13final_50
0	658	Macrophages	FH16final_50
1	605	Fibroblasts	FH16final_50
2	536	Epithelial cells	FH16final_50
3	520	Epithelial cells	FH16final_50
4	446	Epithelial cells	FH16final_50
5	430	Endothelial cells	FH16final_50
6	412	Endothelial cells	FH16final_50
7	341	Epithelial cells	FH16final_50
8	240	Endothelial cells	FH16final_50
9	22	Fibroblasts	FH16final_50
0	658	Epithelial cells	FJ12final_50_2
1	605	Epithelial cells	FJ12final_50_2
2	536	Epithelial cells	FJ12final_50_2
3	520	Fibroblasts	FJ12final_50_2
4	446	B cells	FJ12final_50_2
5	430	Macrophages	FJ12final_50_2
6	412	Macrophages	FJ12final_50_2
0	658	Epithelial cells	FE12final_50
1	605	Epithelial cells	FE12final_50
2	536	Epithelial cells	FE12final_50
3	520	Macrophages	FE12final_50
4	446	Macrophages	FE12final_50
5	430	Epithelial cells	FE12final_50
6	412	Fibroblasts	FE12final_50
7	341	Myofibroblasts	FE12final_50
0	658	Epithelial cells	FH11final_50
1	605	Epithelial cells	FH11final_50
2	536	Epithelial cells	FH11final_50
3	520	Epithelial cells	FH11final_50
4	446	Epithelial cells	FH11final_50
5	430	Epithelial cells	FH11final_50
6	412	Macrophages	FH11final_50
7	341	Epithelial cells	FH11final_50
8	240	Myofibroblasts	FH11final_50
9	22	Fibroblasts	FH11final_50
0	658	Fibroblasts	FH14final_50
1	605	Epithelial cells	FH14final_50
2	536	Myofibroblasts	FH14final_50
3	520	Fibroblasts	FH14final_50
4	446	Macrophages	FH14final_50
5	430	Epithelial cells	FH14final_50', sep='\t')
tt=c(
'FD13final_50',
'FD16final_50',
'FH14final_50',
'FJ13final_50_2')
nn=c(
'FB11final_50',
'FB12final_50',
'FE11final_50',
'FE12final_50',
'FH11final_50',
'FJ12final_50_2')
df=data.frame(k=tt,v='T')
df=rbind(df,data.frame(k=nn,v='N'))
df=rbind(df,data.frame(k='FH16final_50',v='T'))
rownames(df)=df$k
x$V5=df[x$V4,'v']
x$V4=substr(x$V4,1,4)
p=ggplot(x,aes(V4,V2,fill=V3))+geom_col(position='fill')+theme(legend.position='top')+labs(x='',y='')+facet_wrap(~V5,scales='free_x')
ggsave(p,file='cancer_by_patient.pdf',width=6,height=6)
p=ggplot(x,aes(V2,V3,fill=V4))+geom_col(position='fill')+theme(legend.position='top')+labs(x='',y='')
ggsave(p,file='celltypes_by_patient.pdf',width=6,height=6)

