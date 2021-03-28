# cancer

## group-specific markers
```
cut -d, -f5 FJ13markers.csv|head -n11 |tail -n10 > tm1.txt
cut -d, -f17 FJ13markers.csv|head -n11 |tail -n10 > tm2.txt
cut -d, -f2 FJ13markers.csv|head -n11 |tail -n10 > fb1.txt
cut -d, -f14 FJ13markers.csv|head -n11 |tail -n10 > fb2.txt
```
## gene/umi per spot
```
for f in ?ample*;do perl umi_per_spot.pl $f > ${f%.*}.gene 2> ${f%.*}.umi;done
```
## pathway heatmap
```
library(ggplot2)
theme_set(theme_bw())
theme_update(panel.grid.minor = element_line(colour = NA),
panel.grid.major = element_line(colour = NA))

for(i in 1:4) x=rbind(x,data.frame(V1=unique(x$V1),V2=i,V3=1,V4='N'))
p=ggplot(x,aes(V2,V1,fill=V3))+geom_tile()+facet_wrap(~V4,scales='free_x')+scale_fill_gradient2(name='p-value',low='red',high='white',mid='white',midpoint=0.5) + labs(x='Patient',y='')+theme(axis.text.x=element_blank())
ggsave(p,filename='pathways.pdf')
```
