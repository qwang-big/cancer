Sys.setlocale('LC_NUMERIC','C')
options(stringsAsFactors = FALSE)
library(ggplot2)
library(reshape2)

theme_set(theme_bw())
theme_update(panel.grid.minor = element_line(colour = NA),
panel.grid.major = element_line(colour = NA))

ff=dir()
ff=ff[grep('gene$',ff)]
fd=factor(seq_len(length(ff)))
names(fd)=ff
dd=do.call(rbind,lapply(ff[grep('gene$',ff)],function(f){d=read.table(f);d$S=substr(f,7,8);d$T=substr(f,10,10);d$R=fd[f];d}))

p = ggplot(dd, aes(x=R, y=V1, color=T)) + geom_violin() + facet_wrap(~S,scales='free') + labs(title = 'Genes per Spot', x='Patient',y='')+theme(axis.text.x=element_blank())
ggsave(p, filename='gene.pdf')
