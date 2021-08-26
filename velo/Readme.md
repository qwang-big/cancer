## batch script
```sh
ls -1 *.target.bam|perl -ne 'chomp;open(O,">s".(++$i).".sh");print O "/share/app/samtools/1.11/bin/samtools view $_|perl cb.pl > ../w/$_.in 2> ../w/$_.ex\n";close O'
for i in {1..42};do qsub -clear -cwd -l vf=1g,p=1 -binding linear:1 -q st.q -P P20Z10200N0059 s$i.sh;done
for f in $(find -name features.tsv.gz);do echo $f;zgrep -n -P '^A2M\t' $f;done
for f in *.ex;do perl -ne 'chomp;$h{$_}++;END{foreach(keys %h){print "$_\t$h{$_}\n"}}' $f > ${f%.*}.cnt;done
```
## run scVelo
```py
import loompy, os, scipy
import scvelo as scv
s='/outs/filtered_feature_bc_matrix/'

ff=os.listdir()
ff=[f for f in ff if len(f)==16]
loompy.create_from_cellranger(f)
d = scv.read_loom(f+'/'+f+'.loom')
x= scipy.io.mmread(f+s+'spliced.mtx.gz')
d.layers["spliced"] = x.T.tocsr()
d.layers["unspliced"] = d.X - d.layers["spliced"]
scv.tl.velocity(d, mode = "stochastic")
scv.tl.velocity_graph(d)
d.obsm['X_umap'] = pd.read_csv(f+s+'barcodes.tsv.gz',compression='gzip',sep='_',header=None).to_numpy()
scv.pl.velocity_embedding_stream(d, basis='X_umap', save = f+'/'+f+'.pdf')
```
