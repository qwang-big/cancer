## LDSC
```sh
python ~/Downloads/ldsc-master/ldsc.py \
	--h2 PASS_Parkinsons_23andMe_Corces2020.sumstats \
	--ref-ld-chr baseline/baseline. \
	--w-ld-chr weights_hm3_no_hla/weights. \
	--overlap-annot \
	--frqfile-chr 1000G_frq/1000G.mac5eur. \
	--out AD_baseline
	
find -name *.bet|perl -ne '@a=("","Young","Middle","Old","Geriatric");chomp;$s=/\/d/?"Down":"Up";$i=$a[$1] if /(\d+)\.b/;/\/(\S+?)\//;print "cut -f2,4,5 $_|sort -k1,1 -k2,2n > ~/tmp/tt/$1\_$s\_$i.bed\n"' > cp.sh

wget https://hgdownload.soe.ucsc.edu/goldenPath/macFas5/liftOver/macFas5ToHg19.over.chain.gz
for f in *.bed; do /var/bin/liftOver $f macFas5ToHg19.over.chain annot/$f /dev/null;done

for chrom in {1..22};
do
for f in annot/*.bed;
do
python ~/Downloads/ldsc-master/make_annot.py \
--bed-file $f \
--bimfile 1000G_plinkfiles/1000G.mac5eur.${chrom}.bim \
--annot-file $f.${chrom}.annot.gz
done;done

for chrom in {1..22};
do
python ~/Downloads/ldsc-master/ldsc.py \
    --l2 \
    --bfile 1000G_plinkfiles/1000G.mac5eur.${chrom} \
    --ld-wind-cm 1 \
    --annot annot/cortex.${chrom}.annot.gz \
    --out out/cortex.${chrom}
done

python ~/Downloads/ldsc-master/ldsc.py \
	--h2 PASS_Alzheimers_Jansen2019.sumstats \
	--ref-ld-chr out/cortex. \
	--w-ld-chr weights_hm3_no_hla/weights. \
	--overlap-annot \
	--frqfile-chr 1000G_frq/1000G.mac5eur. \
	--out AD_cortex
	
```
