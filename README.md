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
