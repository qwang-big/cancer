# Determine cell number and size in spatialtranscriptome

## Most expressed genes
```
zcat FP200000241TR_F2_web_2.txt.gz |perl -F'\t' -ne 'print if $F[5]>50'|cut -f3|suc > x.xls
```

## Find cell center
```
zcat FP200000241TR_F2_web_2.txt.gz |perl -F'\t' -ne '$d=$F[2];print if $d eq "Ptprd" or $d eq "Kcnip4" or $d eq "Celf2" or $d eq "Rbfox1" or $d eq "Pcdh9" or $d eq "Nrg3" or $d eq "Dlg2" or $d eq "Nrxn3" or $d eq "Csmd1" or $d eq "Nrxn1"' > intron.tsv
```

## Filter centers
```
perl -F'\t' -ne '$x=int($F[0]/3+0.5);$y=int($F[1]/3+0.5);$h{"$x\t$y"}+=$F[3];END{foreach(keys %h){print $h{$_},"\t",$_,"\n"}}' ../intron.tsv  > center.1
```
