reproduce unprocessed files
```
cat all.1|perl -ne 'next if length($_)<88;$s=substr($_,70);chomp $s;$s=~s/ /./g;$s="$s.rds";$h1{$s}=$_;END{@f=`ls -1 *.rds`;foreach(@f){chomp;$h{$_}=1}foreach(keys %h1){print $h1{$_} unless defined $h{$_}}}' |split -a 3 --numeric-suffixes=100 -l 3
```
produce scripts
```
ls -1 *ample*.rda|perl -ne 'chomp;push @a,$_;END{for $i(1..22,"X"){for ($j=0;$j<@a;$j+=2){ $x=$a[$j+1];$y=$a[$j];open(O,">s".(++$k).".sh");print O "/hwfssz5/ST_PRECISION/OCG/wangqi/miniconda3/envs/r4/bin/Rscript run.r $x $y $i\n";close O}}}'
```
