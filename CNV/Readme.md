reproduce unprocessed files
```
cat all.1|perl -ne 'next if length($_)<88;$s=substr($_,70);chomp $s;$s=~s/ /./g;$s="$s.rds";$h1{$s}=$_;END{@f=`ls -1 *.rds`;foreach(@f){chomp;$h{$_}=1}foreach(keys %h1){print $h1{$_} unless defined $h{$_}}}' |split -a 3 --numeric-suffixes=100 -l 3
```
