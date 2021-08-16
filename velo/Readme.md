## batch script
```sh
ls -1 *.target.bam|perl -ne 'chomp;open(O,">s".(++$i).".sh");print O "/share/app/samtools/1.11/bin/samtools view $_|perl cb.pl > ../w/$_.in 2> ../w/$_.ex\n";close O'
for i in {1..42};do qsub -clear -cwd -l vf=1g,p=1 -binding linear:1 -q st.q -P P20Z10200N0059 s$i.sh;done
```
