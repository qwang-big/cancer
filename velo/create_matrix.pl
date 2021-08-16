$b=20;
$f=shift;
$i=0;
$j=0;
$k=0;
open(F,"<$f.in");
while(<F>){chomp;
@t=split(/\t/);
@d=split(/_/,$t[0]);
$x=int($d[0]/$b);
$y=int($d[1]/$b);
$hc{"$x_$y"}=++$i unless defined $hc{"$x_$y"};
$hg{$t[1]} = ++$j unless defined $hg{$t[1]};
$h1{"$x_$y"}->{$t[1]}+=1;
}
close F;
open(F,"<$f.ex");
while(<F>){chomp;
@t=split(/\t/);
@d=split(/_/,$t[0]);
$x=int($d[0]/$b);
$y=int($d[1]/$b);
$hc{"$x_$y"}=++$i unless defined $hc{"$x_$y"};
$hg{$t[1]} = ++$j unless defined $hg{$t[1]};
$h2{"$x_$y"}->{$t[1]}+=1;
}
close F;
open(O,">features.tsv");
foreach (sort {$hg{$a} <=> $hg{$b}} keys %hg){
print O "$_\t$_\tGene Expression\n"
}
close O;
open(O,">barcodes.tsv");
foreach (sort {$hc{$a} <=> $hc{$b}} keys %hc){
print O $_,"\n"
}
close O;
$s='';
open(O1,">spliced.mtx");
open(O2,">unspliced.mtx");
foreach $cell(sort {$hc{$a} <=> $hc{$b}} keys %hc){
$i=$hc{$cell};
foreach $gene(sort {$hg{$a} <=> $hg{$b}} keys %hg){
$j=$hg{$gene};
$n=0;
$n+= $h1{$cell}->{$gene} if defined $h1{$cell}->{$gene};
$n+= $h2{$cell}->{$gene} if defined $h2{$cell}->{$gene};
$k=$n if $k<$n;
$s.= "$j\t$i\t$n\n";
print O1 "$j\t$i\t",$h1{$cell}->{$gene},"\n" if defined $h1{$cell}->{$gene};
print O2 "$j\t$i\t",$h2{$cell}->{$gene},"\n" if defined $h2{$cell}->{$gene};
}
}
close O1;
close O2;
open(O,">matrix.mtx");
print O '%%MatrixMarket matrix coordinate integer general
%metadata_json: {"software_version": "cellranger-4.0.0", "format_version": 2}
',"$j\t$i\t$k\n",$s;
close O;
