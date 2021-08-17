$b=20;
$f=shift;
$i=0;
$j=0;
$k=0;
$sk=0;
use PerlIO::gzip;
use File::Path qw(make_path);
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
$f=substr($f,0,16);
make_path("$f/outs/filtered_feature_bc_matrix");
open(O, '>:gzip', "$f/outs/filtered_feature_bc_matrix/features.tsv.gz");
foreach (sort {$hg{$a} <=> $hg{$b}} keys %hg){
print O "$_\t$_\tGene Expression\n"
}
close O;
open(O, '>:gzip', "$f/outs/filtered_feature_bc_matrix/barcodes.tsv.gz");
foreach (sort {$hc{$a} <=> $hc{$b}} keys %hc){
print O $_,"\n"
}
close O;
$s='';
$ss='';
foreach $cell(sort {$hc{$a} <=> $hc{$b}} keys %hc){
$i=$hc{$cell};
foreach $gene(sort {$hg{$a} <=> $hg{$b}} keys %hg){
$j=$hg{$gene};
$n=0;
$n+= $h1{$cell}->{$gene} if defined $h1{$cell}->{$gene};
$n+= $h2{$cell}->{$gene} if defined $h2{$cell}->{$gene};
next if $n==0;
$k++;
$s .= "$j $i $n\n";
$ss.= "$j $i $n\n" if defined $h2{$cell}->{$gene};
$sk++ if defined $h2{$cell}->{$gene};
}
}
close O1;
close O2;
open(O, '>:gzip', "$f/outs/filtered_feature_bc_matrix/matrix.mtx.gz");
print O '%%MatrixMarket matrix coordinate integer general
%metadata_json: {"software_version": "cellranger-4.0.0", "format_version": 2}
',"$j $i $k\n",$s;
close O;
open(O, '>:gzip', "$f/outs/filtered_feature_bc_matrix/spliced.mtx.gz");
print O '%%MatrixMarket matrix coordinate integer general
%metadata_json: {"software_version": "cellranger-4.0.0", "format_version": 2}
',"$j $i $sk\n",$ss;
close O;
