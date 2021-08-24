use PerlIO::gzip;
while(<DATA>){chomp;
$h{$1}=$_ if /\/(\S+)\/outs/
}
foreach $f(keys %h){
$s='';$i=0;$j=0;$k=0;%p=();
open(F, '<', $h{$f}.'matrix1.mtx');
open(O, ">:gzip", $f.'/barcodes.tsv.gz');
<F>;<F>;
while(<F>){chomp;
@t=split(/ /);
$i=$t[0] if $i<$t[0];
$j=$t[1] if $j<$t[1];
$p{$_}=1;
s/ /_/;
print O $_,"\n"
}
close F;
close O;
open(F, "<:gzip", $h{$f}.'matrix.mtx.gz');
<F>;<F>;
while(<F>){
@t=split(/ /);
if (defined $p{"$t[0] $t[1]"}){
$k++;
$s.=$_
}
}
close F;
open(O, ">:gzip", $f.'/matrix.mtx.gz');
print O '%%MatrixMarket matrix coordinate integer general
%metadata_json: {"software_version": "cellranger-4.0.0", "format_version": 2}
',"$i $j $k\n",$s;
close O;
system("cp $h{$f}features.tsv.gz $f/features.tsv.gz")
}
__DATA__
./FP200000371TR_D3/outs/filtered_feature_bc_matrix/
./FP200000525TR_C5/outs/filtered_feature_bc_matrix/
