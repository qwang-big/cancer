use PerlIO::gzip;
while(<DATA>){chomp;
$h{$1}=$_ if /\/(\S+)\/outs/
}
foreach $f(keys %h){
$s='';$i=0;%p=();%hb=();
open(F, '<', $h{$f}.'matrix1.mtx')||die"$!";
<F>;<F>;
while(<F>){chomp;
s/ /_/;
$p{$_}=1
}
close F;
$i=0;$j=0;
open(F, "<:gzip", "barcodes/$f.tsv.gz")||die"$!";
open(O, ">:gzip", $h{$f}.'barcodes.tsv.gz')||die"$!";
while(<F>){chomp;
$i++;
if (defined $p{$_}){
$j++;
$hb{$i}=$j;
print O $_;
}
}
close F;
close O;
$i=0;$j=0;$k=0;
open(F, "<:gzip", $f.'.mtx.gz')||die"$!";
<F>;<F>;
while(<F>){
@t=split(/ /);
$i=$t[0] if $i<$t[0];
$j=$t[1] if $j<$t[1];
if (defined $hb{$t[1]}){
$s.="$t[0] ".$hb{$t[1]}." $t[2]";
$k++
}
}
close F;
open(O, ">:gzip", $h{$f}.'/matrix.mtx.gz');
print O '%%MatrixMarket matrix coordinate integer general
%metadata_json: {"software_version": "cellranger-4.0.0", "format_version": 2}
',"$i $j $k\n",$s;
close O;
}
__DATA__
./FP200000371TR_D3/outs/filtered_feature_bc_matrix/
./FP200000525TR_C5/outs/filtered_feature_bc_matrix/
./FP200000492TR_A5/outs/filtered_feature_bc_matrix/
./FP200000314BL_E3/outs/filtered_feature_bc_matrix/
./FP200000492TR_B2/outs/filtered_feature_bc_matrix/
./FP200000525TR_D6/outs/filtered_feature_bc_matrix/
./FP200000525TR_A4/outs/filtered_feature_bc_matrix/
./DP8400016859TR_C/outs/filtered_feature_bc_matrix/
./FP200000492TR_B5/outs/filtered_feature_bc_matrix/
./FP200000314BL_C4/outs/filtered_feature_bc_matrix/
./FP200000525TR_E2/outs/filtered_feature_bc_matrix/
./FP200000525TR_D1/outs/filtered_feature_bc_matrix/
./FP200000492TR_B4/outs/filtered_feature_bc_matrix/
./FP200000492TR_B3/outs/filtered_feature_bc_matrix/
./FP200000314BL_E6/outs/filtered_feature_bc_matrix/
./FP200000314BL_D6/outs/filtered_feature_bc_matrix/
./FP200000525TR_D3/outs/filtered_feature_bc_matrix/
./DP8400016859TR_A5/outs/filtered_feature_bc_matrix/
./FP200000525TR_B5/outs/filtered_feature_bc_matrix/
./FP200000492TR_B1/outs/filtered_feature_bc_matrix/
./FP200000525TR_A6/outs/filtered_feature_bc_matrix/
./FP200000314BL_D4/outs/filtered_feature_bc_matrix/
./FP200000525TR_D5/outs/filtered_feature_bc_matrix/
./FP200000525TR_E1/outs/filtered_feature_bc_matrix/
./FP200000492TR_A2/outs/filtered_feature_bc_matrix/
./FP200000492TR_A3/outs/filtered_feature_bc_matrix/
./FP200000525TR_D4/outs/filtered_feature_bc_matrix/
./DP8400016859TR_A3/outs/filtered_feature_bc_matrix/
./FP200000525TR_E3/outs/filtered_feature_bc_matrix/
./FP200000525TR_A1/outs/filtered_feature_bc_matrix/
./FP200000525TR_B6/outs/filtered_feature_bc_matrix/
./FP200000525TR_C6/outs/filtered_feature_bc_matrix/
./FP200000314BL_D3/outs/filtered_feature_bc_matrix/
./FP200000492TR_A4/outs/filtered_feature_bc_matrix/
./FP200000492TR_A6/outs/filtered_feature_bc_matrix/
./FP200000525TR_E4/outs/filtered_feature_bc_matrix/
./FP200000525TR_D2/outs/filtered_feature_bc_matrix/
./FP200000314BL_D5/outs/filtered_feature_bc_matrix/
./FP200000525TR_A5/outs/filtered_feature_bc_matrix/
./FP200000525TR_A2/outs/filtered_feature_bc_matrix/
