use PerlIO::gzip;
while(<DATA>){chomp;
$h{$f}=$_ if /^\d+/;
$f=$_
}
foreach $f(keys %h){
$s='';$i=1;
open(F, "<:gzip", $f.'barcodes.tsv.gz');
while(<F>){chomp;
s/_/ /;
$hi{$i++}=$_
}
close F;
open(F, "<:gzip", $f.'matrix.mtx.gz');
while(<F>){
@t=split(/ /);
$s.="$hi{$t[1]}\t$t[2]" if defined $h{$f} and $h{$f} eq $t[0]
}
close F;
open(O, '>:gzip', $f.'matrix1.mtx.gz');
print O $s;
close O;
}
__DATA__
./FP200000371TR_D3/outs/filtered_feature_bc_matrix/matrix.mtx.gz
20500
./FP200000525TR_C5/outs/filtered_feature_bc_matrix/matrix.mtx.gz
20924
./FP200000492TR_A5/outs/filtered_feature_bc_matrix/matrix.mtx.gz
17945
./DP8400016859TR_A/outs/filtered_feature_bc_matrix/matrix.mtx.gz
21847
./FP200000492TR_A1/outs/filtered_feature_bc_matrix/matrix.mtx.gz
18872
./FP200000314BL_E3/outs/filtered_feature_bc_matrix/matrix.mtx.gz
20040
./FP200000492TR_B2/outs/filtered_feature_bc_matrix/matrix.mtx.gz
17547
./FP200000525TR_D6/outs/filtered_feature_bc_matrix/matrix.mtx.gz
17320
./FP200000525TR_A4/outs/filtered_feature_bc_matrix/matrix.mtx.gz
19947
./DP8400016859TR_C/outs/filtered_feature_bc_matrix/matrix.mtx.gz
20289
./FP200000492TR_B5/outs/filtered_feature_bc_matrix/matrix.mtx.gz
17757
./FP200000314BL_C4/outs/filtered_feature_bc_matrix/matrix.mtx.gz
19184
./FP200000525TR_E2/outs/filtered_feature_bc_matrix/matrix.mtx.gz
19689
./FP200000525TR_D1/outs/filtered_feature_bc_matrix/matrix.mtx.gz
21682
./FP200000492TR_B4/outs/filtered_feature_bc_matrix/matrix.mtx.gz
17595
./FP200000492TR_B3/outs/filtered_feature_bc_matrix/matrix.mtx.gz
18347
./FP200000314BL_E6/outs/filtered_feature_bc_matrix/matrix.mtx.gz
20173
./FP200000314BL_D6/outs/filtered_feature_bc_matrix/matrix.mtx.gz
19226
./FP200000525TR_D3/outs/filtered_feature_bc_matrix/matrix.mtx.gz
21833
./DP8400016859TR_A5/outs/filtered_feature_bc_matrix/matrix.mtx.gz
21847
./FP200000525TR_B5/outs/filtered_feature_bc_matrix/matrix.mtx.gz
19187
./FP200000492TR_B1/outs/filtered_feature_bc_matrix/matrix.mtx.gz
17547
./FP200000525TR_A6/outs/filtered_feature_bc_matrix/matrix.mtx.gz
18917
./FP200000314BL_D4/outs/filtered_feature_bc_matrix/matrix.mtx.gz
19885
./FP200000525TR_D5/outs/filtered_feature_bc_matrix/matrix.mtx.gz
17198
./FP200000525TR_E1/outs/filtered_feature_bc_matrix/matrix.mtx.gz
18494
./FP200000492TR_A2/outs/filtered_feature_bc_matrix/matrix.mtx.gz
19317
./FP200000492TR_A3/outs/filtered_feature_bc_matrix/matrix.mtx.gz
17857
./FP200000525TR_D4/outs/filtered_feature_bc_matrix/matrix.mtx.gz
21019
./DP8400016859TR_A3/outs/filtered_feature_bc_matrix/matrix.mtx.gz
19823
./FP200000525TR_E3/outs/filtered_feature_bc_matrix/matrix.mtx.gz
16784
./FP200000525TR_A1/outs/filtered_feature_bc_matrix/matrix.mtx.gz
18642
./FP200000525TR_B6/outs/filtered_feature_bc_matrix/matrix.mtx.gz
18129
./FP200000525TR_C6/outs/filtered_feature_bc_matrix/matrix.mtx.gz
19719
./FP200000314BL_D3/outs/filtered_feature_bc_matrix/matrix.mtx.gz
19946
./FP200000492TR_A4/outs/filtered_feature_bc_matrix/matrix.mtx.gz
18314
./FP200000492TR_A6/outs/filtered_feature_bc_matrix/matrix.mtx.gz
18434
./FP200000525TR_E4/outs/filtered_feature_bc_matrix/matrix.mtx.gz
16800
./FP200000525TR_D2/outs/filtered_feature_bc_matrix/matrix.mtx.gz
21389
./FP200000314BL_D5/outs/filtered_feature_bc_matrix/matrix.mtx.gz
19903
./FP200000525TR_A5/outs/filtered_feature_bc_matrix/matrix.mtx.gz
19192
./FP200000525TR_A2/outs/filtered_feature_bc_matrix/matrix.mtx.gz
17032
