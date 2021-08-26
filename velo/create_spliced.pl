use PerlIO::gzip;
while(<DATA>){chomp;
$h{$1}=$_ if /\/(\S+)\/outs/
}
foreach $f(keys %h){
$s='';$i=0;%hf=();%hb=();%hf=();
open(F, "<:gzip", $h{$f}.'barcodes.tsv.gz');
while(<F>){chomp;
$i++;
$hb{$_}=$i;
}
close F;
$i=0;
open(F, "<:gzip", $h{$f}.'features.tsv.gz')||die"$!";
while(<F>){chomp;
$i++;
@t=split(/ /);
$hf{$t[0]}=$i;
}
close F;
open(F, "<", $f.'.cnt')||die"$! $f";
while(<F>){
@t=split(/ /);
@d=split(/_/,$t[0]);
$x=int($d[0]/50);
$y=int($d[1]/50);
if (defined $hb{"$x\_$y"}){
$s.="$hf{$t[1]} ".$hb{"$x\_$y"}." $t[2]";
$i = $hf{$t[1]} if $i<$hf{$t[1]};
$j = $hb{"$x\_$y"} if $j<$hb{"$x\_$y"};
$k++
}
}
close F;
open(O, ">:gzip", $h{$f}.'/spliced.mtx.gz');
print O '%%MatrixMarket matrix coordinate integer general
%metadata_json: {"software_version": "cellranger-4.0.0", "format_version": 2}
',"$i $j $k\n",$s;
close O;
}
__DATA__
