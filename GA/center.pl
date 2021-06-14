$N=$1 if $ARGV[0]=~/(\d+)$/;
$n=1;
@rr=(3..15);
open(F,"<$ARGV[0]");
while(<F>){
chomp;
@t=split(/\t/);
foreach $r(@rr){
for ($x=$t[0]-$r;$x<$t[0]+$r;$x++) {
for ($y=$t[1]-$r;$y<$t[1]+$r;$y++) {
$h{"$x\t$y\t$r"}="$n\t$r";
}
}
}
$n++
}
close F;
open(F, '-|', "gzip -dc $ARGV[1]");
while(<F>){
@t=split(/\t/);
$hg{$t[2]}=1;
foreach $r(@rr){
$h1{$t[2]}->{$h{"$t[0]\t$t[1]\t$r"}}+=$t[3] if defined $h{"$t[0]\t$t[1]\t$r"}
}
}
close F;
foreach $r(@rr){
open(O,">$ARGV[0].$r");
foreach $gene(keys %hg){
print O $gene;
foreach $cell(1..$N){
print O "\t",defined $h1{$gene}->{"$cell\t$r"} ? $h1{$gene}->{"$cell\t$r"} : 0
}
print O "\n"
}
close O
}
