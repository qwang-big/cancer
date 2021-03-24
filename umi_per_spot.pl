$bin=50;
while(<>){chomp;
@t=split(/\t/);
$x=int($t[1]/$bin + 0.5);
$y=int($t[2]/$bin + 0.5);
$h{"$x\t$y"}+=$t[3];
$h1{"$x\t$y"}->{$t[0]}=1
}
foreach(keys %h){
$size = keys %{$h1{$_}};
print $size,"\n";
print STDERR $h{$_},"\n"
}
