while(<>){
if(/CB:Z:(\S+)/){
$cb=$1;
if(/GE:Z:(\S+)/){
$ge=$1;
if(/XF:Z:EXONIC/){
print "$cb\t$ge\n"
}elsif(/XF:Z:INTRONIC/){
print STDERR "$cb\t$ge\n"
}}}}
