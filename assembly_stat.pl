#!/usr/bin/perl -w
use strict;

my $in_fasta_stat=shift;
my $in_denovo_path=shift;
my $in_pilon_out=shift;
my $sample=shift;
my $mode=shift;

$sample ||="test";
print "$sample\t";
$mode ||="denovo";
print "$mode\t";

open IN, "$in_fasta_stat" || die "$!\n";
<IN>;
$_=<IN>;
my @info=split(/\s+/, $_);
print "$info[1]\t$info[4]\t$info[2]\t$info[3]\t";
close IN;

open IN4, "$in_pilon_out" || die "$!\n";
while(<IN4>){
	chomp;
	print "$1\t" if (/bases \((\S+)%\)/);
	print "$1\t" if (/Mean total coverage: (\d+)/);
}
close IN4;

open IN2, "$in_denovo_path" ||die "$!\n";
my @info2=split(/\s+/, <IN2>);
#print "$info2[0]\n";
close IN2;

open IN3, "$info2[0]/refer_contigs_denovo.stat" || die "$!\n";
$_=<IN3>;
my @info4=split(/\s+/, $_);
print "$info4[0]\t$info4[1]\t$info4[5]\t$info4[13]\n";
close IN3;
