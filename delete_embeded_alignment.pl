#!/usr/bin/perl -w
use strict;

my $inalign=shift;
my $bestrefer=shift;
my $bestcontig=shift;

my %blocks;
open IN, "$inalign" || die "$!\n";
while(<IN>){
	chomp;
	my @info=split;
	my @info2=split(/\|/, $info[18]);
	next unless (($info[17] eq "$bestrefer") && ($info2[-1] eq "$bestcontig"));
	#print "$_\n";
	my $strand=($info[3]<$info[4]) ? "+" : "-";
	push @{$blocks{"$info[18]\t$strand"}}, [@info];
}
close IN;

open IN2, "$inalign" || die "$!\n";
lable:while(<IN2>){
	chomp;
	my @info=split;
	my @info2=split(/\|/, $info[18]);
	next unless (($info[17] eq "$bestrefer") && ($info2[-1] eq "$bestcontig"));
	my $s=($info[3]<$info[4]) ? "+" : "-";

	my $bid="$info[18]\t$s";
	#print "$bid\n";
	for (my $i=0; $i<@{$blocks{$bid}}; $i++){
		#print "$info[0]\_$info[1]\t$blocks{$bid}[$i][0]\t$blocks{$bid}[$i][1]\n";
		my $olds=$blocks{$bid}[$i][0];
		my $olde=$blocks{$bid}[$i][1];
		if (($info[0]>=$olds && $info[1]<$olde) || ($info[0]>$olds && $info[1]<=$olde)){
			next lable;
		}
	}
	print "$_\n";
}
close IN2;

