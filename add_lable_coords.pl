#!/usr/bin/perl -w
use strict;

my $inalign=shift;
my $lable=shift;

open IN, "$inalign" || die "$!\n";
while(<IN>){
	chomp;
	if (/NODE\_*/){
		print "$_\|$lable\n";
	}else{	
		print "$_\n";
		}
}
close IN;

