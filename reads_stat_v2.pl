#!/usr/bin/perl -w
use strict;

my $incutadapterlog=shift;
my $inbbnormlog=shift;
my $sample=shift;

print "$sample\t";

my ($outputreads, $outputbase);
open IN, "$inbbnormlog" ||die "$!\n";
lable:while(<IN>){
	chomp;
	
	if (/Pass 2/){
		while(<IN>){
			if(/Total reads in:\s+(\d+)\s+(\S+)% Kept/){
				#print "$1\t$2\n";
				#$outputreads=$1;
				$outputreads=int($1*$2/100);


				my $nextline=<IN>;
				if($nextline=~/Total bases in:\s+(\d+)\s+(\S+)% Kept/){
					#print "$1\t$2\n";
					#$outputbase=$1;
					$outputbase=int($1*$2/100);
					last lable;
				}
			}
		}
	}
}
close IN;

my ($inputreads, $inputbase);
open IN2, "$incutadapterlog" || die "$!\n";
while(<IN2>){
	chomp;
	if(/Total read pairs processed:\s+(\S+)/){
		#print "$1\n";
		$inputreads=$1;
	}
	
	if (/Total basepairs processed:\s+(\S+)/){
		#print "$1\n";
		$inputbase=$1;
	}
}
close IN2;


$inputreads=~s/,//g;$inputbase=~s/,//g;
$inputreads *=2;
print "$inputreads\t$inputbase\t$outputreads\t$outputbase";
printf "\t%.2f\n", $outputbase/$inputbase*100;



