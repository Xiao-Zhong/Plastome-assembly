#!/usr/bin/perl -w
use strict;

my $inbbnormlog=shift;
my $incutadapterlog=shift;
my $sample=shift;

print "$sample\t";

my ($rawreads, $rawbase);
open IN, "$inbbnormlog" ||die "$!\n";
lable:while(<IN>){
	chomp;
	if(/Total reads in:\s+(\d+)/){
		#print "$1\t";
		$rawreads=$1;
		my $nextline=<IN>;
		if($nextline=~/Total bases in:\s+(\d+)/){
			#print "$1\t";
			$rawbase=$1;
			last lable;
		}
	}
}
close IN;

my ($inputreads,$inputbase,$inputbase1,$inputbase2,$outputreads,$outputbase,$outputbase1,$outputbase2,$remain_rate);
open IN2, "$incutadapterlog" || die "$!\n";
while(<IN2>){
	chomp;
	if(/Total read pairs processed:\s+(\S+)/){
		#print "$1\n";
		$inputreads=$1;
	}
	
	if (/Pairs written \(passing filters\):\s+(\S+)/){
		#print "$1\n";
		$outputreads=$1;
	}

	if (/Total basepairs processed:\s+(\S+)/){
		#print "$1\n";
		$inputbase=$1;
		my $next1stline=<IN2>;
		if ($next1stline=~/Read 1:\s+(\S+)/){
			#print "$1\n";
			$inputbase1=$1;
		}
		my $next2ndline=<IN2>;
		if ($next2ndline=~/Read 2:\s+(\S+)/){
			#print "$1\n";
			$inputbase2=$1;
		}
	}

	if (/Total written \(filtered\):\s+(\S+) bp \((\S+)%\)/){
		#print "$1\t$2\n";
		$outputbase=$1;
		$remain_rate=$2;
		my $next1stline=<IN2>;
		if ($next1stline=~/Read 1:\s+(\S+)/){
			#print "$1\n";
			$outputbase1=$1;
		}
		my $next2ndline=<IN2>;
		if ($next2ndline=~/Read 2:\s+(\S+)/){
			#print "$1\n";
			$outputbase2=$1;
		}	
	}
}
close IN2;

$inputreads=~s/,//g;$inputbase=~s/,//g;$inputbase1=~s/,//g;$inputbase2=~s/,//g;$outputreads=~s/,//g;$outputbase=~s/,//g;$outputbase1=~s/,//g;$outputbase2=~s/,//g;
print "$rawreads\t$rawbase\t$inputreads\t$inputbase\t$inputbase1\t$inputbase2";
printf "\t%.1f", $inputbase/$rawbase*100;
print "\t$outputreads\t$outputbase\t$outputbase1\t$outputbase2\t$remain_rate\n";


