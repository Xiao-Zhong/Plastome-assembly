#!/usr/bin/perl -w
use strict;

my $id=shift;
my $in_pilon_report=shift;
my $in_samtools_stat=shift;

#print "ID	(raw_read	mapped	mapped_%)	Length(bp)	# of Reads	Mean Coverage	Minimum Coverage	Confirmed(%)"
print "$id\t";

if ($in_samtools_stat){
	open IN, "$in_samtools_stat" || die "$!\n";
	while(<IN>){
		print "$1\t" if (/(\d+).+total/);
		if (/(\d+).+mapped\s+\((\S+)\%/){
			#my $mapped=$1;
			#my $mapped_percent=$2;
			print "$1\t$2\t";
		}
	}
	close IN
}



open IN, "$in_pilon_report" || die "$!\n";
while(<IN>){
	chomp;
	if(/Input genome size: (\d+)/){
		print "$1\t";
	}

	if (/:\s+(\d+)\s+reads.+\s+(\d+)\s+mapped.+\s+(\d+)\+\//){
		print "$1\t$2\t";
		printf "%.2f\t", $2/$1*100;
		print "$3\t";
	}	

	if(/Total Reads: (\d+), Coverage: (\d+), minDepth: (\d+)/){
		print "$2\t$3\t";
	}
	
	if(/Confirmed.+bases \((\S+)\%\)/){
		print "$1\n";
	}
	
}
close IN;

