#!/usr/bin/perl -w
use strict;

my $infile=shift;
my $read_length=shift;

$read_length ||=124;

open IN, "$infile" || die "$!\n";
my %genome;
$/=">";<IN>;$/="\n";
while(<IN>){
	my $scaf=$1 if (/^(\S+)/);
	$/=">";
	my $seq=<IN>;
	chomp($seq);
	$seq=~s/\s+//g;
	$seq=~tr/atcg/ATCG/;
        $/="\n";
	$genome{$scaf}=$seq;
}
close IN;

foreach my $name(keys %genome){
	my $seq=$genome{$name};
	print ">$name\n";
        
	my $start=substr($seq, 0, $read_length);
	#print "$start\n";
	$seq .=$start;
	print "$seq\n";
}
      
