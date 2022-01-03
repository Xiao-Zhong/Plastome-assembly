#!/usr/bin/perl -w
use strict;

my $infile=shift;
my $read_length=shift;

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

	my $seq2=substr($seq, 0, length($seq)-$read_length);
        
	#my $start=substr($seq, 0, 124);
	#print "$start\n";
	#$seq .=$start;
	print "$seq2\n";
}
      
