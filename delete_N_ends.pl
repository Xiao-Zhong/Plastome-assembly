#!/usr/bin/perl -w
use strict;

#die "Usage: $0 <input_seq_file> <output_seq_file> \n" if @ARGV<2;
my $infile=shift;
my $linelen=shift;

$linelen ||=80;

#open IN,$ARGV[0] || die "$!\n";
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
	$seq=~s/^N+//g;
	$seq=~s/N+$//g;
        $/="\n";
	$genome{$scaf}=$seq;
}
close IN;

#open OUT,">$ARGV[1]"|| die "$!\n";
foreach my $name(keys %genome){
	my $seq=$genome{$name};
	print ">$name\n";
        my $len=length($seq);
	for(my $i=0; $i<$len; $i +=$linelen) {
		my $seq2=substr($seq,$i,$linelen);#."\n";
		print "$seq2\n";
	}
}
#close OUT;
      
