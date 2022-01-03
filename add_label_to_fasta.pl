#!/usr/bin/perl -w
use strict;

#die "Usage: $0 <input_seq_file> <output_seq_file> \n" if @ARGV<2;
my $infile=shift;
my $inlable=shift;
my $linelen=shift;

$linelen ||=60;

#open IN,$ARGV[0] || die "$!\n";
open IN, "$infile" || die "$!\n";
my ($scaf,$seq,%genome);
$/=">";<IN>;$/="\n";
while(<IN>){
	$scaf=$1 if (/^(\S+)/);
	$/=">";
	$seq=<IN>;
	chomp($seq);
	$seq=~s/\s+//g;
	$seq=~tr/atcg/ATCG/;
        $/="\n";
	$genome{$scaf}=$seq;
}
close IN;

#open OUT,">$ARGV[1]"|| die "$!\n";
foreach my $name(keys %genome){
	$seq=$genome{$name};
	print ">$name\|$inlable\n";
        my $len=length($seq);
	for(my $i=0; $i<$len; $i +=$linelen) {
		my $seq2=substr($seq,$i,$linelen);#."\n";
		print "$seq2\n";
	}
}
#close OUT;
      
