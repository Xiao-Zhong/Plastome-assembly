#!/usr/bin/perl -w
use strict;
use Getopt::Long;

die "Usage:\n  perl $0 -ref <reference.fa> -pelist <PE_reads_list> -selist <SE_reads_list> -sample <sample> -outdir ./ \n" if @ARGV<2;

my ($inref, $inpe, $inse, $sample, $outdir);
GetOptions(
	"ref:s"=>\$inref,
	"pelist:s"=>\$inpe,
	"selist:s"=>\$inse,
	"sample:s"=>\$sample,
	"outdir:s"=>\$outdir
);

$sample ||="test";
$outdir ||=".";

open OUT, ">$outdir/$sample\_bwa.sh" || die "$!\n";

print OUT "# index reference sequences\n";
if ($inref){
	print OUT "bwa index $inref\n\n";
}

print OUT "# Pair-end reads mapping by BWA\n";
if ($inpe){
	my @inpe=split(/,/, $inpe);
	#for(@inpe){print "$_\n";}
	for (my $i=0; $i<@inpe; $i +=2){
		#print "$inpe[$i]\t$inpe[$i+1]\n";
		#print OUT "bwa mem $inref $inpe[$i] $inpe[$i+1] > $outdir/$sample.pe$i.sam; samtools view -S -b -o $outdir/$sample.pe$i.bam $outdir/$sample.pe$i.sam; samtools sort $outdir/$sample.pe$i.bam $outdir/$sample.pe$i.sort; samtools index $outdir/$sample.pe$i.sort.bam\n";
		print OUT "bwa mem $inref $inpe[$i] $inpe[$i+1] > $outdir/$sample.pe$i.sam";
		#print OUT "rm $outdir/$sample.pe$i.sam $outdir/$sample.pe$i.bam\n";
	}
	print OUT "\n";
}

print OUT "# Single reads mapping by BWA\n";
if ($inse){
	my @inse=split(/,/, $inse);
	for (my $i=0; $i<@inse; $i++){
		#print OUT "bwa mem $inref $inse[$i] > $outdir/$sample.se$i.sam; samtools view -S -b -o $outdir/$sample.se$i.bam $outdir/$sample.se$i.sam; samtools sort $outdir/$sample.se$i.bam $outdir/$sample.se$i.sort; samtools index $outdir/$sample.se$i.sort.bam\n";
		print OUT "bwa mem $inref $inse[$i] > $outdir/$sample.se$i.sam";
		#print OUT "rm $outdir/$sample.se$i.sam $outdir/$sample.se$i.bam\n";
	}
	print OUT "\n";
}
close OUT;

#`sh $sample\_bwa.sh 1>$sample\_bwa.sh.log 2>$sample\_bwa.sh.err`;


