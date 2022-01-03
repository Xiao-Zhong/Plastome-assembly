#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
#use FinBin qw($Bin);
my $Bin="/home/xzhong/working_data_01/01_Assembly/37_pipeline";

die "Usage:\n	perl $0 -refer <input_reference> -query <input_query> -step <1_or_12> -coverlen <e.g. 500> -coverrate <e.g. 5> -outdir ./\n" if @ARGV<2;

my ($refer, $query, $step, $coverlen, $coverrate, $outdir);
GetOptions(
	"refer:s"=>\$refer,
	"query:s"=>\$query,
	"step:s"=>\$step,
	"coverlen:i"=>\$coverlen,
	"coverrate:i"=>\$coverrate,
	"outdir:s"=>\$outdir
);

my $rfname=basename($refer);
my $qfname=basename($query);
#print "$referfilename\n$queryfilename\n";
$step ||='1';
$coverlen ||=0;
$coverrate ||=0;
$outdir ||=".";

if ($step=~/1/){
	`$Bin/nucmer -p $outdir/$rfname\_vs_$qfname $refer $query 1>$outdir/$rfname\_vs_$qfname\_mummer.log 2>$outdir/$rfname\_vs_$qfname\_mummer.err`;
	`$Bin/show-coords -rcl $outdir/$rfname\_vs_$qfname.delta > $outdir/$rfname\_vs_$qfname.coords`;
	`less $outdir/$rfname\_vs_$qfname.coords |grep "|" |grep "% IDY" -v |awk '\$8>=$coverlen && \$16>=$coverrate' > $outdir/$rfname\_vs_$qfname.coords.cp`;
}

if ($step=~/2/){
	unless(-z "$outdir/$rfname\_vs_$qfname.coords.cp"){
	`perl $Bin/fishInWinter.pl -bc 19 -fc 19 $outdir/$rfname\_vs_$qfname.coords.cp $outdir/$rfname\_vs_$qfname.coords > $outdir/$rfname\_vs_$qfname.coords.cp.all.list`;
	`less $outdir/$rfname\_vs_$qfname.coords.cp.all.list |awk '{print \$18"\t"\$1"\t"\$2"\t"\$19"\t"\$4"\t"\$5"\t"\$12"\t"\$13}' |sed 's/|/_/g' > $outdir/$rfname\_vs_$qfname.coords.cp.all.list2`;

	my %chash;
	open IN, "$outdir/$rfname\_vs_$qfname.coords.cp.all.list2" || die "$!\n";
	while(<IN>){
		chomp;
		my @info=split;
		#print "$info[0]\t$info[3]\t$info[6]\t$info[7]\n";
		if (! exists $chash{$info[3]}){
			`perl $Bin/scaffold_pair_svg.pl $outdir/$rfname\_vs_$qfname.coords.cp.all.list2 $outdir/$rfname\_vs_$info[3].svg --chr_1 "$info[0]" --chr_2 "$info[3]" --start_1 1 --end_1 $info[6] --start_2 1 --end_2 $info[7]`;
			$chash{$info[3]}++;
		}		
	}
	close IN;
	}
}

