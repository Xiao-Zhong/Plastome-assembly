#!/usr/bin/perl -w
use strict;

die "Usage:\n  perl $0 <input_contigs> <input_alignment> <refer_id> <sample> <outdir>\n" if @ARGV<2;

my $incontigs=shift;
my $inalignment=shift;
my $referid=shift;
my $sample=shift;
my $outdir=shift;

$referid ||="refer";
$sample ||="sample";
$outdir ||=".";

my %chash;
open IN2, "$incontigs" || die "$!\n";
$/=">";<IN2>;$/="\n";
while(<IN2>){
	my $id=$1 if(/^(\S+)/);
	#print "$id\n";
	$/=">";
	my $seq=<IN2>;
	chomp($seq);
	$seq=~s/\s+//g;
	$seq=~tr/atcg/ATCG/;
	#print "$seq\n";
	$/="\n";
	$chash{$id}=$seq;	
}
close IN2;

my %listhash;my %blockid;my @order;
open IN3, "$inalignment" || die "$!\n";
#<IN3>;<IN3>;<IN3>;<IN3>;<IN3>;
while(<IN3>) {
	chomp;
	#print "$_\n";
	my @info=split;
	next unless ($info[7]>=800 && $info[15]>=6);
	my $strand=($info[3]<$info[4]) ? '+':'-';
	
	if (!exists $blockid{"$info[18]\t$strand"}){
		#print "$info[0]\t$info[18]\t$strand\n";
		push @order, "$info[18]\t$strand";
		$blockid{"$info[18]\t$strand"}++;
	}
	
	push @{$listhash{$info[18]}{$strand}}, [@info];
}
close IN3;

open OUT, ">$outdir/$sample\_linked_cp.log" || die "$!\n";
my %cidhash;my $finalseq;my $conjuction_p;
for(my $i=0; $i<@order; $i++){
	#print "$_\n";
	my $cid; my $s;
	if($order[$i]=~/(\S+)\t(\S+)/){
		$cid=$1;
		$s=$2;
	}
	print OUT ">$cid\t$s\t";
	
	my $num=@{$listhash{$cid}{$s}};
	print OUT "$num\t";
	
	my $seq;
	if(!exists $cidhash{$cid}){
		print OUT "1st_orientation\t";
		$seq=$chash{$cid};
		$cidhash{$cid}++;			
	}else{
		print OUT "2nd_orientation\t";
		#@{$listhash{$cid}{$s}}=sort {$a->[3] <=> $b->[4]} @{$listhash{$cid}{$s}};
		#for (my $i=0; $i<$num; $i++){
		#	print "$listhash{$cid}{$s}[$i][3]\t$listhash{$cid}{$s}[$i][4]\t";
		#}
		#print "\n";
		my $start; my $end;
		if($s eq '+'){
			$start=$listhash{$cid}{$s}[0][3];
			$end=$listhash{$cid}{$s}[-1][4];
		}else{
			$start=$listhash{$cid}{$s}[-1][4];
			$end=$listhash{$cid}{$s}[0][3];
		}
		$seq=substr($chash{$cid}, $start-1, $end-$start+1);
	}
	
	if ($s eq '-'){
		$seq=reverse $seq;
		$seq=~tr/AGCTagct/TCGAtcga/;
	}
	
	#print "$i\n";
	$seq ="nnnnnnnnnn"."$seq" if ($i !=0);
	my $seqlen=length($seq);	
	print OUT "$seqlen\n$seq\n";

	$finalseq .=$seq;
	$conjuction_p .=length($finalseq)."\_";
}
close OUT;

$conjuction_p=~s/\_\d+\_$//g;
open OUT2, ">$outdir/$sample\_linked_cp.fa" || die "$!\n";
print OUT2 ">$sample\_guided_by_$referid\t$conjuction_p\n$finalseq\n";
close OUT2;















