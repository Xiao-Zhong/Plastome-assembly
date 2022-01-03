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

my $cover_highest=0;
my %listhash;my %blockid;my @order;
open IN3, "$inalignment" || die "$!\n";
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
		if ($info[18]=~/cov\_(\S+?)\|/){
			$cover_highest=$1 if ($1>$cover_highest);
		}
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
	#print STDERR  ">$cid\t$s\n";
	
	my $num=@{$listhash{$cid}{$s}};
	#print OUT "$num\t";
	
	my $seq;my $start; my $end;
	if(!exists $cidhash{$cid}){
		print OUT ">$cid\t$s\t$num\t1st_orientation\t";
		#$seq=$chash{$cid};
		$cidhash{$cid}++;

		if(@order>=2){
			my $next_cid; my $next_s;
			if($i<@order-1){
				#print STDERR "$i\t$order[$i]\t$order[$i+1]\n";
				if($order[$i+1]=~/(\S+)\t(\S+)/){
					$next_cid=$1;
					$next_s=$2;
				}
				#print STDERR "$cid\t$s\n";
				#print STDERR "$next_cid\t$next_s\n\n";
			
				if ($s ne $next_s){
					if ($listhash{$cid}{$s}[-1][1] >= $listhash{$next_cid}{$next_s}[0][0] && $listhash{$cid}{$s}[-1][1] <= $listhash{$next_cid}{$next_s}[0][1]){
						if ($s eq "+"){
							$start=1;
							$end=$listhash{$cid}{$s}[-1][4];
						}else{
							$start=$listhash{$cid}{$s}[-1][4];
							$end=$listhash{$cid}{$s}[-1][12];
						}
						#print "$next_cid\t$next_s\t$start\t$end\n";
					}else{
						$start=1;
						$end=$listhash{$cid}{$s}[-1][12];
					}
				}else{
					$start=1;
					$end=$listhash{$cid}{$s}[-1][12];
				}
			}else{
				$start=1;
				$end=$listhash{$cid}{$s}[-1][12];
			}
		}else{
			$start=1;
			$end=$listhash{$cid}{$s}[-1][12];
		}

		$seq=substr($chash{$cid}, $start-1, $end-$start+1);	
	}else{
		#my $tmp=$1 if ($cid=~/cov\_(\S+?)\|/);
		#if ($tmp==$cover_highest){
			print OUT ">$cid\t$s\t$num\t2nd_orientation\t";

			#print "$bnum\t$i\n";
			#if (@order >=4){
				#my $a=$1 if ($order[-4]=~/cov\_(\S+?)\|/);
				#my $b=$1 if ($order[-3]=~/cov\_(\S+?)\|/);
				#my $c=$1 if ($order[-2]=~/cov\_(\S+?)\|/);
				#print "$a\t$b\t$c\t$listhash{$cid}{$s}[0][12]\n";
				if($listhash{$cid}{$s}[0][12]<30000){
					#$seq=$chash{$cid};
					$start=1;
					$end=$listhash{$cid}{$s}[-1][12];
				}else{
					if($s eq '+'){
						$start=$listhash{$cid}{$s}[0][3];
						$end=$listhash{$cid}{$s}[-1][4];
					}else{
						$start=$listhash{$cid}{$s}[-1][4];
						$end=$listhash{$cid}{$s}[0][3];
					}
					#$seq=substr($chash{$cid}, $start-1, $end-$start+1);
				}
				$seq=substr($chash{$cid}, $start-1, $end-$start+1);
			#}else{
			#	if($s eq '+'){
			#		$start=$listhash{$cid}{$s}[0][3];
			#		$end=$listhash{$cid}{$s}[-1][4];
			#	}else{
			#		$start=$listhash{$cid}{$s}[-1][4];
			#		$end=$listhash{$cid}{$s}[0][3];
			#	}
			#	$seq=substr($chash{$cid}, $start-1, $end-$start+1);
			#}
		#}else{
		#	next;
		#}
	}
	
	if ($s eq '-'){
		$seq=reverse $seq;
		$seq=~tr/AGCTagct/TCGAtcga/;
	}
	
	print OUT "$start\t$end\t";

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




