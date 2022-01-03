#!/usr/bin/perl -w
use strict;

die "Usage:\n   perl $0 <alignment_output: *.coords.cp.id> \n" if @ARGV<1;

my $inalign=shift;

my %refer_contigs;
open IN, "$inalign" || die "$!\n";
while(<IN>){
	chomp;
	my @info=split;
	my @info2=split(/\|/, $info[18]);
	#print "$info[17]\t$info2[1]\n";
	push @{$refer_contigs{"$info[17]\t$info2[1]"}}, [@info],
}
close IN;

my %align_stat; my $max_cover=0;
foreach my $rcid(keys %refer_contigs){
	#print "$rcid\t";
	my $refer_cover=0; my $total_c_len=0; my $total_c_aligned_len=0;
	my %chash; my $contig_len_pass=0;my $contiglist;
	for (my $i=0; $i<@{$refer_contigs{$rcid}}; $i++){
		$refer_cover +=$refer_contigs{$rcid}[$i][14];

		my @info=split(/\|/, $refer_contigs{$rcid}[$i][18]);
		unless (exists $chash{$info[0]}){
			$total_c_len +=$refer_contigs{$rcid}[$i][12];
			$contiglist .="$refer_contigs{$rcid}[$i][12],";	
			$chash{$info[0]}++;
		}
		$total_c_aligned_len +=$refer_contigs{$rcid}[$i][7];	
	}

	my $cnum=keys %chash;
	my $aligned_num=@{$refer_contigs{$rcid}};

	my $refer_len=$refer_contigs{$rcid}[0][11];
	#print "$rcid\t$refer_len\n";
	
	#$total_c_len+=23000;# two IR repeats
	#if ($total_c_len > 70000){
		my $contigs_cover=$total_c_aligned_len/$total_c_len*100;
		#print "$rcid\t$refer_cover\t$total_c_len\t$total_c_aligned_len\t$contigs_cover\t$cnum\t$aligned_num\n";
		push @{$align_stat{$rcid}}, ("$refer_cover","$refer_len","$contigs_cover", "$cnum", "$contiglist");
	#}
	
	if ($refer_cover > $max_cover){
		$max_cover=$refer_cover;		
	}
}

foreach my $rcid(keys %align_stat){
	print "$rcid\t$align_stat{$rcid}[0]\t$align_stat{$rcid}[1]\t";
	printf "%.2f", $align_stat{$rcid}[2];
	print "\t$align_stat{$rcid}[3]\t";

	my $rc_score=0;
	$rc_score +=100 if($align_stat{$rcid}[0] >=98);
	$rc_score +=95 if($align_stat{$rcid}[0] >=95 && $align_stat{$rcid}[0] <98);
	$rc_score +=90 if($align_stat{$rcid}[0] >=90 && $align_stat{$rcid}[0] <95);
	$rc_score +=80 if($align_stat{$rcid}[0] >=80 && $align_stat{$rcid}[0] <90);
	$rc_score +=60 if($align_stat{$rcid}[0] >=60 && $align_stat{$rcid}[0] <80);
	$rc_score +=20 if($align_stat{$rcid}[0] >=30 && $align_stat{$rcid}[0] <60);
	print "$rc_score\t";

	my $rl_score=0;
	$rl_score +=50 if ($align_stat{$rcid}[1]>=50);
	print "$rl_score\t";
	
	my $cc_score=0;
	$cc_score +=70 if($align_stat{$rcid}[2] >=90);
	$cc_score +=60 if($align_stat{$rcid}[2] >=80 && $align_stat{$rcid}[2] <90);
	$cc_score +=40 if($align_stat{$rcid}[2] >=60 && $align_stat{$rcid}[2] <80);
	print "$cc_score\t";	

	my $cn_score=0;
	$cn_score +=100 if($align_stat{$rcid}[3] ==3); 
	$cn_score +=50 if($align_stat{$rcid}[3] <=5 && $align_stat{$rcid}[3] !=3);
	#$cn_score +=10 if($align_stat{$rcid}[3] >5 && $align_stat{$rcid}[3] <=10);
	#$cn_score +=10 if($align_stat{$rcid}[3] >10 && $align_stat{$rcid}[3] <=20);
	print "$cn_score\t";

	$align_stat{$rcid}[4]=~s/\,$//g;
	my @info=split (/,/, $align_stat{$rcid}[4]);
	my $cl_score=0;
	for(@info){
		$cl_score +=150 if ($_ >=70000);
		$cl_score +=60 if ($_ >=30000 && $_<70000);
		$cl_score +=20 if ($_ >=10000 && $_<30000);		
	}
	print "$cl_score\t";
	
	#my $an_score=0;
	#$an_score +=30 if($align_stat{$rcid}[3] <=10);
	#$an_score +=20 if($align_stat{$rcid}[3] >10 && $align_stat{$rcid}[3] <=20);
	#$an_score +=5 if($align_stat{$rcid}[3] >20 && $align_stat{$rcid}[3] <=30);
	#print "$an_score\t";
	
	my $rc_bonus=0;
	$rc_bonus=10 if ($align_stat{$rcid}[0] eq $max_cover);
	#print STDERR "$rcid\t$align_stat{$rcid}[0]\t$max_cover\t$rc_bonus\n";
	print "$rc_bonus\t";

	my $total_score=$rc_score+$rl_score+$cc_score+$cn_score+$cl_score+$rc_bonus;
	print "$total_score\t$align_stat{$rcid}[4]\n";
}
