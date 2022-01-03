#!/usr/bin/perl -w
use strict;

die "Usage:\n   perl $0 <alignment_output: *.coords.cp.id> \n" if @ARGV<1;

## add the process "filter overlaps in alignment".

my $inalign=shift;

my %refer_contigs;my %refer_contigs_contig;
open IN, "$inalign" || die "$!\n";
while(<IN>){
	chomp;
	my @info=split;
	my @info2=split(/\|/, $info[18]);
	#print "$info[17]\t$info2[1]\n";
	my $a=$info[3]; my $b=$info[4];
	if ($info[3]>$info[4]){
		$info[3]=$b;
		$info[4]=$a;
	}
	#print join "\t", @info, "\n";
	push @{$refer_contigs{"$info[17]\t$info2[1]"}}, [@info];
	push @{$refer_contigs_contig{"$info[17]\t$info2[1]\t$info2[0]"}}, [@info];
}
close IN;

my %rccalenhash;
foreach my $rccid(keys %refer_contigs_contig){
	#print "$rccid\n";
	@{$refer_contigs_contig{$rccid}}=sort {$a->[3] <=> $b->[3]} @{$refer_contigs_contig{$rccid}};
	my $rccalen;
	for(my $i=0; $i<@{$refer_contigs_contig{$rccid}}; $i++){
		#print join "\t", @{$refer_contigs_contig{$rccid}[$i]}, "\n";
		$rccalen +=$refer_contigs_contig{$rccid}[$i][4]-$refer_contigs_contig{$rccid}[$i][3]+1;
		#print "$rccalen\+\t";
		if ($i<@{$refer_contigs_contig{$rccid}}-1){
			for (my $j=$i+1; $j<@{$refer_contigs_contig{$rccid}}; $j++){
				if ($refer_contigs_contig{$rccid}[$i][4] >= $refer_contigs_contig{$rccid}[$j][3] && $refer_contigs_contig{$rccid}[$i][4] <= $refer_contigs_contig{$rccid}[$j][4]){
					#print "$rccid\t$refer_contigs_contig{$rccid}[$i][3]\ta\t$i\t$j\n";
					$rccalen -=$refer_contigs_contig{$rccid}[$i][4]-$refer_contigs_contig{$rccid}[$j][3]+1;
					#print "$rccalen\-\t";
				}			
				if ($refer_contigs_contig{$rccid}[$i][3] <= $refer_contigs_contig{$rccid}[$j][3] && $refer_contigs_contig{$rccid}[$i][4] > $refer_contigs_contig{$rccid}[$j][4]){
					#print "$rccid\t$refer_contigs_contig{$rccid}[$i][3]\tb\t$i\t$j\n";
					$rccalen -=$refer_contigs_contig{$rccid}[$j][4]-$refer_contigs_contig{$rccid}[$j][3]+1;
					#print "$rccalen\-\t";
				}
			}
		}
	}

	#print "$rccid\t$rccalen\t$refer_contigs_contig{$rccid}[0][12]\t";
	#printf "%.2f", $rccalen/$refer_contigs_contig{$rccid}[0][12]*100;
	#print "\n";		
	my @info=split(/\t/, $rccid);	
	push @{$rccalenhash{"$info[0]\t$info[1]"}}, $rccalen;
}

#for my $key (keys %rccalenhash){
#	for(my $i=0; $i<@{$rccalenhash{$key}};$i++){
		#print "$key\t$rccalenhash{$key}[$i]\n";
#	}
#}

my %align_stat; my $max_cover=0;
foreach my $rcid(keys %refer_contigs){
	#print "$rcid\n";
	my $rcalen;
	for (my $i=0; $i<@{$refer_contigs{$rcid}}; $i++){
		#print join "\t", @{$refer_contigs{$rcid}[$i]},"\n";
		$rcalen +=$refer_contigs{$rcid}[$i][1]-$refer_contigs{$rcid}[$i][0]+1;
		if ($i<@{$refer_contigs{$rcid}}-1){
			for (my $j=$i+1; $j<@{$refer_contigs{$rcid}}; $j++){
				if ($refer_contigs{$rcid}[$i][1] >= $refer_contigs{$rcid}[$j][0] && $refer_contigs{$rcid}[$i][1] <= $refer_contigs{$rcid}[$j][1]){
					#print "$rcid\t$refer_contigs{$rcid}[$i][0]\ta\t$i\t$j\n";
					$rcalen -=$refer_contigs{$rcid}[$i][1]-$refer_contigs{$rcid}[$j][0]+1;
				}
				if ($refer_contigs{$rcid}[$i][0] <= $refer_contigs{$rcid}[$j][0] && $refer_contigs{$rcid}[$i][1] > $refer_contigs{$rcid}[$j][1]){
					#print "$rcid\t$refer_contigs{$rcid}[$i][0]\tb\t$i\t$j\n";
					$rcalen -=$refer_contigs{$rcid}[$j][1]-$refer_contigs{$rcid}[$j][0]+1;
				}
			}
		}
	}

	#print "$rcid\t$rcalen\t$refer_contigs{$rcid}[0][11]\t";
	#printf "%.2f", $rcalen/$refer_contigs{$rcid}[0][11]*100;
	#print "\n";
	
	my $total_c_len=0; my %chash; my $contiglist;
	for (my $i=0; $i<@{$refer_contigs{$rcid}}; $i++){
		#$refer_cover +=$refer_contigs{$rcid}[$i][14];
		my @info=split(/\|/, $refer_contigs{$rcid}[$i][18]);
		unless (exists $chash{$info[0]}){
			$total_c_len +=$refer_contigs{$rcid}[$i][12];
			$contiglist .="$refer_contigs{$rcid}[$i][12],";	
			$chash{$info[0]}++;
		}
		#$total_c_aligned_len +=$refer_contigs{$rcid}[$i][7];
	}

	my $total_c_aligned_len=0;
	for (my $i=0; $i<@{$rccalenhash{$rcid}}; $i++){
		$total_c_aligned_len +=$rccalenhash{$rcid}[$i];	
	}

	my $cnum=keys %chash;
	my $refer_len=$refer_contigs{$rcid}[0][11];
	my $refer_cover=$rcalen/$refer_len*100;
	my $contigs_cover=$total_c_aligned_len/$total_c_len*100;
	#print "$rcid\t$refer_cover\t$refer_len\t$total_c_len\t$total_c_aligned_len\t$contigs_cover\t$cnum\n";
	push @{$align_stat{$rcid}}, ("$refer_cover","$refer_len","$contigs_cover", "$cnum", "$contiglist", "$total_c_len");
	
	if ($refer_cover > $max_cover){
		$max_cover=$refer_cover;		
	}
}

foreach my $rcid(keys %align_stat){
	print "$rcid\t";
	printf "%.2f", $align_stat{$rcid}[0];
	print "\t$align_stat{$rcid}[1]\t";
	printf "%.2f", $align_stat{$rcid}[2];
	print "\t$align_stat{$rcid}[3]\t$align_stat{$rcid}[5]\t";

	my $rc_score=0;
	$rc_score +=$align_stat{$rcid}[0];
	#$rc_score +=100 if($align_stat{$rcid}[0] >=98);
	#$rc_score +=95 if($align_stat{$rcid}[0] >=95 && $align_stat{$rcid}[0] <98);
	#$rc_score +=90 if($align_stat{$rcid}[0] >=90 && $align_stat{$rcid}[0] <95);
	#$rc_score +=80 if($align_stat{$rcid}[0] >=80 && $align_stat{$rcid}[0] <90);
	#$rc_score +=60 if($align_stat{$rcid}[0] >=60 && $align_stat{$rcid}[0] <80);
	#$rc_score +=20 if($align_stat{$rcid}[0] >=30 && $align_stat{$rcid}[0] <60);
	printf "%.2f", "$rc_score";print "\t";

	my $rl_score=0;
	$rl_score +=50 if ($align_stat{$rcid}[1]>=50);
	print "$rl_score\t";
	
	my $cc_score=0;
	$cc_score +=$align_stat{$rcid}[2];
	#$cc_score +=70 if($align_stat{$rcid}[2] >=90);
	#$cc_score +=60 if($align_stat{$rcid}[2] >=80 && $align_stat{$rcid}[2] <90);
	#$cc_score +=40 if($align_stat{$rcid}[2] >=60 && $align_stat{$rcid}[2] <80);
	printf "%.2f", "$cc_score";print "\t";

	my $cn_score=0;
	$cn_score +=100 if($align_stat{$rcid}[3] ==3); 
	$cn_score +=60 if($align_stat{$rcid}[3] <=5 && $align_stat{$rcid}[3] !=3);
	#$cn_score +=10 if($align_stat{$rcid}[3] >5 && $align_stat{$rcid}[3] <=10);
	#$cn_score +=10 if($align_stat{$rcid}[3] >10 && $align_stat{$rcid}[3] <=20);
	print "$cn_score\t";

	$align_stat{$rcid}[4]=~s/\,$//g;
	my @info=split (/,/, $align_stat{$rcid}[4]);
	my $cl_score=0;
	for(@info){
		$cl_score +=100 if ($_ >=70000);
		#$cl_score +=60 if ($_ >=30000 && $_<70000);
		$cl_score +=45 if ($_ >=30000 && $_<70000);
		$cl_score +=20 if ($_ >=10000 && $_<30000);		
	}
	print "$cl_score\t";
	
	my $ctl_score=0;
	#$ctl_score +=40 if ($align_stat{$rcid}[5]>110000);
	$ctl_score +=30 if ($align_stat{$rcid}[5]>100000);
	print "$ctl_score\t";
	
	#my $an_score=0;
	#$an_score +=30 if($align_stat{$rcid}[3] <=10);
	#$an_score +=20 if($align_stat{$rcid}[3] >10 && $align_stat{$rcid}[3] <=20);
	#$an_score +=5 if($align_stat{$rcid}[3] >20 && $align_stat{$rcid}[3] <=30);
	#print "$an_score\t";
	
	my $rc_bonus=0;
	$rc_bonus=10 if ($align_stat{$rcid}[0] eq $max_cover);
	#print STDERR "$rcid\t$align_stat{$rcid}[0]\t$max_cover\t$rc_bonus\n";
	print "$rc_bonus\t";

	my $total_score=$rc_score+$rl_score+$cc_score+$cn_score+$cl_score+$ctl_score+$rc_bonus;
	printf "%.2f", "$total_score"; 
	print "\t$align_stat{$rcid}[4]\n";

}
