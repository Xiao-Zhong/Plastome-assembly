#!/usr/bin/perl -w
use strict;

my $infa=shift;

open IN, "$infa" || die "$!\n";
$/=">";<IN>;$/="\n";
while(<IN>){
	my $id=$1 if(/(\S+)/);
	#print "$id\n";
	$/=">";
	my $seq=<IN>;
	chomp($seq);
	$seq=~s/\s+//g;
	$seq=~tr/atcg/ATCG/;
	#print "$seq\n";
	$/="\n";
	
	my $head_mark=substr($seq, 0, 10);
	#print "$head_mark\n";
	my $mark_num=$seq=~s/$head_mark/$head_mark/g;
	#print "$mark_num\n";

	if ($mark_num>1){
		my $mark_last_s=$-[0]+1;
		my $mark_last_e=$+[0];
		#print "$mark_last_s\t$mark_last_e\n";
		
		my $overlaplen=length($seq)-$mark_last_s+1;
		my $overlap=substr($seq, $mark_last_s-1, $overlaplen);
		#print "$overlap\n";

		if (my $overlap_num=$seq=~s/^$overlap//){
			#print "$overlap_num\n";
			print ">$id\t\"$head_mark\"|$mark_num|$mark_last_s|$mark_last_e|$overlap_num|$overlaplen|\"$overlap\"\n";
			print "$seq\n";
		}else{
			print ">$id\tmultiple_pattern_matchs_but_no_overlap_at_ends\n$seq\n";
			
		}
	}else{
		print ">$id\tno_pattern_match\n$seq\n";
	}
}
close IN;


