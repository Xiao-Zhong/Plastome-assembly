#!/usr/bin/perl -w
use strict;

die "Usage:\n	perl $0 <input_dir> <sample_lable> <cover_list> <kvalue_list> <reference_file_name> \n" if @ARGV<5;

my $indir=shift;
my $sample=shift;
my $clist=shift;
my $klist=shift;
my $refer=shift;

my %align_stat;my $max_cover=0;my $maxcid;
my @cover=split(/\,/, $clist);
my @k=split(/\,/, $klist);
for(my $i=0; $i<@cover; $i++){
	for(my $j=0; $j<@k; $j++){
		#print "$cover[$i]\t$k[$j]\n";
		my $align_file="$indir/$sample\_$cover[$i]\_$k[$j]/$refer\_vs_contigs.fa.coords.cp.id";
		
		if (-e "$align_file"){
		unless (-z "$align_file"){
			#print "$align_file\n";
			my $total_cover; my %idhash;
			open IN, "$align_file" || die "$!\n";
			while(<IN>){
				chomp;
				my @info=split;
				$total_cover +=$info[14];
				$idhash{$info[18]}++ if (! exists $idhash{$info[18]});
			}
			close IN;
			#print "$sample\_$cover[$i]\_$k[$j]\t$total_cover\t";
			my $idnum=keys %idhash; #print "$idnum\n";
			#print scalar keys %idhash, "\n";		
			push @{$align_stat{"$sample\_$cover[$i]\_$k[$j]"}}, ("$total_cover", "$idnum");
			
			if ($total_cover >$max_cover){
				$max_cover=$total_cover;
				$maxcid="$sample\_$cover[$i]\_$k[$j]";
			}
		}
		}
	}
}
#print "$max_cover\n";

my $maxsid;my $max_score=0;
foreach my $dirid(keys %align_stat){
	print "$dirid\t$align_stat{$dirid}[0]\t$align_stat{$dirid}[1]\t";
	my $score=0;
	$score=100 if ($align_stat{$dirid}[0]>=99);
	$score=90 if ($align_stat{$dirid}[0]>=98 && $align_stat{$dirid}[0]<99);
	$score=75 if ($align_stat{$dirid}[0]>=97 && $align_stat{$dirid}[0]<98);
	$score=55 if ($align_stat{$dirid}[0]>=96 && $align_stat{$dirid}[0]<97);
	$score=30 if ($align_stat{$dirid}[0]>=95 && $align_stat{$dirid}[0]<96);
		
	$score +=50 if ($align_stat{$dirid}[1]<=5);
	$score +=40 if ($align_stat{$dirid}[1]>5 && $align_stat{$dirid}[1]<=10);
	$score +=20 if ($align_stat{$dirid}[1]>10 && $align_stat{$dirid}[1]<=20);
	
	$score +=10 if ($dirid eq "$maxcid");
	print "$score\n";

	if ($score >$max_score){
		$max_score=$score;
		$maxsid=$dirid;
	}
}
print "Best_candidate: $maxsid\n";











