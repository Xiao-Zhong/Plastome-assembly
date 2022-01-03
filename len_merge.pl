#!/usr/bin/perl -w
use strict;

my $list=shift;

open IN, "$list" || die "$!\n";
while(<IN>){
	chomp;
	my $id=$_;
	print "$id\t";
	if(-e "$id\_denovo_Chloe.fa.len"){
	if (-z "$id\_denovo_Chloe.fa.len"){
		print "Chloe_failed\t";
	}else{
		open IN2, "$id\_denovo_Chloe.fa.len" || die "$!\n";
		while(<IN2>){
			chomp;
			my @info=split;
			print "$info[1]\t";
		}
		close IN2;
	}
	
	if (-z "$id\_after_fill_gap_check.fa.len"){
		print "NA\tNA\t";
	}else{
		open IN3, "$id\_after_fill_gap_check.fa.len" || die "$!\n";
		while(<IN3>){
			chomp;
			my @info2=split;
			my $fillgap=$info2[1]-124;
			print "$fillgap\t$info2[1]\t";
		}
		close IN3;
	}
	
	if (-z "$id\.pilon.final.fasta.len"){
		print "NA\tNA\t";
	}else{
		open IN4, "$id\.pilon.final.fasta.len"|| die "$!\n";
		while(<IN4>){
			chomp;
			my @info3=split;
			my $pilon=$info3[1]+124;
			print "$pilon\t$info3[1]\t";
		}
		close IN4;
	}
	}else{
		print "denovo_failed\tNA\tNA\tNA\tNA\t";
	}

	if (-z "$id\_guided_1st_pilon.fasta.len"){
		print "NA\t";
	}else{
		open IN5, "$id\_guided_1st_pilon.fasta.len" || die "$!\n";
		while(<IN5>){
			chomp;
			my @info4=split;
			print "$info4[1]\t";
		}
		close IN5;
	}
	
	if (-z "$id\_after_fill_gap_check_guided.fa.len"){
		print "NA\tNA\t";
	}else{
		open IN6, "$id\_after_fill_gap_check_guided.fa.len" || die "$!\n";
		while(<IN6>){
			chomp;
			my @info5=split;
			my $fillgap=$info5[1]-124;
			print "$fillgap\t$info5[1]\t";
		}
		close IN6;
	}
	
	if (-z "$id\.pilon.guided.final.fasta.len"){
		print "NA\tNA\t";
	}else{
		open IN7, "$id\.pilon.guided.final.fasta.len"|| die "$!\n";
		while(<IN7>){
			chomp;
			my @info6=split;
			my $pilon=$info6[1]+124;
			print "$pilon\t$info6[1]\t";
		}
		close IN7;
	}

	if (-z "$id\.chlo.fa.len"){
		print "NA\t";
	}else{
		my $tmp;
		open IN8, "$id\.chlo.fa.len"|| die "$!\n";
		while(<IN8>){
			chomp;
			my @info7=split;
			$tmp .= "$info7[1]\|";
		}
		close IN8;
		$tmp=~s/\|$//g;
		print "$tmp\t";
	}

	if (-z "$id\_org_after_fill_gap_check.fa.len"){
		print "NA\tNA\t";
	}else{
		my $tmp; my $tmp2;
		open IN9, "$id\_org_after_fill_gap_check.fa.len"|| die "$!\n";
		while(<IN9>){
			chomp;
			my @info8=split;
			my $fillgap=$info8[1]-124;
			$tmp .="$fillgap\|";
			$tmp2 .="$info8[1]\|";
			#print "$fillgap\t$info8[1]\t";
		}
		close IN9;
		$tmp=~s/\|$//g;
		$tmp2=~s/\|$//g;
		print "$tmp\t$tmp2\t";
	}
	
	if (-z "$id\.org.pilon.final.fasta.len"){
		print "NA\tNA\t";
	}else{
		my $tmp; my $tmp2;
		open IN10, "$id\.org.pilon.final.fasta.len"|| die "$!\n";
		while(<IN10>){
			chomp;
			my @info9=split;
			my $pilon=$info9[1]+124;
			$tmp .="$pilon\|";
			$tmp2 .="$info9[1]\|";
			#print "$pilon\t$info9[1]\n";
		}
		close IN10;
		$tmp=~s/\|$//g;
		$tmp2=~s/\|$//g;
		print "$tmp\t$tmp2\t";
	}

	if (-z "$id\.chlo.fa.raw.len"){
		print "NA\n";
	}else{
		my $tmp;
		open IN11, "$id\.chlo.fa.raw.len"|| die "$!\n";
		while(<IN11>){
			chomp;
			my @info10=split;
			$tmp .= "$info10[1]\|";
		}
		close IN11;
		$tmp=~s/\|$//g;
		print "$tmp\n";
	}
	
}
close IN;
