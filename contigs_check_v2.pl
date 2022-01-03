#!/usr/bin/perl -w
use strict;

my $innucmer=shift;
my $len_cutoff=shift;

$len_cutoff ||="74,92,10,22,18,29";# xxKB, LSC/SSC/IR "74,92,10,22,17,29"
my @info=split(/,/,$len_cutoff);

my %bhash;
my %chash;
my %cidhash;
open IN, "$innucmer" || die "$!\n";
while(<IN>){
	chomp;
	next unless (/\s+NODE\_/);
	#print "$_\n";
	my @info=split;
	$bhash{$info[18]}=$info[12];
	#my $cover=$1 if ($info[18]=~/\_cov\_(\S+)/);
	my $cover=$1 if ($info[18]=~/\_cov\_(\S+)\|/);
	#print "$info[18]\t$cover\n";
	$chash{$info[18]}=$cover;

	my $s=($info[3]<$info[4])?("+"):("-");
	$cidhash{$info[18]}{$s}++ if (not exists $cidhash{$info[18]}{$s});
}
close IN;

my($x, $y, $z);
my $cnumber=scalar keys %bhash;
#print "$cnumber\n";
if($cnumber==3){
	my $num=1;
	foreach my $id(sort {$chash{$b}<=>$chash{$a}} keys %chash){
		$x="LSC-YES" if ($num!=1 && $bhash{$id}>=$info[0]*1000 && $bhash{$id}<=$info[1]*1000);
		$y="SSC-YES" if ($num!=1 && $bhash{$id}>=$info[2]*1000 && $bhash{$id}<=$info[3]*1000);
		if (exists $cidhash{$id}{"+"} && exists $cidhash{$id}{"-"}){
			$z="IR-YES" if ($num==1 && $bhash{$id}>=$info[4]*1000 && $bhash{$id}<=$info[5]*1000);
		}
		$num++;
	}
}

if ((defined $x && $x eq "LSC-YES") && (defined $y && $y eq "SSC-YES") && (defined $z && $z eq "IR-YES")){
	print "PASS";
}else{
	print "FAIL";
}


