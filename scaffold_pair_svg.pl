#!/usr/bin/perl-w
=head1 Program Description

This program is designed to draw list match situation in assigned region of two chromosomes

=head1 Contact & Version

Author: Zengpeng, zengpeng@genomics.org.cn
Version: 1.0,  Date: 2011-12-8

	=head1 Command-line Option
	--chr_1 <str>           chrome_up (needed)
	--start_1 <num>         start_cutoff_up (bp) (needed)
	--end_1  <num>          end_cutoff_up (needed)
	--chr_2 <str>           chrome_down (needed)
	--start_2 <num>         start_cutoff_down (bp) (needed)
	--end_2  <num>          end_cutoff_down (needed)	
	--stand <str>	          match direction, minus stands for reversr direction ,default=plus		


	=head1 Usage Exmples

	perl gene_pair_svg.pl list out_svg --chr_1 chr1 --start_1 2000 --end_1 5000 --chr_2 chr2 --start_2 7000 --end_2 10000 

=cut




	use strict;
	use lib "/home/xzhong/common_bin/Draw_scripts";
	use SVG;
	use FontSize;
	use Data::Dumper;
	use Getopt::Long;

	my ($chr1,$st1,$en1,$chr2,$st2,$en2,$stand);
	my ($Verbose,$Help);
	GetOptions(
			"chr_1:s"=>\$chr1,
			"start_1:i"=>\$st1,
			"end_1:i"=>\$en1,
			"chr_2:s"=>\$chr2,
			"start_2:i"=>\$st2,
			"end_2:i"=>\$en2,        
			"stand:s"=>\$stand,
			"verbose"=>\$Verbose,
			"help"=>\$Help
		  );
##$File_in_dir_num ||= 1000;

	die `pod2text $0` if ($Help);

	my ($list,$out_svg)=@ARGV;
#my ($gff,$blast,$chr1,$st1,$en1,$chr2,$st2,$en2,$out_svg)=@ARGV;

##chr06	PGSC0003DMP200000010_SOLTU	41322150	41327569	+
	$stand ||="plus";
	my @gene_up;
	my @gene_down;
	my %match;

	open BLAST,$list or die "$!";
	while(<BLAST>){
		chomp;
		my @f=split;
		#next if($_[1]<$st1 || $_[1]>$en1);
		#next if($_[2]<$st1 || $_[2]>$en1);
		#next if($_[5]<$st2 || $_[5]>$en2);	
		#next if($_[4]<$st2 || $_[4]>$en2);
		push @{$match{"$f[0]&$f[3]"}},[$f[1],$f[2],$f[4],$f[5]];
	}
close BLAST;

my ($start_gene_aa,$end_gene_aa,$start_gene_bb,$end_gene_bb)=($st1,$en1,$st2,$en2);
my $divider=($end_gene_aa-$start_gene_aa)/1500;


&draw_svg($out_svg);


sub draw_svg{
	my $out=shift;
	my $len_aa=$end_gene_aa-$start_gene_aa;
	my $len_bb=$end_gene_bb-$start_gene_bb;
	print "$len_aa\t$len_bb\n";
	my $width=2000;
	#($len_aa>$len_bb)?$len_aa:$len_bb;
#	$width=$width/$divider;
	$width+=200;
	my $height=500;
	my $svg = SVG->new('width',$width,'height',$height);
	$svg->line(x1=>100,y1=>100,x2=>100+$len_aa/$divider,y2=>100, 'stroke-width'=>5,stroke=>'black',opacity=>0.4);
	$svg->text('x',20,'y',100,'-cdata',$chr1,'font-size',15);
	$svg->line(x1=>100,y1=>290,x2=>100+$len_bb/$divider,y2=>290, 'stroke-width'=>5,stroke=>'black',opacity=>0.4);
	$svg->text('x',20,'y',290,'-cdata',$chr2,'font-size',15);

	my $st_1=substr($start_gene_aa/1000000,0,7);
	my $st_2=substr($start_gene_bb/1000000,0,7);
	$svg->text('x',20,'y',75,'-cdata',"$st_1 MB",'font-size',15);      #line up start
		my $en_1=substr($end_gene_aa/1000000,0,7);
	my $en_2=substr($end_gene_bb/1000000,0,7);
	$svg->text('x',$width-100,'y',75,'-cdata',"$en_1 MB",'font-size',15);   #line up end
		if($stand eq 'plus'){
			$svg->text('x',20,'y',310,'-cdata',"$st_2 MB",'font-size',15);
			$svg->text('x',$width-100,'y',310,'-cdata',"$en_2 MB",'font-size',15);
		}else{
			$svg->text('x',20,'y',310,'-cdata',"$en_2 MB",'font-size',15);
			$svg->text('x',$width-100,'y',310,'-cdata',"$st_2 MB",'font-size',15);
		}
	$svg->text('x',20,'y',200,'-cdata',$stand,'font-size',15);

	my $pair="$chr1&$chr2";
	my @pair=@{$match{$pair}};
	for my $aa(@pair){	
		my @aa=@$aa;
		my $st_a=100+($aa[0]-$start_gene_aa)/$divider;my $en_a=100+($aa[1]-$start_gene_aa)/$divider;
		my $st_b=100+($aa[2]-$start_gene_bb)/$divider;my $en_b=100+($aa[3]-$start_gene_bb)/$divider;
		$st_b=(200-$st_b+$len_bb/$divider) if($stand eq 'minus');
		$en_b=(200-$en_b+$len_bb/$divider) if($stand eq 'minus');

		$svg->polygon('points',[$st_a,105,$en_a,105,$en_b,285,$st_b,285],'fill','blue','stroke','blue',opacity=>0.5);
	}





	open  OUT, ">$out_svg" || die "fail create $!";
	print OUT $svg->xmlify();
	close OUT;	
}

