#!/usr/bin/perl

=head1 Usage

perl auto_assembly.pl [options] <sample_list> <reads_path>

=head1 Example
    perl auto_assembly.pl -outdir output -step 0123456789 sample.list /home/xzhong/working_data_01/01_Assembly/37_pipeline
##### each step description #####
    step 0: bbnorm.sh;
    step 1: cutadapt.sh;
    step 2: spades.sh;
    step 3: flash.sh or bbmerge.sh[default];
    step 4: velveth.sh & velvetg.sh & nucmer.sh & contigs_check.pl;
    step 5: chloe.sh for denovo assembly, or carry on reference-guided assembly;
    step 6: fill gaps by gapfiller;
    step 7: bwa;
    step 8: pilon;
    step 9: statistics.
#################################
=cut

use strict;
use Getopt::Long;
use File::Basename;
#use FinBin qw($Bin);
my $Bin="/home/xzhong/working_data_01/01_Assembly/37_pipeline";

my ($Step, $Outdir, $Help);
my ($Maxinputreads, $Maxoutputreads);
my ($Adapter1, $Adapter2);
my $Merge;
my ($Mink, $Maxk, $Stepk, $Coverlist,$Refer_DB,$Insert_size,$Insert_size_SD,$Exp_cover,$Contigs_len_cutoff);
my $Clean;

GetOptions(
	"outdir:s"=>\$Outdir,
	"step:s"=>\$Step,

	"maxinputreads:i"=>\$Maxinputreads,
	"maxoutputreads:i"=>\$Maxoutputreads,

	"adapter1:s"=>\$Adapter1, #adapter1, Sequence of an adapter that was ligated to the 3' end.
	"adapter2:s"=>\$Adapter2, #adapter2, 3' adapter to be removed from second read in a pair.

	"merge:s"=>\$Merge,

	"mink:i"=>\$Mink,
	"maxk:i"=>\$Maxk,
	"stepk:i"=>\$Stepk, #Velvet will then hash from k=m to k=M with a step of s.
	"coverlist:s"=>\$Coverlist, #Velvet's -cov_cutoff.
	"insertsize:i"=>\$Insert_size,
	"insertsizesd:i"=>\$Insert_size_SD,
	"expectedcover:i"=>\$Exp_cover,
	"referdb:s"=>\$Refer_DB,
	"contigslen:s"=>\$Contigs_len_cutoff, # LSG,SSG and IR 's length cutoff(Kb).
	
	"clean:s"=>\$Clean,	
	"help!"=>\$Help
);

$Outdir ||=".";
#$Step ||='0123456789';
$Maxinputreads ||=-1;
$Maxoutputreads ||=-1;
$Adapter1 ||="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
$Adapter2 ||="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";
$Merge ||="bbmerge";
$Mink ||=51;
$Maxk ||=115;
$Stepk ||=20;
$Refer_DB ||="$Bin/ncbi_cp.fa";
$Insert_size ||=350;
$Insert_size_SD ||=$Insert_size*1/10;
$Exp_cover ||=30;
$Coverlist ||="10,15,7,20";# print "$Coverlist\n";
my @Cover=split(/\,/,$Coverlist);
$Contigs_len_cutoff ||="74,92,10,22,18,29";# "Kb", LSC/SSC/IR
$Clean ||="no";
die `pod2text $0` if (@ARGV == 0 || $Help);

mkdir($Outdir) unless(-d $Outdir);
my $rundir=`pwd`;
chomp $rundir;#print "$rundir";
$Outdir="$rundir/$Outdir"; #print "$Outdir\n";
#`cp $Refer_DB $Outdir/Cp_DB_updated.fa`;

my $sample_list=shift;
my $reads_dir=shift;

`rm "$Outdir/denovo_pass_id.log"` if (-e "$Outdir/denovo_pass_id.log");
`rm "$Outdir/guided_pass_id.log"` if (-e "$Outdir/guided_pass_id.log");
`rm "$Outdir/both_fail_id.log"` if (-e "$Outdir/both_fail_id.log");

open IN, "$sample_list" || die "fail $sample_list";
while(<IN>){
	chomp;
	my $sample="$_";
	#print "$sample\n";
	mkdir("$Outdir/$sample") unless(-d "$Outdir/$sample");
	
	## run bbnorm.sh ## 
	if ($Step=~/0/){
		`ls $reads_dir/$sample\_*\_R1*.fastq.gz $reads_dir/$sample\_*\_R2*.fastq.gz > tmp`;
		open IN1, "tmp" || die "$!\n";
		my $input1=<IN1>; chomp $input1;
		my $input2=<IN1>; chomp $input2;
		#print "$input1\t$input2\n";
		close IN1;
		
		my $bbnorm_shell="$Outdir/$sample/$sample.bbnorm.shell";
		open OUT, ">$bbnorm_shell" || die "fail $bbnorm_shell";
		print OUT "$Bin/bbnorm.sh in=$input1 in2=$input2 tablereads=$Maxinputreads reads=$Maxoutputreads zerobin=t prefilter=t maxdepth=1000 lowbindepth=10 highbindepth=500 hist=$Outdir/$sample/$sample\_bbnorm.hist out=$Outdir/$sample/$sample\_bbnorm.1.fastq.gz out2=$Outdir/$sample/$sample\_bbnorm.2.fastq.gz outt=$Outdir/$sample/$sample\_bbnorm_excluded.fastq.gz\n";
		close OUT;
		`rm tmp`;
		`sh $bbnorm_shell 1>$bbnorm_shell.log 2>$bbnorm_shell.err`;
	}

	## run cutadapt.sh ##	
	if ($Step=~/1/){
		my $cutadapter_shell="$Outdir/$sample/$sample.cutadapter.shell";
		open OUT1, ">$cutadapter_shell" || die "fail $cutadapter_shell";
		print OUT1 "$Bin/cutadapt -a $Adapter1 -A $Adapter2 -o $Outdir/$sample/$sample\_trimmed.1.fastq.gz -p $Outdir/$sample/$sample\_trimmed.2.fastq.gz $Outdir/$sample/$sample\_bbnorm.1.fastq.gz $Outdir/$sample/$sample\_bbnorm.2.fastq.gz\n";
		close OUT1;
		`sh $cutadapter_shell 1>$cutadapter_shell.log 2>$cutadapter_shell.err`;
	}
	
	## run spades.sh ##
	if ($Step=~/2/){
		my $spades_shell="$Outdir/$sample/$sample.spades.shell";
		open OUT1x, ">$spades_shell" || die "fail $spades_shell";
		print OUT1x "$Bin/python $Bin/spades.py --only-error-correction --careful --pe1-1 $Outdir/$sample/$sample\_trimmed.1.fastq.gz --pe1-2 $Outdir/$sample/$sample\_trimmed.2.fastq.gz -o $Outdir/$sample/spades_outdir\n";
		close OUT1x;
		`sh $spades_shell 1>$spades_shell.log 2>$spades_shell.err`;		
	}

	## run flash.sh ##
	if($Step=~/3/ && ($Merge ne "bbmerge")){
		my $flash_shell="$Outdir/$sample/$sample.flash.shell";
		open OUT2, ">$flash_shell" || die "fail $flash_shell";
		print OUT2 "$Bin/flash -q -m 10 -x 0.2 --max-overlap 125 -z -d $Outdir/$sample -o $sample $Outdir/$sample/spades_outdir/corrected/$sample\_trimmed.1.fastq.00.0_0.cor.fastq.gz $Outdir/$sample/spades_outdir/corrected/$sample\_trimmed.2.fastq.00.0_0.cor.fastq.gz\n";
		close OUT2;
		`sh $flash_shell 1>$flash_shell.log 2>$flash_shell.err`;
	}	
	## run bbmerge.sh [default] ##
	if($Step=~/3/ && ($Merge eq "bbmerge")){
		my $bbmerge_shell="$Outdir/$sample/$sample.bbmerge.shell";
		open OUT2x, ">$bbmerge_shell" || die "fail $bbmerge_shell";
		print OUT2x "$Bin/bbmerge.sh in1=$Outdir/$sample/spades_outdir/corrected/$sample\_trimmed.1.fastq.00.0_0.cor.fastq.gz in2=$Outdir/$sample/spades_outdir/corrected/$sample\_trimmed.2.fastq.00.0_0.cor.fastq.gz out=$Outdir/$sample/$sample.extendedFrags.fastq.gz outu1=$Outdir/$sample/$sample.notCombined_1.fastq.gz outu2=$Outdir/$sample/$sample.notCombined_2.fastq.gz\n";
		close OUT2x;
		`sh $bbmerge_shell 1>$bbmerge_shell.log 2>$bbmerge_shell.err`;
	}

	## run velveth.sh & velvetg.sh & nucmer.sh & contigs_check.pl ##
	if($Step=~/4/){
		open OUTa, ">$Outdir/$sample/cover_kvalue.list" || die "$!\n";
		lable: for(my $j=0; $j<@Cover; $j++){
			#print "$Cover[$j]\n";
			my $velveth_shell="$Outdir/$sample/$sample\_$Cover[$j]\_velveth_shell";
			open OUT3, ">$velveth_shell" || die "fail $velveth_shell";
			print OUT3 "$Bin/velveth $Outdir/$sample/$sample\_$Cover[$j] $Mink,$Maxk,$Stepk -shortPaired -fastq.gz -separate $Outdir/$sample/$sample.notCombined_1.fastq.gz $Outdir/$sample/$sample.notCombined_2.fastq.gz -long $Outdir/$sample/$sample.extendedFrags.fastq.gz -short $Outdir/$sample/spades_outdir/corrected/$sample\_trimmed.00.0_0.cor.fastq.gz\n";
			close OUT3;
			`sh $velveth_shell 1>$velveth_shell.log 2>$velveth_shell.err`;

			my $addnum=int(($Maxk-$Mink)/$Stepk);# caculate No. of times to increase k values.
			my $kvalue=$Mink;
			for(my $i=0; $i<=$addnum; $i++){
				print OUTa "$sample\_$Cover[$j]\_$kvalue\n";
				my $velvetg_shell="$Outdir/$sample/$sample\_$Cover[$j]\_$kvalue\_velvetg_shell";
				open OUT4, ">$velvetg_shell" || die "fail $velvetg_shell";
				print OUT4 "$Bin/velvetg $Outdir/$sample/$sample\_$Cover[$j]\_$kvalue -ins_length $Insert_size -ins_length_sd $Insert_size_SD -cov_cutoff $Cover[$j] -exp_cov $Exp_cover -min_contig_lgth 2000 -scaffolding yes -clean yes\n";
				close OUT4;
				`sh $velvetg_shell 1>$velvetg_shell.log 2>$velvetg_shell.err`;
			
				my $contigs_dir="$Outdir/$sample/$sample\_$Cover[$j]\_$kvalue";	
				unless (-z "$contigs_dir/contigs.fa"){
					`perl $Bin/fasta_shorter.pl $contigs_dir/contigs.fa "$sample\_$Cover[$j]\_$kvalue" > $contigs_dir/contigs.lable.fa`;
					`perl $Bin/align_draw.pl -refer $Outdir/Cp_DB_updated.fa -query $contigs_dir/contigs.lable.fa -fig no -coverlen 800 -coverrate 6 -outdir $contigs_dir`;
					unless (-z "$contigs_dir/Cp_DB_updated.fa_vs_contigs.lable.fa.coords.cp.id"){
						`perl $Bin/choose_refer_contigs_v2.pl $contigs_dir/Cp_DB_updated.fa_vs_contigs.lable.fa.coords.cp.id |$Bin/msort -k nr15 >$contigs_dir/refer_contigs_denovo.stat`;
                          
						open INc, "$contigs_dir/refer_contigs_denovo.stat" || die "$!\n";
                         			my $bestr=<INc>;
                         			chomp($bestr);
                         			close INc;

                         			my @info=split(/\s+/, $bestr);
						my $bestrefer=$info[0]; my $bestcontigs=$info[1];
                         			#`less $contigs_dir/Cp_DB_updated.fa_vs_contigs.lable.fa.coords.cp.id |grep "$bestrefer" >$contigs_dir/refer_contigs_denovo_best.align`;
						`perl $Bin/delete_embeded_alignment.pl $contigs_dir/Cp_DB_updated.fa_vs_contigs.lable.fa.coords.cp.id "$bestrefer" "$bestcontigs" > $contigs_dir/refer_contigs_denovo_best.align`;
						`perl $Bin/contigs_check.pl $contigs_dir/refer_contigs_denovo_best.align $Contigs_len_cutoff > tmp`;

						open IN2, "tmp" || die "$!\n";
						$_=<IN2>;
						`rm tmp`;
						if(/PASS/){
							my $pass_path="$Outdir\/$sample\/denovo_pass_path.log";
							open OUT8, ">$pass_path" || die "$!\n";
							print OUT8  "$Outdir/$sample/$sample\_$Cover[$j]\_$kvalue\t$bestrefer\n";
							close OUT8;
							last lable;
						}
						close IN2;
					}
				}
				$kvalue +=$Stepk;
			}
		}
		close OUTa;
	}

	## run chloe.sh ##
	if($Step=~/5/){
		if (-e "$Outdir\/$sample\/denovo_pass_path.log"){
			`perl $Bin/fishInWinter.pl -bc 2 $Outdir/$sample/denovo_pass_path.log  $Outdir/Cp_DB_updated.fa > $Outdir/$sample/best_refer.fa`;

			open IN3, "$Outdir/$sample/denovo_pass_path.log" || die "$!\n";
			$_=<IN3>;
			chomp;
			my @info=split(/\s+/, $_);	
			my $chloe_shell="$info[0]\_chloe_shell";
			open OUT9, ">$chloe_shell" || die "fail $!\n";
			print OUT9 "$Bin/java -jar $Bin/Chloe.jar $info[0]/contigs.lable.fa $info[0]/Cp_DB_updated.fa_vs_contigs.lable.fa.coords $Outdir/$sample/$sample\_denovo_Chloe.fa $Outdir/$sample/best_refer.fa\n";
			close OUT9;
			`sh $chloe_shell 1>$chloe_shell.log 2>$chloe_shell.err`;
			
			my $bestcontigs=basename($_);
			open OUT7, ">>$Outdir/denovo_pass_id.log" || die "$!\n";
			print OUT7 "$bestcontigs\tPASS\n";
			close OUT7;
			close IN3;
		}else{
			`rm $Outdir/$sample/Cp_DB_updated.fa_vs_all.contigs.lable.fa.coords.cp.id` if (-e "$Outdir/$sample/Cp_DB_updated.fa_vs_all.contigs.lable.fa.coords.cp.id");
			open INa, "$Outdir/$sample/cover_kvalue.list" || die "$!\n";
			while(<INa>){
				chomp;
				if (-e "$Outdir/$sample/$_/Cp_DB_updated.fa_vs_contigs.lable.fa.coords.cp.id"){
					unless (-z "$Outdir/$sample/$_/Cp_DB_updated.fa_vs_contigs.lable.fa.coords.cp.id"){
					`cat $Outdir/$sample/$_/Cp_DB_updated.fa_vs_contigs.lable.fa.coords.cp.id >>$Outdir/$sample/Cp_DB_updated.fa_vs_all.contigs.lable.fa.coords.cp.id`;
					}
				}	
			}
			close INa;
			
			if (-e "$Outdir/$sample/Cp_DB_updated.fa_vs_all.contigs.lable.fa.coords.cp.id"){
				`perl $Bin/choose_refer_contigs_v2.pl $Outdir/$sample/Cp_DB_updated.fa_vs_all.contigs.lable.fa.coords.cp.id |$Bin/msort -k nr15 >$Outdir/$sample/refer_contigs_guided.stat`;
			
				open INb, "$Outdir/$sample/refer_contigs_guided.stat" || die "$!\n";
				my $bestrc=<INb>;
				chomp($bestrc);
				#print "$bestrc\n";
				close INb;
	
				my @info=split(/\s+/, $bestrc);
				my $bestrefer=$info[0]; my $bestcontigs=$info[1];
				#print "$bestrefer\t$bestcontigs\n";
				#`less $Outdir/$sample/Cp_DB_updated.fa\_vs_all.contigs.lable.fa.coords.cp.id |grep "$bestrefer" |grep "$bestcontigs" >$Outdir/$sample/refer_contigs_guided_best.align`;
				`perl $Bin/delete_embeded_alignment.pl $Outdir/$sample/Cp_DB_updated.fa\_vs_all.contigs.lable.fa.coords.cp.id "$bestrefer" "$bestcontigs" >$Outdir/$sample/refer_contigs_guided_best.align`;
				`perl $Bin/link_contigs_v3.pl $Outdir/$sample/$bestcontigs/contigs.lable.fa $Outdir/$sample/refer_contigs_guided_best.align Cp_DB $bestcontigs $Outdir/$sample`;
				`perl $Bin/bwa_perl.pl -ref $Outdir/$sample/$bestcontigs\_linked_cp.fa -pelist $Outdir/$sample/$sample.notCombined_1.fastq.gz,$Outdir/$sample/$sample.notCombined_2.fastq.gz -selist $Outdir/$sample/$sample.extendedFrags.fastq.gz,$Outdir/$sample/spades_outdir/corrected/$sample\_trimmed.00.0_0.cor.fastq.gz -sample $bestcontigs\_1st -outdir $Outdir/$sample`;
				`sh $Outdir/$sample/$bestcontigs\_1st_bwa.sh 1>$Outdir/$sample/$bestcontigs\_1st_bwa.sh.log 2>$Outdir/$sample/$bestcontigs\_1st_bwa.sh.err`;	
				`java -jar $Bin/pilon-1.16.jar --genome $Outdir/$sample/$bestcontigs\_linked_cp.fa --frags $Outdir/$sample/$bestcontigs\_1st.pe0.sort.bam --unpaired $Outdir/$sample/$bestcontigs\_1st.se0.sort.bam --unpaired $Outdir/$sample/$bestcontigs\_1st.se1.sort.bam --output $Outdir/$sample/$sample\_guided_1st_pilon 1>$Outdir/$sample/$sample\_guided_1st_pilon.log 2>$Outdir/$sample/$sample\_guided_1st_pilon.err`; 
      
                        	open OUT6, ">>$Outdir/guided_pass_id.log" || die "$!\n";
                        	print OUT6 "$bestcontigs\t$bestrefer\tPASS\n";
                        	close OUT6;
			}else{
				open OUT5, ">>$Outdir/both_fail_id.log" || die "$!\n";
				print OUT5 "$sample\tFAIL\n";
				close OUT5;			
			}
		}
	}
	
	## run gapfiller ##
	if ($Step=~/6/){
		my $final_fa; my $after_fill_gap;
		$final_fa="$Outdir/$sample/$sample\_denovo_Chloe.fa" if(-e "$Outdir/$sample/$sample\_denovo_Chloe.fa");
		$final_fa="$Outdir/$sample/$sample\_guided_1st_pilon.fasta" if (-e "$Outdir\/$sample\/$sample\_guided_1st_pilon.fasta");
		unless (-z "$Outdir/$sample/$sample\_denovo_Chloe.fa"){
		if (-e $final_fa){
			`perl $Bin/delete_N_ends.pl $final_fa > $Outdir/$sample/fasta_no_N_end.fa`;
			$final_fa="$Outdir/$sample/fasta_no_N_end.fa";
			`perl $Bin/fa_quality.pl -Head -len -gap -N -gc $final_fa >$final_fa.stat`;
			open INd,"$final_fa.stat" || die "$!\n";
			<INd>;
			$_=<INd>;
			my @info=split(/\s+/, $_);
			close INd;
			
			if ($info[2] !=0){			
				open OUT12,">$Outdir/$sample/lib.txt" || die "$!\n";
				print OUT12 "$sample\_PE bwa $Outdir/$sample/$sample\_bbnorm.1.fastq.00.0_0.cor.fastq.gz $Outdir/$sample/$sample\_bbnorm.2.fastq.00.0_0.cor.fastq.gz $Insert_size 0.25 FR\n";
 	        		close OUT12;
                		open OUT13,">$Outdir/$sample/gapfiller.sh" || die "$!\n";
                		print OUT13 "cd $Outdir/$sample;\nperl /home/xzhong/Software/GapFiller/GapFiller.pl -l $Outdir/$sample/lib.txt -s $final_fa -b gap_filling;\ncd $rundir;\n";
                		close OUT13;

                		`sh $Outdir/$sample/gapfiller.sh 1>$Outdir/$sample/gapfiller.sh.log 2>$Outdir/$sample/gapfiller.sh.err`;
				$after_fill_gap="$Outdir/$sample/gap_filling/gap_filling.gapfilled.final.fa";
			}else{
				$after_fill_gap=$final_fa;
			}
                `perl $Bin/end_check.pl $after_fill_gap >$Outdir/$sample/after_fill_gap_check.fa`;
		}
		}
	}

	## run bwa ##
	if($Step=~/7/){
		if (-e "$Outdir/$sample/$sample\_denovo_Chloe.fa" || -e "$Outdir\/$sample\/$sample\_guided_1st_pilon.fasta"){
		unless (-z "$Outdir/$sample/$sample\_denovo_Chloe.fa"){
		#if (-e "$Outdir/$sample/$sample\_denovo_Chloe.fa"){
		my $final_fa="$Outdir/$sample/after_fill_gap_check.fa";

		`perl $Bin/bwa_perl.pl -ref $final_fa -pelist $Outdir/$sample/$sample\_bbnorm.1.fastq.00.0_0.cor.fastq.gz,$Outdir/$sample/$sample\_bbnorm.2.fastq.00.0_0.cor.fastq.gz -sample $sample -outdir $Outdir/$sample`;
		`sh $Outdir/$sample/$sample\_bwa.sh 1>$Outdir/$sample/$sample\_bwa.sh.log 2>$Outdir/$sample/$sample\_bwa.sh.err`;
		}
		}
	}	

	## run pilon ##	
	if($Step=~/8/){
		if (-e "$Outdir/$sample/$sample\_denovo_Chloe.fa" || -e "$Outdir\/$sample\/$sample\_guided_1st_pilon.fasta"){
		unless (-z "$Outdir/$sample/$sample\_denovo_Chloe.fa"){
		#if (-e "$Outdir/$sample/$sample\_denovo_Chloe.fa"){
		my $final_fa="$Outdir/$sample/after_fill_gap_check.fa";

		open OUT11, ">$Outdir/$sample/$sample\_pilon.sh" || die "$!\n";
		print OUT11 " $Bin/java -jar $Bin/pilon-1.16.jar --genome $final_fa --frags $Outdir/$sample/$sample.pe0.sort.bam --output $Outdir/$sample/$sample\_pilon\n";
		close OUT11;
		`sh $Outdir/$sample/$sample\_pilon.sh 1>$Outdir/$sample/$sample\_pilon.sh.log 2>$Outdir/$sample/$sample\_pilon.sh.err`;

		`perl $Bin/end_restore.pl $Outdir/$sample/$sample\_pilon.fasta >$Outdir/$sample/$sample\_pilon\_final.fasta`;
		#open OUTcp, ">$Outdir/$sample/$sample\_permutation.sh" || die "$!\n";
		#print OUTcp "$Bin/java -jar $Bin/CircularPermutation.jar $Outdir/$sample/$sample\_pilon\_restore.fasta $Bin/NC_030504_refs.fa $Outdir/$sample/$sample\_pilon\_permutation.fasta\n";
		#close OUTcp;
		#`sh $Outdir/$sample/$sample\_permutation.sh 1>$Outdir/$sample/$sample\_permutation.sh.log 2>$Outdir/$sample/$sample\_permutation.sh.err`;

		#if (-e "$Outdir/$sample/$sample\_pilon\_permutation.fasta"){
		#	`ln -s $Outdir/$sample/$sample\_pilon\_permutation.fasta $Outdir/$sample/$sample\_pilon\_final.fasta`;
		#}else{
		#	`ln -s $Outdir/$sample/$sample\_pilon\_restore.fasta $Outdir/$sample/$sample\_pilon\_final.fasta`;
		#}

		`perl $Bin/fa_quality.pl -Head -len -gap -N -gc $Outdir/$sample/$sample\_pilon\_final.fasta > $Outdir/$sample/$sample\_pilon\_final\_fasta.stat\n`;
		}
		
		#if (-e "$Outdir/$sample/$sample\_denovo_Chloe.fa"){
		#	`cat $Outdir/Cp_DB_updated.fa $Outdir/$sample/$sample.pilon.final.fasta >$Outdir/Cp_DB_updated.fa.tmp`;
		#	`mv $Outdir/Cp_DB_updated.fa.tmp $Outdir/Cp_DB_updated.fa`;
		#}
		}
	}
	
	## run statistics ##		
	if ($Step=~/9/){
		#`rm $Outdir/reads.stat` if (-e "$Outdir/reads.stat");
		`perl $Bin/reads_stat.pl $Outdir/$sample/$sample.bbnorm.shell.err $Outdir/$sample/$sample.cutadapter.shell.log $sample >>$Outdir/reads.stat`;
		
		if (-e "$Outdir/$sample/$sample\_denovo_Chloe.fa"){
			`perl $Bin/assembly_stat.pl $Outdir/$sample/$sample.pilon.final.fasta.stat $Outdir/$sample/denovo_pass_path.log $sample\_pilon.sh.log $sample denovo >> $Outdir/denovo_assembly.stat`;
		}

		if (-e "$Outdir\/$sample\/$sample\_guided_1st_pilon.fasta"){
			`perl $Bin/assembly_stat2.pl $Outdir/$sample/$sample.pilon.final.fasta.stat $Outdir/$sample/refer_contigs_guided.stat $sample\_pilon.sh.log $sample reference-guided >> $Outdir/refer_assembly.stat`;
		}
	}
	
	## remove the temporary ##	
	if ($Clean eq 'yes'){
		`rm $Outdir/$sample/$sample\_bbnorm*.gz`;
		`rm $Outdir/$sample/$sample\_trimmed*.gz`;
		`rm $Outdir/$sample/spades_outdir/corrected/$sample\_trimmed.*.fastq.*.fastq.gz`;
		`rm $Outdir/$sample/$sample\_*/Sequences`;
		`rm $Outdir/$sample/$sample\_*/Roadmaps`;
		`rm $Outdir/$sample/$sample\_*/Graph2`;
		`rm $Outdir/$sample/$sample\_*/PreGraph`;
		`rm -rf $Outdir/$sample/gap_filling/reads`;
		`rm -rf $Outdir/$sample/gap_filling/intermediate_results`;	
		`rm $Outdir/$sample/*_1st.*.sort.bam`;	
	}
}
close IN;
	

