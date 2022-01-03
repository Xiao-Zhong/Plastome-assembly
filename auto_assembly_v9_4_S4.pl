#!/usr/bin/perl

=head1 Usage:

perl auto_assembly.pl [options] <samples_list> <reads_path>

=head1 Example:
    perl auto_assembly.pl -outdir output -step abcdefghi samples_list /home/xzhong/working_data_01/01_Assembly
########## description ##########
    step a: remove adapter sequences using cutadapt.sh;
    step b: normalize read depth based on kmer counts using bbnorm.sh;
    step c: correct read errors using spades.sh;
    step d: merge paired-end reads into single reads by overlap detection using flash.sh or bbmerge.sh[default];
    step e: assemble using Velvet;
    setp f: link assembled contigs to be a single long contig using Chloe2.jar or *.pl;
    step g: fill gaps using Gapfiller;
    step h: map reads back to the single contig  using BWA;
    step i: check, improve and report the assembly quality using Pilon;
    step j: statistics;
    step k: remove temporary files.
#################################
=cut

use strict;
use Getopt::Long;
use File::Basename;
#use FinBin qw($Bin);
my $Bin="/home/xzhong/working_data_01/01_Assembly/37_pipeline";

my ($Step, $Outdir, $Help);
my ($Maxinputreads, $Maxoutputreads, $BBnorm_mem);
my ($Read_len, $Adapter1, $Adapter2);
my $Merge;
my ($Mink, $Maxk, $Stepk, $Coverlist,$Refer_DB,$Insert_size,$Insert_size_SD,$Exp_cover,$Contigs_len_range);
my $Clean;

GetOptions(
	"outdir:s"=>\$Outdir,
	"step:s"=>\$Step,

	"maxinputreads:i"=>\$Maxinputreads,
	"maxoutputreads:i"=>\$Maxoutputreads,
	"bbnorm_mem:i"=>\$BBnorm_mem,

	"readlength:i"=>\$Read_len,
	"adapter1:s"=>\$Adapter1, # adapter1: Sequence of an adapter that was ligated to the 3' end.
	"adapter2:s"=>\$Adapter2, # adapter2: 3' adapter to be removed from second read in a pair.

	"merge:s"=>\$Merge,

	"mink:i"=>\$Mink,
	"maxk:i"=>\$Maxk,
	"stepk:i"=>\$Stepk, # Velvet will hash from k=m to k=M with a step of s.
	"coverlist:s"=>\$Coverlist, # Velvet's -cov_cutoff.
	"insertsize:i"=>\$Insert_size,
	"insertsizesd:i"=>\$Insert_size_SD,
	"expectedcover:i"=>\$Exp_cover,
	"referdb:s"=>\$Refer_DB,
	"contigslen:s"=>\$Contigs_len_range, # LSG,SSG and IR 's length range(Kb).
	
	#"clean:s"=>\$Clean,	
	"help!"=>\$Help
);

$Outdir ||=".";
$Step ||='abcdefghijk';
$Maxinputreads ||=-1;
$Maxoutputreads ||=-1;
$BBnorm_mem ||=40; # (Gb)
$Read_len ||=124;
$Adapter1 ||="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
$Adapter2 ||="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";
$Merge ||="bbmerge";
$Mink ||=51;
$Maxk ||=115;
$Stepk ||=20;
$Refer_DB ||="$Bin/Cp_DB.fa";
$Insert_size ||=350;
$Insert_size_SD ||=$Insert_size*1/10;
$Exp_cover ||=30;
$Coverlist ||="10,15,7,20";# print "$Coverlist\n";
my @Cover=split(/\,/,$Coverlist);
$Contigs_len_range ||="74,92,10,22,18,29";# "Kb", LSC/SSC/IR
#$Clean ||="no";
die `pod2text $0` if (@ARGV == 0 || $Help);

mkdir($Outdir) unless(-d $Outdir);
my $rundir=`pwd`;
chomp $rundir;#print "$rundir";
$Outdir="$rundir/$Outdir"; #print "$Outdir\n";
`cp $Refer_DB $Outdir/Cp_DB_updated.fa`;

my $sample_list=shift;
my $reads_dir=shift;
#`cp $Refer_DB $Outdir/$sample_list\_Cp_DB_updated.fa`;

open IN, "$sample_list" || die "fail $sample_list";
while(<IN>){
	chomp;
	my $sample="$_";
	#print "$sample\n";
	mkdir("$Outdir/$sample") unless(-d "$Outdir/$sample");
	
	## run cutadapt.sh ##
	if ($Step=~/a/){
		`ls $reads_dir/$sample*\_R1*.fastq.gz $reads_dir/$sample*\_R2*.fastq.gz > $sample\_tmp`;
		open IN_reads, "$sample\_tmp" || die "$!\n";
		my $input1=<IN_reads>; chomp $input1;
		my $input2=<IN_reads>; chomp $input2;
		#print "$input1\t$input2\n";
		close IN_reads;
	
		my $cutadapter_shell="$Outdir/$sample/$sample\_cutadapter.sh";
		open OUTca, ">$cutadapter_shell" || die "fail $cutadapter_shell";
		print OUTca "/dd_groupdata/xzhong/tools/cutadapt-1.9.1/bin/cutadapt -a $Adapter1 -A $Adapter2 -o $Outdir/$sample/$sample\_trimmed.1.fastq.gz -p $Outdir/$sample/$sample\_trimmed.2.fastq.gz $input1 $input2\n";
		close OUTca;

		`rm $sample\_tmp`;
		`sh $cutadapter_shell 1>$cutadapter_shell.log 2>$cutadapter_shell.err`;
	}

	## run bbnorm.sh ##
	 if ($Step=~/b/){		
		my $bbnorm_shell="$Outdir/$sample/$sample\_bbnorm.sh";
		open OUTbn, ">$bbnorm_shell" || die "fail $bbnorm_shell";
		print OUTbn "$Bin/bbnorm.sh -Xmx$BBnorm_mem\g in=$Outdir/$sample/$sample\_trimmed.1.fastq.gz in2=$Outdir/$sample/$sample\_trimmed.2.fastq.gz tablereads=$Maxinputreads reads=$Maxoutputreads zerobin=t prefilter=t maxdepth=1000 lowbindepth=10 highbindepth=500 hist=$Outdir/$sample/$sample\_bbnorm.hist out=$Outdir/$sample/$sample\_bbnorm.1.fastq.gz out2=$Outdir/$sample/$sample\_bbnorm.2.fastq.gz\n";
		close OUTbn;
		`sh $bbnorm_shell 1>$bbnorm_shell.log 2>$bbnorm_shell.err`;
	}

	## run spades.sh ##
	if ($Step=~/c/){
		my $spades_shell="$Outdir/$sample/$sample\_spades.sh";
		open OUTsd, ">$spades_shell" || die "fail $spades_shell";
		print OUTsd "$Bin/python $Bin/spades.py --only-error-correction --careful --pe1-1 $Outdir/$sample/$sample\_bbnorm.1.fastq.gz --pe1-2 $Outdir/$sample/$sample\_bbnorm.2.fastq.gz -o $Outdir/$sample/spades_outdir\n";
		close OUTsd;
		`sh $spades_shell 1>$spades_shell.log 2>$spades_shell.err`;		
	}

	## run flash.sh or bbmerge.sh[default] ##
	if($Step=~/d/ && ($Merge ne "bbmerge")){
		my $flash_shell="$Outdir/$sample/$sample\_flash.sh";
		open OUTfs, ">$flash_shell" || die "fail $flash_shell";
		print OUTfs "$Bin/flash -q -m 10 -x 0.2 --max-overlap 125 -z -d $Outdir/$sample -o $sample $Outdir/$sample/spades_outdir/corrected/$sample\_bbnorm.1.fastq.00.0_0.cor.fastq.gz $Outdir/$sample/spades_outdir/corrected/$sample\_bbnorm.2.fastq.00.0_0.cor.fastq.gz\n";
		close OUTfs;
		`sh $flash_shell 1>$flash_shell.log 2>$flash_shell.err`;
	}	
	if($Step=~/d/ && ($Merge eq "bbmerge")){
		my $bbmerge_shell="$Outdir/$sample/$sample\_bbmerge.sh";
		open OUTbm, ">$bbmerge_shell" || die "fail $bbmerge_shell";
		print OUTbm "$Bin/bbmerge.sh in1=$Outdir/$sample/spades_outdir/corrected/$sample\_bbnorm.1.fastq.00.0_0.cor.fastq.gz in2=$Outdir/$sample/spades_outdir/corrected/$sample\_bbnorm.2.fastq.00.0_0.cor.fastq.gz out=$Outdir/$sample/$sample.extendedFrags.fastq.gz outu1=$Outdir/$sample/$sample.notCombined_1.fastq.gz outu2=$Outdir/$sample/$sample.notCombined_2.fastq.gz\n";
		close OUTbm;
		`sh $bbmerge_shell 1>$bbmerge_shell.log 2>$bbmerge_shell.err`;
	}

	## run velveth.sh & velvetg.sh & nucmer.sh & contigs_check.pl ##
	if($Step=~/e/){
		open OUTck, ">$Outdir/$sample/cover_kvalue.list" || die "$!\n";
		lable: for(my $j=0; $j<@Cover; $j++){
			#print "$Cover[$j]\n";
			my $velveth_shell="$Outdir/$sample/$sample\_$Cover[$j]\_velveth_shell";
			open OUTvh, ">$velveth_shell" || die "fail $velveth_shell";
			print OUTvh "$Bin/velveth $Outdir/$sample/$sample\_$Cover[$j] $Mink,$Maxk,$Stepk -shortPaired -fastq.gz -separate $Outdir/$sample/$sample.notCombined_1.fastq.gz $Outdir/$sample/$sample.notCombined_2.fastq.gz -long $Outdir/$sample/$sample.extendedFrags.fastq.gz -short $Outdir/$sample/spades_outdir/corrected/$sample\_bbnorm.00.0_0.cor.fastq.gz\n";
			close OUTvh;
			`sh $velveth_shell 1>$velveth_shell.log 2>$velveth_shell.err`;

			my $addnum=int(($Maxk-$Mink)/$Stepk);# caculate No. of times to increase k values.
			my $kvalue=$Mink;
			for(my $i=0; $i<=$addnum; $i++){
				print OUTck "$sample\_$Cover[$j]\_$kvalue\n";
				my $velvetg_shell="$Outdir/$sample/$sample\_$Cover[$j]\_$kvalue\_velvetg_shell";
				open OUTvg, ">$velvetg_shell" || die "fail $velvetg_shell";
				print OUTvg "$Bin/velvetg $Outdir/$sample/$sample\_$Cover[$j]\_$kvalue -ins_length $Insert_size -ins_length_sd $Insert_size_SD -cov_cutoff $Cover[$j] -exp_cov $Exp_cover -min_contig_lgth 2000 -scaffolding yes -clean yes\n";
				close OUTvg;
				`sh $velvetg_shell 1>$velvetg_shell.log 2>$velvetg_shell.err`;
			
				my $contigs_dir="$Outdir/$sample/$sample\_$Cover[$j]\_$kvalue";	
				unless (-z "$contigs_dir/contigs.fa"){
					`perl $Bin/add_label_to_fasta.pl $contigs_dir/contigs.fa "$sample\_$Cover[$j]\_$kvalue" > $contigs_dir/contigs.lable.fa`;
					`perl $Bin/align_draw_for_pipe.pl -refer $Outdir/Cp_DB_updated.fa -query $contigs_dir/contigs.lable.fa -coverlen 800 -coverrate 6 -outdir $contigs_dir`;
					unless (-z "$contigs_dir/Cp_DB_updated.fa_vs_contigs.lable.fa.coords.cp"){
						`perl $Bin/choose_refer_contigs_v2.pl $contigs_dir/Cp_DB_updated.fa_vs_contigs.lable.fa.coords.cp |$Bin/msort -k nr15 >$contigs_dir/refer_contigs_denovo.stat`;
                          
						open INds, "$contigs_dir/refer_contigs_denovo.stat" || die "$!\n";
                         			my $bestr=<INds>;
                         			chomp($bestr);
                         			close INds;

                         			my @info=split(/\s+/, $bestr);
						my $bestrefer=$info[0]; my $bestcontigs=$info[1];
						`perl $Bin/delete_embeded_alignment.pl $contigs_dir/Cp_DB_updated.fa_vs_contigs.lable.fa.coords.cp "$bestrefer" "$bestcontigs" > $contigs_dir/refer_contigs_denovo_best.align`;
						`perl $Bin/contigs_check.pl $contigs_dir/refer_contigs_denovo_best.align $Contigs_len_range > tmp`;

						open INtmp, "tmp" || die "$!\n";
						$_=<INtmp>;
						`rm tmp`;
						if(/PASS/){
							my $pass_path="$Outdir\/$sample\/denovo_pass_path.log";
							open OUTdp, ">$pass_path" || die "$!\n";
							print OUTdp  "$Outdir/$sample/$sample\_$Cover[$j]\_$kvalue\t$bestrefer\n";
							close OUTdp;
							last lable;
						}
						close INtmp;
					}
				}
				$kvalue +=$Stepk;
			}
		}
		close OUTck;
	}

	## run chloe.sh ##
	if($Step=~/f/){
		if (-e "$Outdir\/$sample\/denovo_pass_path.log"){
			#`perl $Bin/fishInWinter.pl -bc 2 $Outdir/$sample/denovo_pass_path.log  $Outdir/Cp_DB_updated.fa > $Outdir/$sample/best_refer.fa`;

			open INdp, "$Outdir/$sample/denovo_pass_path.log" || die "$!\n";
			$_=<INdp>;
			chomp;
			my @info=split(/\s+/, $_);	
			my $chloe_shell="$info[0]\_chloe_shell";
			open OUTce, ">$chloe_shell" || die "fail $!\n";
			print OUTce "$Bin/java -jar $Bin/Chloe2.jar $info[0]/contigs.lable.fa $info[0]/refer_contigs_denovo_best.align $Outdir/$sample/$sample\_denovo_Chloe.fa\n";
			close OUTce;
			`sh $chloe_shell 1>$chloe_shell.log 2>$chloe_shell.err`;
			
			my $bestcontigs=basename($_);
			open OUTdi, ">>$Outdir/denovo_pass_id.log" || die "$!\n";
			print OUTdi "$bestcontigs\tPASS\n";
			close OUTdi;
			close INdp;
		}else{
			
			if (-e "$Outdir/$sample/Cp_DB_updated.fa_vs_all.contigs.lable.fa.coords.cp"){
				`perl $Bin/choose_refer_contigs_v2.pl $Outdir/$sample/Cp_DB_updated.fa_vs_all.contigs.lable.fa.coords.cp |$Bin/msort -k nr15 >$Outdir/$sample/refer_contigs_guided.stat`;
			
				open INgs, "$Outdir/$sample/refer_contigs_guided.stat" || die "$!\n";
				my $bestrc=<INgs>;
				chomp($bestrc);
				#print "$bestrc\n";
				close INgs;
	
				my @info=split(/\s+/, $bestrc);
				my $bestrefer=$info[0]; my $bestcontigs=$info[1];
				`perl $Bin/delete_embeded_alignment.pl $Outdir/$sample/Cp_DB_updated.fa\_vs_all.contigs.lable.fa.coords.cp "$bestrefer" "$bestcontigs" >$Outdir/$sample/refer_contigs_guided_best.align`;
				`perl $Bin/link_contigs_v4_tmp.pl $Outdir/$sample/$bestcontigs/contigs.lable.fa $Outdir/$sample/refer_contigs_guided_best.align Cp_DB $bestcontigs $Outdir/$sample`;
				`perl $Bin/delete_overlap.pl $Outdir/$sample/$bestcontigs\_linked_cp.fa >$Outdir/$sample/$bestcontigs\_linked_cp_delete.fa`;
				
				`perl $Bin/bwa_perl.pl -ref $Outdir/$sample/$bestcontigs\_linked_cp_delete.fa -pelist $Outdir/$sample/$sample.notCombined_1.fastq.gz,$Outdir/$sample/$sample.notCombined_2.fastq.gz -selist $Outdir/$sample/$sample.extendedFrags.fastq.gz,$Outdir/$sample/spades_outdir/corrected/$sample\_bbnorm.00.0_0.cor.fastq.gz -sample $bestcontigs\_1st -outdir $Outdir/$sample`;
				#`perl $Bin/bwa_perl.pl -ref $Outdir/$sample/$bestcontigs\_linked_cp_delete.fa -pelist $Outdir/$sample/spades_outdir/corrected/$sample\_bbnorm.1.fastq.00.0_0.cor.fastq.gz,$Outdir/$sample/spades_outdir/corrected/$sample\_bbnorm.2.fastq.00.0_0.cor.fastq.gz -sample $bestcontigs\_1st -outdir $Outdir/$sample`;
				`sh $Outdir/$sample/$bestcontigs\_1st_bwa.sh 1>$Outdir/$sample/$bestcontigs\_1st_bwa.sh.log 2>$Outdir/$sample/$bestcontigs\_1st_bwa.sh.err`;

				open OUTgp, ">$Outdir/$sample/$sample\_guided_1st_pilon.sh" || die "$!\n";
				print OUTgp "$Bin/java -jar $Bin/pilon-1.22.jar --genome $Outdir/$sample/$bestcontigs\_linked_cp_delete.fa --frags $Outdir/$sample/$bestcontigs\_1st.pe0.sort.bam --unpaired $Outdir/$sample/$bestcontigs\_1st.se0.sort.bam --unpaired $Outdir/$sample/$bestcontigs\_1st.se1.sort.bam --output $Outdir/$sample/$sample\_guided_1st_pilon\n";
				#print OUTgp "$Bin/java -jar $Bin/pilon-1.22.jar --genome $Outdir/$sample/$bestcontigs\_linked_cp_delete.fa --frags $Outdir/$sample/$bestcontigs\_1st.pe0.sort.bam --output $Outdir/$sample/$sample\_guided_1st_pilon\n"; 
      				close OUTgp;
				`sh $Outdir/$sample/$sample\_guided_1st_pilon.sh 1>$Outdir/$sample/$sample\_guided_1st_pilon.sh.log 2>$Outdir/$sample/$sample\_guided_1st_pilon.sh.err`;

                        	open OUTgi, ">>$Outdir/guided_pass_id.log" || die "$!\n";
                        	print OUTgi "$bestcontigs\t$bestrefer\tPASS\n";
                        	close OUTgi;
			}else{
				open OUTbi, ">>$Outdir/both_fail_id.log" || die "$!\n";
				print OUTbi "$sample\tFAIL\n";
				close OUTbi;			
			}
		}
	}
	
	## run gapfiller ##
	if ($Step=~/g/){
		my $final_fa; my $after_fill_gap;
		if(-e "$Outdir/$sample/$sample\_denovo_Chloe.fa"){
			$final_fa="$Outdir/$sample/$sample\_denovo_Chloe.fa";
		}elsif (-e "$Outdir\/$sample\/$sample\_guided_1st_pilon.fasta"){
			$final_fa="$Outdir/$sample/$sample\_guided_1st_pilon.fasta";
		}

		if (-e $final_fa){
			`perl $Bin/fa_quality.pl -Head -len -gap -N -gc $final_fa >$final_fa.stat`;
			open INfs,"$final_fa.stat" || die "$!\n";
			<INfs>;
			$_=<INfs>;
			my @info=split(/\s+/, $_);
			close INfs;
			
			if ($info[2] !=0){			
				open OUTlt,">$Outdir/$sample/lib.txt" || die "$!\n";
				print OUTlt "$sample\_PE bwa $Outdir/$sample/4_S4_L001_R1_001.fastq.gz $Outdir/$sample/4_S4_L001_R2_001.fastq.gz $Insert_size 0.25 FR\n";
 	        		close OUTlt;

                		open OUTgs,">$Outdir/$sample/gapfiller.sh" || die "$!\n";
                		print OUTgs "cd $Outdir/$sample;\n/usr/bin/perl /home/xzhong/Software/GapFiller/GapFiller.pl -l $Outdir/$sample/lib.txt -s $final_fa -b gap_filling;\ncd $rundir;\n";
                		close OUTgs;

                		`sh $Outdir/$sample/gapfiller.sh 1>$Outdir/$sample/gapfiller.sh.log 2>$Outdir/$sample/gapfiller.sh.err`;
				$after_fill_gap="$Outdir/$sample/gap_filling/gap_filling.gapfilled.final.fa";
			}else{
				$after_fill_gap=$final_fa;
			}
                `perl $Bin/end_check.pl $after_fill_gap $Read_len >$Outdir/$sample/after_fill_gap_check.fa`;
		}
	}

	## run bwa ##
	if($Step=~/h/){
		#if (-e "$Outdir/$sample/$sample\_denovo_Chloe.fa" || -e "$Outdir\/$sample\/$sample\_guided_1st_pilon.fasta"){
		if (-e "$Outdir/$sample/after_fill_gap_check.fa"){
			my $final_fa="$Outdir/$sample/after_fill_gap_check.fa";

			`perl $Bin/bwa_perl.pl -ref $final_fa -pelist $Outdir/$sample/$sample.notCombined_1.fastq.gz,$Outdir/$sample/$sample.notCombined_2.fastq.gz -selist $Outdir/$sample/$sample.extendedFrags.fastq.gz,$Outdir/$sample/spades_outdir/corrected/$sample\_bbnorm.00.0_0.cor.fastq.gz -sample $sample -outdir $Outdir/$sample`;
			#`perl $Bin/bwa_perl.pl -ref $final_fa -pelist $Outdir/$sample/spades_outdir/corrected/$sample\_bbnorm.1.fastq.00.0_0.cor.fastq.gz,$Outdir/$sample/spades_outdir/corrected/$sample\_bbnorm.2.fastq.00.0_0.cor.fastq.gz -sample $sample -outdir $Outdir/$sample`;
			`sh $Outdir/$sample/$sample\_bwa.sh 1>$Outdir/$sample/$sample\_bwa.sh.log 2>$Outdir/$sample/$sample\_bwa.sh.err`;
		}
	}	

	## run pilon ##	
	if($Step=~/i/){
		#if (-e "$Outdir/$sample/$sample\_denovo_Chloe.fa" || -e "$Outdir\/$sample\/$sample\_guided_1st_pilon.fasta"){
		if (-e "$Outdir/$sample/after_fill_gap_check.fa"){
			my $final_fa="$Outdir/$sample/after_fill_gap_check.fa";

			open OUTpl, ">$Outdir/$sample/$sample\_pilon.sh" || die "$!\n";
			print OUTpl "$Bin/java -jar $Bin/pilon-1.22.jar --genome $final_fa --frags $Outdir/$sample/$sample.pe0.sort.bam --unpaired $Outdir/$sample/$sample.se0.sort.bam --unpaired $Outdir/$sample/$sample.se1.sort.bam --output $Outdir/$sample/$sample\_pilon\n";
			#print OUTpl "$Bin/java -jar $Bin/pilon-1.22.jar --genome $final_fa --frags $Outdir/$sample/$sample.pe0.sort.bam --output $Outdir/$sample/$sample\_pilon\n";
			close OUTpl;
			`sh $Outdir/$sample/$sample\_pilon.sh 1>$Outdir/$sample/$sample\_pilon.sh.log 2>$Outdir/$sample/$sample\_pilon.sh.err`;

			`perl $Bin/end_restore.pl $Outdir/$sample/$sample\_pilon.fasta >$Outdir/$sample/$sample\_pilon\_restore.fasta`;

			`rm $Outdir/$sample/$sample\_pilon\_final.fasta` if (-e "$Outdir/$sample/$sample\_pilon\_final.fasta");
			`rm $Outdir/$sample/$sample\_pilon\_permutation.fasta ` if (-e "$Outdir/$sample/$sample\_pilon\_permutation.fasta");
			open OUTcp, ">$Outdir/$sample/$sample\_permutation.sh" || die "$!\n";
			print OUTcp "$Bin/java -jar $Bin/CircularPermutation.jar $Outdir/$sample/$sample\_pilon\_restore.fasta $Bin/NC_030504_refs.fa $Outdir/$sample/$sample\_pilon\_permutation.fasta\n";
			close OUTcp;
			`sh $Outdir/$sample/$sample\_permutation.sh 1>$Outdir/$sample/$sample\_permutation.sh.log 2>$Outdir/$sample/$sample\_permutation.sh.err`;
			
			#`rm $Outdir/$sample/$sample\_pilon\_final.fasta` if (-e "$Outdir/$sample/$sample\_pilon\_final.fasta");			
			if (-e "$Outdir/$sample/$sample\_pilon\_permutation.fasta"){
				`ln -s $Outdir/$sample/$sample\_pilon\_permutation.fasta $Outdir/$sample/$sample\_pilon\_final.fasta`;
			}else{
				`ln -s $Outdir/$sample/$sample\_pilon\_restore.fasta $Outdir/$sample/$sample\_pilon\_final.fasta`;
			}

			`perl $Bin/fa_quality.pl -Head -len -gap -N -gc $Outdir/$sample/$sample\_pilon\_final.fasta > $Outdir/$sample/$sample\_pilon\_final\_fasta.stat\n`;
			}
		
		if (-e "$Outdir/$sample/$sample\_denovo_Chloe.fa"){
			open INff, "$Outdir/$sample/$sample\_pilon_final_fasta.stat" || die "$!\n";
			<INff>;
			my $tmp=<INff>;
			my @info=split(/\t/, $tmp);
			if ($info[1]<200000){
				`cat $Outdir/Cp_DB_updated.fa $Outdir/$sample/$sample\_pilon\_final.fasta >$Outdir/Cp_DB_updated.fa.tmp`;
				`mv $Outdir/Cp_DB_updated.fa.tmp $Outdir/Cp_DB_updated.fa`;
			}
			close INff;
		}
	}
	
	## run statistics ##		
	if ($Step=~/j/){
		#`rm $Outdir/reads.stat` if (-e "$Outdir/reads.stat");
		`perl $Bin/reads_stat.pl $Outdir/$sample/$sample.bbnorm.shell.err $Outdir/$sample/$sample.cutadapter.shell.log $sample >>$Outdir/reads.stat`;
		
		if (-e "$Outdir/$sample/$sample\_denovo_Chloe.fa"){
			`perl $Bin/assembly_stat.pl $Outdir/$sample/$sample\_pilon\_final\_fasta.stat $Outdir/$sample/denovo_pass_path.log $sample\_pilon.sh.log $sample denovo >> $Outdir/denovo_assembly.stat`;
		}

		if (-e "$Outdir\/$sample\/$sample\_guided_1st_pilon.fasta"){
			`perl $Bin/assembly_stat2.pl $Outdir/$sample/$sample\_pilon\_final\_fasta.stat $Outdir/$sample/refer_contigs_guided.stat $sample\_pilon.sh.log $sample reference-guided >> $Outdir/refer_assembly.stat`;
		}
	}
	
	## remove the temporary ##	
	if ($Step=~/k/){
		#`rm $Outdir/$sample/$sample\_bbnorm*.gz`;
		`rm $Outdir/$sample/$sample\_trimmed*.gz`;
		#`rm $Outdir/$sample/spades_outdir/corrected/$sample\_trimmed.*.fastq.*.fastq.gz`;
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
	

