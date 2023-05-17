# Plastome-assembly
The pipeline can assemble a (near) complete plastome genome from short reads automatically. It condains both denovo and reference-based assembly, and is useful for large-scale chloroplast genome sequencing projects. For a single or few samples, please try NOVOPlasty, GetOrganelle or others.

Usage:
    perl auto_assembly.pl [options] <samples_list> <reads_path>

Example:
    perl auto_assembly.pl -outdir output -step abcdefghijk samples_list /home/xzhong/working_data_01/01_Assembly
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

The workflow is lack of maintenance nowadays...
