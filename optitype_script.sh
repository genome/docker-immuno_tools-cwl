#!/bin/bash

# $1 is the directory where all intermediate files will be written
# $2 is the prefix that will be added to all files
# $3 is the location of the input (normal) cram file
# $4 is the directory final outputs should be written to
#bsub -R 'select[mem>64000] rusage[mem=64000]' -e <error_file_name> -o <output_file_name> -q research-hpc -a 'docker(johnegarza/immuno-testing:latest)' /bin/bash /usr/bin/optitype_script.sh <intermediate files directory> <final results directory> <output file prefix>  <cram path>

# Optitype DNA reference file
dnaref="/ref_data/optitype_ref/hla_reference_dna.fasta";

TEMPDIR=$1
outdir=$2;
name=$3;
cram=$4;
mkdir -p $outdir

echo Converting cram to bam
/opt/samtools/bin/samtools view -b $cram > $TEMPDIR/$name.unsorted.bam

echo Sorting bam
sambamba sort --tmpdir $TEMPDIR -n -t 4 -m 8G -o $TEMPDIR/$name.qsorted.bam $TEMPDIR/$name.unsorted.bam  ## 4-threaded replacement sorting with sambamba:
rm -f $TEMPDIR/$name.unsorted.bam

echo Running bedtools bamtofastq
/usr/bin/bedtools bamtofastq -fq $TEMPDIR/$name.q.fwd.fastq -fq2 $TEMPDIR/$name.q.rev.fastq -i $TEMPDIR/$name.qsorted.bam 2>/dev/null;
rm -f $TEMPDIR/$name.qsorted.bam;

#echo step 0
#0 index the Optitype reference file:
#/usr/local/bin/bwa index $dnaref

echo Aligning forward reads to reference HLA locus sequence
/usr/local/bin/bwa mem -t 4 $dnaref $TEMPDIR/$name.q.fwd.fastq > $TEMPDIR/$name.aln.fwd.sam # use bwa mem, store output IN TEMP, and skip samse step
rm -f $TEMPDIR/$name.q.fwd.fastq;

echo Aligning reverse reads to reference HLA locus sequence
/usr/local/bin/bwa mem -t 4 $dnaref $TEMPDIR/$name.q.rev.fastq > $TEMPDIR/$name.aln.rev.sam # use bwa mem, store output IN TEMP, and skip samse step
rm -f $TEMPDIR/$name.q.rev.fastq

echo Select only the mapped reads from the sam files:
/opt/samtools/bin/samtools view -S -F 4 $TEMPDIR/$name.aln.fwd.sam > $TEMPDIR/$name.aln.map.fwd.sam
/opt/samtools/bin/samtools view -S -F 4 $TEMPDIR/$name.aln.rev.sam > $TEMPDIR/$name.aln.map.rev.sam
rm -f $TEMPDIR/$name.aln.fwd.sam
rm -f $TEMPDIR/$name.aln.rev.sam

echo Convert sam files to fastq files, also stored in temp dir
cat $TEMPDIR/$name.aln.map.fwd.sam | grep -v ^@ | /usr/bin/awk '{print "@"$1"\n"$10"\n+\n"$11}' > $outdir/$name.hla.fwd.fastq
cat $TEMPDIR/$name.aln.map.rev.sam | grep -v ^@ | /usr/bin/awk '{print "@"$1"\n"$10"\n+\n"$11}' > $outdir/$name.hla.rev.fastq
rm -f $TEMPDIR/$name.aln.map.fwd.sam
rm -f $TEMPDIR/$name.aln.map.rev.sam

echo step 5: run Optitype
# run optitype 
/usr/bin/python /usr/local/bin/OptiType/OptiTypePipeline.py -i $outdir/$name.hla.fwd.fastq $outdir/$name.hla.rev.fastq --dna -v -p $name -o $outdir
