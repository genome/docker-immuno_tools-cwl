#!/bin/bash

# this script accepts as input a sample name, a bam file, and an output directory.  The output directory will contain a new subdirectory based on the "name"
# it should be called with a command such as the following:
# bsub -oo $stdoutfile -q long -M 32000000 -R 'select[type==LINUX64 && gtmp>100] rusage[mem=32000,gtmp=100]' -J $jobname "bash optitype.161208.sh $name $cram $dir"

#TEMPDIR=${TMPDIR}/${LSB_JOBID}.tmpdir
TEMPDIR=$1
echo $TEMPDIR

# Optitype DNA reference file

#dnaref="/gscmnt/gc6134/cancer-genomics/medseq/shared/Software/Optitype/Data/hla_reference_dna.fasta";
dnaref="/ref_data/hla_reference_dna.fasta";

# Directory into which to put the fastq files:

#newdir="/gscmnt/gc6134/cancer-genomics/medseq/shared/Software/Optitype/Output";

name=$2;
cram=$3;
dir=$4;
echo $name
echo $cram

# Directory into which Optitype results should be put:
outdir=$dir/$name;
mkdir -p $outdir

echo converting cram to bam
# convert cram to bam:
/opt/samtools/bin/samtools view -b $cram > $TEMPDIR/$name.unsorted.bam

echo sorting bam
# Sort the bam:
sambamba sort --tmpdir $TEMPDIR -n -t 4 -m 8G -o $TEMPDIR/$name.qsorted.bam $TEMPDIR/$name.unsorted.bam  ## 4-threaded replacement sorting with sambamba:
#rm -f $TEMPDIR/$name.unsorted.bam

# convert bam to two fastq files:

echo running bedtools bamtofastq
/usr/bin/bedtools bamtofastq -fq $TEMPDIR/$name.q.fwd.fastq -fq2 $TEMPDIR/$name.q.rev.fastq -i $TEMPDIR/$name.qsorted.bam 2>/dev/null;
#ls -l $TEMPDIR/$name.q.fwd.fastq
#ls -l $TEMPDIR/$name.q.rev.fastq
#rm -f $TEMPDIR/$name.qsorted.bam;

#echo step 0
#0 index the Optitype reference file:
#/usr/local/bin/bwa index $dnaref

echo step 1
#1 Align forward reads to reference HLA locus sequence
/usr/local/bin/bwa mem -t 4 $dnaref $TEMPDIR/$name.q.fwd.fastq > $TEMPDIR/$name.aln.fwd.sam # use bwa mem, store output IN TEMP, and skip samse step
#rm -f $TEMPDIR/$name.q.fwd.fastq;

echo step 2
#2 Align reverse reads to reference HLA locus sequence
/usr/local/bin/bwa mem -t 4 $dnaref $TEMPDIR/$name.q.rev.fastq > $TEMPDIR/$name.aln.rev.sam # use bwa mem, store output IN TEMP, and skip samse step
#rm -f $TEMPDIR/$name.q.rev.fastq

echo step 3
# select only the mapped reads from the sam files:
/opt/samtools/bin/samtools view -S -F 4 $TEMPDIR/$name.aln.fwd.sam > $TEMPDIR/$name.aln.map.fwd.sam
/opt/samtools/bin/samtools view -S -F 4 $TEMPDIR/$name.aln.rev.sam > $TEMPDIR/$name.aln.map.rev.sam
#echo output step 3:
#head $TEMPDIR/$name.aln.fwd.sam
#ls -l $TEMPDIR/$name.aln.map.fwd.sam
#ls -l $TEMPDIR/$name.aln.map.rev.sam

#rm -f $TEMPDIR/$name.aln.fwd.sam
#rm -f $TEMPDIR/$name.aln.rev.sam

echo step 4
# convert sam files to fastq files, also stored in temp dir
#touch $outdir/$name.hla.fwd.fastq
#touch $outdir/$name.hla.rev.fastq
cat $TEMPDIR/$name.aln.map.fwd.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > $outdir/$name.hla.fwd.fastq
cat $TEMPDIR/$name.aln.map.rev.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > $outdir/$name.hla.rev.fastq
#rm -f $TEMPDIR/$name.aln.map.fwd.sam
#rm -f $TEMPDIR/$name.aln.map.rev.sam

echo step 5: run Optitype
# run optitype 
#PATH=/usr/local/bin/OptiType:/gscuser/apetti/.Renv/shims:/gscuser/apetti/.Renv/bin:/gscuser/apetti/bin:/gapp/noarch/bin:/gapp/x64linux/bin:/gapp/ia32linux/bin:/gsc/scripts/opt/genome/current/user/bin:/gsc/scripts/bin:/gsc/scripts/opt/lims/snapshots/current/bin:/gsc/bin:/gsc/java/bin:/opt/lsf9/9.1/linux2.6-glibc2.3-x86_64/etc:/opt/lsf9/9.1/linux2.6-glibc2.3-x86_64/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/gsc/scripts/gsc/info:/gsc/scripts/gsc/medseq:/gsc/pkg/oracle/10gR2/db_1/bin:/gscuser/apetti/meme/bin:/gscuser/apetti/usr/bin

#LSB_SUB_ADDITIONAL='docker(fred2/optitype)' bsub -oo $outdir/out.err -q research-hpc -J Optitype -M 32000000 -R 'select[type==LINUX64 && mem>32000] rusage[mem=32000]' -- -i $outdir/$name.hla.fwd.fastq $outdir/$name.hla.rev.fastq --dna -v --outdir $outdir # using space in network for fastq files
python /usr/local/bin/OptiType/OptiTypePipeline.py -i $outdir/$name.hla.fwd.fastq $outdir/$name.hla.rev.fastq --dna -v --outdir $outdir
