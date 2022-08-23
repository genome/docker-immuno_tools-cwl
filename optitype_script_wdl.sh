#!/bin/bash

# $1 is the directory where all intermediate files will be written
# $2 is the directory final outputs should be written to
# $3 is the prefix that will be added to all files
# $4 is the location of the input (normal) cram file
# $5 is the reference file used to create the cram
# $6 is the number of cores that will be used by processes
# $7 is the total memory in the node
#bsub -R 'select[mem>64000] rusage[mem=64000]' -e <error_file_name> -o <output_file_name> -q research-hpc -a 'docker(johnegarza/immuno-testing:latest)' /bin/bash /usr/bin/optitype_script.sh <intermediate files directory> <final results directory> <output file prefix>  <cram path>

set -e -o pipefail

# Optitype DNA reference file
dnaref="/ref_data/optitype_ref/hla_reference_dna.fasta";

DEFAULT_THREADS=4
DEFAULT_MEM=8
MEM_UTIL=80 # percent

TEMPDIR="$1";
outdir="$2";
name="$3";
cram="$4";
reference="$5";
THREADS=${6:-$DEFAULT_THREADS}
if [ ${7:-$DEFAULT_MEM} -eq 1 ]; then MEM=1G; else MEM="$((${7:-$DEFAULT_MEM}*$MEM_UTIL/100))"G; fi

mkdir -p $TEMPDIR
mkdir -p $outdir

echo "Step 1: extracting hla region from chr6, reads to alternate HLA sequences, and all unmapped reads ..."

# filter out reads directly from an existing CRAM file of alignments, only those reads that align to this region: chr6:29836259-33148325
echo "[INFO] /opt/samtools/bin/samtools view -h -T $reference $cram chr6:29836259-33148325 >$TEMPDIR/reads.sam"
time /opt/samtools/bin/samtools view -h -T $reference $cram chr6:29836259-33148325 >$TEMPDIR/reads.sam

# pull out only the *header* lines from the CRAM with the -H parameter. then get the sequence names and for those that match the string "HLA" do the following
echo "[INFO] /opt/samtools/bin/samtools view -H -T $reference $cram | grep "^@SQ" | cut -f 2 | cut -f 2- -d : | grep HLA | while read chr;do ..."
time /opt/samtools/bin/samtools view -H -T $reference $cram | grep "^@SQ" | cut -f 2 | cut -f 2- -d : | grep HLA | while read chr;do 
# echo "checking $chr:1-9999999"
# grab all the reads that align to each alternate "HLA" sequence
/opt/samtools/bin/samtools view -T $reference $cram "$chr:1-9999999" >>$TEMPDIR/reads.sam
done

# grab all the reads that are unaligned
echo "[INFO]  /opt/samtools/bin/samtools view -f 4 -T $reference $cram >>$TEMPDIR/reads.sam"
time /opt/samtools/bin/samtools view -f 4 -T $reference $cram >>$TEMPDIR/reads.sam
# covert from .sam to .bam format
echo "[INFO] /opt/samtools/bin/samtools view -Sb -o $TEMPDIR/reads.bam $TEMPDIR/reads.sam"
time /opt/samtools/bin/samtools view -Sb -o $TEMPDIR/reads.bam $TEMPDIR/reads.sam 

echo "Step 2: running picard ..."
echo "[INFO] /usr/bin/java -Xmx6g -jar /usr/picard/picard.jar SamToFastq VALIDATION_STRINGENCY=LENIENT F=$TEMPDIR/$name.q.fwd.fastq.gz F2=$TEMPDIR/$name.q.rev.fastq.gz I=$TEMPDIR/reads.bam R=$reference FU=$TEMPDIR/unpaired.fastq.gz"
time /usr/bin/java -Xmx6g -jar /usr/picard/picard.jar SamToFastq VALIDATION_STRINGENCY=LENIENT F=$TEMPDIR/$name.q.fwd.fastq.gz F2=$TEMPDIR/$name.q.rev.fastq.gz I=$TEMPDIR/reads.bam R=$reference FU=$TEMPDIR/unpaired.fastq.gz


#echo step 0
#0 index the Optitype reference file:
#/usr/local/bin/bwa index $dnaref

echo "Step 3: Aligning forward reads to reference HLA locus sequence"
echo "[INFO] /usr/local/bin/bwa mem -t $THREADS $dnaref $TEMPDIR/$name.q.fwd.fastq.gz > $TEMPDIR/$name.aln.fwd.sam"
time /usr/local/bin/bwa mem -t $THREADS $dnaref $TEMPDIR/$name.q.fwd.fastq.gz > $TEMPDIR/$name.aln.fwd.sam # use bwa mem, store output IN TEMP, and skip samse step
rm -f $TEMPDIR/$name.q.fwd.fastq.gz;

echo "Step 4: Aligning reverse reads to reference HLA locus sequence"
echo "[INFO] /usr/local/bin/bwa mem -t $THREADS $dnaref $TEMPDIR/$name.q.rev.fastq.gz > $TEMPDIR/$name.aln.rev.sam"
time /usr/local/bin/bwa mem -t $THREADS $dnaref $TEMPDIR/$name.q.rev.fastq.gz > $TEMPDIR/$name.aln.rev.sam # use bwa mem, store output IN TEMP, and skip samse step
rm -f $TEMPDIR/$name.q.rev.fastq.gz

echo "Step 5: Select only the mapped reads from the sam files:"
echo "[INFO] /opt/samtools/bin/samtools view -S -F 4 $TEMPDIR/$name.aln.fwd.sam > $TEMPDIR/$name.aln.map.fwd.sam"
time /opt/samtools/bin/samtools view -S -F 4 $TEMPDIR/$name.aln.fwd.sam > $TEMPDIR/$name.aln.map.fwd.sam
echo "[INFO] /opt/samtools/bin/samtools view -S -F 4 $TEMPDIR/$name.aln.rev.sam > $TEMPDIR/$name.aln.map.rev.sam"
time /opt/samtools/bin/samtools view -S -F 4 $TEMPDIR/$name.aln.rev.sam > $TEMPDIR/$name.aln.map.rev.sam
rm -f $TEMPDIR/$name.aln.fwd.sam
rm -f $TEMPDIR/$name.aln.rev.sam

echo "Step 6: Convert sam files to fastq files, also stored in temp dir"
time cat $TEMPDIR/$name.aln.map.fwd.sam | grep -v ^@ | /usr/bin/awk '{print "@"$1"\n"$10"\n+\n"$11}' > $outdir/$name.hla.fwd.fastq
time cat $TEMPDIR/$name.aln.map.rev.sam | grep -v ^@ | /usr/bin/awk '{print "@"$1"\n"$10"\n+\n"$11}' > $outdir/$name.hla.rev.fastq
rm -f $TEMPDIR/$name.aln.map.fwd.sam
rm -f $TEMPDIR/$name.aln.map.rev.sam

echo "Step 7: run Optitype"
# run optitype
echo "[INFO] /usr/bin/python /usr/local/bin/OptiType/OptiTypePipeline.py -i $outdir/$name.hla.fwd.fastq $outdir/$name.hla.rev.fastq --dna -v -p $name -o $outdir" 
time /usr/bin/python /usr/local/bin/OptiType/OptiTypePipeline.py -i $outdir/$name.hla.fwd.fastq $outdir/$name.hla.rev.fastq --dna -v -p $name -o $outdir
