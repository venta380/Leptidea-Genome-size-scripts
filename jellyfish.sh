#!/bin/bash -l

mkdir $TMPDIR/jelly
mkdir "/proj/b2014034/nobackup/jellyfish_genome_size_predict/genome_size_17/"
outdir=$TMPDIR/jelly
cd /proj/b2014034/nobackup/jellyfish_genome_size_predict/


for i in {1..50}
do
var=Sample_$i
file=./paths.txt
if grep -q -w -F $var $file;
then
read1="/proj/b2014034/private/reseq_analysis/trimmed_fastq/"$i"_*R1_001_val_1.fq.gz"
read2="/proj/b2014034/private/reseq_analysis/trimmed_fastq/"$i"_*R2_001_val_2.fq.gz"
gunzip -c $read1  > $TMPDIR/jelly/reads1.$i.fq
gunzip -c $read2  > $TMPDIR/jelly/reads2.$i.fq
jellyfish count -t 16 -C -m 17 -s 5G  -o $outdir/genome_size.$i $TMPDIR/jelly/reads1.$i.fq $TMPDIR/jelly/reads2.$i.fq
jellyfish  histo -t 16 -o $outdir/genome_size.$i.histo $outdir/genome_size.$i
cp $outdir/genome_size.$i.histo /proj/b2014034/nobackup/jellyfish_genome_size_predict/genome_size_17/
rm $TMPDIR/jelly/reads1.$i.fq
rm $TMPDIR/jelly/reads2.$i.fq
fi
done

rm -r $TMPDIR/jelly
