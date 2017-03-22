#!/bin/bash -l

#SBATCH -A b2013146
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 6:00:00
#SBATCH -J read_details_old_data

python_script="/proj/b2014034/nobackup/jellyfish_genome_size_predict/avg.py"

mkdir $TMPDIR/jelly/

cd /proj/b2014034/nobackup/jellyfish_genome_size_predict/

for i in 1 2 3 6 7 8 9 10 11 12 13 14 15 17 18 20 21 22 23 24 26 27 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 46 47

do
var=Sample_$i
read1="/proj/b2014034/private/raw_data/Population_resequencing_data/140221_SN7001335_0100_AC3MG8ACXX/$var/*R1_001.fastq.gz"

ls $read1
gunzip -c $read1  > $TMPDIR/jelly/reads1.$i.fq

awk '{if(NR%4==2) print length($0)}' $TMPDIR/jelly/reads1.$i.fq > $TMPDIR/jelly/reads1.$i.readlengths

echo $i  $(cat $TMPDIR/jelly/reads1.$i.fq | wc -l) $(python2.7 $python_script  $TMPDIR/jelly/reads1.$i.readlengths ) >> genome_size/reads.txt

rm $TMPDIR/jelly/reads1.$i.fq
rm $TMPDIR/jelly/reads1.$i.readlengths


done

rm -r $TMPDIR/jelly/
