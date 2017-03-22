module load bioinfo-tools
module load jellyfish
module load seqtk

mkdir $TMPDIR/jelly
output="/proj/b2014034/nobackup/jellyfish_genome_size_predict/genome_size_newdata_trimmed_K17_subsample_2/"
mkdir $output
subsample="/home/venkat/bin/sample-master/sample"
outdir=$TMPDIR/jelly
cd /proj/b2014034/nobackup/jellyfish_genome_size_predict/



re='^[0-9]+$'
for i in $(cat /proj/b2014034/nobackup/jellyfish_genome_size_predict/genome_size_all_data/subsample_data.txt | cut -f4)
do
if ! [[ $i =~ $re ]]
then
        N=$(awk '$4 == "'$i'"' /proj/b2014034/nobackup/jellyfish_genome_size_predict/reads_set2.txt | awk '{print  $2}')
        var=$i
        read1="/proj/b2014034/private/reseq_analysis/trimmed_fastq/"$i"_*R1_001_val_1.fq.gz"
        read2="/proj/b2014034/private/reseq_analysis/trimmed_fastq/"$i"_*R2_001_val_2.fq.gz"

        gunzip -c $read1  > $TMPDIR/jelly/reads1.$i.fq 
        gunzip -c $read2  > $TMPDIR/jelly/reads2.$i.fq 


        seqtk sample -s100 $TMPDIR/jelly/reads1.$i.fq $N > $TMPDIR/jelly/subsampled1.$i.fq
		seqtk sample -s100 $TMPDIR/jelly/reads2.$i.fq $N > TMPDIR/jelly/subsampled2.$i.fq
		rm $TMPDIR/jelly/reads1.$i.fq  $TMPDIR/jelly/reads2.$i.fq
else
        N=$(awk '$4 == '$i'' /proj/b2014034/nobackup/jellyfish_genome_size_predict/reads_set2.txt | awk '{print  $2}')
        var=Sample_$i
        read1=$(ls /proj/b2014034/private/reseq_analysis/trimmed_fastq/"$i"_*R1_001_val_1.fq.gz | head -1)
        read2=$(ls /proj/b2014034/private/reseq_analysis/trimmed_fastq/"$i"_*R2_001_val_2.fq.gz | head -1)
        read3=$(ls /proj/b2014034/private/reseq_analysis/trimmed_fastq/"$i"_*R1_001_val_1.fq.gz | tail -1)
        read4=$(ls /proj/b2014034/private/reseq_analysis/trimmed_fastq/"$i"_*R2_001_val_2.fq.gz | tail -1)

        gunzip -c $read1  > $TMPDIR/jelly/reads1.$i.fq
        gunzip -c $read2  > $TMPDIR/jelly/reads2.$i.fq
        gunzip -c $read3  >> $TMPDIR/jelly/reads1.$i.fq
        gunzip -c $read4  >> $TMPDIR/jelly/reads2.$i.fq

		seqtk sample -s100 $TMPDIR/jelly/reads1.$i.fq $N > $TMPDIR/jelly/subsampled1.$i.fq
		seqtk sample -s100 $TMPDIR/jelly/reads2.$i.fq $N > TMPDIR/jelly/subsampled2.$i.fq
		rm $TMPDIR/jelly/reads1.$i.fq  $TMPDIR/jelly/reads2.$i.fq

        
fi


jellyfish count -t 16 -C -m 17 -s 5G  -o $outdir/genome_size.$i $TMPDIR/jelly/subsampled1.$i.fq $TMPDIR/jelly/subsampled2.$i.fq
jellyfish  histo -t 16 -o $outdir/genome_size.$i.histo $outdir/genome_size.$i
cp $outdir/genome_size.$i.histo $output


done

rm -r $TMPDIR/jelly

