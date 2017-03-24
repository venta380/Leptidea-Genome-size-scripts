#!/bin/bash -l

module load bioinfo-tools
module load jellyfish

mkdir $TMPDIR/jelly
output="/proj/b2014034/nobackup/jellyfish_genome_size_predict/genome_size_newdata_trimmed_K17_subsample_1/"
mkdir $output
subsample="/home/venkat/bin/sample-master/sample"
outdir=$TMPDIR/jelly
cd /proj/b2014034/nobackup/jellyfish_genome_size_predict/



re='^[0-9]+$'
for i in $(cat /proj/b2014034/nobackup/jellyfish_genome_size_predict/genome_size_all_data/subsample_data.txt | cut -f4)
do
if ! [[ $i =~ $re ]]
then
        N=$(awk '$4 == "'$i'"' /proj/b2014034/nobackup/jellyfish_genome_size_predict/genome_size_all_data/subsample_data2.txt | awk '{print  $2}')
        var=$i
        read1="/proj/b2014034/nobackup/private/annaj/trimmed_fastq/"$i"_*R1_001_val_1.fq.gz"
        read2="/proj/b2014034/nobackup/private/annaj/trimmed_fastq/"$i"_*R2_001_val_2.fq.gz"
        ls $read1
        ls $read2

        gunzip -c $read1  > $TMPDIR/jelly/reads1.$i.fq
        gunzip -c $read2  > $TMPDIR/jelly/reads2.$i.fq

        paste  -d "\t " $TMPDIR/jelly/reads1.$i.fq  $TMPDIR/jelly/reads2.$i.fq > $TMPDIR/jelly/reads.$i.fq
        rm $TMPDIR/jelly/reads1.$i.fq  $TMPDIR/jelly/reads2.$i.fq
        $subsample --lines-per-offset=4 --sample-size=$N $TMPDIR/jelly/reads.$i.fq > $TMPDIR/jelly/subsampled.$i.fq
        cut -f1 $TMPDIR/jelly/subsampled.$i.fq > $TMPDIR/jelly/subsampled1.$i.fq
        cut -f2 $TMPDIR/jelly/subsampled.$i.fq > $TMPDIR/jelly/subsampled2.$i.fq
        rm $TMPDIR/jelly/subsampled.$i.fq

else
        N=$(awk '$4 == '$i'' /proj/b2014034/nobackup/jellyfish_genome_size_predict/genome_size_all_data/subsample_data.txt | awk '{print  $2}')
        var=Sample_$i
        read1=$(ls /proj/b2014034/nobackup/private/annaj/trimmed_fastq/"$i"_*_*_R1_001_val_1.fq.gz )
        read2=$(ls /proj/b2014034/nobackup/private/annaj/trimmed_fastq/"$i"_*_*_R2_001_val_2.fq.gz )
        read3=$(ls /proj/b2014034/nobackup/private/annaj/trimmed_fastq/"$i"_*_*_R1_001_val_1.fq.gz )
        read4=$(ls /proj/b2014034/nobackup/private/annaj/trimmed_fastq/"$i"_*_*_R2_001_val_2.fq.gz )
        
        
        echo $(ls $read1| cut -f 1 -d " " | head -1)
        echo $(ls $read1| cut -f 1 -d " " | tail -1)
        echo $(ls $read2| cut -f 1 -d " " | head -1)
        echo $(ls $read2| cut -f 1 -d " " | tail -1)
        
        gunzip -c $(ls $read1| cut -f 1 -d " " | head -1)  > $TMPDIR/jelly/reads1.$i.fq
        gunzip -c $(ls $read1| cut -f 1 -d " " | tail -1)  >> $TMPDIR/jelly/reads1.$i.fq
        gunzip -c $(ls $read2| cut -f 1 -d " " | head -1)  > $TMPDIR/jelly/reads2.$i.fq
        gunzip -c $(ls $read2| cut -f 1 -d " " | tail -1)  >> $TMPDIR/jelly/reads2.$i.fq

        paste  -d "\t " $TMPDIR/jelly/reads1.$i.fq  $TMPDIR/jelly/reads2.$i.fq > $TMPDIR/jelly/reads.$i.fq
        rm $TMPDIR/jelly/reads1.$i.fq  $TMPDIR/jelly/reads2.$i.fq
        $subsample --lines-per-offset=4 --sample-size=$N $TMPDIR/jelly/reads.$i.fq > $TMPDIR/jelly/subsampled.$i.fq
        cut -f1 $TMPDIR/jelly/subsampled.$i.fq > $TMPDIR/jelly/subsampled1.$i.fq
        cut -f2 $TMPDIR/jelly/subsampled.$i.fq > $TMPDIR/jelly/subsampled2.$i.fq
        rm $TMPDIR/jelly/subsampled.$i.fq


fi


jellyfish count -t 16 -C -m 17 -s 5G  -o $outdir/genome_size.$i $TMPDIR/jelly/subsampled1.$i.fq $TMPDIR/jelly/subsampled2.$i.fq
jellyfish  histo -t 16 -o $outdir/genome_size.$i.histo $outdir/genome_size.$i
cp $outdir/genome_size.$i.histo $output


done

rm -r $TMPDIR/jelly

