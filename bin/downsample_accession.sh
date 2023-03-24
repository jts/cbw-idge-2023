#! /bin/bash

#ACC=$1
#DATA=$2 # Module3, current dir
#REPO=$3 # covid-19-signal

while getopts "a:r:" option; do
	case "${option}" in
		a) ACC=$OPTARG;;
		r) REPO=$OPTARG;;
	esac
done

mkdir -p downsample_tmp
bwa index $REPO/resources/primer_schemes/artic_v3/nCoV-2019.reference.fasta
bwa mem -t 16 $REPO/resources/primer_schemes/artic_v3/nCoV-2019.reference.fasta raw_data/${ACC}_R1.fastq.gz raw_data/${ACC}_R2.fastq.gz | samtools sort -T tmp - > downsample_tmp/$ACC.sorted.bam
samtools index downsample_tmp/$ACC.sorted.bam

python $(realpath $(dirname $0))/downsample_amplicons.py --bed $REPO/resources/primer_schemes/artic_v3/nCoV-2019.bed downsample_tmp/$ACC.sorted.bam | samtools sort -T tmp - > downsample_tmp/$ACC.downsampled.sorted.bam

mkdir -p raw_reads_downsampled
samtools sort -n downsample_tmp/$ACC.downsampled.sorted.bam | samtools fastq -1 raw_reads_downsampled/${ACC}_R1.fastq -2 raw_reads_downsampled/${ACC}_R2.fastq -0 /dev/null -n -
gzip raw_reads_downsampled/${ACC}_R1.fastq
gzip raw_reads_downsampled/${ACC}_R2.fastq

# cleanup
rm -r downsample_tmp
