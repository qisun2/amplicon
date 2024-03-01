#!/bin/bash

fastqDir=$1
refGenome=$2
numJobs=$3

if [ ! -f $refGenome ]
then
	echo "Reference genome fasta file $refGenome not found.\n"
	exit
fi

if [ ! -f ${refGenome}.bwt ]
then
        echo "BWA index for $refGenome not found. There should be a ${refGenome}.bwt file.\n"
	exit
fi

if [ ! -d $fastqDir ]
then
        echo "Input fastq directory $fastqDir not found.\n"
        exit
fi

if [ -n $numJobs ];
then
    echo "Parellel jobs number is set to $numJobs.\n"
else
    echo "Parellel jobs number is not set. Default to 4.\n"
    numJobs=4
fi 

mkdir ./bam

ls $fastqDir/*fastq | xargs -s 1000 -I {} basename {} | sed "s/.fastq//" | xargs -s 2000 -I {} echo "bwa mem -t 2 -R \"@RG\tID:{}\tSM:{}\tLB:{}\tPL:ILLUMINA\" $refGenome $fastqDir/{}.fastq | samtools sort -O bam -o bam/{}.bam" > tasks_align

parallel -j $numJobs < tasks_align


ls -1 bam/*.bam > bamlist

bcftools mpileup -b bamlist --fasta-ref $refGenome --threads $numJobs  | bcftools call -mv  -O v -o ./unfiltered.vcf

