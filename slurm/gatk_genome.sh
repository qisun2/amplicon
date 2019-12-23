#!/bin/bash

# index reference genome for bwa, create fasta indexes (fai and dict)
INPUT=$1
OUTPUTDIR=$2

TMP=/workdir
GATKDIR=/programs/gatk-4.1.4.0
export PATH=$GATKDIR:$PATH

[ -d $OUTPUTDIR ] || mkdir $OUTPUTDIR
cp $INPUT  $OUTPUTDIR/genome.fasta

MYCWD==$(pwd)
cd $OUTPUTDIR

# Genome summary files needed and by GATK tools
gatk CreateSequenceDictionary -R genome.fasta -O genome.dict
samtools faidx genome.fasta

# index for BWA alignment
bwa index genome.fasta

# index image file needed by some Spark-based tools (if used)
#gatk --java-options "-Djava.io.tmpdir=$TMP" BwaMemIndexImageCreator \
#     -I genome.fasta \
#     -O genome.fasta.img



